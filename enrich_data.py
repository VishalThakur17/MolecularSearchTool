"""
Enrich existing proteins and binders with more detailed data.

Run this after your normal ingestion:
    python enrich_existing_data.py

What it does:
1. Re-fetches richer UniProt fields for every gene already in your database.
2. Searches ChEMBL for existing binder names and fills missing status, modality,
   molecular weight, SMILES, and mechanism where available.
3. Applies a small curated demo override table for common therapeutics so your
   UI does not show avoidable "Not available" values for well-known examples.
"""

from __future__ import annotations

import re
import requests
import psycopg2.extras

from ingest.db import get_connection
from ingest.ingest_uniprot import upsert_gene_and_protein
from ingest.ingest_binders import (
    fetch_molecule_details,
    infer_binder_type,
    infer_clinical_status,
    infer_modality,
    extract_sequence,
    get_or_create_modality,
)

CHEMBL_MOLECULE_SEARCH_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"

# These overrides are intentionally small and transparent. They are useful for
# demo quality and for entries where public APIs return sparse commercial data.
CURATED_BINDER_OVERRIDES = {
    "osimertinib": {
        "binder_type": "Small Molecule",
        "modality_name": "Small Molecule",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "AstraZeneca",
        "mechanism_of_action": "Third-generation EGFR tyrosine kinase inhibitor used for EGFR-mutated non-small cell lung cancer, including common activating mutations and T790M resistance mutations.",
        "binder_description": "Osimertinib is a small-molecule EGFR inhibitor. It is retained as therapeutic context, but it is separated from the sponsor protein-binder taxonomy.",
        "sequence": None,
    },
    "trastuzumab": {
        "binder_type": "IgG",
        "modality_name": "Antibody",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Genentech / Roche",
        "mechanism_of_action": "Humanized monoclonal antibody targeting HER2/ERBB2.",
        "binder_description": "Trastuzumab is a therapeutic monoclonal antibody that binds HER2/ERBB2 and is used in HER2-positive cancers.",
    },
    "trastuzumab emtansine": {
        "binder_type": "ADC",
        "modality_name": "Antibody-Drug Conjugate",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Genentech / Roche",
        "mechanism_of_action": "HER2-directed antibody-drug conjugate that combines trastuzumab with the cytotoxic payload DM1.",
        "binder_description": "Trastuzumab emtansine is an antibody-drug conjugate used for HER2-positive cancers.",
    },
    "trastuzumab deruxtecan": {
        "binder_type": "ADC",
        "modality_name": "Antibody-Drug Conjugate",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Daiichi Sankyo / AstraZeneca",
        "mechanism_of_action": "HER2-directed antibody-drug conjugate with a topoisomerase I inhibitor payload.",
        "binder_description": "Trastuzumab deruxtecan is a HER2-directed ADC used across HER2-positive and HER2-low tumor settings.",
    },
    "amivantamab": {
        "binder_type": "Bispecific Antibody",
        "modality_name": "Bispecific Antibody",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Janssen",
        "mechanism_of_action": "Bispecific antibody targeting EGFR and MET.",
        "binder_description": "Amivantamab is a bispecific antibody that binds EGFR and MET and is used in selected lung cancer settings.",
    },
    "blinatumomab": {
        "binder_type": "Bispecific Antibody",
        "modality_name": "Bispecific Antibody",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Amgen",
        "mechanism_of_action": "Bispecific T-cell engager targeting CD19 and CD3.",
        "binder_description": "Blinatumomab is a bispecific T-cell engager that links CD19-positive cells to CD3-positive T cells.",
    },
    "caplacizumab": {
        "binder_type": "VHH",
        "modality_name": "Nanobody / VHH",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Sanofi",
        "mechanism_of_action": "Nanobody targeting von Willebrand factor.",
        "binder_description": "Caplacizumab is a VHH/nanobody therapeutic that targets von Willebrand factor.",
    },
    "aflibercept": {
        "binder_type": "Fc Fusion",
        "modality_name": "Fc Fusion",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Regeneron / Bayer",
        "mechanism_of_action": "VEGF trap fusion protein that binds VEGF ligands.",
        "binder_description": "Aflibercept is an Fc fusion protein that acts as a soluble decoy receptor for VEGF pathway ligands.",
    },
    "romiplostim": {
        "binder_type": "Fc Fusion",
        "modality_name": "Fc Fusion",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Amgen",
        "mechanism_of_action": "Peptibody/Fc fusion that activates the thrombopoietin receptor pathway.",
        "binder_description": "Romiplostim is an Fc-containing peptibody therapeutic.",
    },
    "cilengitide": {
        "binder_type": "Peptide",
        "modality_name": "Peptide",
        "clinical_status": "Clinical / Investigational",
        "approval_status": "Investigational",
        "developer_company": "Merck KGaA",
        "mechanism_of_action": "Cyclic RGD peptide inhibitor of selected integrins.",
        "binder_description": "Cilengitide is a cyclic peptide therapeutic candidate that binds integrins.",
    },
    "envafolimab": {
        "binder_type": "VHH",
        "modality_name": "Nanobody / VHH",
        "clinical_status": "Approved / Clinical",
        "approval_status": "Approved or clinical use depending on region/indication",
        "developer_company": "Alphamab Oncology / 3D Medicines",
        "mechanism_of_action": "Single-domain antibody targeting PD-L1.",
        "binder_description": "Envafolimab is a single-domain antibody/nanobody-style therapeutic targeting PD-L1.",
    },
    "ozoralizumab": {
        "binder_type": "VHH",
        "modality_name": "Nanobody / VHH",
        "clinical_status": "Approved / Clinical",
        "approval_status": "Approved or clinical use depending on region/indication",
        "developer_company": "Taisho / Ablynx",
        "mechanism_of_action": "Trivalent NANOBODY therapeutic targeting TNF-alpha and albumin.",
        "binder_description": "Ozoralizumab is a VHH/nanobody therapeutic included to improve representation of the VHH class.",
    },
    "zanidatamab": {
        "binder_type": "Bispecific Antibody",
        "modality_name": "Bispecific Antibody",
        "clinical_status": "Approved / Clinical",
        "approval_status": "Approved or late-stage clinical depending on region/indication",
        "developer_company": "Jazz Pharmaceuticals / Zymeworks",
        "mechanism_of_action": "Bispecific HER2-directed antibody.",
        "binder_description": "Zanidatamab is a HER2-directed bispecific antibody.",
    },
    "faricimab": {
        "binder_type": "Bispecific Antibody",
        "modality_name": "Bispecific Antibody",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Roche / Genentech",
        "mechanism_of_action": "Bispecific antibody targeting VEGF-A and Angiopoietin-2.",
        "binder_description": "Faricimab is a bispecific antibody included as a clean example of the bispecific class.",
    },
    "etanercept": {
        "binder_type": "Fc Fusion",
        "modality_name": "Fc Fusion",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Amgen / Pfizer",
        "mechanism_of_action": "TNF receptor Fc fusion protein that binds TNF.",
        "binder_description": "Etanercept is a receptor-Fc fusion therapeutic.",
    },
    "abatacept": {
        "binder_type": "Fc Fusion",
        "modality_name": "Fc Fusion",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Bristol Myers Squibb",
        "mechanism_of_action": "CTLA-4 Fc fusion protein that modulates T-cell costimulation.",
        "binder_description": "Abatacept is an Fc fusion protein included to improve Fc fusion coverage.",
    },
    "leuprolide": {
        "binder_type": "Peptide",
        "modality_name": "Peptide",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Multiple manufacturers",
        "mechanism_of_action": "Peptide GnRH receptor agonist.",
        "binder_description": "Leuprolide is a peptide therapeutic included as a clear peptide-class example.",
    },
    "goserelin": {
        "binder_type": "Peptide",
        "modality_name": "Peptide",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "AstraZeneca",
        "mechanism_of_action": "Peptide GnRH receptor agonist.",
        "binder_description": "Goserelin is a peptide therapeutic included as a clear peptide-class example.",
    },
}


def phase_to_approval(max_phase):
    if max_phase is None:
        return None
    try:
        phase = int(float(max_phase))
    except Exception:
        return f"Max Phase {max_phase}"
    if phase >= 4:
        return "Approved / Max Phase 4"
    return f"Max Phase {phase}"


def clean_key(name: str | None) -> str:
    return re.sub(r"\s+", " ", (name or "").strip().lower())

SMALL_MOLECULE_NAME_HINTS = [
    "inib", "tinib", "rafenib", "metinib", "parib", "sotorasib", "adagrasib",
    "osimertinib", "gefitinib", "erlotinib", "afatinib", "lapatinib",
    "neratinib", "tucatinib", "crizotinib", "alectinib", "brigatinib",
    "lorlatinib", "capmatinib", "tepotinib", "selpercatinib", "pralsetinib",
    "fulvestrant", "tamoxifen", "anastrozole", "letrozole", "exemestane",
]


def infer_type_from_name_only(name: str) -> dict:
    """Fallback for records where ChEMBL lookup is sparse or fails."""
    key = clean_key(name)
    updates = {}

    if any(hint in key for hint in SMALL_MOLECULE_NAME_HINTS):
        updates.update({
            "binder_type": "Small Molecule",
            "modality_name": "Small Molecule",
            "sequence": None,
            "binder_description": f"{name} appears to be a small-molecule therapeutic based on its name pattern. It is stored as therapeutic context rather than a sponsor protein-binder class.",
        })

    return updates



def find_chembl_molecule_by_name(name: str) -> dict | None:
    params = {"pref_name__iexact": name, "limit": 1}
    response = requests.get(CHEMBL_MOLECULE_SEARCH_URL, params=params, timeout=45)
    response.raise_for_status()
    molecules = response.json().get("molecules", [])
    if not molecules:
        params = {"molecule_synonyms__molecule_synonym__iexact": name, "limit": 1}
        response = requests.get(CHEMBL_MOLECULE_SEARCH_URL, params=params, timeout=45)
        response.raise_for_status()
        molecules = response.json().get("molecules", [])
    if not molecules:
        return None
    chembl_id = molecules[0].get("molecule_chembl_id")
    if not chembl_id:
        return molecules[0]
    return fetch_molecule_details(chembl_id)


def get_mechanism_for_molecule(chembl_id: str | None) -> str | None:
    if not chembl_id:
        return None
    url = "https://www.ebi.ac.uk/chembl/api/data/mechanism.json"
    params = {"molecule_chembl_id": chembl_id, "limit": 5}
    response = requests.get(url, params=params, timeout=45)
    response.raise_for_status()
    mechanisms = response.json().get("mechanisms", [])
    values = []
    for item in mechanisms:
        moa = item.get("mechanism_of_action")
        action = item.get("action_type")
        target = item.get("target_name")
        pieces = [p for p in [moa, action, target] if p]
        if pieces:
            values.append("; ".join(pieces))
    return " | ".join(values) if values else None


def enrich_proteins():
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT DISTINCT gene_symbol FROM genes WHERE gene_symbol IS NOT NULL ORDER BY gene_symbol;")
            genes = [row[0] for row in cur.fetchall()]

    print(f"\n=== Enriching proteins from UniProt ({len(genes)} genes) ===")
    for gene in genes:
        try:
            upsert_gene_and_protein(gene)
        except Exception as exc:
            print(f"[Protein Enrichment] Failed for {gene}: {exc}")


def enrich_binders():
    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("SELECT binder_id, binder_name FROM binders ORDER BY binder_name;")
            binders = cur.fetchall()

            print(f"\n=== Enriching binders from ChEMBL/curated overrides ({len(binders)} binders) ===")

            for binder in binders:
                binder_id = binder["binder_id"]
                binder_name = binder["binder_name"]
                key = clean_key(binder_name)

                updates = {}
                molecule = None

                try:
                    molecule = find_chembl_molecule_by_name(binder_name)
                except Exception as exc:
                    print(f"[Binder Enrichment] ChEMBL lookup failed for {binder_name}: {exc}")

                if molecule:
                    max_phase = molecule.get("max_phase")
                    molecule_structures = molecule.get("molecule_structures") or {}
                    molecule_properties = molecule.get("molecule_properties") or {}

                    updates["binder_type"] = infer_binder_type(molecule)
                    updates["clinical_status"] = infer_clinical_status(molecule)
                    updates["approval_status"] = phase_to_approval(max_phase)
                    updates["sequence"] = extract_sequence(molecule)
                    updates["smiles"] = molecule_structures.get("canonical_smiles")
                    updates["molecular_weight"] = molecule_properties.get("full_mwt")

                    try:
                        updates["mechanism_of_action"] = get_mechanism_for_molecule(molecule.get("molecule_chembl_id"))
                    except Exception:
                        updates["mechanism_of_action"] = None

                    updates["binder_description"] = (
                        molecule.get("pref_name") and f"Imported/enriched from ChEMBL molecule record for {molecule.get('pref_name')}."
                    )

                    modality_name = infer_modality(molecule)
                    updates["modality_id"] = get_or_create_modality(cur, modality_name)

                # If ChEMBL did not return a useful result, still catch obvious
                # small-molecule drug names so they do not remain as "Other".
                if not updates:
                    updates.update(infer_type_from_name_only(binder_name))

                override = CURATED_BINDER_OVERRIDES.get(key)
                if override:
                    override = override.copy()
                    modality_name = override.pop("modality_name", None)
                    if modality_name:
                        override["modality_id"] = get_or_create_modality(cur, modality_name)
                    # Curated values win when provided because they are more readable for demo/report use.
                    updates.update({k: v for k, v in override.items() if v is not None})

                # Convert fallback modality_name to modality_id if needed.
                fallback_modality = updates.pop("modality_name", None)
                if fallback_modality and "modality_id" not in updates:
                    updates["modality_id"] = get_or_create_modality(cur, fallback_modality)

                if not updates:
                    continue

                allowed_cols = [
                    "binder_type", "sequence", "clinical_status", "binder_description",
                    "mechanism_of_action", "approval_status", "developer_company",
                    "smiles", "molecular_weight", "modality_id",
                ]
                assignments = []
                values = []
                for col in allowed_cols:
                    if col in updates and updates[col] is not None:
                        assignments.append(f"{col} = %s")
                        values.append(updates[col])

                if not assignments:
                    continue

                values.append(binder_id)
                cur.execute(
                    f"UPDATE binders SET {', '.join(assignments)} WHERE binder_id = %s;",
                    values,
                )
                print(f"[Binder Enrichment] Updated {binder_name}")

        conn.commit()


def normalize_remaining_other_records():
    print("\n=== Normalizing remaining 'Other' binders ===")

    conn = get_connection()
    cur = conn.cursor()

    cur.execute("""
        SELECT binder_id, binder_name
        FROM binders
        WHERE binder_type = 'Other';
    """)

    rows = cur.fetchall()

    if not rows:
        print("No 'Other' binders to normalize.")
        cur.close()
        conn.close()
        return

    updated_count = 0

    for row in rows:
        try:
            # Handle dict cursor
            if isinstance(row, dict):
                binder_id = row.get("binder_id")
                binder_name = row.get("binder_name", "")
            else:
                # Handle tuple safely
                if not row or len(row) < 2:
                    continue
                binder_id = row[0]
                binder_name = row[1]

            if not binder_id or not binder_name:
                continue

            name_lower = binder_name.lower()

            # Detect small molecules (robust rules)
            if any(keyword in name_lower for keyword in [
                "tinib", "statin", "azole", "ib", "nib", "zole",
                "acid", "ate", "ol", "one"
            ]):
                cur.execute("""
                    UPDATE binders
                    SET binder_type = 'Small Molecule'
                    WHERE binder_id = %s;
                """, (binder_id,))
                updated_count += 1

        except Exception as e:
            print(f"[WARNING] Skipping row due to error: {e}")
            continue

    conn.commit()
    cur.close()
    conn.close()

    print(f"Normalized {updated_count} binders from 'Other' → 'Small Molecule'")


def main():
    enrich_proteins()
    enrich_binders()
    normalize_remaining_other_records()
    print("\n=== ENRICHMENT COMPLETE ===")


if __name__ == "__main__":
    main()
