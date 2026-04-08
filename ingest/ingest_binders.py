import requests
from ingest.db import get_connection

CHEMBL_MECHANISM_URL = "https://www.ebi.ac.uk/chembl/api/data/mechanism.json"
CHEMBL_MOLECULE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule"


def fetch_mechanisms_for_gene(gene_symbol: str, limit: int = 15):
    """
    Fetch ChEMBL mechanism records filtered by target gene symbol.
    """
    params = {
        "target_components__accession": gene_symbol,
        "limit": limit,
    }

    response = requests.get(CHEMBL_MECHANISM_URL, params=params, timeout=45)
    response.raise_for_status()
    data = response.json()
    return data.get("mechanisms", [])


def fetch_molecule_details(chembl_molecule_id: str):
    url = f"{CHEMBL_MOLECULE_URL}/{chembl_molecule_id}.json"
    response = requests.get(url, timeout=45)
    response.raise_for_status()
    return response.json()


def infer_modality(molecule_record: dict) -> str:
    """
    Existing simple modality inference for alpha compatibility.
    """
    molecule_type = (molecule_record.get("molecule_type") or "").lower()

    if "antibody" in molecule_type:
        return "Antibody"
    if "protein" in molecule_type or "oligo" in molecule_type or "peptide" in molecule_type:
        return "Peptide"
    return "Small Molecule"


def infer_binder_type(molecule_record: dict) -> str:
    """
    New sponsor-aligned binder type classification.
    MVP target types:
    - IgG
    - VHH
    - Peptide

    For now, keep the rules simple and safe.
    """
    molecule_type = (molecule_record.get("molecule_type") or "").lower()
    pref_name = (molecule_record.get("pref_name") or "").lower()

    if "nanobody" in molecule_type or "vhh" in molecule_type or "nanobody" in pref_name or "vhh" in pref_name:
        return "VHH"
    if "peptide" in molecule_type or "oligopeptide" in molecule_type:
        return "Peptide"
    if "antibody" in molecule_type or "protein" in molecule_type:
        return "IgG"

    return "Other"


def infer_clinical_status(molecule_record: dict) -> str | None:
    """
    Convert ChEMBL max_phase into a human-readable clinical status.
    """
    max_phase = molecule_record.get("max_phase")

    phase_map = {
        0: "Preclinical",
        1: "Phase 1",
        2: "Phase 2",
        3: "Phase 3",
        4: "Approved",
    }

    if max_phase is None:
        return None

    return phase_map.get(max_phase, f"Max Phase {max_phase}")


def extract_sequence(molecule_record: dict) -> str | None:
    """
    Placeholder for sequence extraction.

    ChEMBL small molecules usually will not have a FASTA sequence.
    Some biologics may require a different endpoint or richer parsing.
    For task 1, we safely store NULL when sequence is unavailable.
    """
    sequence = molecule_record.get("sequence")
    if sequence:
        return sequence.strip()

    return None


def upsert_source(cur, source_name: str, source_url: str):
    cur.execute(
        """
        INSERT INTO sources (source_name, source_url)
        VALUES (%s, %s)
        ON CONFLICT (source_name)
        DO UPDATE SET source_url = EXCLUDED.source_url
        RETURNING source_id;
        """,
        (source_name, source_url),
    )
    return cur.fetchone()[0]


def get_or_create_modality(cur, modality_name: str):
    cur.execute(
        """
        INSERT INTO binder_modalities (modality_name)
        VALUES (%s)
        ON CONFLICT (modality_name)
        DO UPDATE SET modality_name = EXCLUDED.modality_name
        RETURNING modality_id;
        """,
        (modality_name,),
    )
    return cur.fetchone()[0]


def get_protein_id_by_gene(cur, gene_symbol: str):
    cur.execute(
        """
        SELECT p.protein_id
        FROM proteins p
        JOIN genes g ON p.gene_id = g.gene_id
        WHERE g.gene_symbol = %s
        LIMIT 1;
        """,
        (gene_symbol,),
    )
    row = cur.fetchone()
    return row[0] if row else None


def upsert_binders_for_gene(gene_symbol: str, limit: int = 15):
    """
    Pull therapeutic mechanism records from ChEMBL and map them into:
    - binders
    - binder_modalities
    - protein_binders
    """
    mechanisms = fetch_mechanisms_for_gene(gene_symbol, limit=limit)

    if not mechanisms:
        print(f"[Binders] No ChEMBL mechanism records found for {gene_symbol}")
        return

    with get_connection() as conn:
        with conn.cursor() as cur:
            source_id = upsert_source(cur, "ChEMBL", "https://www.ebi.ac.uk/chembl/")
            protein_id = get_protein_id_by_gene(cur, gene_symbol)

            if not protein_id:
                print(f"[Binders] No protein found in DB for gene {gene_symbol}. Skipping binder linking.")
                return

            inserted = 0

            for mech in mechanisms:
                molecule_chembl_id = mech.get("molecule_chembl_id")
                if not molecule_chembl_id:
                    continue

                try:
                    molecule = fetch_molecule_details(molecule_chembl_id)
                except Exception as e:
                    print(f"[Binders] Failed molecule fetch for {molecule_chembl_id}: {e}")
                    continue

                binder_name = molecule.get("pref_name") or molecule_chembl_id
                mechanism_text = mech.get("mechanism_of_action")
                action_type = mech.get("action_type")
                developer_company = None

                modality_name = infer_modality(molecule)
                modality_id = get_or_create_modality(cur, modality_name)

                binder_type = infer_binder_type(molecule)
                clinical_status = infer_clinical_status(molecule)
                sequence = extract_sequence(molecule)

                max_phase = molecule.get("max_phase")
                approval_status = f"Max Phase {max_phase}" if max_phase is not None else None

                structures = molecule.get("molecule_structures") or {}
                properties = molecule.get("molecule_properties") or {}

                smiles = structures.get("canonical_smiles")
                mw_raw = properties.get("full_mwt")

                try:
                    molecular_weight = float(mw_raw) if mw_raw is not None else None
                except Exception:
                    molecular_weight = None

                cur.execute(
                    """
                    INSERT INTO binders (
                        binder_name,
                        modality_id,
                        binder_type,
                        sequence,
                        clinical_status,
                        binder_description,
                        mechanism_of_action,
                        smiles,
                        molecular_weight,
                        developer_company,
                        approval_status,
                        source_id
                    )
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                    ON CONFLICT (binder_name)
                    DO UPDATE SET
                        modality_id = EXCLUDED.modality_id,
                        binder_type = EXCLUDED.binder_type,
                        sequence = COALESCE(EXCLUDED.sequence, binders.sequence),
                        clinical_status = EXCLUDED.clinical_status,
                        binder_description = EXCLUDED.binder_description,
                        mechanism_of_action = EXCLUDED.mechanism_of_action,
                        smiles = EXCLUDED.smiles,
                        molecular_weight = EXCLUDED.molecular_weight,
                        developer_company = EXCLUDED.developer_company,
                        approval_status = EXCLUDED.approval_status,
                        source_id = EXCLUDED.source_id
                    RETURNING binder_id;
                    """,
                    (
                        binder_name,
                        modality_id,
                        binder_type,
                        sequence,
                        clinical_status,
                        f"Imported from ChEMBL for {gene_symbol}",
                        mechanism_text,
                        smiles,
                        molecular_weight,
                        developer_company,
                        approval_status,
                        source_id,
                    ),
                )
                binder_id = cur.fetchone()[0]

                cur.execute(
                    """
                    INSERT INTO protein_binders (
                        protein_id,
                        binder_id,
                        interaction_type,
                        evidence_summary,
                        source_id
                    )
                    VALUES (%s, %s, %s, %s, %s)
                    ON CONFLICT (protein_id, binder_id)
                    DO UPDATE SET
                        interaction_type = EXCLUDED.interaction_type,
                        evidence_summary = EXCLUDED.evidence_summary,
                        source_id = EXCLUDED.source_id;
                    """,
                    (
                        protein_id,
                        binder_id,
                        action_type,
                        f"Imported from ChEMBL mechanism data for {gene_symbol}",
                        source_id,
                    ),
                )

                inserted += 1

        conn.commit()

    print(f"[Binders] Upserted {inserted} binders for {gene_symbol}")