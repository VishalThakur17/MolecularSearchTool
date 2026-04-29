import re
import requests
from ingest.db import get_connection

CHEMBL_MECHANISM_URL = "https://www.ebi.ac.uk/chembl/api/data/mechanism.json"
CHEMBL_MOLECULE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule"
CHEMBL_TARGET_URL = "https://www.ebi.ac.uk/chembl/api/data/target.json"

VALID_AA_CHARS = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")


def fetch_target_chembl_ids_for_gene(gene_symbol: str, limit: int = 25):
    params = {
        "target_synonym__icontains": gene_symbol,
        "limit": limit,
    }

    response = requests.get(CHEMBL_TARGET_URL, params=params, timeout=45)
    response.raise_for_status()
    data = response.json()
    targets = data.get("targets", [])

    chembl_ids = []
    seen = set()

    for target in targets:
        tid = target.get("target_chembl_id")
        if not tid or tid in seen:
            continue

        target_type = (target.get("target_type") or "").lower()
        if target_type and "single protein" not in target_type:
            continue

        chembl_ids.append(tid)
        seen.add(tid)

    return chembl_ids


def fetch_mechanisms_for_target_chembl_id(target_chembl_id: str, limit: int = 50):
    params = {
        "target_chembl_id": target_chembl_id,
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


def _molecule_text(molecule_record: dict) -> str:
    parts = [
        molecule_record.get("pref_name"),
        molecule_record.get("molecule_type"),
        molecule_record.get("molecule_chembl_id"),
    ]

    for synonym in molecule_record.get("molecule_synonyms", []) or []:
        if isinstance(synonym, dict):
            parts.append(synonym.get("molecule_synonym"))
            parts.append(synonym.get("syn_type"))

    biotherapeutic = molecule_record.get("biotherapeutic") or {}
    if isinstance(biotherapeutic, dict):
        parts.extend(str(v) for v in biotherapeutic.values() if v)

    return " ".join(str(p) for p in parts if p).lower()


def infer_modality(molecule_record: dict) -> str:
    text = _molecule_text(molecule_record)
    molecule_type = (molecule_record.get("molecule_type") or "").lower()

    if any(term in text for term in [
        "antibody-drug conjugate", "antibody drug conjugate", " adc ",
        "emtansine", "deruxtecan", "vedotin", "ozogamicin",
        "tesirine", "mafodotin", "duocarmazine"
    ]):
        return "Antibody-Drug Conjugate"

    if any(term in text for term in [
        "bispecific", "bi-specific", "bi specific", "dual-specific",
        "dual specific", "bite", "t-cell engager", "t cell engager"
    ]):
        return "Bispecific Antibody"

    if any(term in text for term in [
        "nanobody", "vhh", "single-domain antibody", "single domain antibody",
        "camelid", "sdab"
    ]):
        return "Nanobody / VHH"

    if any(term in text for term in [
        "fc fusion", "fc-fusion", "fusion protein", "receptor-fc",
        "receptor fc", "fc region", "trap"
    ]):
        return "Fc Fusion"

    if "peptide" in text or "oligopeptide" in molecule_type:
        return "Peptide"

    if any(term in molecule_type for term in ["antibody", "protein", "monoclonal"]):
        return "Antibody"

    if any(term in text for term in ["mab", "monoclonal"]):
        return "Antibody"

    return "Small Molecule"


def _clean_sequence(raw_value) -> str | None:
    if raw_value is None:
        return None

    seq = str(raw_value).strip().upper()
    if not seq:
        return None

    seq = re.sub(r"[\s\-_*.:;,\|/\\]+", "", seq)

    if len(seq) < 5:
        return None

    if not set(seq).issubset(VALID_AA_CHARS):
        return None

    return seq


def _collect_candidate_sequences(value, candidates: list[str]):
    if value is None:
        return

    if isinstance(value, str):
        cleaned = _clean_sequence(value)
        if cleaned:
            candidates.append(cleaned)
        return

    if isinstance(value, list):
        for item in value:
            _collect_candidate_sequences(item, candidates)
        return

    if isinstance(value, dict):
        preferred_keys = [
            "sequence",
            "full_sequence",
            "sequence_text",
            "protein_sequence",
            "peptide_sequence",
            "chain_sequence",
            "component_sequence",
        ]

        for key in preferred_keys:
            if key in value:
                _collect_candidate_sequences(value.get(key), candidates)

        for nested_value in value.values():
            if isinstance(nested_value, (dict, list)):
                _collect_candidate_sequences(nested_value, candidates)


def extract_sequence(molecule_record: dict) -> str | None:
    candidates: list[str] = []

    top_level_fields = [
        molecule_record.get("sequence"),
        molecule_record.get("full_sequence"),
        molecule_record.get("protein_sequence"),
        molecule_record.get("peptide_sequence"),
    ]
    for value in top_level_fields:
        cleaned = _clean_sequence(value)
        if cleaned:
            candidates.append(cleaned)

    nested_sections = [
        molecule_record.get("molecule_properties"),
        molecule_record.get("molecule_structures"),
        molecule_record.get("molecule_hierarchy"),
        molecule_record.get("biotherapeutic"),
        molecule_record.get("biotherapeutic_properties"),
        molecule_record.get("molecule_components"),
        molecule_record.get("cross_references"),
    ]
    for section in nested_sections:
        _collect_candidate_sequences(section, candidates)

    if not candidates:
        return None

    return max(candidates, key=len)


def infer_binder_type(molecule_record: dict) -> str:
    text = _molecule_text(molecule_record)
    molecule_type = (molecule_record.get("molecule_type") or "").lower().strip()
    molecule_structures = molecule_record.get("molecule_structures") or {}
    molecule_properties = molecule_record.get("molecule_properties") or {}
    sequence = (extract_sequence(molecule_record) or "").strip()

    # Order matters. ADCs and bispecifics usually include antibody-like terms,
    # so they must be checked before general IgG / monoclonal antibody.
    if any(term in text for term in [
        "antibody-drug conjugate", "antibody drug conjugate", " adc ",
        "emtansine", "deruxtecan", "vedotin", "ozogamicin",
        "tesirine", "mafodotin", "duocarmazine", "maytansinoid"
    ]):
        return "ADC"

    if any(term in text for term in [
        "bispecific", "bi-specific", "bi specific", "dual-specific",
        "dual specific", "bite", "t-cell engager", "t cell engager"
    ]):
        return "Bispecific Antibody"

    if any(term in text for term in [
        "nanobody", "vhh", "single-domain antibody", "single domain antibody",
        "camelid", "sdab"
    ]):
        return "VHH"

    if any(term in text for term in [
        "fc fusion", "fc-fusion", "fusion protein", "receptor-fc",
        "receptor fc", "fc region", "trap"
    ]):
        return "Fc Fusion"

    if "peptide" in text or "oligopeptide" in molecule_type:
        return "Peptide"

    if sequence:
        seq_len = len(sequence)
        if 2 <= seq_len <= 80:
            return "Peptide"

    if "antibody" in molecule_type or "monoclonal" in molecule_type:
        return "IgG"

    pref_name = (molecule_record.get("pref_name") or "").lower().strip()
    if pref_name.endswith("mab") or " monoclonal " in text:
        return "IgG"

    # ChEMBL records with SMILES, chemistry properties, or molecule_type = molecule
    # are usually small molecules. Classify them explicitly instead of vague Other.
    has_smiles = bool(molecule_structures.get("canonical_smiles"))
    has_chem_props = bool(molecule_properties.get("full_mwt") or molecule_properties.get("alogp"))
    if (
        molecule_type in {"small molecule", "synthetic small molecule", "molecule", ""}
        or "small molecule" in text
        or has_smiles
        or has_chem_props
    ):
        return "Small Molecule"

    return "Other"


def infer_clinical_status(molecule_record: dict) -> str | None:
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
    if not modality_name:
        modality_name = "Unknown"

    # This project already uses binder_modalities as the referenced table for
    # binders.modality_id. Do not use a separate modalities table here.
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
    row = cur.fetchone()
    if isinstance(row, dict):
        return row["modality_id"]
    return row[0]


def get_protein_id_by_gene(cur, gene_symbol: str):
    cur.execute(
        """
        SELECT p.protein_id
        FROM proteins p
        JOIN genes g ON p.gene_id = g.gene_id
        WHERE LOWER(g.gene_symbol) = LOWER(%s)
        LIMIT 1;
        """,
        (gene_symbol,),
    )
    row = cur.fetchone()
    return row[0] if row else None


def binder_already_linked(cur, protein_id: int, binder_id: int) -> bool:
    cur.execute(
        """
        SELECT 1
        FROM protein_binders
        WHERE protein_id = %s AND binder_id = %s
        LIMIT 1;
        """,
        (protein_id, binder_id),
    )
    return cur.fetchone() is not None


def upsert_binder(cur, molecule: dict, gene_symbol: str, source_id: int, mechanism_text: str | None):
    binder_name = molecule.get("pref_name") or molecule.get("molecule_chembl_id")
    if not binder_name:
        return None

    modality_name = infer_modality(molecule)
    modality_id = get_or_create_modality(cur, modality_name)

    sequence = extract_sequence(molecule)
    binder_type = infer_binder_type(molecule)
    clinical_status = infer_clinical_status(molecule)

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
            clinical_status = COALESCE(EXCLUDED.clinical_status, binders.clinical_status),
            binder_description = EXCLUDED.binder_description,
            mechanism_of_action = COALESCE(EXCLUDED.mechanism_of_action, binders.mechanism_of_action),
            smiles = COALESCE(EXCLUDED.smiles, binders.smiles),
            molecular_weight = COALESCE(EXCLUDED.molecular_weight, binders.molecular_weight),
            developer_company = COALESCE(EXCLUDED.developer_company, binders.developer_company),
            approval_status = COALESCE(EXCLUDED.approval_status, binders.approval_status),
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
            None,
            approval_status,
            source_id,
        ),
    )
    return cur.fetchone()[0]


def link_binder_to_protein(cur, protein_id: int, binder_id: int, action_type: str | None, source_id: int, gene_symbol: str):
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
            interaction_type = COALESCE(EXCLUDED.interaction_type, protein_binders.interaction_type),
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


def upsert_binders_for_gene(
    gene_symbol: str,
    limit_per_target: int = 50,
    target_search_limit: int = 25,
):
    with get_connection() as conn:
        with conn.cursor() as cur:
            protein_id = get_protein_id_by_gene(cur, gene_symbol)
            if not protein_id:
                print(f"[Binders] No protein found in DB for gene {gene_symbol}.")
                return

            source_id = upsert_source(cur, "ChEMBL", "https://www.ebi.ac.uk/chembl/")

            target_chembl_ids = fetch_target_chembl_ids_for_gene(
                gene_symbol,
                limit=target_search_limit,
            )
            if not target_chembl_ids:
                print(f"[Binders] No ChEMBL target IDs found for gene {gene_symbol}.")
                return

            total_mechanisms = 0
            inserted_binders = 0
            linked_pairs = 0
            sequences_found = 0

            seen_molecule_ids = set()

            for target_chembl_id in target_chembl_ids:
                try:
                    mechanisms = fetch_mechanisms_for_target_chembl_id(
                        target_chembl_id,
                        limit=limit_per_target,
                    )
                except Exception as e:
                    print(f"[Binders] Failed mechanism fetch for target {target_chembl_id}: {e}")
                    continue

                total_mechanisms += len(mechanisms)

                for mech in mechanisms:
                    molecule_chembl_id = mech.get("molecule_chembl_id")
                    if not molecule_chembl_id or molecule_chembl_id in seen_molecule_ids:
                        continue

                    seen_molecule_ids.add(molecule_chembl_id)

                    try:
                        molecule = fetch_molecule_details(molecule_chembl_id)
                    except Exception as e:
                        print(f"[Binders] Failed molecule fetch for {molecule_chembl_id}: {e}")
                        continue

                    if extract_sequence(molecule):
                        sequences_found += 1

                    binder_id = upsert_binder(
                        cur,
                        molecule,
                        gene_symbol,
                        source_id,
                        mech.get("mechanism_of_action"),
                    )
                    if not binder_id:
                        continue

                    inserted_binders += 1

                    action_type = mech.get("action_type")
                    was_linked = binder_already_linked(cur, protein_id, binder_id)

                    link_binder_to_protein(
                        cur,
                        protein_id,
                        binder_id,
                        action_type,
                        source_id,
                        gene_symbol,
                    )

                    if not was_linked:
                        linked_pairs += 1

            conn.commit()

    print(
        f"[Binders] Gene={gene_symbol} | "
        f"Targets found={len(target_chembl_ids)} | "
        f"Mechanisms read={total_mechanisms} | "
        f"Binders upserted={inserted_binders} | "
        f"New protein-binder links={linked_pairs} | "
        f"Molecules with extracted sequence={sequences_found}"
    )