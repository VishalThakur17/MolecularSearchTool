"""
Add a small curated set of binder examples for sponsor-taxonomy coverage.

Run after enrich_data.py if your audit still shows missing/low categories:

    python seed_curated_binder_examples.py

This does NOT replace your real ingestion. It adds or updates a small set of
well-known examples so the UI can demonstrate IgG, ADC, bispecific antibody,
VHH/nanobody, Fc fusion, peptide, and small molecule/context separation.
"""

from ingest.db import get_connection
from ingest.ingest_binders import get_or_create_modality


CURATED_BINDERS = [
    {
        "binder_name": "Trastuzumab",
        "binder_type": "IgG",
        "modality_name": "Antibody",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Genentech / Roche",
        "mechanism_of_action": "Humanized monoclonal antibody targeting HER2/ERBB2.",
        "binder_description": "Therapeutic monoclonal antibody that binds HER2/ERBB2.",
        "gene_symbol": "ERBB2",
    },
    {
        "binder_name": "Cetuximab",
        "binder_type": "IgG",
        "modality_name": "Antibody",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Eli Lilly / Bristol Myers Squibb",
        "mechanism_of_action": "Monoclonal antibody targeting EGFR.",
        "binder_description": "Therapeutic monoclonal antibody that binds EGFR.",
        "gene_symbol": "EGFR",
    },
    {
        "binder_name": "Trastuzumab emtansine",
        "binder_type": "ADC",
        "modality_name": "Antibody-Drug Conjugate",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Genentech / Roche",
        "mechanism_of_action": "HER2-directed antibody-drug conjugate using trastuzumab linked to DM1.",
        "binder_description": "ADC example for HER2-positive cancers.",
        "gene_symbol": "ERBB2",
    },
    {
        "binder_name": "Trastuzumab deruxtecan",
        "binder_type": "ADC",
        "modality_name": "Antibody-Drug Conjugate",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Daiichi Sankyo / AstraZeneca",
        "mechanism_of_action": "HER2-directed antibody-drug conjugate with a topoisomerase I inhibitor payload.",
        "binder_description": "ADC example for HER2-positive and HER2-low tumor settings.",
        "gene_symbol": "ERBB2",
    },
    {
        "binder_name": "Amivantamab",
        "binder_type": "Bispecific Antibody",
        "modality_name": "Bispecific Antibody",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Janssen",
        "mechanism_of_action": "Bispecific antibody targeting EGFR and MET.",
        "binder_description": "Bispecific antibody example that binds EGFR and MET.",
        "gene_symbol": "EGFR",
    },
    {
        "binder_name": "Zanidatamab",
        "binder_type": "Bispecific Antibody",
        "modality_name": "Bispecific Antibody",
        "clinical_status": "Approved / Clinical",
        "approval_status": "Approved or late-stage clinical depending on indication/region",
        "developer_company": "Jazz Pharmaceuticals / Zymeworks",
        "mechanism_of_action": "HER2-directed bispecific antibody.",
        "binder_description": "Bispecific antibody example for HER2/ERBB2.",
        "gene_symbol": "ERBB2",
    },
    {
        "binder_name": "Caplacizumab",
        "binder_type": "VHH",
        "modality_name": "Nanobody / VHH",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Sanofi",
        "mechanism_of_action": "Nanobody targeting von Willebrand factor.",
        "binder_description": "VHH/nanobody therapeutic example.",
        "gene_symbol": None,
    },
    {
        "binder_name": "Envafolimab",
        "binder_type": "VHH",
        "modality_name": "Nanobody / VHH",
        "clinical_status": "Approved / Clinical",
        "approval_status": "Approved or clinical depending on indication/region",
        "developer_company": "Alphamab Oncology / 3D Medicines",
        "mechanism_of_action": "Single-domain antibody targeting PD-L1.",
        "binder_description": "VHH/nanobody-style therapeutic example.",
        "gene_symbol": None,
    },
    {
        "binder_name": "Aflibercept",
        "binder_type": "Fc Fusion",
        "modality_name": "Fc Fusion",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Regeneron / Bayer",
        "mechanism_of_action": "VEGF trap Fc fusion protein.",
        "binder_description": "Fc fusion therapeutic example.",
        "gene_symbol": None,
    },
    {
        "binder_name": "Etanercept",
        "binder_type": "Fc Fusion",
        "modality_name": "Fc Fusion",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "Amgen / Pfizer",
        "mechanism_of_action": "TNF receptor Fc fusion protein that binds TNF.",
        "binder_description": "Fc fusion therapeutic example.",
        "gene_symbol": None,
    },
    {
        "binder_name": "Cilengitide",
        "binder_type": "Peptide",
        "modality_name": "Peptide",
        "clinical_status": "Clinical / Investigational",
        "approval_status": "Investigational",
        "developer_company": "Merck KGaA",
        "mechanism_of_action": "Cyclic RGD peptide inhibitor of selected integrins.",
        "binder_description": "Peptide binder example.",
        "gene_symbol": None,
    },
    {
        "binder_name": "Goserelin",
        "binder_type": "Peptide",
        "modality_name": "Peptide",
        "clinical_status": "Approved",
        "approval_status": "Approved / Max Phase 4",
        "developer_company": "AstraZeneca",
        "mechanism_of_action": "Peptide GnRH receptor agonist.",
        "binder_description": "Peptide therapeutic example.",
        "gene_symbol": None,
    },
]


def scalar(row, key_or_index):
    if row is None:
        return None
    if isinstance(row, dict):
        return row.get(key_or_index)
    return row[key_or_index]


def get_or_create_source(cur):
    cur.execute(
        """
        INSERT INTO sources (source_name, source_url)
        VALUES (%s, %s)
        ON CONFLICT (source_name)
        DO UPDATE SET source_url = EXCLUDED.source_url
        RETURNING source_id;
        """,
        ("Curated Demo Set", "Internal curated taxonomy examples"),
    )
    row = cur.fetchone()
    return scalar(row, "source_id") if isinstance(row, dict) else row[0]


def find_protein_id(cur, gene_symbol):
    if not gene_symbol:
        return None

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
    if not row:
        return None
    return scalar(row, "protein_id") if isinstance(row, dict) else row[0]


def upsert_curated_binder(cur, item, source_id):
    modality_id = get_or_create_modality(cur, item["modality_name"])

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
            developer_company,
            approval_status,
            source_id
        )
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT (binder_name)
        DO UPDATE SET
            modality_id = EXCLUDED.modality_id,
            binder_type = EXCLUDED.binder_type,
            clinical_status = COALESCE(EXCLUDED.clinical_status, binders.clinical_status),
            binder_description = COALESCE(EXCLUDED.binder_description, binders.binder_description),
            mechanism_of_action = COALESCE(EXCLUDED.mechanism_of_action, binders.mechanism_of_action),
            developer_company = COALESCE(EXCLUDED.developer_company, binders.developer_company),
            approval_status = COALESCE(EXCLUDED.approval_status, binders.approval_status),
            source_id = EXCLUDED.source_id
        RETURNING binder_id;
        """,
        (
            item["binder_name"],
            modality_id,
            item["binder_type"],
            item.get("sequence"),
            item.get("clinical_status"),
            item.get("binder_description"),
            item.get("mechanism_of_action"),
            item.get("developer_company"),
            item.get("approval_status"),
            source_id,
        ),
    )
    row = cur.fetchone()
    return scalar(row, "binder_id") if isinstance(row, dict) else row[0]


def link_to_protein_if_possible(cur, protein_id, binder_id, source_id, item):
    if not protein_id:
        return False

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
            evidence_summary = COALESCE(EXCLUDED.evidence_summary, protein_binders.evidence_summary),
            source_id = EXCLUDED.source_id;
        """,
        (
            protein_id,
            binder_id,
            "Curated target association",
            f"Curated sponsor-taxonomy example for {item['binder_type']}.",
            source_id,
        ),
    )
    return True


def main():
    with get_connection() as conn:
        with conn.cursor() as cur:
            source_id = get_or_create_source(cur)

            print("\n=== Seeding curated sponsor-taxonomy binder examples ===")
            for item in CURATED_BINDERS:
                binder_id = upsert_curated_binder(cur, item, source_id)
                protein_id = find_protein_id(cur, item.get("gene_symbol"))
                linked = link_to_protein_if_possible(cur, protein_id, binder_id, source_id, item)

                link_msg = f" linked to {item['gene_symbol']}" if linked else " no protein link"
                print(f"[Curated] {item['binder_name']} -> {item['binder_type']} |{link_msg}")

        conn.commit()

    print("\n=== CURATED BINDER SEED COMPLETE ===")


if __name__ == "__main__":
    main()
