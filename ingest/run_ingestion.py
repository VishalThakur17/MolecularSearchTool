from ingest.catalog import GENES_TO_INGEST, TRIAL_TARGETS
from ingest.ingest_uniprot import upsert_gene_and_protein
from ingest.ingest_binders import upsert_binders_for_gene
from ingest.ingest_trials import upsert_trials_for_target
from ingest.ingest_pdb import upsert_structures_for_gene


def main():
    print("\n=== STEP 1: PROTEINS / GENES ===")
    for gene in GENES_TO_INGEST:
        try:
            upsert_gene_and_protein(gene)
        except Exception as e:
            print(f"[UniProt] Failed for {gene}: {e}")

    print("\n=== STEP 2: BINDERS ===")
    for gene in GENES_TO_INGEST:
        try:
            upsert_binders_for_gene(
                gene_symbol=gene,
                limit_per_target=100,
                target_search_limit=50,
            )
        except Exception as e:
            print(f"[Binders] Failed for {gene}: {e}")

    print("\n=== STEP 3: CLINICAL TRIALS ===")
    for item in TRIAL_TARGETS:
        try:
            upsert_trials_for_target(
                gene_symbol=item["gene"],
                disease_name=item["disease"],
                max_trials=100,
                max_pages=5,
            )
        except Exception as e:
            print(f"[Trials] Failed for {item['gene']} / {item['disease']}: {e}")

    print("\n=== STEP 4: STRUCTURES / PDB ===")
    for gene in GENES_TO_INGEST:
        try:
            upsert_structures_for_gene(
                gene_symbol=gene,
                max_results=50,
            )
        except Exception as e:
            print(f"[PDB] Failed for {gene}: {e}")

    print("\n=== INGESTION COMPLETE ===")
    print("Now run: python finalize_relationships.py")


if __name__ == "__main__":
    main()