from ingest.catalog import GENES_TO_INGEST, TRIAL_TARGETS
from ingest.ingest_uniprot import upsert_gene_and_protein
from ingest.ingest_binders import upsert_binders_for_gene
from ingest.ingest_trials import upsert_trials_for_target


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
                limit_per_target=50,
                target_search_limit=25,
            )
        except Exception as e:
            print(f"[Binders] Failed for {gene}: {e}")

    print("\n=== STEP 3: CLINICAL TRIALS ===")
    for item in TRIAL_TARGETS:
        gene = item["gene"]
        disease = item["disease"]

        try:
            upsert_trials_for_target(
                gene_symbol=gene,
                disease_name=disease,
                max_trials=50,
                max_pages=3,
            )
        except Exception as e:
            print(f"[Trials] Failed for {gene} / {disease}: {e}")

    print("\n=== INGESTION COMPLETE ===")


if __name__ == "__main__":
    main()