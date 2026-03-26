from ingest.ingest_uniprot import upsert_gene_and_protein
from ingest.ingest_trials import upsert_trials_for_target
from ingest.ingest_binders import upsert_binders_for_gene
from ingest.ingest_pdb import upsert_structures_for_gene

TARGETS = [
    {"gene": "EGFR", "disease": "Lung Cancer"},
    {"gene": "ERBB2", "disease": "Breast Cancer"},
    {"gene": "KRAS", "disease": "Lung Cancer"},
    {"gene": "ALK", "disease": "Lung Cancer"},
]


def main():
    print("Starting ingestion...")

    for item in TARGETS:
        gene = item["gene"]
        disease = item["disease"]

        print(f"\n--- Processing {gene} / {disease} ---")

        upsert_gene_and_protein(gene)
        upsert_trials_for_target(gene, disease, max_trials=8)
        upsert_binders_for_gene(gene, limit=10)
        upsert_structures_for_gene(gene, max_results=8)

    print("\nIngestion complete.")


if __name__ == "__main__":
    main()