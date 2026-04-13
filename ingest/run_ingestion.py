from ingest.ingest_binders import upsert_binders_for_gene


def main():
    genes = ["EGFR", "ERBB2", "KRAS", "ALK"]

    for gene in genes:
        print(f"\n=== Processing binders for {gene} ===")
        try:
            upsert_binders_for_gene(gene, limit_per_target=25)
        except Exception as e:
            print(f"Failed for {gene}: {e}")


if __name__ == "__main__":
    main()