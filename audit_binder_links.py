from ingest.db import get_connection


def audit_target_binder_links():
    sql = """
    SELECT
        g.gene_symbol,
        p.protein_name,
        COUNT(DISTINCT pb.binder_id) AS binder_count
    FROM proteins p
    JOIN genes g ON p.gene_id = g.gene_id
    LEFT JOIN protein_binders pb ON p.protein_id = pb.protein_id
    GROUP BY g.gene_symbol, p.protein_name
    ORDER BY binder_count DESC, g.gene_symbol;
    """

    detail_sql = """
    SELECT
        g.gene_symbol,
        b.binder_name,
        b.binder_type,
        b.clinical_status,
        pb.interaction_type
    FROM protein_binders pb
    JOIN proteins p ON pb.protein_id = p.protein_id
    JOIN genes g ON p.gene_id = g.gene_id
    JOIN binders b ON pb.binder_id = b.binder_id
    WHERE LOWER(g.gene_symbol) = LOWER(%s)
    ORDER BY b.binder_name;
    """

    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute(sql)
            rows = cur.fetchall()

            print("\n=== TARGET → BINDER COVERAGE ===")
            for row in rows:
                gene_symbol, protein_name, binder_count = row
                print(f"{gene_symbol:<10} | binders={binder_count:<3} | {protein_name}")

            print("\n=== DETAILED CHECKS ===")
            for gene in ["EGFR", "ERBB2", "KRAS", "ALK"]:
                cur.execute(detail_sql, (gene,))
                detail_rows = cur.fetchall()

                print(f"\n--- {gene} ---")
                if not detail_rows:
                    print("No binders linked.")
                    continue

                for d in detail_rows[:20]:
                    gene_symbol, binder_name, binder_type, clinical_status, interaction_type = d
                    print(
                        f"{binder_name} | type={binder_type or 'N/A'} | "
                        f"status={clinical_status or 'N/A'} | interaction={interaction_type or 'N/A'}"
                    )


if __name__ == "__main__":
    audit_target_binder_links()