from ingest.db import get_connection


def audit_binder_types():
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute("""
                SELECT binder_type, COUNT(*)
                FROM binders
                GROUP BY binder_type
                ORDER BY COUNT(*) DESC, binder_type;
            """)
            rows = cur.fetchall()

            print("\n=== BINDER TYPE DISTRIBUTION ===")
            for binder_type, count in rows:
                print(f"{(binder_type or 'NULL'):<20} | {count}")

            print("\n=== SAMPLE BINDERS BY TYPE ===")
            for target_type in ["IgG", "VHH", "Peptide", "Small Molecule", "Other"]:
                print(f"\n--- {target_type} ---")
                cur.execute("""
                    SELECT binder_name, clinical_status
                    FROM binders
                    WHERE binder_type = %s
                    ORDER BY binder_name
                    LIMIT 10;
                """, (target_type,))
                samples = cur.fetchall()

                if not samples:
                    print("No entries.")
                    continue

                for binder_name, clinical_status in samples:
                    print(f"{binder_name} | {clinical_status or 'N/A'}")


if __name__ == "__main__":
    audit_binder_types()