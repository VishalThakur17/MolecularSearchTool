from ingest.db import get_connection


def pct(value, total):
    return 0 if total == 0 else round((value / total) * 100, 1)


def main():
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute("SELECT COUNT(*) FROM binders;")
            total = cur.fetchone()[0]

            checks = {
                "Has binder type": "COALESCE(NULLIF(binder_type, ''), '') <> ''",
                "Has sequence": "COALESCE(NULLIF(sequence, ''), '') <> ''",
                "Has clinical status": "COALESCE(NULLIF(clinical_status, ''), NULLIF(approval_status, ''), '') <> ''",
                "Has linked target": "EXISTS (SELECT 1 FROM protein_binders pb WHERE pb.binder_id = binders.binder_id)",
                "Has disease tag": "EXISTS (SELECT 1 FROM binder_diseases bd WHERE bd.binder_id = binders.binder_id)",
                "Has clinical trial evidence": "EXISTS (SELECT 1 FROM binder_trials bt WHERE bt.binder_id = binders.binder_id)",
                "Has direct structure": "EXISTS (SELECT 1 FROM binder_structures bs WHERE bs.binder_id = binders.binder_id)",
                "Has direct or target fallback structure": "EXISTS (SELECT 1 FROM binder_structures bs WHERE bs.binder_id = binders.binder_id) OR EXISTS (SELECT 1 FROM protein_binders pb JOIN protein_structures ps ON ps.protein_id = pb.protein_id WHERE pb.binder_id = binders.binder_id)",
            }

            print("\nBINDER COMPLETENESS AUDIT")
            print("=" * 72)
            print(f"Total binders: {total}")
            print("-" * 72)

            for label, condition in checks.items():
                cur.execute(f"SELECT COUNT(*) FROM binders WHERE {condition};")
                count = cur.fetchone()[0]
                print(f"{label:<40} {count:>6}/{total:<6} {pct(count, total):>6}%")

            print("-" * 72)
            cur.execute("""
                SELECT
                    binder_id,
                    binder_name,
                    binder_type,
                    clinical_status,
                    data_completeness_score
                FROM binders
                ORDER BY data_completeness_score ASC NULLS FIRST, binder_name
                LIMIT 20;
            """)
            low_rows = cur.fetchall()

            print("\nLowest-completeness binders:")
            for row in low_rows:
                binder_id, name, binder_type, status, score = row
                print(f"{binder_id:<5} {str(name)[:32]:<32} {str(binder_type or 'Missing')[:18]:<18} {str(status or 'Missing')[:14]:<14} score={score}")

            print("\nTip: In pgAdmin, run: SELECT * FROM binder_completeness_audit ORDER BY data_completeness_score ASC;")


if __name__ == "__main__":
    main()
