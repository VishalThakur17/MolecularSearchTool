"""
Audit binder classification coverage for the full sponsor taxonomy.

Run:
    python audit_binder_types.py
"""

from ingest.db import get_connection

FULL_BINDER_TYPES = [
    "IgG",
    "Bispecific Antibody",
    "ADC",
    "VHH",
    "Fc Fusion",
    "Peptide",
    "Other",
]


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
            counts = {binder_type or "NULL": count for binder_type, count in rows}

            print("\n=== BINDER TYPE DISTRIBUTION ===")
            for binder_type in FULL_BINDER_TYPES:
                print(f"{binder_type:<24} | {counts.get(binder_type, 0)}")

            extras = [(binder_type, count) for binder_type, count in rows if binder_type not in FULL_BINDER_TYPES]
            if extras:
                print("\n=== NON-STANDARD TYPES THAT SHOULD BE NORMALIZED ===")
                for binder_type, count in extras:
                    print(f"{(binder_type or 'NULL'):<24} | {count}")

            print("\n=== SAMPLE BINDERS BY TYPE ===")
            for target_type in FULL_BINDER_TYPES:
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
                    print("No entries yet. Add/import binders for this category if it is required for your demo.")
                    continue

                for binder_name, clinical_status in samples:
                    print(f"{binder_name} | {clinical_status or 'N/A'}")


if __name__ == "__main__":
    audit_binder_types()
