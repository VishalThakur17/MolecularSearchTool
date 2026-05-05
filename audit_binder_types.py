"""
Audit binder classification coverage for the full sponsor taxonomy.

Run:
    python audit_binder_types.py
"""

from ingest.db import get_connection

FULL_BINDER_TYPES = [
    "Small Molecule",
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
                SELECT COUNT(*)
                FROM binders;
            """)
            total = cur.fetchone()[0] or 0

            cur.execute("""
                SELECT binder_type, COUNT(*)
                FROM binders
                GROUP BY binder_type
                ORDER BY COUNT(*) DESC, binder_type;
            """)
            rows = cur.fetchall()
            counts = {binder_type or "NULL": count for binder_type, count in rows}

            print("\n=== BINDER TYPE DISTRIBUTION ===")
            print(f"Total binders: {total}")

            for binder_type in FULL_BINDER_TYPES:
                count = counts.get(binder_type, 0)
                pct = (count / total * 100) if total else 0
                print(f"{binder_type:<24} | {count:>5} | {pct:>6.1f}%")

            extras = [
                (binder_type, count)
                for binder_type, count in rows
                if (binder_type or "NULL") not in FULL_BINDER_TYPES
            ]

            if extras:
                print("\n=== NON-STANDARD TYPES THAT SHOULD BE NORMALIZED ===")
                for binder_type, count in extras:
                    print(f"{(binder_type or 'NULL'):<24} | {count}")
            else:
                print("\nNo non-standard binder types found.")

            sponsor_required_types = [
                "IgG",
                "Bispecific Antibody",
                "ADC",
                "VHH",
                "Fc Fusion",
                "Peptide",
            ]

            missing = [
                binder_type
                for binder_type in sponsor_required_types
                if counts.get(binder_type, 0) == 0
            ]

            if missing:
                print("\n=== TAXONOMY GAPS ===")
                for binder_type in missing:
                    print(f"No examples found for: {binder_type}")
                print("Add curated examples or improve ingestion keywords before calling the taxonomy fully demonstrated.")
            else:
                print("\nEvery sponsor-required taxonomy category has at least one binder example.")

            other_count = counts.get("Other", 0)
            other_pct = (other_count / total * 100) if total else 0

            print("\n=== CLASSIFICATION QUALITY CHECK ===")
            print(f"Other count: {other_count} / {total} ({other_pct:.1f}%)")

            if other_pct > 25:
                print("Warning: A large percentage of binders are still classified as Other.")
                print("This usually means more normalization keywords or better source data are needed.")
            else:
                print("Good: Other is reasonably low after Small Molecule separation.")

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
                    print("No entries yet.")
                    continue

                for binder_name, clinical_status in samples:
                    print(f"{binder_name} | {clinical_status or 'N/A'}")


if __name__ == "__main__":
    audit_binder_types()