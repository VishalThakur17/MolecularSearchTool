from ingest.db import get_connection


def normalize_binder_types():
    updates = [
        # --- IgG rules ---
        """
        UPDATE binders
        SET binder_type = 'IgG'
        WHERE
            LOWER(COALESCE(binder_name, '')) LIKE '%mab'
            OR LOWER(COALESCE(binder_type, '')) IN ('antibody', 'monoclonal antibody', 'igg antibody', 'mab', 'igg');
        """,

        # --- VHH rules ---
        """
        UPDATE binders
        SET binder_type = 'VHH'
        WHERE
            LOWER(COALESCE(binder_name, '')) LIKE '%vhh%'
            OR LOWER(COALESCE(binder_name, '')) LIKE '%nanobody%'
            OR LOWER(COALESCE(binder_type, '')) IN ('nanobody', 'vhh', 'single domain antibody');
        """,

        # --- Peptide rules ---
        """
        UPDATE binders
        SET binder_type = 'Peptide'
        WHERE
            LOWER(COALESCE(binder_type, '')) IN ('peptide', 'oligopeptide', 'protein peptide');
        """,

        # --- Small Molecule rules ---
        """
        UPDATE binders
        SET binder_type = 'Small Molecule'
        WHERE
            LOWER(COALESCE(binder_type, '')) IN (
                'small molecule',
                'small-molecule',
                'small_molecule',
                'drug',
                'compound'
            );
        """,

        # --- Convert noisy "Other" entries with no biologic pattern into Small Molecule ---
        """
        UPDATE binders
        SET binder_type = 'Small Molecule'
        WHERE
            COALESCE(binder_type, '') = 'Other'
            AND LOWER(COALESCE(binder_name, '')) NOT LIKE '%mab'
            AND LOWER(COALESCE(binder_name, '')) NOT LIKE '%nanobody%'
            AND LOWER(COALESCE(binder_name, '')) NOT LIKE '%vhh%';
        """,

        # --- Clean empty/null values ---
        """
        UPDATE binders
        SET binder_type = 'Other'
        WHERE binder_type IS NULL OR TRIM(binder_type) = '';
        """
    ]

    with get_connection() as conn:
        with conn.cursor() as cur:
            for sql in updates:
                cur.execute(sql)
        conn.commit()

    print("✅ Binder types normalized successfully")


if __name__ == "__main__":
    normalize_binder_types()