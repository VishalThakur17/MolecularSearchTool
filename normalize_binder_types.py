"""
Normalize binders into the full sponsor-requested binder taxonomy:
- IgG
- Bispecific Antibody
- ADC
- VHH
- Fc Fusion
- Peptide
- Other

Run:
    python normalize_binder_types.py
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


def normalize_binder_types():
    updates = [
        # Clean empty/null values first.
        """
        UPDATE binders
        SET binder_type = 'Other'
        WHERE binder_type IS NULL OR TRIM(binder_type) = '';
        """,

        # Antibody-drug conjugates must be checked before IgG because many ADCs end in -mab.
        """
        UPDATE binders
        SET binder_type = 'ADC'
        WHERE
            LOWER(COALESCE(binder_name, '')) LIKE '%emtansine%'
            OR LOWER(COALESCE(binder_name, '')) LIKE '%deruxtecan%'
            OR LOWER(COALESCE(binder_name, '')) LIKE '%vedotin%'
            OR LOWER(COALESCE(binder_name, '')) LIKE '%ozogamicin%'
            OR LOWER(COALESCE(binder_name, '')) LIKE '%tesirine%'
            OR LOWER(COALESCE(binder_name, '')) LIKE '%mafodotin%'
            OR LOWER(COALESCE(binder_name, '')) LIKE '%duocarmazine%'
            OR LOWER(COALESCE(binder_type, '')) IN (
                'adc',
                'antibody-drug conjugate',
                'antibody drug conjugate',
                'drug conjugate'
            )
            OR LOWER(COALESCE(binder_description, '')) LIKE '%antibody-drug conjugate%'
            OR LOWER(COALESCE(binder_description, '')) LIKE '%antibody drug conjugate%'
            OR LOWER(COALESCE(mechanism_of_action, '')) LIKE '%antibody-drug conjugate%'
            OR LOWER(COALESCE(mechanism_of_action, '')) LIKE '%antibody drug conjugate%';
        """,

        # Bispecific antibodies / T-cell engagers.
        """
        UPDATE binders
        SET binder_type = 'Bispecific Antibody'
        WHERE
            LOWER(COALESCE(binder_name, '')) LIKE '%bispecific%'
            OR LOWER(COALESCE(binder_name, '')) LIKE '%bi-specific%'
            OR LOWER(COALESCE(binder_name, '')) LIKE '%bite%'
            OR LOWER(COALESCE(binder_type, '')) IN (
                'bispecific',
                'bispecific antibody',
                'bi-specific antibody',
                'dual-specific antibody',
                't-cell engager',
                't cell engager',
                'bite'
            )
            OR LOWER(COALESCE(binder_description, '')) LIKE '%bispecific%'
            OR LOWER(COALESCE(binder_description, '')) LIKE '%t-cell engager%'
            OR LOWER(COALESCE(binder_description, '')) LIKE '%t cell engager%'
            OR LOWER(COALESCE(mechanism_of_action, '')) LIKE '%bispecific%'
            OR LOWER(COALESCE(mechanism_of_action, '')) LIKE '%t-cell engager%'
            OR LOWER(COALESCE(mechanism_of_action, '')) LIKE '%t cell engager%';
        """,

        # Nanobody / VHH.
        """
        UPDATE binders
        SET binder_type = 'VHH'
        WHERE
            LOWER(COALESCE(binder_name, '')) LIKE '%nanobody%'
            OR LOWER(COALESCE(binder_name, '')) LIKE '%vhh%'
            OR LOWER(COALESCE(binder_type, '')) IN (
                'nanobody',
                'vhh',
                'single domain antibody',
                'single-domain antibody',
                'sdab'
            )
            OR LOWER(COALESCE(binder_description, '')) LIKE '%nanobody%'
            OR LOWER(COALESCE(binder_description, '')) LIKE '%single-domain antibody%'
            OR LOWER(COALESCE(binder_description, '')) LIKE '%single domain antibody%';
        """,

        # Fc fusion proteins.
        """
        UPDATE binders
        SET binder_type = 'Fc Fusion'
        WHERE
            LOWER(COALESCE(binder_name, '')) LIKE '%fc fusion%'
            OR LOWER(COALESCE(binder_name, '')) LIKE '%fc-fusion%'
            OR LOWER(COALESCE(binder_type, '')) IN (
                'fc fusion',
                'fc-fusion',
                'fusion protein',
                'receptor-fc',
                'receptor fc'
            )
            OR LOWER(COALESCE(binder_description, '')) LIKE '%fc fusion%'
            OR LOWER(COALESCE(binder_description, '')) LIKE '%fc-fusion%'
            OR LOWER(COALESCE(binder_description, '')) LIKE '%fusion protein%';
        """,

        # Peptides.
        """
        UPDATE binders
        SET binder_type = 'Peptide'
        WHERE
            LOWER(COALESCE(binder_type, '')) IN (
                'peptide',
                'oligopeptide',
                'protein peptide',
                'peptidomimetic'
            )
            OR LOWER(COALESCE(binder_name, '')) LIKE '%peptide%'
            OR LOWER(COALESCE(binder_description, '')) LIKE '%peptide%'
            OR (
                sequence IS NOT NULL
                AND LENGTH(REGEXP_REPLACE(sequence, '[^A-Za-z]', '', 'g')) BETWEEN 2 AND 80
                AND COALESCE(binder_type, '') NOT IN ('ADC', 'Bispecific Antibody', 'VHH', 'Fc Fusion')
            );
        """,

        # General monoclonal antibodies / IgG checked after the more specific antibody categories.
        """
        UPDATE binders
        SET binder_type = 'IgG'
        WHERE
            binder_type NOT IN ('ADC', 'Bispecific Antibody', 'VHH', 'Fc Fusion', 'Peptide')
            AND (
                LOWER(COALESCE(binder_name, '')) LIKE '%mab'
                OR LOWER(COALESCE(binder_type, '')) IN (
                    'antibody',
                    'monoclonal antibody',
                    'igg antibody',
                    'mab',
                    'igg'
                )
                OR LOWER(COALESCE(binder_description, '')) LIKE '%monoclonal antibody%'
            );
        """,

        # Any legacy Small Molecule values are outside the sponsor binder taxonomy.
        """
        UPDATE binders
        SET binder_type = 'Other'
        WHERE LOWER(COALESCE(binder_type, '')) IN (
            'small molecule',
            'small-molecule',
            'small_molecule',
            'drug',
            'compound'
        );
        """,

        # Final safety pass.
        """
        UPDATE binders
        SET binder_type = 'Other'
        WHERE binder_type NOT IN ('IgG', 'Bispecific Antibody', 'ADC', 'VHH', 'Fc Fusion', 'Peptide', 'Other');
        """,
    ]

    with get_connection() as conn:
        with conn.cursor() as cur:
            for sql in updates:
                cur.execute(sql)
        conn.commit()

    print("✅ Binder types normalized to full sponsor taxonomy:")
    for item in FULL_BINDER_TYPES:
        print(f"   - {item}")


if __name__ == "__main__":
    normalize_binder_types()
