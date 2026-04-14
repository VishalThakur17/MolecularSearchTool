from ingest.db import get_connection

sql_statements = [
    # protein_diseases from existing disease_proteins
    """
    INSERT INTO protein_diseases (protein_id, disease_id, source_id, tag_reason)
    SELECT DISTINCT
        dp.protein_id,
        dp.disease_id,
        NULL::INTEGER,
        'Backfilled from disease_proteins'
    FROM disease_proteins dp
    ON CONFLICT (protein_id, disease_id) DO NOTHING;
    """,

    # binder_diseases from protein_binders + protein_diseases
    """
    INSERT INTO binder_diseases (binder_id, disease_id, source_id, tag_reason)
    SELECT DISTINCT
        pb.binder_id,
        pd.disease_id,
        NULL::INTEGER,
        'Backfilled from protein_diseases through protein_binders'
    FROM protein_binders pb
    JOIN protein_diseases pd
        ON pb.protein_id = pd.protein_id
    ON CONFLICT (binder_id, disease_id) DO NOTHING;
    """,

    # trial_diseases from existing disease_trials
    """
    INSERT INTO trial_diseases (trial_id, disease_id, source_id, tag_reason)
    SELECT DISTINCT
        dt.trial_id,
        dt.disease_id,
        NULL::INTEGER,
        'Backfilled from disease_trials'
    FROM disease_trials dt
    ON CONFLICT (trial_id, disease_id) DO NOTHING;
    """
]

with get_connection() as conn:
    with conn.cursor() as cur:
        for sql in sql_statements:
            cur.execute(sql)
    conn.commit()

print("✅ Disease tags backfilled successfully")