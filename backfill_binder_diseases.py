from ingest.db import get_connection

sql = """
INSERT INTO binder_diseases (binder_id, disease_id, source_id, tag_reason)
SELECT DISTINCT
    pb.binder_id,
    dp.disease_id,
    NULL::INTEGER,
    'Inherited from protein-disease relationship during Task 2 backfill'
FROM protein_binders pb
JOIN disease_proteins dp
    ON pb.protein_id = dp.protein_id
ON CONFLICT (binder_id, disease_id) DO NOTHING;
"""

with get_connection() as conn:
    with conn.cursor() as cur:
        cur.execute(sql)
    conn.commit()

print("✅ Binder disease tags backfilled successfully")