from ingest.db import get_connection

binder_id = 1
disease_id = 1

with get_connection() as conn:
    with conn.cursor() as cur:
        cur.execute("""
            INSERT INTO binder_diseases (binder_id, disease_id, tag_reason)
            VALUES (%s, %s, %s)
            ON CONFLICT (binder_id, disease_id)
            DO UPDATE SET tag_reason = EXCLUDED.tag_reason;
        """, (binder_id, disease_id, "Manual sponsor-aligned disease tag"))
    conn.commit()

print("✅ Binder disease tag inserted")