from ingest.db import get_connection


def main():
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute("""
                INSERT INTO binder_diseases (binder_id, disease_id, tag_reason)
                SELECT DISTINCT
                    bt.binder_id,
                    td.disease_id,
                    'Linked through clinical trial disease association'
                FROM binder_trials bt
                JOIN trial_diseases td
                    ON bt.trial_id = td.trial_id
                ON CONFLICT (binder_id, disease_id) DO NOTHING;
            """)

        conn.commit()

    print("Binder disease tags backfilled from linked clinical trials.")


if __name__ == "__main__":
    main()