import os
from dotenv import load_dotenv
import psycopg2

load_dotenv()
DB_HOST = os.getenv("DB_HOST")
DB_PORT = os.getenv("DB_PORT", "5432")
DB_NAME = os.getenv("DB_NAME", "molecular_search_db")
DB_USER = os.getenv("DB_USER", "postgres")
DB_PASSWORD = os.getenv("DB_PASSWORD")


def get_connection():
    return psycopg2.connect(
        host=DB_HOST,
        port=DB_PORT,
        dbname=DB_NAME,
        user=DB_USER,
        password=DB_PASSWORD,
        sslmode="require",
    )


def main():
    with get_connection() as conn:
        with conn.cursor() as cur:
            cur.execute("""
                CREATE TABLE IF NOT EXISTS binder_binding_sites (
                    binding_site_id SERIAL PRIMARY KEY,
                    binder_id INTEGER NOT NULL REFERENCES binders(binder_id) ON DELETE CASCADE,
                    protein_id INTEGER NOT NULL REFERENCES proteins(protein_id) ON DELETE CASCADE,
                    region_start INTEGER,
                    region_end INTEGER,
                    region_label VARCHAR(255),
                    evidence_type VARCHAR(120),
                    source_note TEXT,
                    color_hex VARCHAR(20),
                    created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                    UNIQUE (binder_id, protein_id, region_start, region_end, region_label)
                );
            """)
            cur.execute("ALTER TABLE binder_binding_sites ADD COLUMN IF NOT EXISTS updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP;")
            cur.execute("CREATE INDEX IF NOT EXISTS idx_binder_binding_sites_binder ON binder_binding_sites(binder_id);")
            cur.execute("CREATE INDEX IF NOT EXISTS idx_binder_binding_sites_protein ON binder_binding_sites(protein_id);")
            cur.execute("CREATE INDEX IF NOT EXISTS idx_binder_binding_sites_region ON binder_binding_sites(protein_id, region_start, region_end);")
            cur.execute("""
                CREATE OR REPLACE FUNCTION set_binder_binding_sites_updated_at()
                RETURNS TRIGGER AS $$
                BEGIN
                    NEW.updated_at = CURRENT_TIMESTAMP;
                    RETURN NEW;
                END;
                $$ LANGUAGE plpgsql;
            """)
            cur.execute("DROP TRIGGER IF EXISTS trg_binder_binding_sites_updated_at ON binder_binding_sites;")
            cur.execute("""
                CREATE TRIGGER trg_binder_binding_sites_updated_at
                BEFORE UPDATE ON binder_binding_sites
                FOR EACH ROW
                EXECUTE FUNCTION set_binder_binding_sites_updated_at();
            """)
        conn.commit()
    print("✅ binder_binding_sites table is ready for final visualization support.")
    print("Use this table for real or curated residue-level target-binder binding regions.")


if __name__ == "__main__":
    main()
