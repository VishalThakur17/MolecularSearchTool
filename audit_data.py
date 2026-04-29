import os
from collections import Counter
from dotenv import load_dotenv
import psycopg2
import psycopg2.extras

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
    )


def pct(part, total):
    return round((part / total) * 100, 1) if total else 0


def print_counter(title, counter):
    print(f"\n=== {title} ===")
    if not counter:
        print("No data")
        return
    for label, count in counter.most_common():
        print(f"{label or 'Not available':30} | {count}")


def main():
    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""
                SELECT
                    COUNT(*) AS total,
                    COUNT(nct_id) AS with_nct,
                    COUNT(trial_title) AS with_title,
                    COUNT(condition_name) AS with_condition,
                    COUNT(phase) AS with_phase,
                    COUNT(recruitment_status) AS with_status,
                    COUNT(study_type) AS with_study_type,
                    COUNT(sponsor_name) AS with_sponsor,
                    COUNT(brief_summary) AS with_summary,
                    COUNT(trial_url) AS with_url
                FROM clinical_trials;
            """)
            row = cur.fetchone()
            total = row["total"] or 0

            print("\n=== CLINICAL TRIAL DATA COMPLETENESS ===")
            print(f"Total trials:       {total}")
            for key, label in [
                ("with_nct", "NCT ID"),
                ("with_title", "Title"),
                ("with_condition", "Condition"),
                ("with_phase", "Phase"),
                ("with_status", "Recruitment"),
                ("with_study_type", "Study type"),
                ("with_sponsor", "Sponsor"),
                ("with_summary", "Brief summary"),
                ("with_url", "External URL"),
            ]:
                value = row[key] or 0
                print(f"{label:17}: {value} / {total} ({pct(value, total)}%)")

            cur.execute("SELECT COALESCE(phase, 'Not available') AS phase, COUNT(*) AS count FROM clinical_trials GROUP BY COALESCE(phase, 'Not available') ORDER BY count DESC;")
            print_counter("TRIAL PHASE DISTRIBUTION", Counter({r["phase"]: r["count"] for r in cur.fetchall()}))

            cur.execute("SELECT COALESCE(recruitment_status, 'Not available') AS status, COUNT(*) AS count FROM clinical_trials GROUP BY COALESCE(recruitment_status, 'Not available') ORDER BY count DESC;")
            print_counter("TRIAL STATUS DISTRIBUTION", Counter({r["status"]: r["count"] for r in cur.fetchall()}))

            cur.execute("""
                SELECT
                    (SELECT COUNT(*) FROM binder_trials) AS binder_links,
                    (SELECT COUNT(*) FROM protein_trials) AS protein_links,
                    (SELECT COUNT(*) FROM trial_diseases) AS disease_links;
            """)
            links = cur.fetchone()
            print("\n=== CLINICAL EVIDENCE LINKS ===")
            print(f"Binder ↔ Trial links:  {links['binder_links']}")
            print(f"Protein ↔ Trial links: {links['protein_links']}")
            print(f"Disease ↔ Trial tags:  {links['disease_links']}")

            cur.execute("""
                SELECT ct.trial_id, ct.nct_id, ct.trial_title,
                       COUNT(DISTINCT bt.binder_id) AS binder_count,
                       COUNT(DISTINCT pt.protein_id) AS protein_count,
                       COUNT(DISTINCT td.disease_id) AS disease_count
                FROM clinical_trials ct
                LEFT JOIN binder_trials bt ON ct.trial_id = bt.trial_id
                LEFT JOIN protein_trials pt ON ct.trial_id = pt.trial_id
                LEFT JOIN trial_diseases td ON ct.trial_id = td.trial_id
                GROUP BY ct.trial_id, ct.nct_id, ct.trial_title
                ORDER BY binder_count DESC, protein_count DESC, disease_count DESC, ct.trial_title
                LIMIT 10;
            """)
            print("\n=== TOP LINKED TRIALS ===")
            for r in cur.fetchall():
                name = r["nct_id"] or r["trial_title"] or "Unknown trial"
                print(f"{name:18} | binders={r['binder_count']} | proteins={r['protein_count']} | disease tags={r['disease_count']}")


if __name__ == "__main__":
    main()
