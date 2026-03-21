import os
from flask import Flask, render_template, request, jsonify
from dotenv import load_dotenv
import psycopg2
import psycopg2.extras

load_dotenv()

app = Flask(__name__)

# --------------------------------------------------
# DATABASE CONFIG
# Set these as environment variables or replace them
# directly for quick testing.
# --------------------------------------------------
DB_HOST = os.getenv("DB_HOST", "YOUR_RDS_ENDPOINT")
DB_PORT = os.getenv("DB_PORT", "5432")
DB_NAME = os.getenv("DB_NAME", "molecular_search_db")
DB_USER = os.getenv("DB_USER", "postgres")
DB_PASSWORD = os.getenv("DB_PASSWORD", "YOUR_PASSWORD")


def get_connection():
    return psycopg2.connect(
        host=DB_HOST,
        port=DB_PORT,
        dbname=DB_NAME,
        user=DB_USER,
        password=DB_PASSWORD
    )


def search_database(search_term: str):
    """
    Search across diseases, proteins, binders, and clinical trials.
    Returns a dictionary grouped by entity type.
    """
    like_term = f"%{search_term}%"

    results = {
        "diseases": [],
        "proteins": [],
        "binders": [],
        "trials": []
    }

    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:

            # ----------------------------------------
            # Diseases
            # ----------------------------------------
            cur.execute("""
                SELECT
                    d.disease_id,
                    d.disease_name,
                    d.disease_category,
                    d.description
                FROM diseases d
                WHERE d.disease_name ILIKE %s
                   OR d.description ILIKE %s
                ORDER BY d.disease_name
                LIMIT 20;
            """, (like_term, like_term))
            results["diseases"] = cur.fetchall()

            # ----------------------------------------
            # Proteins
            # ----------------------------------------
            cur.execute("""
                SELECT
                    p.protein_id,
                    p.protein_name,
                    p.uniprot_accession,
                    p.organism_name,
                    p.functional_description,
                    g.gene_symbol,
                    g.gene_name
                FROM proteins p
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                WHERE p.protein_name ILIKE %s
                   OR p.uniprot_accession ILIKE %s
                   OR p.functional_description ILIKE %s
                   OR g.gene_symbol ILIKE %s
                   OR g.gene_name ILIKE %s
                ORDER BY p.protein_name
                LIMIT 20;
            """, (like_term, like_term, like_term, like_term, like_term))
            results["proteins"] = cur.fetchall()

            # ----------------------------------------
            # Binders
            # ----------------------------------------
            cur.execute("""
                SELECT
                    b.binder_id,
                    b.binder_name,
                    b.mechanism_of_action,
                    b.approval_status,
                    bm.modality_name
                FROM binders b
                LEFT JOIN binder_modalities bm ON b.modality_id = bm.modality_id
                WHERE b.binder_name ILIKE %s
                   OR b.mechanism_of_action ILIKE %s
                   OR b.binder_description ILIKE %s
                ORDER BY b.binder_name
                LIMIT 20;
            """, (like_term, like_term, like_term))
            results["binders"] = cur.fetchall()

            # ----------------------------------------
            # Clinical Trials
            # ----------------------------------------
            cur.execute("""
                SELECT
                    trial_id,
                    nct_id,
                    trial_title,
                    condition_name,
                    phase,
                    recruitment_status
                FROM clinical_trials
                WHERE trial_title ILIKE %s
                   OR condition_name ILIKE %s
                   OR nct_id ILIKE %s
                   OR brief_summary ILIKE %s
                ORDER BY trial_title
                LIMIT 20;
            """, (like_term, like_term, like_term, like_term))
            results["trials"] = cur.fetchall()

    return results


@app.route("/", methods=["GET"])
def index():
    query = request.args.get("q", "").strip()
    results = None
    error = None

    if query:
        try:
            results = search_database(query)
        except Exception as e:
            error = str(e)

    return render_template("index.html", query=query, results=results, error=error)


@app.route("/api/search", methods=["GET"])
def api_search():
    query = request.args.get("q", "").strip()

    if not query:
        return jsonify({"error": "Missing search query"}), 400

    try:
        results = search_database(query)
        return jsonify(results)
    except Exception as e:
        return jsonify({"error": str(e)}), 500


if __name__ == "__main__":
    app.run(debug=True)