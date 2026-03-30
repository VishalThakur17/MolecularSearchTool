import os
from flask import Flask, render_template, request, jsonify
from dotenv import load_dotenv
import psycopg2
import psycopg2.extras

load_dotenv()

app = Flask(__name__)

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
        sslmode="require"
    )


def search_database(search_term: str):
    like_term = f"%{search_term}%"

    results = {
        "diseases": [],
        "proteins": [],
        "binders": [],
        "trials": []
    }

    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:

            cur.execute("""
                SELECT DISTINCT
                    d.disease_id,
                    d.disease_name,
                    d.disease_category,
                    d.description
                FROM diseases d
                LEFT JOIN disease_proteins dp ON d.disease_id = dp.disease_id
                LEFT JOIN proteins p ON dp.protein_id = p.protein_id
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                LEFT JOIN protein_binders pb ON p.protein_id = pb.protein_id
                LEFT JOIN binders b ON pb.binder_id = b.binder_id
                LEFT JOIN disease_trials dt ON d.disease_id = dt.disease_id
                LEFT JOIN clinical_trials ct ON dt.trial_id = ct.trial_id
                WHERE d.disease_name ILIKE %s
                   OR d.description ILIKE %s
                   OR p.protein_name ILIKE %s
                   OR g.gene_symbol ILIKE %s
                   OR b.binder_name ILIKE %s
                   OR ct.trial_title ILIKE %s
                   OR ct.condition_name ILIKE %s
                ORDER BY d.disease_name
                LIMIT 20;
            """, (like_term, like_term, like_term, like_term, like_term, like_term, like_term))
            results["diseases"] = cur.fetchall()

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


def build_free_summary(query: str, results: dict) -> str:
    """
    Free lightweight RAG-style summary generated from database results only.
    No LLM required.
    """
    diseases = results.get("diseases", [])
    proteins = results.get("proteins", [])
    binders = results.get("binders", [])
    trials = results.get("trials", [])

    disease_count = len(diseases)
    protein_count = len(proteins)
    binder_count = len(binders)
    trial_count = len(trials)

    if disease_count == 0 and protein_count == 0 and binder_count == 0 and trial_count == 0:
        return (
            f'No matching records were found for "{query}" in the current dataset. '
            "Try searching a different disease, protein target, therapeutic, or trial keyword."
        )

    parts = []

    parts.append(
        f'The search for "{query}" returned '
        f'{disease_count} disease record{"s" if disease_count != 1 else ""}, '
        f'{protein_count} protein target{"s" if protein_count != 1 else ""}, '
        f'{binder_count} therapeutic binder{"s" if binder_count != 1 else ""}, and '
        f'{trial_count} clinical trial{"s" if trial_count != 1 else ""}.'
    )

    if proteins:
        top_proteins = ", ".join(
            [
                p.get("gene_symbol") or p.get("protein_name") or "Unknown protein"
                for p in proteins[:3]
            ]
        )
        parts.append(
            f'The strongest protein-related matches include {top_proteins}.'
        )

    if diseases:
        top_diseases = ", ".join(
            [d.get("disease_name") or "Unknown disease" for d in diseases[:3]]
        )
        parts.append(
            f'Related disease context in the current dataset includes {top_diseases}.'
        )

    if binders:
        top_binders = ", ".join(
            [b.get("binder_name") or "Unknown therapeutic" for b in binders[:3]]
        )
        parts.append(
            f'Linked therapeutic or binder matches include {top_binders}.'
        )

    if trials:
        top_trials = ", ".join(
            [
                t.get("nct_id") or t.get("trial_title") or "Unknown trial"
                for t in trials[:2]
            ]
        )
        parts.append(
            f'Relevant clinical study records include {top_trials}.'
        )

    parts.append(
        "This summary is generated directly from retrieved database results, so it reflects only the records currently stored in the platform."
    )

    return " ".join(parts)


@app.route("/", methods=["GET"])
def index():
    query = request.args.get("q", "").strip()
    results = None
    error = None
    summary = None

    if query:
        try:
            results = search_database(query)
            summary = build_free_summary(query, results)
        except Exception as e:
            error = str(e)

    return render_template(
        "index.html",
        query=query,
        results=results,
        error=error,
        summary=summary
    )


@app.route("/api/search", methods=["GET"])
def api_search():
    query = request.args.get("q", "").strip()

    if not query:
        return jsonify({"error": "Missing search query"}), 400

    try:
        results = search_database(query)
        summary = build_free_summary(query, results)
        return jsonify({
            "summary": summary,
            "results": results
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/protein/<int:protein_id>")
def protein_detail(protein_id):
    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""
                SELECT
                    p.protein_id,
                    p.protein_name,
                    p.uniprot_accession,
                    p.organism_name,
                    p.sequence_length,
                    p.functional_description,
                    p.subcellular_location,
                    g.gene_symbol,
                    g.gene_name
                FROM proteins p
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                WHERE p.protein_id = %s;
            """, (protein_id,))
            protein = cur.fetchone()

            if not protein:
                return "Protein not found", 404

            cur.execute("""
                SELECT d.disease_id, d.disease_name, d.disease_category
                FROM disease_proteins dp
                JOIN diseases d ON dp.disease_id = d.disease_id
                WHERE dp.protein_id = %s
                ORDER BY d.disease_name;
            """, (protein_id,))
            diseases = cur.fetchall()

            cur.execute("""
                SELECT
                    b.binder_id,
                    b.binder_name,
                    b.mechanism_of_action,
                    b.approval_status,
                    bm.modality_name,
                    pb.interaction_type
                FROM protein_binders pb
                JOIN binders b ON pb.binder_id = b.binder_id
                LEFT JOIN binder_modalities bm ON b.modality_id = bm.modality_id
                WHERE pb.protein_id = %s
                ORDER BY b.binder_name;
            """, (protein_id,))
            binders = cur.fetchall()

            cur.execute("""
                SELECT
                    ct.trial_id,
                    ct.nct_id,
                    ct.trial_title,
                    ct.phase,
                    ct.recruitment_status,
                    ct.condition_name
                FROM protein_trials pt
                JOIN clinical_trials ct ON pt.trial_id = ct.trial_id
                WHERE pt.protein_id = %s
                ORDER BY ct.trial_title;
            """, (protein_id,))
            trials = cur.fetchall()

            cur.execute("""
                SELECT
                    s.structure_id,
                    s.pdb_id,
                    s.structure_title,
                    s.experimental_method,
                    s.resolution
                FROM protein_structures ps
                JOIN structures s ON ps.structure_id = s.structure_id
                WHERE ps.protein_id = %s
                ORDER BY s.pdb_id;
            """, (protein_id,))
            structures = cur.fetchall()

    return render_template(
        "protein_detail.html",
        protein=protein,
        diseases=diseases,
        binders=binders,
        trials=trials,
        structures=structures
    )


@app.route("/binder/<int:binder_id>")
def binder_detail(binder_id):
    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""
                SELECT
                    b.binder_id,
                    b.binder_name,
                    b.binder_description,
                    b.mechanism_of_action,
                    b.approval_status,
                    b.developer_company,
                    bm.modality_name
                FROM binders b
                LEFT JOIN binder_modalities bm ON b.modality_id = bm.modality_id
                WHERE b.binder_id = %s;
            """, (binder_id,))
            binder = cur.fetchone()

            if not binder:
                return "Binder not found", 404

            cur.execute("""
                SELECT
                    p.protein_id,
                    p.protein_name,
                    g.gene_symbol,
                    pb.interaction_type
                FROM protein_binders pb
                JOIN proteins p ON pb.protein_id = p.protein_id
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                WHERE pb.binder_id = %s
                ORDER BY p.protein_name;
            """, (binder_id,))
            proteins = cur.fetchall()

            cur.execute("""
                SELECT
                    ct.trial_id,
                    ct.nct_id,
                    ct.trial_title,
                    ct.phase,
                    ct.recruitment_status
                FROM binder_trials bt
                JOIN clinical_trials ct ON bt.trial_id = ct.trial_id
                WHERE bt.binder_id = %s
                ORDER BY ct.trial_title;
            """, (binder_id,))
            trials = cur.fetchall()

            cur.execute("""
                SELECT
                    s.structure_id,
                    s.pdb_id,
                    s.structure_title,
                    s.experimental_method,
                    s.resolution
                FROM binder_structures bs
                JOIN structures s ON bs.structure_id = s.structure_id
                WHERE bs.binder_id = %s
                ORDER BY s.pdb_id;
            """, (binder_id,))
            structures = cur.fetchall()

    return render_template(
        "binder_detail.html",
        binder=binder,
        proteins=proteins,
        trials=trials,
        structures=structures
    )


@app.route("/trial/<int:trial_id>")
def trial_detail(trial_id):
    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""
                SELECT
                    trial_id,
                    nct_id,
                    trial_title,
                    condition_name,
                    phase,
                    recruitment_status,
                    study_type,
                    sponsor_name,
                    brief_summary,
                    trial_url,
                    start_date,
                    completion_date
                FROM clinical_trials
                WHERE trial_id = %s;
            """, (trial_id,))
            trial = cur.fetchone()

            if not trial:
                return "Trial not found", 404

            cur.execute("""
                SELECT d.disease_id, d.disease_name, d.disease_category
                FROM disease_trials dt
                JOIN diseases d ON dt.disease_id = d.disease_id
                WHERE dt.trial_id = %s
                ORDER BY d.disease_name;
            """, (trial_id,))
            diseases = cur.fetchall()

            cur.execute("""
                SELECT b.binder_id, b.binder_name, bm.modality_name
                FROM binder_trials bt
                JOIN binders b ON bt.binder_id = b.binder_id
                LEFT JOIN binder_modalities bm ON b.modality_id = bm.modality_id
                WHERE bt.trial_id = %s
                ORDER BY b.binder_name;
            """, (trial_id,))
            binders = cur.fetchall()

            cur.execute("""
                SELECT p.protein_id, p.protein_name, g.gene_symbol
                FROM protein_trials pt
                JOIN proteins p ON pt.protein_id = p.protein_id
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                WHERE pt.trial_id = %s
                ORDER BY p.protein_name;
            """, (trial_id,))
            proteins = cur.fetchall()

    return render_template(
        "trial_detail.html",
        trial=trial,
        diseases=diseases,
        binders=binders,
        proteins=proteins
    )


@app.route("/disease/<int:disease_id>")
def disease_detail(disease_id):
    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""
                SELECT
                    disease_id,
                    disease_name,
                    disease_category,
                    description,
                    mesh_id,
                    doid,
                    ncit_code
                FROM diseases
                WHERE disease_id = %s;
            """, (disease_id,))
            disease = cur.fetchone()

            if not disease:
                return "Disease not found", 404

            cur.execute("""
                SELECT
                    p.protein_id,
                    p.protein_name,
                    g.gene_symbol
                FROM disease_proteins dp
                JOIN proteins p ON dp.protein_id = p.protein_id
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                WHERE dp.disease_id = %s
                ORDER BY p.protein_name;
            """, (disease_id,))
            proteins = cur.fetchall()

            cur.execute("""
                SELECT
                    ct.trial_id,
                    ct.nct_id,
                    ct.trial_title,
                    ct.phase,
                    ct.recruitment_status
                FROM disease_trials dt
                JOIN clinical_trials ct ON dt.trial_id = ct.trial_id
                WHERE dt.disease_id = %s
                ORDER BY ct.trial_title;
            """, (disease_id,))
            trials = cur.fetchall()

            cur.execute("""
                SELECT DISTINCT
                    b.binder_id,
                    b.binder_name,
                    b.mechanism_of_action,
                    b.approval_status,
                    bm.modality_name
                FROM disease_proteins dp
                JOIN protein_binders pb ON dp.protein_id = pb.protein_id
                JOIN binders b ON pb.binder_id = b.binder_id
                LEFT JOIN binder_modalities bm ON b.modality_id = bm.modality_id
                WHERE dp.disease_id = %s
                ORDER BY b.binder_name;
            """, (disease_id,))
            binders = cur.fetchall()

    return render_template(
        "disease_detail.html",
        disease=disease,
        proteins=proteins,
        trials=trials,
        binders=binders
    )


if __name__ == "__main__":
    app.run(debug=True)