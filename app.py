import os
import re
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


def get_filter_options():
    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""
                SELECT DISTINCT binder_type
                FROM binders
                WHERE binder_type IS NOT NULL AND binder_type <> ''
                ORDER BY binder_type;
            """)
            binder_types = [row["binder_type"] for row in cur.fetchall()]

            cur.execute("""
                SELECT DISTINCT clinical_status
                FROM binders
                WHERE clinical_status IS NOT NULL AND clinical_status <> ''
                ORDER BY clinical_status;
            """)
            clinical_statuses = [row["clinical_status"] for row in cur.fetchall()]

            cur.execute("""
                SELECT disease_name
                FROM diseases
                ORDER BY disease_name;
            """)
            diseases = [row["disease_name"] for row in cur.fetchall()]

    return {
        "binder_types": binder_types,
        "clinical_statuses": clinical_statuses,
        "diseases": diseases
    }


def is_fasta_like(query: str) -> bool:
    """
    Detects simple FASTA / sequence-like input.
    """
    if not query:
        return False

    q = query.strip()

    if q.startswith(">"):
        return True

    compact = re.sub(r"\s+", "", q).upper()

    if len(compact) < 15:
        return False

    allowed = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*-")
    return all(char in allowed for char in compact)


def detect_query_intent(query: str) -> str:
    """
    Returns one of:
    - binder
    - protein
    - disease
    - sequence
    - general
    """
    q = (query or "").strip().lower()

    if not q:
        return "general"

    if is_fasta_like(query):
        return "sequence"

    binder_terms = [
        "binder", "binders", "antibody", "antibodies",
        "peptide", "peptides", "nanobody", "nanobodies", "vhh"
    ]
    disease_terms = [
        "cancer", "tumor", "tumour", "carcinoma", "disease", "oncology"
    ]

    if any(term in q for term in binder_terms):
        return "binder"

    if any(term in q for term in disease_terms):
        return "disease"

    # Short all-caps or gene-like queries often indicate protein/target intent
    compact = q.replace("-", "").replace("_", "")
    if len(q.split()) == 1 and len(compact) <= 10:
        return "protein"

    return "general"


def reorder_results_by_intent(results: dict, intent: str) -> dict:
    """
    Adds an ordering hint used by the UI.
    """
    if intent == "binder":
        results["primary_section"] = "binders"
    elif intent == "protein":
        results["primary_section"] = "proteins"
    elif intent == "disease":
        results["primary_section"] = "diseases"
    elif intent == "sequence":
        results["primary_section"] = "binders"
    else:
        results["primary_section"] = "binders"

    return results


def search_database(search_term: str, binder_type=None, clinical_status=None, disease_name=None):
    like_term = f"%{search_term}%"
    intent = detect_query_intent(search_term)

    results = {
        "diseases": [],
        "proteins": [],
        "binders": [],
        "trials": [],
        "intent": intent,
        "primary_section": "binders"
    }

    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:

            # --- DISEASES ---
            disease_sql = """
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
            """
            disease_params = [like_term, like_term, like_term, like_term, like_term, like_term, like_term]

            if disease_name:
                disease_sql += " AND d.disease_name = %s"
                disease_params.append(disease_name)

            disease_sql += " ORDER BY d.disease_name LIMIT 20;"
            cur.execute(disease_sql, tuple(disease_params))
            results["diseases"] = cur.fetchall()

            # --- PROTEINS ---
            protein_sql = """
                SELECT
                    p.protein_id,
                    p.protein_name,
                    p.uniprot_accession,
                    p.organism_name,
                    p.functional_description,
                    g.gene_symbol,
                    g.gene_name,
                    CASE
                        WHEN LOWER(g.gene_symbol) = LOWER(%s) THEN 1
                        WHEN LOWER(p.protein_name) = LOWER(%s) THEN 2
                        ELSE 3
                    END AS rank_score
                FROM proteins p
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                WHERE p.protein_name ILIKE %s
                   OR p.uniprot_accession ILIKE %s
                   OR p.functional_description ILIKE %s
                   OR g.gene_symbol ILIKE %s
                   OR g.gene_name ILIKE %s
                ORDER BY rank_score, p.protein_name
                LIMIT 20;
            """
            cur.execute(
                protein_sql,
                (
                    search_term, search_term,
                    like_term, like_term, like_term, like_term, like_term
                )
            )
            results["proteins"] = cur.fetchall()

            # --- BINDERS ---
            binder_sql = """
                SELECT DISTINCT
                    b.binder_id,
                    b.binder_name,
                    b.binder_type,
                    b.sequence,
                    b.clinical_status,
                    b.mechanism_of_action,
                    b.approval_status,
                    bm.modality_name,
                    CASE
                        WHEN LOWER(b.binder_name) = LOWER(%s) THEN 1
                        WHEN LOWER(b.binder_type) = LOWER(%s) THEN 2
                        ELSE 3
                    END AS rank_score
                FROM binders b
                LEFT JOIN binder_modalities bm ON b.modality_id = bm.modality_id
                LEFT JOIN binder_diseases bd ON b.binder_id = bd.binder_id
                LEFT JOIN diseases d ON bd.disease_id = d.disease_id
                WHERE (
                    b.binder_name ILIKE %s
                    OR b.binder_type ILIKE %s
                    OR b.clinical_status ILIKE %s
                    OR b.mechanism_of_action ILIKE %s
                    OR b.binder_description ILIKE %s
                    OR d.disease_name ILIKE %s
                )
            """
            binder_params = [
                search_term, search_term,
                like_term, like_term, like_term, like_term, like_term, like_term
            ]

            if binder_type:
                binder_sql += " AND b.binder_type = %s"
                binder_params.append(binder_type)

            if clinical_status:
                binder_sql += " AND b.clinical_status = %s"
                binder_params.append(clinical_status)

            if disease_name:
                binder_sql += " AND d.disease_name = %s"
                binder_params.append(disease_name)

            binder_sql += " ORDER BY rank_score, b.binder_name LIMIT 50;"
            cur.execute(binder_sql, tuple(binder_params))
            results["binders"] = cur.fetchall()

            # --- TRIALS ---
            trial_sql = """
                SELECT
                    trial_id,
                    nct_id,
                    trial_title,
                    condition_name,
                    phase,
                    recruitment_status,
                    CASE
                        WHEN LOWER(nct_id) = LOWER(%s) THEN 1
                        WHEN LOWER(trial_title) = LOWER(%s) THEN 2
                        ELSE 3
                    END AS rank_score
                FROM clinical_trials
                WHERE trial_title ILIKE %s
                   OR condition_name ILIKE %s
                   OR nct_id ILIKE %s
                   OR brief_summary ILIKE %s
                ORDER BY rank_score, trial_title
                LIMIT 20;
            """
            cur.execute(
                trial_sql,
                (search_term, search_term, like_term, like_term, like_term, like_term)
            )
            results["trials"] = cur.fetchall()

    return reorder_results_by_intent(results, intent)


def build_free_summary(query: str, results: dict, binder_type=None, clinical_status=None, disease_name=None) -> str:
    diseases = results.get("diseases", [])
    proteins = results.get("proteins", [])
    binders = results.get("binders", [])
    trials = results.get("trials", [])
    intent = results.get("intent", "general")

    disease_count = len(diseases)
    protein_count = len(proteins)
    binder_count = len(binders)
    trial_count = len(trials)

    if disease_count == 0 and protein_count == 0 and binder_count == 0 and trial_count == 0:
        return (
            f'No matching records were found for "{query}" in the current dataset. '
            "Try a different keyword or adjust the binder filters."
        )

    intent_labels = {
        "binder": "binder-focused",
        "protein": "target/protein-focused",
        "disease": "disease-focused",
        "sequence": "sequence-style",
        "general": "general"
    }

    filter_parts = []
    if binder_type:
        filter_parts.append(f"binder type = {binder_type}")
    if clinical_status:
        filter_parts.append(f"clinical status = {clinical_status}")
    if disease_name:
        filter_parts.append(f"disease = {disease_name}")

    filter_text = ""
    if filter_parts:
        filter_text = " Active binder filters: " + ", ".join(filter_parts) + "."

    parts = []

    parts.append(
        f'The search for "{query}" was interpreted as a {intent_labels.get(intent, "general")} query.'
    )

    if intent == "sequence":
        parts.append(
            "This looks like FASTA or sequence-like input, so the current MVP is prioritizing binder-style matches as a placeholder for future sequence similarity search."
        )

    parts.append(
        f'It returned {disease_count} disease record{"s" if disease_count != 1 else ""}, '
        f'{protein_count} protein target{"s" if protein_count != 1 else ""}, '
        f'{binder_count} therapeutic binder{"s" if binder_count != 1 else ""}, and '
        f'{trial_count} clinical trial{"s" if trial_count != 1 else ""}.'
        f'{filter_text}'
    )

    if intent == "protein" and proteins:
        top = ", ".join([p.get("gene_symbol") or p.get("protein_name") or "Unknown protein" for p in proteins[:3]])
        parts.append(f"Top target matches include {top}.")

    elif intent == "binder" and binders:
        top = ", ".join([b.get("binder_name") or "Unknown therapeutic" for b in binders[:3]])
        parts.append(f"Top binder matches include {top}.")

    elif intent == "disease" and diseases:
        top = ", ".join([d.get("disease_name") or "Unknown disease" for d in diseases[:3]])
        parts.append(f"Top disease matches include {top}.")

    else:
        if binders:
            top_binders = ", ".join([b.get("binder_name") or "Unknown therapeutic" for b in binders[:3]])
            parts.append(f"Top binder matches include {top_binders}.")
        if proteins:
            top_proteins = ", ".join([p.get("gene_symbol") or p.get("protein_name") or "Unknown protein" for p in proteins[:3]])
            parts.append(f"Protein-related matches include {top_proteins}.")

    if trials:
        top_trials = ", ".join([t.get("nct_id") or t.get("trial_title") or "Unknown trial" for t in trials[:2]])
        parts.append(f"Relevant trial matches include {top_trials}.")

    parts.append(
        "This summary is generated directly from retrieved database results currently stored in the platform."
    )

    return " ".join(parts)


@app.route("/", methods=["GET"])
def index():
    query = request.args.get("q", "").strip()
    binder_type = request.args.get("binder_type", "").strip() or None
    clinical_status = request.args.get("clinical_status", "").strip() or None
    disease_name = request.args.get("disease_name", "").strip() or None

    results = None
    error = None
    summary = None
    filter_options = get_filter_options()

    if query:
        try:
            results = search_database(
                query,
                binder_type=binder_type,
                clinical_status=clinical_status,
                disease_name=disease_name
            )
            summary = build_free_summary(
                query,
                results,
                binder_type=binder_type,
                clinical_status=clinical_status,
                disease_name=disease_name
            )
        except Exception as e:
            error = str(e)

    return render_template(
        "index.html",
        query=query,
        results=results,
        error=error,
        summary=summary,
        filter_options=filter_options,
        selected_binder_type=binder_type,
        selected_clinical_status=clinical_status,
        selected_disease_name=disease_name
    )


@app.route("/api/search", methods=["GET"])
def api_search():
    query = request.args.get("q", "").strip()
    binder_type = request.args.get("binder_type", "").strip() or None
    clinical_status = request.args.get("clinical_status", "").strip() or None
    disease_name = request.args.get("disease_name", "").strip() or None

    if not query:
        return jsonify({"error": "Missing search query"}), 400

    try:
        results = search_database(
            query,
            binder_type=binder_type,
            clinical_status=clinical_status,
            disease_name=disease_name
        )
        summary = build_free_summary(
            query,
            results,
            binder_type=binder_type,
            clinical_status=clinical_status,
            disease_name=disease_name
        )
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
                SELECT
                    d.disease_id,
                    d.disease_name,
                    d.disease_category
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
                    b.binder_type,
                    b.clinical_status,
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
                    b.binder_type,
                    b.sequence,
                    b.clinical_status,
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

            cur.execute("""
                SELECT
                    d.disease_id,
                    d.disease_name,
                    d.disease_category,
                    bd.tag_reason
                FROM binder_diseases bd
                JOIN diseases d ON bd.disease_id = d.disease_id
                WHERE bd.binder_id = %s
                ORDER BY d.disease_name;
            """, (binder_id,))
            diseases = cur.fetchall()

    return render_template(
        "binder_detail.html",
        binder=binder,
        proteins=proteins,
        trials=trials,
        structures=structures,
        diseases=diseases
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
                SELECT
                    d.disease_id,
                    d.disease_name,
                    d.disease_category
                FROM disease_trials dt
                JOIN diseases d ON dt.disease_id = d.disease_id
                WHERE dt.trial_id = %s
                ORDER BY d.disease_name;
            """, (trial_id,))
            diseases = cur.fetchall()

            cur.execute("""
                SELECT
                    b.binder_id,
                    b.binder_name,
                    b.binder_type,
                    b.clinical_status,
                    bm.modality_name
                FROM binder_trials bt
                JOIN binders b ON bt.binder_id = b.binder_id
                LEFT JOIN binder_modalities bm ON b.modality_id = bm.modality_id
                WHERE bt.trial_id = %s
                ORDER BY b.binder_name;
            """, (trial_id,))
            binders = cur.fetchall()

            cur.execute("""
                SELECT
                    p.protein_id,
                    p.protein_name,
                    g.gene_symbol
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
                    b.binder_type,
                    b.clinical_status,
                    b.mechanism_of_action,
                    b.approval_status,
                    bm.modality_name
                FROM binder_diseases bd
                JOIN binders b ON bd.binder_id = b.binder_id
                LEFT JOIN binder_modalities bm ON b.modality_id = bm.modality_id
                WHERE bd.disease_id = %s
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