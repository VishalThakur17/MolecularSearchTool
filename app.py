import os
import re
from collections import Counter
from flask import Flask, render_template, request, jsonify, redirect, url_for
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
                WHERE binder_type IS NOT NULL
                  AND TRIM(binder_type) <> ''
                  AND binder_type IN ('IgG', 'VHH', 'Peptide', 'Small Molecule', 'Other')
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


def normalize_text(value: str) -> str:
    value = (value or "").strip().lower()
    value = re.sub(r"[_\-]+", " ", value)
    value = re.sub(r"\s+", " ", value)
    return value


def safe_string(value):
    return (value or "").strip()


def summarize_counter(counter_obj):
    return [
        {"label": label, "count": count}
        for label, count in counter_obj.most_common()
    ]


def build_molstar_pdb_viewer_url(pdb_id: str) -> str:
    pdb_id = safe_string(pdb_id)
    if not pdb_id:
        return ""
    return f"https://molstar.org/viewer/?pdb={pdb_id.lower()}&hide-controls=0&collapse-left-panel=1&pdb-provider=rcsb"


def build_molstar_afdb_viewer_url(uniprot_accession: str) -> str:
    accession = safe_string(uniprot_accession)
    if not accession:
        return ""
    return f"https://molstar.org/viewer/?afdb={accession}&hide-controls=0&collapse-left-panel=1"


def dedupe_structure_candidates(candidates):
    seen = set()
    deduped = []
    for item in candidates:
        key = (item.get("source_type"), item.get("source_id"))
        if key in seen:
            continue
        seen.add(key)
        deduped.append(item)
    return deduped


def build_structure_candidates(structures, fallback_uniprot=None, fallback_title=None):
    candidates = []

    for structure in structures:
        pdb_id = safe_string(structure.get("pdb_id"))
        if not pdb_id:
            continue

        title = safe_string(structure.get("structure_title")) or f"PDB {pdb_id.upper()}"
        experimental_method = safe_string(structure.get("experimental_method")) or "Experimental structure"
        resolution = structure.get("resolution")

        resolution_text = ""
        if resolution is not None:
            resolution_text = f" · Resolution: {resolution}"

        candidates.append({
            "source_type": "pdb",
            "source_id": pdb_id.upper(),
            "label": f"PDB {pdb_id.upper()}",
            "title": title,
            "subtitle": f"{experimental_method}{resolution_text}",
            "viewer_url": build_molstar_pdb_viewer_url(pdb_id),
            "external_url": f"https://www.rcsb.org/structure/{pdb_id.upper()}",
            "origin": "Linked structure"
        })

    if not candidates and safe_string(fallback_uniprot):
        accession = safe_string(fallback_uniprot)
        candidates.append({
            "source_type": "afdb",
            "source_id": accession,
            "label": f"AlphaFold {accession}",
            "title": safe_string(fallback_title) or f"Predicted model for {accession}",
            "subtitle": "Computed structure model",
            "viewer_url": build_molstar_afdb_viewer_url(accession),
            "external_url": f"https://alphafold.ebi.ac.uk/entry/{accession}",
            "origin": "Fallback model"
        })

    return dedupe_structure_candidates(candidates)


def detect_query_intent(query: str) -> str:
    q = normalize_text(query)

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
    trial_terms = [
        "trial", "trials", "clinical trial", "clinical trials", "nct"
    ]

    if any(term in q for term in trial_terms):
        return "trial"

    if any(term in q for term in binder_terms):
        return "binder"

    if any(term in q for term in disease_terms):
        return "disease"

    compact = q.replace(" ", "")
    if len(q.split()) == 1 and len(compact) <= 12:
        return "protein"

    return "general"


def extract_focus_query(query: str, intent: str) -> str:
    q = normalize_text(query)
    if not q:
        return query

    stopwords = {
        "show", "find", "search", "for", "all", "the", "a", "an", "of",
        "with", "related", "about", "associated", "linked", "known", "current"
    }

    intent_words = {
        "binder": {"binder", "binders", "antibody", "antibodies", "peptide", "peptides", "nanobody", "nanobodies", "vhh"},
        "protein": {"protein", "proteins", "target", "targets", "gene", "genes"},
        "trial": {"trial", "trials", "clinical", "nct", "study", "studies"},
        "disease": {"disease", "diseases"}
    }

    words = []
    for token in q.split():
        if token in stopwords:
            continue
        if token in intent_words.get(intent, set()):
            continue
        words.append(token)

    cleaned = " ".join(words).strip()
    return cleaned if cleaned else query.strip()


def get_exact_protein_match(query: str):
    q = normalize_text(query)
    if not q:
        return None

    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""
                SELECT
                    p.protein_id,
                    p.protein_name,
                    g.gene_symbol
                FROM proteins p
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                WHERE LOWER(TRIM(COALESCE(g.gene_symbol, ''))) = %s
                   OR LOWER(TRIM(COALESCE(p.protein_name, ''))) = %s
                   OR LOWER(TRIM(COALESCE(p.uniprot_accession, ''))) = %s
                ORDER BY
                    CASE
                        WHEN LOWER(TRIM(COALESCE(g.gene_symbol, ''))) = %s THEN 1
                        WHEN LOWER(TRIM(COALESCE(p.protein_name, ''))) = %s THEN 2
                        ELSE 3
                    END
                LIMIT 1;
            """, (q, q, q, q, q))
            return cur.fetchone()


def get_exact_binder_match(query: str):
    q = normalize_text(query)
    if not q:
        return None

    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""
                SELECT
                    binder_id,
                    binder_name,
                    binder_type
                FROM binders
                WHERE LOWER(TRIM(COALESCE(binder_name, ''))) = %s
                ORDER BY binder_name
                LIMIT 1;
            """, (q,))
            return cur.fetchone()


def get_exact_disease_match(query: str):
    q = normalize_text(query)
    if not q:
        return None

    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""
                SELECT
                    disease_id,
                    disease_name,
                    disease_category
                FROM diseases
                WHERE LOWER(TRIM(COALESCE(disease_name, ''))) = %s
                ORDER BY disease_name
                LIMIT 1;
            """, (q,))
            return cur.fetchone()


def get_exact_trial_match(query: str):
    q = normalize_text(query)
    if not q:
        return None

    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""
                SELECT
                    trial_id,
                    nct_id,
                    trial_title
                FROM clinical_trials
                WHERE LOWER(TRIM(COALESCE(nct_id, ''))) = %s
                   OR LOWER(TRIM(COALESCE(trial_title, ''))) = %s
                ORDER BY
                    CASE
                        WHEN LOWER(TRIM(COALESCE(nct_id, ''))) = %s THEN 1
                        ELSE 2
                    END
                LIMIT 1;
            """, (q, q, q))
            return cur.fetchone()


def reorder_results_by_intent(results: dict, intent: str) -> dict:
    if intent == "binder":
        results["primary_section"] = "binders"
        results["section_order"] = ["binders", "proteins", "trials", "diseases"]
    elif intent == "protein":
        results["primary_section"] = "proteins"
        results["section_order"] = ["proteins", "binders", "trials", "diseases"]
    elif intent == "trial":
        results["primary_section"] = "trials"
        results["section_order"] = ["trials", "binders", "proteins", "diseases"]
    elif intent == "disease":
        results["primary_section"] = "binders"
        results["section_order"] = ["binders", "proteins", "trials", "diseases"]
    elif intent == "sequence":
        results["primary_section"] = "binders"
        results["section_order"] = ["binders", "proteins", "trials", "diseases"]
    else:
        results["primary_section"] = "binders"
        results["section_order"] = ["binders", "proteins", "trials", "diseases"]

    return results


def decide_search_route(query: str, binder_type=None, clinical_status=None, disease_name=None):
    q = (query or "").strip()
    if not q:
        return {
            "mode": "results",
            "intent": "general",
            "effective_query": q,
            "route_hint": None,
            "auto_disease_name": disease_name
        }

    intent = detect_query_intent(q)
    focus_query = extract_focus_query(q, intent)
    filters_active = any([binder_type, clinical_status, disease_name])

    if intent == "sequence":
        return {
            "mode": "results",
            "intent": intent,
            "effective_query": q,
            "route_hint": "Sequence-like input detected. Showing binder-first matches as a placeholder for future sequence search.",
            "auto_disease_name": disease_name
        }

    if intent == "disease":
        exact_disease = get_exact_disease_match(q) or get_exact_disease_match(focus_query)
        if exact_disease and not disease_name:
            return {
                "mode": "results",
                "intent": intent,
                "effective_query": focus_query if focus_query else exact_disease["disease_name"],
                "route_hint": f'Exact disease match found for "{exact_disease["disease_name"]}". Results are filtered using that disease tag.',
                "auto_disease_name": exact_disease["disease_name"]
            }

        return {
            "mode": "results",
            "intent": intent,
            "effective_query": focus_query if focus_query else q,
            "route_hint": "Disease-style query detected. Showing binder-first filtered results.",
            "auto_disease_name": disease_name
        }

    if intent == "binder":
        return {
            "mode": "results",
            "intent": intent,
            "effective_query": focus_query if focus_query else q,
            "route_hint": "Binder-focused query detected. Showing binder-first results.",
            "auto_disease_name": disease_name
        }

    if not filters_active:
        exact_trial = get_exact_trial_match(q)
        if exact_trial:
            return {
                "mode": "redirect",
                "intent": "trial",
                "endpoint": "trial_detail",
                "values": {"trial_id": exact_trial["trial_id"]},
                "route_hint": None,
                "effective_query": q,
                "auto_disease_name": disease_name
            }

        exact_binder = get_exact_binder_match(q)
        if exact_binder:
            return {
                "mode": "redirect",
                "intent": "binder",
                "endpoint": "binder_detail",
                "values": {"binder_id": exact_binder["binder_id"]},
                "route_hint": None,
                "effective_query": q,
                "auto_disease_name": disease_name
            }

        exact_protein = get_exact_protein_match(q)
        if exact_protein:
            return {
                "mode": "redirect",
                "intent": "protein",
                "endpoint": "protein_detail",
                "values": {"protein_id": exact_protein["protein_id"]},
                "route_hint": None,
                "effective_query": q,
                "auto_disease_name": disease_name
            }

    return {
        "mode": "results",
        "intent": intent,
        "effective_query": focus_query if focus_query else q,
        "route_hint": None,
        "auto_disease_name": disease_name
    }


def enrich_search_results(results: dict):
    """Add exploration metadata to cards on the search results page."""
    binders = results.get("binders", [])
    proteins = results.get("proteins", [])
    trials = results.get("trials", [])
    diseases = results.get("diseases", [])

    binder_ids = [b["binder_id"] for b in binders if b.get("binder_id")]
    protein_ids = [p["protein_id"] for p in proteins if p.get("protein_id")]
    trial_ids = [t["trial_id"] for t in trials if t.get("trial_id")]
    disease_ids = [d["disease_id"] for d in diseases if d.get("disease_id")]

    if not any([binder_ids, protein_ids, trial_ids, disease_ids]):
        return results

    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            binder_targets_map = {}
            binder_diseases_map = {}
            protein_binder_count_map = {}
            protein_trial_count_map = {}
            trial_binder_count_map = {}
            trial_disease_names_map = {}
            disease_binder_count_map = {}
            disease_trial_count_map = {}

            if binder_ids:
                cur.execute("""
                    SELECT
                        pb.binder_id,
                        p.protein_id,
                        p.protein_name,
                        g.gene_symbol
                    FROM protein_binders pb
                    JOIN proteins p ON pb.protein_id = p.protein_id
                    LEFT JOIN genes g ON p.gene_id = g.gene_id
                    WHERE pb.binder_id = ANY(%s)
                    ORDER BY p.protein_name;
                """, (binder_ids,))
                for row in cur.fetchall():
                    binder_targets_map.setdefault(row["binder_id"], []).append({
                        "protein_id": row["protein_id"],
                        "label": row["gene_symbol"] or row["protein_name"] or "Unknown target"
                    })

                cur.execute("""
                    SELECT
                        bd.binder_id,
                        d.disease_name
                    FROM binder_diseases bd
                    JOIN diseases d ON bd.disease_id = d.disease_id
                    WHERE bd.binder_id = ANY(%s)
                    ORDER BY d.disease_name;
                """, (binder_ids,))
                for row in cur.fetchall():
                    binder_diseases_map.setdefault(row["binder_id"], []).append(row["disease_name"])

            if protein_ids:
                cur.execute("""
                    SELECT
                        pb.protein_id,
                        COUNT(DISTINCT pb.binder_id) AS binder_count
                    FROM protein_binders pb
                    WHERE pb.protein_id = ANY(%s)
                    GROUP BY pb.protein_id;
                """, (protein_ids,))
                for row in cur.fetchall():
                    protein_binder_count_map[row["protein_id"]] = row["binder_count"]

                cur.execute("""
                    SELECT
                        pt.protein_id,
                        COUNT(DISTINCT pt.trial_id) AS trial_count
                    FROM protein_trials pt
                    WHERE pt.protein_id = ANY(%s)
                    GROUP BY pt.protein_id;
                """, (protein_ids,))
                for row in cur.fetchall():
                    protein_trial_count_map[row["protein_id"]] = row["trial_count"]

            if trial_ids:
                cur.execute("""
                    SELECT
                        bt.trial_id,
                        COUNT(DISTINCT bt.binder_id) AS binder_count
                    FROM binder_trials bt
                    WHERE bt.trial_id = ANY(%s)
                    GROUP BY bt.trial_id;
                """, (trial_ids,))
                for row in cur.fetchall():
                    trial_binder_count_map[row["trial_id"]] = row["binder_count"]

                cur.execute("""
                    SELECT
                        td.trial_id,
                        d.disease_name
                    FROM trial_diseases td
                    JOIN diseases d ON td.disease_id = d.disease_id
                    WHERE td.trial_id = ANY(%s)
                    ORDER BY d.disease_name;
                """, (trial_ids,))
                for row in cur.fetchall():
                    trial_disease_names_map.setdefault(row["trial_id"], []).append(row["disease_name"])

            if disease_ids:
                cur.execute("""
                    SELECT
                        bd.disease_id,
                        COUNT(DISTINCT bd.binder_id) AS binder_count
                    FROM binder_diseases bd
                    WHERE bd.disease_id = ANY(%s)
                    GROUP BY bd.disease_id;
                """, (disease_ids,))
                for row in cur.fetchall():
                    disease_binder_count_map[row["disease_id"]] = row["binder_count"]

                cur.execute("""
                    SELECT
                        td.disease_id,
                        COUNT(DISTINCT td.trial_id) AS trial_count
                    FROM trial_diseases td
                    WHERE td.disease_id = ANY(%s)
                    GROUP BY td.disease_id;
                """, (disease_ids,))
                for row in cur.fetchall():
                    disease_trial_count_map[row["disease_id"]] = row["trial_count"]

    for binder in binders:
        binder_id = binder.get("binder_id")
        targets = binder_targets_map.get(binder_id, [])
        diseases_for_binder = binder_diseases_map.get(binder_id, [])
        binder["linked_targets_preview"] = targets[:3]
        binder["linked_target_count"] = len(targets)
        binder["linked_diseases_preview"] = diseases_for_binder[:2]

    for protein in proteins:
        protein_id = protein.get("protein_id")
        protein["linked_binder_count"] = protein_binder_count_map.get(protein_id, 0)
        protein["linked_trial_count"] = protein_trial_count_map.get(protein_id, 0)

    for trial in trials:
        trial_id = trial.get("trial_id")
        trial["linked_binder_count"] = trial_binder_count_map.get(trial_id, 0)
        trial["linked_diseases_preview"] = trial_disease_names_map.get(trial_id, [])[:2]

    for disease in diseases:
        disease_id = disease.get("disease_id")
        disease["linked_binder_count"] = disease_binder_count_map.get(disease_id, 0)
        disease["linked_trial_count"] = disease_trial_count_map.get(disease_id, 0)

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
        "primary_section": "binders",
        "section_order": ["binders", "proteins", "trials", "diseases"]
    }

    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            disease_sql = """
                SELECT DISTINCT
                    d.disease_id,
                    d.disease_name,
                    d.disease_category,
                    d.description
                FROM diseases d
                LEFT JOIN protein_diseases pd ON d.disease_id = pd.disease_id
                LEFT JOIN proteins p ON pd.protein_id = p.protein_id
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                LEFT JOIN binder_diseases bd ON d.disease_id = bd.disease_id
                LEFT JOIN binders b ON bd.binder_id = b.binder_id
                LEFT JOIN trial_diseases td ON d.disease_id = td.disease_id
                LEFT JOIN clinical_trials ct ON td.trial_id = ct.trial_id
                WHERE (
                    d.disease_name ILIKE %s
                    OR d.description ILIKE %s
                    OR p.protein_name ILIKE %s
                    OR g.gene_symbol ILIKE %s
                    OR b.binder_name ILIKE %s
                    OR ct.trial_title ILIKE %s
                    OR ct.condition_name ILIKE %s
                )
            """
            disease_params = [like_term, like_term, like_term, like_term, like_term, like_term, like_term]

            if disease_name:
                disease_sql += " AND d.disease_name = %s"
                disease_params.append(disease_name)

            disease_sql += " ORDER BY d.disease_name LIMIT 20;"
            cur.execute(disease_sql, tuple(disease_params))
            results["diseases"] = cur.fetchall()

            protein_sql = """
                SELECT DISTINCT
                    p.protein_id,
                    p.protein_name,
                    p.uniprot_accession,
                    p.organism_name,
                    p.functional_description,
                    g.gene_symbol,
                    g.gene_name,
                    CASE
                        WHEN LOWER(COALESCE(g.gene_symbol, '')) = LOWER(%s) THEN 1
                        WHEN LOWER(COALESCE(p.protein_name, '')) = LOWER(%s) THEN 2
                        ELSE 3
                    END AS rank_score
                FROM proteins p
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                LEFT JOIN protein_diseases pd ON p.protein_id = pd.protein_id
                LEFT JOIN diseases d ON pd.disease_id = d.disease_id
                WHERE (
                    p.protein_name ILIKE %s
                    OR p.uniprot_accession ILIKE %s
                    OR p.functional_description ILIKE %s
                    OR g.gene_symbol ILIKE %s
                    OR g.gene_name ILIKE %s
                    OR d.disease_name ILIKE %s
                )
            """
            protein_params = [
                search_term, search_term,
                like_term, like_term, like_term, like_term, like_term, like_term
            ]

            if disease_name:
                protein_sql += " AND d.disease_name = %s"
                protein_params.append(disease_name)

            protein_sql += " ORDER BY rank_score, p.protein_name LIMIT 20;"
            cur.execute(protein_sql, tuple(protein_params))
            results["proteins"] = cur.fetchall()

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
                        WHEN LOWER(COALESCE(b.binder_name, '')) = LOWER(%s) THEN 1
                        WHEN LOWER(COALESCE(b.binder_type, '')) = LOWER(%s) THEN 2
                        ELSE 3
                    END AS rank_score
                FROM binders b
                LEFT JOIN binder_modalities bm ON b.modality_id = bm.modality_id
                LEFT JOIN binder_diseases bd ON b.binder_id = bd.binder_id
                LEFT JOIN diseases d ON bd.disease_id = d.disease_id
                LEFT JOIN protein_binders pb ON b.binder_id = pb.binder_id
                LEFT JOIN proteins p ON pb.protein_id = p.protein_id
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                WHERE (
                    b.binder_name ILIKE %s
                    OR b.binder_type ILIKE %s
                    OR b.clinical_status ILIKE %s
                    OR b.mechanism_of_action ILIKE %s
                    OR b.binder_description ILIKE %s
                    OR d.disease_name ILIKE %s
                    OR p.protein_name ILIKE %s
                    OR g.gene_symbol ILIKE %s
                )
            """
            binder_params = [
                search_term, search_term,
                like_term, like_term, like_term, like_term, like_term, like_term, like_term, like_term
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

            trial_sql = """
                SELECT DISTINCT
                    ct.trial_id,
                    ct.nct_id,
                    ct.trial_title,
                    ct.condition_name,
                    ct.phase,
                    ct.recruitment_status,
                    CASE
                        WHEN LOWER(COALESCE(ct.nct_id, '')) = LOWER(%s) THEN 1
                        WHEN LOWER(COALESCE(ct.trial_title, '')) = LOWER(%s) THEN 2
                        ELSE 3
                    END AS rank_score
                FROM clinical_trials ct
                LEFT JOIN trial_diseases td ON ct.trial_id = td.trial_id
                LEFT JOIN diseases d ON td.disease_id = d.disease_id
                LEFT JOIN binder_trials bt ON ct.trial_id = bt.trial_id
                LEFT JOIN binders b ON bt.binder_id = b.binder_id
                LEFT JOIN protein_trials pt ON ct.trial_id = pt.trial_id
                LEFT JOIN proteins p ON pt.protein_id = p.protein_id
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                WHERE (
                    ct.trial_title ILIKE %s
                    OR ct.condition_name ILIKE %s
                    OR ct.nct_id ILIKE %s
                    OR ct.brief_summary ILIKE %s
                    OR d.disease_name ILIKE %s
                    OR b.binder_name ILIKE %s
                    OR p.protein_name ILIKE %s
                    OR g.gene_symbol ILIKE %s
                )
            """
            trial_params = [
                search_term, search_term,
                like_term, like_term, like_term, like_term, like_term, like_term, like_term, like_term
            ]

            if disease_name:
                trial_sql += " AND d.disease_name = %s"
                trial_params.append(disease_name)

            trial_sql += " ORDER BY rank_score, ct.trial_title LIMIT 20;"
            cur.execute(trial_sql, tuple(trial_params))
            results["trials"] = cur.fetchall()

    results = reorder_results_by_intent(results, intent)
    results = enrich_search_results(results)
    return results


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
            "Try a different keyword or adjust the filters."
        )

    intent_labels = {
        "binder": "binder-focused",
        "protein": "target/protein-focused",
        "disease": "disease-filtered",
        "trial": "clinical-trial-focused",
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
        filter_text = " Active filters: " + ", ".join(filter_parts) + "."

    parts = []
    parts.append(
        f'The search for "{query}" was interpreted as a {intent_labels.get(intent, "general")} query.'
    )

    if intent == "sequence":
        parts.append(
            "This looks like sequence-like input, so the current MVP is prioritizing binder-style matches as a placeholder for future sequence similarity search."
        )
    elif intent == "binder":
        parts.append(
            "Results are ordered to keep binders first, followed by linked targets, trials, and disease context."
        )
    elif intent == "disease":
        parts.append(
            "This disease-style query is shown as a filtered exploration view rather than redirecting to a single disease page."
        )
    elif intent == "protein":
        parts.append(
            "Exact target-style searches may redirect directly to a protein page when there are no extra filters."
        )

    parts.append(
        f'It returned {disease_count} disease record{"s" if disease_count != 1 else ""}, '
        f'{protein_count} protein target{"s" if protein_count != 1 else ""}, '
        f'{binder_count} therapeutic binder{"s" if binder_count != 1 else ""}, and '
        f'{trial_count} clinical trial{"s" if trial_count != 1 else ""}.'
        f'{filter_text}'
    )

    if binders:
        top_binders = ", ".join([b.get("binder_name") or "Unknown binder" for b in binders[:3]])
        parts.append(f"Top binder matches include {top_binders}.")
    if proteins:
        top_proteins = ", ".join([p.get("gene_symbol") or p.get("protein_name") or "Unknown protein" for p in proteins[:3]])
        parts.append(f"Protein-related matches include {top_proteins}.")
    if diseases:
        top_diseases = ", ".join([d.get("disease_name") or "Unknown disease" for d in diseases[:3]])
        parts.append(f"Disease-related matches include {top_diseases}.")
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
    route_hint = None
    effective_query = query
    filter_options = get_filter_options()

    if query:
        try:
            route_decision = decide_search_route(
                query,
                binder_type=binder_type,
                clinical_status=clinical_status,
                disease_name=disease_name
            )

            if route_decision["mode"] == "redirect":
                return redirect(url_for(route_decision["endpoint"], **route_decision["values"]))

            effective_query = route_decision["effective_query"]
            disease_name = route_decision["auto_disease_name"]
            route_hint = route_decision["route_hint"]

            results = search_database(
                effective_query,
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
        effective_query=effective_query,
        results=results,
        error=error,
        summary=summary,
        route_hint=route_hint,
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
        route_decision = decide_search_route(
            query,
            binder_type=binder_type,
            clinical_status=clinical_status,
            disease_name=disease_name
        )

        if route_decision["mode"] == "redirect":
            return jsonify({
                "redirect": True,
                "endpoint": route_decision["endpoint"],
                "values": route_decision["values"]
            })

        effective_query = route_decision["effective_query"]
        disease_name = route_decision["auto_disease_name"]

        results = search_database(
            effective_query,
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
            "redirect": False,
            "route_hint": route_decision["route_hint"],
            "effective_query": effective_query,
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
                    d.disease_category,
                    pd.tag_reason
                FROM protein_diseases pd
                JOIN diseases d ON pd.disease_id = d.disease_id
                WHERE pd.protein_id = %s
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
                ORDER BY
                    CASE
                        WHEN b.clinical_status ILIKE 'approved%%' THEN 1
                        WHEN b.clinical_status ILIKE 'phase 4%%' THEN 2
                        WHEN b.clinical_status ILIKE 'phase 3%%' THEN 3
                        WHEN b.clinical_status ILIKE 'phase 2%%' THEN 4
                        WHEN b.clinical_status ILIKE 'phase 1%%' THEN 5
                        ELSE 6
                    END,
                    b.binder_name;
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

    binder_type_breakdown = summarize_counter(
        Counter((b.get("binder_type") or "Unclassified") for b in binders)
    )
    binder_status_breakdown = summarize_counter(
        Counter((b.get("clinical_status") or b.get("approval_status") or "Unknown") for b in binders)
    )

    structure_candidates = build_structure_candidates(
        structures,
        fallback_uniprot=protein.get("uniprot_accession"),
        fallback_title=protein.get("protein_name")
    )
    default_structure = structure_candidates[0] if structure_candidates else None

    return render_template(
        "protein_detail.html",
        protein=protein,
        diseases=diseases,
        binders=binders,
        trials=trials,
        structures=structures,
        binder_type_breakdown=binder_type_breakdown,
        binder_status_breakdown=binder_status_breakdown,
        structure_candidates=structure_candidates,
        default_structure=default_structure
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
                    p.uniprot_accession,
                    g.gene_symbol,
                    g.gene_name,
                    p.functional_description,
                    pb.interaction_type,
                    COUNT(DISTINCT pt.trial_id) AS linked_trial_count
                FROM protein_binders pb
                JOIN proteins p ON pb.protein_id = p.protein_id
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                LEFT JOIN protein_trials pt ON pt.protein_id = p.protein_id
                WHERE pb.binder_id = %s
                GROUP BY
                    p.protein_id,
                    p.protein_name,
                    p.uniprot_accession,
                    g.gene_symbol,
                    g.gene_name,
                    p.functional_description,
                    pb.interaction_type
                ORDER BY p.protein_name;
            """, (binder_id,))
            proteins = cur.fetchall()

            cur.execute("""
                SELECT
                    ct.trial_id,
                    ct.nct_id,
                    ct.trial_title,
                    ct.phase,
                    ct.recruitment_status,
                    ct.condition_name,
                    ct.sponsor_name
                FROM binder_trials bt
                JOIN clinical_trials ct ON bt.trial_id = ct.trial_id
                WHERE bt.binder_id = %s
                ORDER BY
                    CASE
                        WHEN ct.phase ILIKE 'phase 4%%' THEN 1
                        WHEN ct.phase ILIKE 'phase 3%%' THEN 2
                        WHEN ct.phase ILIKE 'phase 2%%' THEN 3
                        WHEN ct.phase ILIKE 'phase 1%%' THEN 4
                        ELSE 5
                    END,
                    ct.trial_title;
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

            cur.execute("""
                SELECT DISTINCT
                    b2.binder_id,
                    b2.binder_name,
                    b2.binder_type,
                    b2.clinical_status,
                    bm2.modality_name
                FROM binders b2
                LEFT JOIN binder_modalities bm2 ON b2.modality_id = bm2.modality_id
                WHERE b2.binder_id <> %s
                  AND (
                    EXISTS (
                        SELECT 1
                        FROM protein_binders pb1
                        JOIN protein_binders pb2 ON pb1.protein_id = pb2.protein_id
                        WHERE pb1.binder_id = %s
                          AND pb2.binder_id = b2.binder_id
                    )
                    OR EXISTS (
                        SELECT 1
                        FROM binder_diseases bd1
                        JOIN binder_diseases bd2 ON bd1.disease_id = bd2.disease_id
                        WHERE bd1.binder_id = %s
                          AND bd2.binder_id = b2.binder_id
                    )
                  )
                ORDER BY b2.binder_name
                LIMIT 8;
            """, (binder_id, binder_id, binder_id))
            related_binders = cur.fetchall()

            cur.execute("""
                SELECT DISTINCT
                    p.protein_id,
                    p.protein_name,
                    p.uniprot_accession,
                    g.gene_symbol,
                    s.structure_id,
                    s.pdb_id,
                    s.structure_title,
                    s.experimental_method,
                    s.resolution
                FROM protein_binders pb
                JOIN proteins p ON pb.protein_id = p.protein_id
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                JOIN protein_structures ps ON p.protein_id = ps.protein_id
                JOIN structures s ON ps.structure_id = s.structure_id
                WHERE pb.binder_id = %s
                ORDER BY p.protein_name, s.pdb_id;
            """, (binder_id,))
            target_structures = cur.fetchall()

    target_count = len(proteins)
    trial_count = len(trials)
    disease_count = len(diseases)
    structure_count = len(structures)
    sequence_length = len((binder.get("sequence") or "").replace("\n", "").replace(" ", "")) if binder.get("sequence") else 0

    trial_phase_breakdown = summarize_counter(
        Counter((t.get("phase") or "Unknown") for t in trials)
    )
    disease_category_breakdown = summarize_counter(
        Counter((d.get("disease_category") or "Uncategorized") for d in diseases)
    )

    binder_story_parts = []
    binder_story_parts.append(
        f'{binder.get("binder_name", "This binder")} is currently classified as '
        f'{binder.get("binder_type") or "an unclassified binder"}.'
    )
    binder_story_parts.append(
        f'It is linked to {target_count} target{"s" if target_count != 1 else ""}, '
        f'{trial_count} clinical trial{"s" if trial_count != 1 else ""}, '
        f'{disease_count} disease tag{"s" if disease_count != 1 else ""}, and '
        f'{structure_count} directly linked structure{"s" if structure_count != 1 else ""}.'
    )
    if binder.get("clinical_status") or binder.get("approval_status"):
        binder_story_parts.append(
            f'Clinical status: {binder.get("clinical_status") or binder.get("approval_status")}.'
        )
    if proteins:
        protein_names = ", ".join([(p.get("gene_symbol") or p.get("protein_name") or "Unknown target") for p in proteins[:3]])
        binder_story_parts.append(f'Primary linked targets include {protein_names}.')
    binder_story = " ".join(binder_story_parts)

    target_structure_candidates = []
    for row in target_structures:
        pdb_id = safe_string(row.get("pdb_id"))
        if not pdb_id:
            continue

        gene_or_name = safe_string(row.get("gene_symbol")) or safe_string(row.get("protein_name")) or "Target"
        method = safe_string(row.get("experimental_method")) or "Experimental structure"
        resolution = row.get("resolution")
        resolution_text = f" · Resolution: {resolution}" if resolution is not None else ""

        target_structure_candidates.append({
            "source_type": "pdb",
            "source_id": pdb_id.upper(),
            "label": f"{gene_or_name} · PDB {pdb_id.upper()}",
            "title": safe_string(row.get("structure_title")) or f"{gene_or_name} structure",
            "subtitle": f"{method}{resolution_text}",
            "viewer_url": build_molstar_pdb_viewer_url(pdb_id),
            "external_url": f"https://www.rcsb.org/structure/{pdb_id.upper()}",
            "origin": f"Target-linked structure ({gene_or_name})"
        })

    fallback_uniprot = None
    fallback_title = None
    if proteins:
        first_target = proteins[0]
        fallback_uniprot = first_target.get("uniprot_accession")
        fallback_title = first_target.get("protein_name")

    direct_candidates = build_structure_candidates(structures)
    structure_candidates = dedupe_structure_candidates(direct_candidates + target_structure_candidates)

    if not structure_candidates and fallback_uniprot:
        structure_candidates = build_structure_candidates(
            [],
            fallback_uniprot=fallback_uniprot,
            fallback_title=fallback_title
        )

    default_structure = structure_candidates[0] if structure_candidates else None

    return render_template(
        "binder_detail.html",
        binder=binder,
        proteins=proteins,
        trials=trials,
        structures=structures,
        diseases=diseases,
        related_binders=related_binders,
        target_count=target_count,
        trial_count=trial_count,
        disease_count=disease_count,
        structure_count=structure_count,
        sequence_length=sequence_length,
        trial_phase_breakdown=trial_phase_breakdown,
        disease_category_breakdown=disease_category_breakdown,
        binder_story=binder_story,
        structure_candidates=structure_candidates,
        default_structure=default_structure
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
                    d.disease_category,
                    td.tag_reason
                FROM trial_diseases td
                JOIN diseases d ON td.disease_id = d.disease_id
                WHERE td.trial_id = %s
                ORDER BY d.disease_name;
            """, (trial_id,))
            diseases = cur.fetchall()

            cur.execute("""
                SELECT
                    b.binder_id,
                    b.binder_name,
                    b.binder_type,
                    b.clinical_status,
                    b.approval_status,
                    b.mechanism_of_action,
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
                    p.uniprot_accession,
                    g.gene_symbol,
                    g.gene_name
                FROM protein_trials pt
                JOIN proteins p ON pt.protein_id = p.protein_id
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                WHERE pt.trial_id = %s
                ORDER BY p.protein_name;
            """, (trial_id,))
            proteins = cur.fetchall()

            cur.execute("""
                SELECT DISTINCT
                    ct2.trial_id,
                    ct2.nct_id,
                    ct2.trial_title,
                    ct2.phase,
                    ct2.recruitment_status
                FROM clinical_trials ct2
                LEFT JOIN trial_diseases td2 ON ct2.trial_id = td2.trial_id
                LEFT JOIN binder_trials bt2 ON ct2.trial_id = bt2.trial_id
                LEFT JOIN protein_trials pt2 ON ct2.trial_id = pt2.trial_id
                WHERE ct2.trial_id <> %s
                  AND (
                    EXISTS (
                        SELECT 1
                        FROM trial_diseases td_self
                        WHERE td_self.trial_id = %s
                          AND td_self.disease_id = td2.disease_id
                    )
                    OR EXISTS (
                        SELECT 1
                        FROM binder_trials bt_self
                        WHERE bt_self.trial_id = %s
                          AND bt_self.binder_id = bt2.binder_id
                    )
                    OR EXISTS (
                        SELECT 1
                        FROM protein_trials pt_self
                        WHERE pt_self.trial_id = %s
                          AND pt_self.protein_id = pt2.protein_id
                    )
                  )
                ORDER BY ct2.trial_title
                LIMIT 8;
            """, (trial_id, trial_id, trial_id, trial_id))
            related_trials = cur.fetchall()

    binder_type_breakdown = summarize_counter(
        Counter((b.get("binder_type") or "Unclassified") for b in binders)
    )
    trial_context_counts = {
        "binders": len(binders),
        "proteins": len(proteins),
        "diseases": len(diseases),
        "related_trials": len(related_trials)
    }

    return render_template(
        "trial_detail.html",
        trial=trial,
        diseases=diseases,
        binders=binders,
        proteins=proteins,
        related_trials=related_trials,
        binder_type_breakdown=binder_type_breakdown,
        trial_context_counts=trial_context_counts
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
                    p.uniprot_accession,
                    g.gene_symbol,
                    g.gene_name
                FROM protein_diseases pd
                JOIN proteins p ON pd.protein_id = p.protein_id
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                WHERE pd.disease_id = %s
                ORDER BY p.protein_name;
            """, (disease_id,))
            proteins = cur.fetchall()

            cur.execute("""
                SELECT
                    ct.trial_id,
                    ct.nct_id,
                    ct.trial_title,
                    ct.phase,
                    ct.recruitment_status,
                    ct.condition_name
                FROM trial_diseases td
                JOIN clinical_trials ct ON td.trial_id = ct.trial_id
                WHERE td.disease_id = %s
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

            cur.execute("""
                SELECT DISTINCT
                    d2.disease_id,
                    d2.disease_name,
                    d2.disease_category
                FROM diseases d2
                LEFT JOIN protein_diseases pd2 ON d2.disease_id = pd2.disease_id
                LEFT JOIN binder_diseases bd2 ON d2.disease_id = bd2.disease_id
                LEFT JOIN trial_diseases td2 ON d2.disease_id = td2.disease_id
                WHERE d2.disease_id <> %s
                  AND (
                    EXISTS (
                        SELECT 1
                        FROM protein_diseases pd_self
                        WHERE pd_self.disease_id = %s
                          AND pd_self.protein_id = pd2.protein_id
                    )
                    OR EXISTS (
                        SELECT 1
                        FROM binder_diseases bd_self
                        WHERE bd_self.disease_id = %s
                          AND bd_self.binder_id = bd2.binder_id
                    )
                    OR EXISTS (
                        SELECT 1
                        FROM trial_diseases td_self
                        WHERE td_self.disease_id = %s
                          AND td_self.trial_id = td2.trial_id
                    )
                  )
                ORDER BY d2.disease_name
                LIMIT 8;
            """, (disease_id, disease_id, disease_id, disease_id))
            related_diseases = cur.fetchall()

    binder_type_breakdown = summarize_counter(
        Counter((b.get("binder_type") or "Unclassified") for b in binders)
    )
    protein_gene_breakdown = summarize_counter(
        Counter((p.get("gene_symbol") or p.get("protein_name") or "Unknown target") for p in proteins[:12])
    )

    return render_template(
        "disease_detail.html",
        disease=disease,
        proteins=proteins,
        trials=trials,
        binders=binders,
        related_diseases=related_diseases,
        binder_type_breakdown=binder_type_breakdown,
        protein_gene_breakdown=protein_gene_breakdown
    )


if __name__ == "__main__":
    app.run(debug=True)