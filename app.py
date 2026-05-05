import os
import re
import math
import shlex
from functools import lru_cache
from collections import Counter
from difflib import SequenceMatcher

try:
    from Bio import pairwise2
    from Bio.Align import substitution_matrices
    BIOPYTHON_AVAILABLE = True
except Exception:
    pairwise2 = None
    substitution_matrices = None
    BIOPYTHON_AVAILABLE = False
from flask import Flask, render_template, request, jsonify, redirect, url_for
from dotenv import load_dotenv
import requests
import psycopg2
import psycopg2.extras

PAGE_SIZE = 5

load_dotenv()

app = Flask(__name__)

DB_HOST = os.getenv("DB_HOST")
DB_PORT = os.getenv("DB_PORT", "5432")
DB_NAME = os.getenv("DB_NAME", "molecular_search_db")
DB_USER = os.getenv("DB_USER", "postgres")
DB_PASSWORD = os.getenv("DB_PASSWORD")


FULL_BINDER_TYPES = [
    "IgG",
    "Bispecific Antibody",
    "ADC",
    "VHH",
    "Fc Fusion",
    "Peptide",
]

# Full sponsor taxonomy. Non-protein therapeutics and unclear records are kept
# in Other so the UI/filtering matches the sponsor-requested binder classes.
EXTENDED_BINDER_TYPES = FULL_BINDER_TYPES + ["Other"]
PRIMARY_BINDER_TYPES = FULL_BINDER_TYPES

BINDER_TYPE_LABELS = {
    "IgG": "Monoclonal Antibody (IgG)",
    "Bispecific Antibody": "Bispecific Antibody",
    "ADC": "Antibody-Drug Conjugate (ADC)",
    "VHH": "Nanobody (VHH)",
    "Fc Fusion": "Fc Fusion",
    "Peptide": "Peptide",
    "Other": "Other / Unclassified",
}

BINDER_TYPE_DEFINITIONS = {
    "IgG": "Standard monoclonal antibody / IgG therapeutic.",
    "Bispecific Antibody": "Engineered antibody that can bind two targets or epitopes.",
    "ADC": "Antibody-drug conjugate with an antibody linked to a payload.",
    "VHH": "Nanobody or VHH single-domain antibody binder.",
    "Fc Fusion": "Protein or receptor domain fused to an Fc region.",
    "Peptide": "Short protein/peptide binder sequence.",
    "Other": "Unclassified, non-protein, or unsupported binder type.",
}


def canonicalize_binder_type(value: str) -> str:
    raw = safe_string(value).lower()
    if not raw:
        return "Other"

    normalized = re.sub(r"[_\-]+", " ", raw)
    normalized = re.sub(r"\s+", " ", normalized).strip()

    # Order matters. ADCs and bispecific antibodies are also antibodies,
    # so classify these before the general IgG / monoclonal antibody rule.
    if any(term in normalized for term in [
        "antibody drug conjugate", "antibody drug", "drug conjugate",
        "adc", "emtansine", "deruxtecan", "vedotin", "ozogamicin",
        "tesirine", "mafodotin", "duocarmazine", "maytansinoid"
    ]):
        return "ADC"

    if any(term in normalized for term in [
        "bispecific", "bi specific", "dual specific", "dual-specific", "multispecific",
        "trispecific", "bite", "t cell engager", "t-cell engager", "engager",
        "crossmab", "duobody", "tandab", "dart"
    ]):
        return "Bispecific Antibody"

    if any(term in normalized for term in [
        "nanobody", "vhh", "single domain antibody", "single-domain antibody",
        "camelid", "sdab", "single domain", "single-domain"
    ]):
        return "VHH"

    if any(term in normalized for term in [
        "fc fusion", "fc-fusion", "fusion protein", "receptor fc",
        "receptor-fc", "fc region", "trap", "immunoadhesin",
        "etanercept", "aflibercept", "abatacept", "belatacept", "romiplostim"
    ]):
        return "Fc Fusion"

    if any(term in normalized for term in [
        "peptide", "oligopeptide", "protein peptide", "peptidomimetic"
    ]):
        return "Peptide"

    if (
        normalized in {"igg", "igg antibody", "antibody", "monoclonal antibody", "mab"}
        or normalized.endswith("mab")
        or "monoclonal" in normalized
    ):
        return "IgG"

    # Anything outside the sponsor taxonomy stays in Other instead of becoming
    # a separate filter category.
    return "Other"


def binder_classification_family(value: str) -> str:
    binder_type = canonicalize_binder_type(value)
    if binder_type in FULL_BINDER_TYPES:
        return "Full sponsor taxonomy"
    return "Unresolved / other"


def decorate_binder_record(record: dict) -> dict:
    if not record:
        return record
    binder_type = canonicalize_binder_type(record.get("binder_type") or record.get("modality_name"))
    record["binder_type"] = binder_type
    record["binder_type_label"] = BINDER_TYPE_LABELS.get(binder_type, binder_type)
    record["binder_class_family"] = binder_classification_family(binder_type)
    record["is_primary_binder_class"] = binder_type in FULL_BINDER_TYPES
    record["disease_tags_label"] = "Disease tags"
    return record


def decorate_binder_records(records):
    return [decorate_binder_record(r) for r in (records or [])]


def summarize_binder_classes(records):
    counter = Counter(canonicalize_binder_type((r or {}).get("binder_type") or (r or {}).get("modality_name")) for r in (records or []))
    ordered = []
    for binder_type in EXTENDED_BINDER_TYPES:
        count = counter.get(binder_type, 0)
        if count:
            ordered.append({"label": binder_type, "count": count})
    return ordered


def build_binder_classification(binder: dict, proteins=None, trials=None, diseases=None, structures=None) -> dict:
    binder = binder or {}
    proteins = proteins or []
    trials = trials or []
    diseases = diseases or []
    structures = structures or []

    binder_type = canonicalize_binder_type(binder.get("binder_type") or binder.get("modality_name"))
    return {
        "binder_type": binder_type,
        "binder_type_label": BINDER_TYPE_LABELS.get(binder_type, binder_type),
        "binder_type_definition": BINDER_TYPE_DEFINITIONS.get(binder_type, ""),
        "class_family": binder.get("binder_class_family") or binder_classification_family(binder_type),
        "is_primary_class": binder_type in FULL_BINDER_TYPES,
        "linked_target_count": len(proteins),
        "linked_trial_count": len(trials),
        "linked_disease_count": len(diseases),
        "linked_structure_count": len(structures),
    }


def build_binder_visualization(binder: dict, proteins=None, diseases=None, trials=None, related_binders=None) -> dict:
    proteins = proteins or []
    diseases = diseases or []
    trials = trials or []
    related_binders = related_binders or []

    node_specs = []
    for row in proteins[:3]:
        label = safe_string(row.get("gene_symbol")) or safe_string(row.get("protein_name")) or "Target"
        node_specs.append({
            "category": "protein",
            "label": label[:16],
            "full_label": label,
            "subtitle": "Target",
            "href": url_for("protein_detail", protein_id=row.get("protein_id")),
        })
    for row in diseases[:2]:
        label = safe_string(row.get("disease_name")) or "Disease"
        node_specs.append({
            "category": "disease",
            "label": label[:16],
            "full_label": label,
            "subtitle": "Disease",
            "href": url_for("disease_detail", disease_id=row.get("disease_id")),
        })
    for row in trials[:2]:
        label = safe_string(row.get("nct_id")) or safe_string(row.get("trial_title")) or "Trial"
        node_specs.append({
            "category": "trial",
            "label": label[:16],
            "full_label": safe_string(row.get("trial_title")) or label,
            "subtitle": safe_string(row.get("phase")) or "Trial",
            "href": url_for("trial_detail", trial_id=row.get("trial_id")),
        })
    for row in related_binders[:1]:
        label = safe_string(row.get("binder_name")) or "Binder"
        node_specs.append({
            "category": "binder",
            "label": label[:16],
            "full_label": label,
            "subtitle": safe_string(row.get("binder_type")) or "Related",
            "href": url_for("binder_detail", binder_id=row.get("binder_id")),
        })

    positions = [
        (250, 62), (405, 110), (438, 240), (330, 314),
        (170, 314), (62, 240), (95, 110), (250, 338),
    ]

    nodes = []
    for index, spec in enumerate(node_specs[:len(positions)]):
        x, y = positions[index]
        nodes.append({
            **spec,
            "x": x,
            "y": y,
            "line_x1": 250,
            "line_y1": 180,
            "line_x2": x,
            "line_y2": y,
        })

    legend = [
        {"label": "Targets", "count": len(proteins)},
        {"label": "Disease Tags", "count": len(diseases)},
        {"label": "Trials", "count": len(trials)},
        {"label": "Related Binders", "count": len(related_binders)},
    ]

    return {
        "center_label": safe_string((binder or {}).get("binder_name"))[:22] or "Binder",
        "center_type": safe_string((binder or {}).get("binder_type")) or "Binder",
        "nodes": nodes,
        "legend": legend,
        "note": "The center node is the selected binder. Outer nodes summarize the strongest linked targets, diseases, trials, and neighboring binders.",
    }




def display_value(value, missing_text="Data not available from current source"):
    if value is None:
        return missing_text
    if isinstance(value, str) and not value.strip():
        return missing_text
    return value


def build_record_completeness(record: dict, fields: list[tuple[str, str]]) -> dict:
    total = len(fields)
    present_items = []
    missing_items = []

    for field_name, label in fields:
        value = (record or {}).get(field_name)
        is_present = value is not None and (not isinstance(value, str) or bool(value.strip()))
        if is_present:
            present_items.append(label)
        else:
            missing_items.append(label)

    score = round((len(present_items) / total) * 100) if total else 0
    return {
        "score": score,
        "present_count": len(present_items),
        "total_count": total,
        "present_items": present_items,
        "missing_items": missing_items,
    }


def binder_sequence_message(binder: dict) -> str:
    modality = safe_string((binder or {}).get("modality_name")).lower()
    binder_type = safe_string((binder or {}).get("binder_type")).lower()
    name = safe_string((binder or {}).get("binder_name"))

    if "small molecule" in modality or binder_type == "other":
        return f"Sequence is not applicable or not expected for {name or 'this record'} because it is currently stored as a small-molecule or non-protein therapeutic context item."

    return "Sequence data is not available from the current source. For protein therapeutics, this can be added from a curated FASTA record or a biologics-specific source."


def execute_search(query: str, binder_type=None, clinical_status=None, disease_name=None):
    route_decision = decide_search_route(
        query,
        binder_type=binder_type,
        clinical_status=clinical_status,
        disease_name=disease_name,
    )

    if route_decision["mode"] == "redirect":
        return {
            "route_decision": route_decision,
            "effective_query": route_decision.get("effective_query", query),
            "results": None,
            "summary": None,
        }

    effective_query = route_decision["effective_query"]
    resolved_disease_name = route_decision["auto_disease_name"]

    if route_decision.get("intent") == "sequence":
        results = fetch_sequence_similarity_results(
            effective_query,
            binder_type=binder_type,
            clinical_status=clinical_status,
            disease_name=resolved_disease_name,
        )
    else:
        results = search_database(
            effective_query,
            binder_type=binder_type,
            clinical_status=clinical_status,
            disease_name=resolved_disease_name,
            intent_override=route_decision.get("intent"),
        )

    summary = build_free_summary(
        query,
        results,
        binder_type=binder_type,
        clinical_status=clinical_status,
        disease_name=resolved_disease_name,
    )

    return {
        "route_decision": route_decision,
        "effective_query": effective_query,
        "resolved_disease_name": resolved_disease_name,
        "results": results,
        "summary": summary,
    }



def build_free_summary(query, results, binder_type=None, clinical_status=None, disease_name=None):
    """Create a short user-facing summary for the search page and API response."""
    results = results or {}
    binder_count = len(results.get("binders") or [])
    protein_count = len(results.get("proteins") or [])
    trial_count = len(results.get("trials") or [])
    disease_count = len(results.get("diseases") or [])

    active_filters = []
    if binder_type:
        active_filters.append(f"binder type: {binder_type}")
    if clinical_status:
        active_filters.append(f"clinical status: {clinical_status}")
    if disease_name:
        active_filters.append(f"disease tag: {disease_name}")

    filter_text = ""
    if active_filters:
        filter_text = " Filters applied: " + "; ".join(active_filters) + "."

    return (
        f'Search results for "{query}": showing {binder_count} binder(s), '
        f'{protein_count} protein target(s), {trial_count} clinical trial(s), '
        f'and {disease_count} disease tag(s).' + filter_text
    )

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
                ORDER BY binder_type;
            """)
            binder_types = [canonicalize_binder_type(row["binder_type"]) for row in cur.fetchall()]
            binder_types = [item for item in EXTENDED_BINDER_TYPES if item in set(binder_types)]
            # Always expose the full sponsor taxonomy in the filter, even if one class has zero current records.
            for item in EXTENDED_BINDER_TYPES:
                if item not in binder_types:
                    binder_types.append(item)

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


def get_filter_options_safe():
    try:
        return get_filter_options()
    except Exception:
        return {
            "binder_types": [],
            "clinical_statuses": [],
            "diseases": []
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


AMINO_ACID_CHARS = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
ALLOWED_STRUCTURE_EXTENSIONS = {".pdb", ".cif", ".mmcif"}

def normalize_biological_sequence(sequence: str) -> str:
    if not sequence:
        return ""
    cleaned_lines = []
    for raw_line in sequence.splitlines():
        line = raw_line.strip()
        if not line or line.startswith(">"): continue
        cleaned_lines.append(line)
    compact = "".join(cleaned_lines) if cleaned_lines else sequence
    compact = compact.upper().replace("*", "").replace("-", "")
    compact = re.sub(r"[^A-Z]", "", compact)
    return "".join(ch for ch in compact if ch in AMINO_ACID_CHARS)

def sequence_preview(sequence: str, preview_len: int = 24) -> str:
    seq = normalize_biological_sequence(sequence)
    return seq if len(seq) <= preview_len else f"{seq[:preview_len]}…"

def build_kmer_set(sequence: str, k: int = 3):
    seq = normalize_biological_sequence(sequence)
    if not seq: return set()
    if len(seq) < k: return {seq}
    return {seq[i:i+k] for i in range(len(seq)-k+1)}

def calculate_sequence_similarity_details(query_sequence: str, candidate_sequence: str) -> dict:
    """
    Compare a user FASTA/amino-acid query against a stored binder sequence.

    Primary method: Biopython local alignment with BLOSUM62. This is better
    than plain text matching because it handles substitutions, gaps, and
    partial FASTA fragments. If Biopython is not installed, the function falls
    back to the older k-mer/text similarity method so the app still runs.
    """
    query = normalize_biological_sequence(query_sequence)
    candidate = normalize_biological_sequence(candidate_sequence)

    empty_result = {
        "similarity_score": 0.0,
        "percent_identity": 0.0,
        "query_coverage": 0.0,
        "alignment_score": 0.0,
        "method": "none",
    }

    if not query or not candidate:
        return empty_result

    if query == candidate:
        return {
            "similarity_score": 1.0,
            "percent_identity": 100.0,
            "query_coverage": 100.0,
            "alignment_score": float(len(query)),
            "method": "exact",
        }

    if BIOPYTHON_AVAILABLE:
        try:
            matrix = substitution_matrices.load("BLOSUM62")
            alignments = pairwise2.align.localds(
                query,
                candidate,
                matrix,
                -10,
                -0.5,
                one_alignment_only=True,
            )

            if not alignments:
                return empty_result

            aligned_query, aligned_candidate, raw_score, _start, _end = alignments[0]

            aligned_pairs = [
                (q, c)
                for q, c in zip(aligned_query, aligned_candidate)
                if q != "-" and c != "-"
            ]
            aligned_query_residues = sum(1 for q in aligned_query if q != "-")

            if not aligned_pairs:
                return empty_result

            matches = sum(1 for q, c in aligned_pairs if q == c)
            percent_identity = matches / len(aligned_pairs)
            query_coverage = aligned_query_residues / max(len(query), 1)
            length_ratio = min(len(query), len(candidate)) / max(len(query), len(candidate), 1)

            # Composite score used for ranking. Identity is most important,
            # coverage prevents tiny local matches from ranking too high, and
            # length ratio slightly favors comparable sequence lengths.
            similarity_score = (0.65 * percent_identity) + (0.25 * query_coverage) + (0.10 * length_ratio)

            return {
                "similarity_score": round(min(similarity_score, 1.0), 6),
                "percent_identity": round(percent_identity * 100, 2),
                "query_coverage": round(query_coverage * 100, 2),
                "alignment_score": round(float(raw_score), 2),
                "method": "Biopython local alignment (BLOSUM62)",
            }
        except Exception:
            # Fall through to safe fallback below.
            pass

    query_kmers = build_kmer_set(query, 3)
    candidate_kmers = build_kmer_set(candidate, 3)
    union = query_kmers | candidate_kmers
    jaccard = (len(query_kmers & candidate_kmers) / len(union)) if union else 0.0
    sequence_ratio = SequenceMatcher(None, query, candidate).ratio()
    length_ratio = min(len(query), len(candidate)) / max(len(query), len(candidate), 1)
    containment_bonus = 0.08 if query in candidate or candidate in query else 0.0
    score = min((0.50 * sequence_ratio) + (0.35 * jaccard) + (0.15 * length_ratio) + containment_bonus, 1.0)

    return {
        "similarity_score": round(score, 6),
        "percent_identity": round(sequence_ratio * 100, 2),
        "query_coverage": round(length_ratio * 100, 2),
        "alignment_score": round(score * 100, 2),
        "method": "Fallback k-mer/text similarity",
    }

def calculate_sequence_similarity(query_sequence: str, candidate_sequence: str) -> float:
    return calculate_sequence_similarity_details(query_sequence, candidate_sequence)["similarity_score"]

def allowed_structure_filename(filename: str) -> bool:
    filename = (filename or "").lower()
    return any(filename.endswith(ext) for ext in ALLOWED_STRUCTURE_EXTENSIONS)

def safe_float(value):
    try:
        if value in (None, "", ".", "?"): return None
        return float(value)
    except (TypeError, ValueError):
        return None

@lru_cache(maxsize=128)
def fetch_structure_text(url: str, timeout: int = 8) -> str:
    response = requests.get(url, timeout=timeout)
    response.raise_for_status()
    return response.text

def parse_pdb_coordinates(text: str):
    atoms=[]; residues=set(); chains=set()
    for line in text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")): continue
        try:
            x=float(line[30:38].strip()); y=float(line[38:46].strip()); z=float(line[46:54].strip())
        except ValueError:
            continue
        chain=line[21:22].strip() or "?"; resseq=line[22:26].strip() or "?"; resname=line[17:20].strip() or "UNK"
        atoms.append((x,y,z)); residues.add((chain,resseq,resname)); chains.add(chain)
    return atoms,residues,chains

def tokenize_cif_line(line: str):
    try: return shlex.split(line, posix=False)
    except Exception: return line.split()

def parse_mmcif_coordinates(text: str):
    atoms=[]; residues=set(); chains=set(); lines=text.splitlines(); i=0
    while i < len(lines):
        if lines[i].strip() != 'loop_': i += 1; continue
        i += 1; headers=[]
        while i < len(lines) and lines[i].strip().startswith('_'):
            headers.append(lines[i].strip()); i += 1
        if not headers or not any(h.startswith('_atom_site.') for h in headers): continue
        hi={h: idx for idx,h in enumerate(headers)}
        x_idx=hi.get('_atom_site.Cartn_x'); y_idx=hi.get('_atom_site.Cartn_y'); z_idx=hi.get('_atom_site.Cartn_z')
        chain_idx=hi.get('_atom_site.label_asym_id') or hi.get('_atom_site.auth_asym_id')
        comp_idx=hi.get('_atom_site.label_comp_id'); seq_idx=hi.get('_atom_site.label_seq_id') or hi.get('_atom_site.auth_seq_id')
        if x_idx is None or y_idx is None or z_idx is None: continue
        while i < len(lines):
            row=lines[i].strip()
            if not row or row=='#' or row=='loop_' or row.startswith('data_') or row.startswith('_'): break
            tokens=tokenize_cif_line(row)
            if len(tokens) >= len(headers):
                x=safe_float(tokens[x_idx]); y=safe_float(tokens[y_idx]); z=safe_float(tokens[z_idx])
                if x is not None and y is not None and z is not None:
                    chain=tokens[chain_idx] if chain_idx is not None and chain_idx < len(tokens) else '?'
                    comp=tokens[comp_idx] if comp_idx is not None and comp_idx < len(tokens) else 'UNK'
                    seq=tokens[seq_idx] if seq_idx is not None and seq_idx < len(tokens) else '?'
                    atoms.append((x,y,z)); residues.add((chain,seq,comp)); chains.add(chain)
            i += 1
        continue
    return atoms,residues,chains

def compute_structure_signature_from_text(text: str, filename: str = 'uploaded_structure') -> dict:
    filename=(filename or 'uploaded_structure').lower()
    if filename.endswith('.pdb'): atoms,residues,chains=parse_pdb_coordinates(text)
    else:
        atoms,residues,chains=parse_mmcif_coordinates(text)
        if not atoms: atoms,residues,chains=parse_pdb_coordinates(text)
    if not atoms: raise ValueError('No atomic coordinates could be extracted from the uploaded structure. Please upload a valid PDB or mmCIF file.')
    xs=[a[0] for a in atoms]; ys=[a[1] for a in atoms]; zs=[a[2] for a in atoms]
    centroid=(sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))
    rg=math.sqrt(sum((x-centroid[0])**2 + (y-centroid[1])**2 + (z-centroid[2])**2 for x,y,z in atoms)/len(atoms))
    span_x=max(xs)-min(xs); span_y=max(ys)-min(ys); span_z=max(zs)-min(zs); max_span=max(span_x,span_y,span_z,1.0)
    return {'atom_count':len(atoms),'residue_count':len(residues),'chain_count':len(chains) or 1,'radius_gyration':round(rg,4),'span_x':round(span_x,4),'span_y':round(span_y,4),'span_z':round(span_z,4),'volume_hint':round(span_x*span_y*span_z,4),'compactness':round(rg/max_span,6)}

def compare_structure_signatures(query_sig: dict, candidate_sig: dict) -> dict:
    rel=lambda a,b: abs(a-b)/max(abs(a),abs(b),1.0)
    atom_diff=rel(query_sig['atom_count'],candidate_sig['atom_count']); residue_diff=rel(query_sig['residue_count'],candidate_sig['residue_count'])
    chain_diff=abs(query_sig['chain_count']-candidate_sig['chain_count']); rg_diff=rel(query_sig['radius_gyration'],candidate_sig['radius_gyration'])
    span_diff=(rel(query_sig['span_x'],candidate_sig['span_x'])+rel(query_sig['span_y'],candidate_sig['span_y'])+rel(query_sig['span_z'],candidate_sig['span_z']))/3.0
    compactness_diff=rel(query_sig['compactness'],candidate_sig['compactness'])
    distance=(0.26*atom_diff + 0.20*residue_diff + 0.10*min(chain_diff,5)/5.0 + 0.20*rg_diff + 0.18*span_diff + 0.06*compactness_diff)
    similarity=max(0.0, 1.0-min(distance,1.0))
    return {'distance':round(distance,4),'similarity_score':round(similarity*100.0,2),'score_breakdown':{'atom_diff':round(atom_diff,4),'residue_diff':round(residue_diff,4),'chain_diff':int(chain_diff),'rg_diff':round(rg_diff,4),'span_diff':round(span_diff,4),'compactness_diff':round(compactness_diff,4)}}

def load_structure_similarity_candidates(limit: int = 24):
    sql="""
        SELECT s.structure_id,s.pdb_id,s.structure_title,s.experimental_method,s.resolution,s.structure_file_url,
               p.protein_id,p.protein_name,p.uniprot_accession,g.gene_symbol
        FROM structures s
        LEFT JOIN protein_structures ps ON s.structure_id = ps.structure_id
        LEFT JOIN proteins p ON ps.protein_id = p.protein_id
        LEFT JOIN genes g ON p.gene_id = g.gene_id
        WHERE s.pdb_id IS NOT NULL
        ORDER BY s.deposition_date DESC NULLS LAST, s.pdb_id
        LIMIT %s;
    """
    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(sql,(limit,)); return cur.fetchall()

def run_structure_similarity_search(uploaded_filename: str, uploaded_bytes: bytes, top_k: int = 10):
    if not uploaded_bytes:
        raise ValueError('No structure file was uploaded.')

    decoded = uploaded_bytes.decode('utf-8', errors='ignore')
    query_signature = compute_structure_signature_from_text(decoded, uploaded_filename)
    candidates = load_structure_similarity_candidates(limit=max(top_k * 2, 24))
    matches = []
    failures = []

    for row in candidates:
        pdb_id = safe_string(row.get('pdb_id'))
        structure_url = safe_string(row.get('structure_file_url'))
        if not structure_url and pdb_id:
            structure_url = f'https://files.rcsb.org/download/{pdb_id.upper()}.cif'
        if not structure_url:
            continue

        try:
            candidate_text = fetch_structure_text(structure_url)
            candidate_signature = compute_structure_signature_from_text(candidate_text, structure_url)
            comparison = compare_structure_signatures(query_signature, candidate_signature)
            matches.append({
                'structure_id': row.get('structure_id'),
                'pdb_id': pdb_id.upper() if pdb_id else None,
                'structure_title': row.get('structure_title'),
                'experimental_method': row.get('experimental_method'),
                'resolution': row.get('resolution'),
                'protein_id': row.get('protein_id'),
                'protein_name': row.get('protein_name'),
                'gene_symbol': row.get('gene_symbol'),
                'uniprot_accession': row.get('uniprot_accession'),
                'structure_file_url': structure_url,
                'viewer_url': build_molstar_pdb_viewer_url(pdb_id) if pdb_id else '',
                'external_url': f'https://www.rcsb.org/structure/{pdb_id.upper()}' if pdb_id else structure_url,
                **comparison,
            })
        except Exception as exc:
            failures.append(f"{pdb_id or structure_url}: {exc}")

    matches.sort(key=lambda item: item['similarity_score'], reverse=True)
    top_matches = matches[:top_k]
    summary = (
        f'The uploaded structure "{uploaded_filename}" was compared against {len(matches)} linked structure candidates using an MVP structural signature based on atom count, residue count, chain count, overall span, and radius of gyration. The top match scored {top_matches[0]["similarity_score"]:.2f}% similarity.'
        if top_matches else
        'The structure similarity workflow ran, but no comparable structures were available in the current dataset.'
    )
    return {
        'query_filename': uploaded_filename,
        'query_signature': query_signature,
        'matches': top_matches,
        'summary': summary,
        'failure_count': len(failures),
        'failures': failures[:5],
    }

def fetch_sequence_similarity_results(query_sequence: str, binder_type=None, clinical_status=None, disease_name=None, limit: int = 25):
    normalized_query=normalize_biological_sequence(query_sequence)
    results={'diseases':[],'proteins':[],'binders':[],'trials':[],'intent':'sequence','primary_section':'binders','section_order':['binders','proteins','trials','diseases'],'sequence_query':normalized_query,'sequence_query_length':len(normalized_query),'sequence_strategy':'Biopython local alignment with BLOSUM62 across stored binder FASTA sequences' if BIOPYTHON_AVAILABLE else 'Fallback k-mer/text similarity across stored binder FASTA sequences'}
    if len(normalized_query) < 12: return results
    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            binder_sql="""
                SELECT DISTINCT b.binder_id,b.binder_name,b.binder_type,b.sequence,b.clinical_status,b.mechanism_of_action,b.approval_status,bm.modality_name
                FROM binders b
                LEFT JOIN binder_modalities bm ON b.modality_id = bm.modality_id
                LEFT JOIN binder_diseases bd ON b.binder_id = bd.binder_id
                LEFT JOIN diseases d ON bd.disease_id = d.disease_id
                WHERE b.sequence IS NOT NULL AND TRIM(b.sequence) <> ''
            """
            params=[]
            if binder_type: binder_sql += ' AND b.binder_type = %s'; params.append(binder_type)
            if clinical_status: binder_sql += ' AND b.clinical_status = %s'; params.append(clinical_status)
            if disease_name: binder_sql += ' AND d.disease_name = %s'; params.append(disease_name)
            binder_sql += ' ORDER BY b.binder_name;'; cur.execute(binder_sql, tuple(params)); binder_rows=cur.fetchall()
            scored=[]
            for row in binder_rows:
                candidate_sequence=normalize_biological_sequence(row.get('sequence') or '')
                if len(candidate_sequence) < 12: continue
                similarity_details = calculate_sequence_similarity_details(normalized_query, candidate_sequence)
                similarity = similarity_details['similarity_score']
                if similarity < 0.20:
                    continue
                row['sequence_length'] = len(candidate_sequence)
                row['sequence_similarity_score'] = round(similarity, 4)
                row['sequence_similarity_percent'] = round(similarity * 100, 1)
                row['percent_identity'] = similarity_details['percent_identity']
                row['query_coverage'] = similarity_details['query_coverage']
                row['alignment_score'] = similarity_details['alignment_score']
                row['sequence_similarity_method'] = similarity_details['method']
                row['query_sequence_preview'] = sequence_preview(normalized_query)
                row['matched_sequence_preview'] = sequence_preview(candidate_sequence)
                scored.append(row)
            scored.sort(key=lambda item:(-item['sequence_similarity_score'], item.get('binder_name') or ''))
            top_binders=scored[:limit]; results['binders']=top_binders
            binder_ids=[row['binder_id'] for row in top_binders if row.get('binder_id')]
            if not binder_ids: return results
            sim_map={row['binder_id']: row['sequence_similarity_score'] for row in top_binders if row.get('binder_id')}
            cur.execute("""SELECT DISTINCT p.protein_id,p.protein_name,p.uniprot_accession,p.organism_name,p.functional_description,g.gene_symbol,g.gene_name,pb.binder_id FROM protein_binders pb JOIN proteins p ON pb.protein_id = p.protein_id LEFT JOIN genes g ON p.gene_id = g.gene_id WHERE pb.binder_id = ANY(%s) ORDER BY p.protein_name;""", (binder_ids,)); protein_rows=cur.fetchall(); pmap={}
            for row in protein_rows:
                pid=row['protein_id']; best=sim_map.get(row['binder_id'],0.0); existing=pmap.get(pid)
                if not existing or best > existing['best_sequence_similarity']:
                    cleaned=dict(row); cleaned.pop('binder_id',None); cleaned['best_sequence_similarity']=round(best,4); cleaned['best_sequence_similarity_percent']=round(best*100,1); pmap[pid]=cleaned
            results['proteins']=sorted(pmap.values(), key=lambda item:(-item.get('best_sequence_similarity',0), item.get('protein_name') or ''))[:20]
            cur.execute("""SELECT DISTINCT ct.trial_id,ct.nct_id,ct.trial_title,ct.condition_name,ct.phase,ct.recruitment_status,bt.binder_id FROM binder_trials bt JOIN clinical_trials ct ON bt.trial_id = ct.trial_id WHERE bt.binder_id = ANY(%s) ORDER BY ct.trial_title;""", (binder_ids,)); trial_rows=cur.fetchall(); tmap={}
            for row in trial_rows:
                tid=row['trial_id']; best=sim_map.get(row['binder_id'],0.0); existing=tmap.get(tid)
                if not existing or best > existing['best_sequence_similarity']:
                    cleaned=dict(row); cleaned.pop('binder_id',None); cleaned['best_sequence_similarity']=round(best,4); cleaned['best_sequence_similarity_percent']=round(best*100,1); tmap[tid]=cleaned
            results['trials']=sorted(tmap.values(), key=lambda item:(-item.get('best_sequence_similarity',0), item.get('trial_title') or ''))[:20]
            cur.execute("""SELECT DISTINCT d.disease_id,d.disease_name,d.disease_category,d.description,bd.binder_id FROM binder_diseases bd JOIN diseases d ON bd.disease_id = d.disease_id WHERE bd.binder_id = ANY(%s) ORDER BY d.disease_name;""", (binder_ids,)); disease_rows=cur.fetchall(); dmap={}
            for row in disease_rows:
                did=row['disease_id']; best=sim_map.get(row['binder_id'],0.0); existing=dmap.get(did)
                if not existing or best > existing['best_sequence_similarity']:
                    cleaned=dict(row); cleaned.pop('binder_id',None); cleaned['best_sequence_similarity']=round(best,4); cleaned['best_sequence_similarity_percent']=round(best*100,1); dmap[did]=cleaned
            results['diseases']=sorted(dmap.values(), key=lambda item:(-item.get('best_sequence_similarity',0), item.get('disease_name') or ''))[:20]
    results['binders']=decorate_binder_records(results.get('binders',[])); results=enrich_search_results(results); results['binder_class_breakdown']=summarize_binder_classes(results.get('binders',[])); return results

def apply_result_sorting(results: dict, binder_sort='relevance', protein_sort='relevance', trial_sort='relevance', disease_sort='relevance'):
    s=lambda v: safe_string(v).lower()
    binders=results.get('binders') or []; proteins=results.get('proteins') or []; trials=results.get('trials') or []; diseases=results.get('diseases') or []
    if binder_sort == 'name_asc': binders.sort(key=lambda x: s(x.get('binder_name')))
    elif binder_sort == 'name_desc': binders.sort(key=lambda x: s(x.get('binder_name')), reverse=True)
    elif binder_sort == 'binder_type': binders.sort(key=lambda x: (s(x.get('binder_type')), s(x.get('binder_name'))))
    elif binder_sort == 'clinical_status': binders.sort(key=lambda x: (s(x.get('clinical_status') or x.get('approval_status')), s(x.get('binder_name'))))
    if protein_sort == 'gene_symbol': proteins.sort(key=lambda x: (s(x.get('gene_symbol')), s(x.get('protein_name'))))
    elif protein_sort == 'protein_name': proteins.sort(key=lambda x: s(x.get('protein_name')))
    elif protein_sort == 'uniprot_accession': proteins.sort(key=lambda x: s(x.get('uniprot_accession')))
    if trial_sort == 'nct_id': trials.sort(key=lambda x: s(x.get('nct_id')))
    elif trial_sort == 'phase': trials.sort(key=lambda x: (s(x.get('phase')), s(x.get('trial_title'))))
    elif trial_sort == 'recruitment_status': trials.sort(key=lambda x: (s(x.get('recruitment_status')), s(x.get('trial_title'))))
    elif trial_sort == 'title': trials.sort(key=lambda x: s(x.get('trial_title')))
    if disease_sort == 'disease_name': diseases.sort(key=lambda x: s(x.get('disease_name')))
    elif disease_sort == 'disease_category': diseases.sort(key=lambda x: (s(x.get('disease_category')), s(x.get('disease_name'))))
    results.update({'binders':binders,'proteins':proteins,'trials':trials,'diseases':diseases,'sorts':{'binder_sort':binder_sort,'protein_sort':protein_sort,'trial_sort':trial_sort,'disease_sort':disease_sort}})
    return results



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


BINDING_REGION_COLORS = [
    "#dc2626", "#2563eb", "#16a34a", "#d97706", "#7c3aed",
    "#0891b2", "#db2777", "#4f46e5", "#65a30d", "#ea580c"
]


def table_exists(cur, table_name: str) -> bool:
    cur.execute("""
        SELECT EXISTS (
            SELECT 1
            FROM information_schema.tables
            WHERE table_schema = 'public' AND table_name = %s
        ) AS exists;
    """, (table_name,))
    row = cur.fetchone()
    if isinstance(row, dict):
        return bool(row.get('exists'))
    return bool(row[0]) if row else False


def binding_region_color(index: int) -> str:
    return BINDING_REGION_COLORS[index % len(BINDING_REGION_COLORS)]


def _safe_int(value, default=None):
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _build_annotation_entry(row: dict, protein_length: int, fallback_index: int = 0, inferred: bool = False) -> dict:
    protein_length = _safe_int(protein_length, 0) or 0
    if protein_length <= 0:
        protein_length = 1000

    start = _safe_int(row.get('region_start'))
    end = _safe_int(row.get('region_end'))
    if start is None or end is None or end < start:
        # fallback conceptual band if a record exists but lacks coordinates
        width = max(20, protein_length // 10)
        start = min(protein_length, 1 + fallback_index * max(25, protein_length // 12))
        end = min(protein_length, start + width)
        inferred = True

    start = max(1, min(start, protein_length))
    end = max(start, min(end, protein_length))
    span = max(1, end - start + 1)
    left_percent = ((start - 1) / protein_length) * 100.0
    width_percent = max(1.5, (span / protein_length) * 100.0)
    center_percent = min(100.0, left_percent + (width_percent / 2.0))
    color = row.get('color_hex') or binding_region_color(fallback_index)

    evidence = safe_string(row.get('evidence_type')) or ('Inferred layout' if inferred else 'Annotation')
    label = safe_string(row.get('region_label')) or ('Conceptual interaction zone' if inferred else 'Binding region')

    return {
        'binder_id': row.get('binder_id'),
        'binder_name': row.get('binder_name') or 'Unknown binder',
        'protein_id': row.get('protein_id'),
        'protein_name': row.get('protein_name') or 'Unknown protein',
        'gene_symbol': row.get('gene_symbol'),
        'region_start': start,
        'region_end': end,
        'region_span': span,
        'region_label': label,
        'evidence_type': evidence,
        'source_note': row.get('source_note') or '',
        'color_hex': color,
        'left_percent': round(left_percent, 3),
        'width_percent': round(width_percent, 3),
        'center_percent': round(center_percent, 3),
        'inferred': inferred,
    }


def infer_binding_annotations_for_protein(protein: dict, binders: list) -> list:
    binders = binders or []
    if not binders:
        return []
    protein_length = _safe_int(protein.get('sequence_length'), 0) or max(300, 80 * len(binders))
    spacing = max(24, protein_length // max(len(binders) + 1, 2))
    zone_width = max(18, min(70, protein_length // max(len(binders), 4)))
    annotations = []
    for idx, binder in enumerate(binders):
        start = 1 + idx * spacing
        end = min(protein_length, start + zone_width)
        annotations.append(_build_annotation_entry({
            'binder_id': binder.get('binder_id'),
            'binder_name': binder.get('binder_name'),
            'protein_id': protein.get('protein_id'),
            'protein_name': protein.get('protein_name'),
            'gene_symbol': protein.get('gene_symbol'),
            'region_start': start,
            'region_end': end,
            'region_label': f"Conceptual site for {binder.get('binder_name') or 'binder'}",
            'evidence_type': 'Inferred layout',
            'source_note': 'Generated from binder-target links because no residue-level annotation is stored yet.',
            'color_hex': binding_region_color(idx),
        }, protein_length, fallback_index=idx, inferred=True))
    return annotations


def load_binding_annotations_for_protein(cur, protein: dict, binders: list) -> tuple[list, str]:
    """Load binding annotations only for the binders currently visible on the page.

    This is an important performance fix for target pages like EGFR. The old
    version loaded every binding-site row for the protein, even when the page
    only displayed 5 binders. For common targets, that hidden visualization
    query can make the page feel slow.
    """
    binders = binders or []
    visible_binder_ids = [b.get('binder_id') for b in binders if b.get('binder_id') is not None]

    if not visible_binder_ids:
        return [], 'none'

    if not table_exists(cur, 'binder_binding_sites'):
        inferred = infer_binding_annotations_for_protein(protein, binders)
        return inferred, ('inferred' if inferred else 'none')

    cur.execute("""
        SELECT
            bbs.binding_site_id,
            bbs.binder_id,
            bbs.protein_id,
            bbs.region_start,
            bbs.region_end,
            bbs.region_label,
            bbs.evidence_type,
            bbs.source_note,
            bbs.color_hex,
            b.binder_name,
            p.protein_name,
            g.gene_symbol
        FROM binder_binding_sites bbs
        JOIN binders b ON bbs.binder_id = b.binder_id
        JOIN proteins p ON bbs.protein_id = p.protein_id
        LEFT JOIN genes g ON p.gene_id = g.gene_id
        WHERE bbs.protein_id = %s
          AND bbs.binder_id = ANY(%s)
        ORDER BY COALESCE(bbs.region_start, 999999), b.binder_name
        LIMIT %s;
    """, (protein.get('protein_id'), visible_binder_ids, max(len(visible_binder_ids) * 3, PAGE_SIZE + 1)))
    rows = cur.fetchall() or []
    if rows:
        protein_length = _safe_int(protein.get('sequence_length'), 0) or 1000
        annotations = [
            _build_annotation_entry(row, protein_length, fallback_index=i, inferred=False)
            for i, row in enumerate(rows)
        ]
        return annotations, 'annotated'

    inferred = infer_binding_annotations_for_protein(protein, binders)
    return inferred, ('inferred' if inferred else 'none')


def build_binder_binding_maps(cur, binder: dict, proteins: list) -> tuple[list, str]:
    proteins = proteins or []
    if not proteins:
        return [], 'none'

    annotations_by_protein = {}
    mode = 'none'
    if table_exists(cur, 'binder_binding_sites'):
        cur.execute("""
            SELECT
                bbs.binding_site_id,
                bbs.binder_id,
                bbs.protein_id,
                bbs.region_start,
                bbs.region_end,
                bbs.region_label,
                bbs.evidence_type,
                bbs.source_note,
                bbs.color_hex,
                b.binder_name,
                p.protein_name,
                g.gene_symbol
            FROM binder_binding_sites bbs
            JOIN binders b ON bbs.binder_id = b.binder_id
            JOIN proteins p ON bbs.protein_id = p.protein_id
            LEFT JOIN genes g ON p.gene_id = g.gene_id
            WHERE bbs.binder_id = %s
            ORDER BY p.protein_name, COALESCE(bbs.region_start, 999999);
        """, (binder.get('binder_id'),))
        rows = cur.fetchall() or []
        if rows:
            mode = 'annotated'
            for i, row in enumerate(rows):
                annotations_by_protein.setdefault(row.get('protein_id'), []).append(row)

    maps = []
    for idx, protein in enumerate(proteins):
        protein_length = _safe_int(protein.get('sequence_length'), 0) or 1000
        if mode == 'annotated' and annotations_by_protein.get(protein.get('protein_id')):
            annotations = [
                _build_annotation_entry(row, protein_length, fallback_index=j, inferred=False)
                for j, row in enumerate(annotations_by_protein[protein.get('protein_id')])
            ]
        else:
            mode = 'inferred'
            annotations = [
                _build_annotation_entry({
                    'binder_id': binder.get('binder_id'),
                    'binder_name': binder.get('binder_name'),
                    'protein_id': protein.get('protein_id'),
                    'protein_name': protein.get('protein_name'),
                    'gene_symbol': protein.get('gene_symbol'),
                    'region_start': 1 + idx * 75,
                    'region_end': min(protein_length, (1 + idx * 75) + max(24, protein_length // 8)),
                    'region_label': f"Conceptual site on {protein.get('gene_symbol') or protein.get('protein_name')}",
                    'evidence_type': 'Inferred layout',
                    'source_note': 'Generated from the binder-target relationship because no explicit binding-site coordinates are stored yet.',
                    'color_hex': binding_region_color(idx),
                }, protein_length, fallback_index=idx, inferred=True)
            ]
        maps.append({
            'protein_id': protein.get('protein_id'),
            'protein_label': protein.get('gene_symbol') or protein.get('protein_name') or 'Target',
            'protein_name': protein.get('protein_name') or 'Unknown protein',
            'gene_symbol': protein.get('gene_symbol') or '',
            'uniprot_accession': protein.get('uniprot_accession') or '',
            'sequence_length': protein_length,
            'interaction_type': protein.get('interaction_type') or 'Interaction not available',
            'annotations': annotations,
        })
    return maps, mode


def build_binding_region_summary(annotations: list, mode: str) -> str:
    if not annotations:
        return 'No binding-region map is available yet for this page.'
    if mode == 'annotated':
        return 'Highlighted bands represent stored binding-site annotations tied to the current target-binder relationships.'
    return 'Highlighted bands are conceptual interaction zones inferred from current target-binder links. They improve the visual story without claiming residue-level experimental certainty.'


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
            "route_hint": "Sequence-like input detected. Running binder-sequence similarity search across stored FASTA sequences.",
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


def search_database(
    search_term: str,
    binder_type=None,
    clinical_status=None,
    disease_name=None,
    binder_limit=5,
    section_limits=None,
    intent_override=None,
):
    """Fast home-page search.

    This version fixes the real performance issue from searches such as
    "EGFR binders": the route stripped the word "binders" and passed only
    "EGFR" into this function, so the function re-detected the intent as a
    target/general search and still ran every large search query.

    It also avoids broad DISTINCT/ILIKE joins when the query resolves to an
    exact disease or target. In those cases it starts from the small link table
    using IDs, which is much faster.
    """
    search_term = (search_term or "").strip()
    like_term = f"%{search_term}%"
    intent = intent_override or detect_query_intent(search_term)
    section_limits = section_limits or {}

    def safe_limit(section_name, default=PAGE_SIZE):
        try:
            return max(int(section_limits.get(section_name, default) or default), PAGE_SIZE)
        except (TypeError, ValueError):
            return default

    try:
        binder_limit = max(int(binder_limit or PAGE_SIZE), PAGE_SIZE)
    except (TypeError, ValueError):
        binder_limit = PAGE_SIZE

    limits = {
        "binders": binder_limit,
        "proteins": safe_limit("proteins"),
        "trials": safe_limit("trials"),
        "diseases": safe_limit("diseases"),
    }

    if intent == "binder":
        sections_to_run = {"binders"}
    else:
        sections_to_run = {"binders", "proteins", "trials", "diseases"}

    results = {
        "diseases": [],
        "proteins": [],
        "binders": [],
        "trials": [],
        "intent": intent,
        "primary_section": "binders",
        "section_order": ["binders", "proteins", "trials", "diseases"],
    }

    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            exact_protein = get_exact_protein_match(search_term)
            exact_disease = None
            if disease_name:
                exact_disease = get_exact_disease_match(disease_name)
            if not exact_disease and intent == "disease":
                exact_disease = get_exact_disease_match(search_term)

            # -----------------------------
            # Binder results
            # -----------------------------
            if "binders" in sections_to_run:
                if exact_protein:
                    binder_sql = """
                        SELECT
                            b.binder_id,
                            b.binder_name,
                            b.binder_type,
                            b.sequence,
                            b.clinical_status,
                            b.mechanism_of_action,
                            b.approval_status,
                            bm.modality_name,
                            1 AS rank_score
                        FROM protein_binders pb
                        JOIN binders b ON pb.binder_id = b.binder_id
                        LEFT JOIN binder_modalities bm ON b.modality_id = bm.modality_id
                        WHERE pb.protein_id = %s
                    """
                    binder_params = [exact_protein["protein_id"]]
                    if binder_type:
                        binder_sql += " AND b.binder_type = %s"
                        binder_params.append(binder_type)
                    if clinical_status:
                        binder_sql += " AND b.clinical_status = %s"
                        binder_params.append(clinical_status)
                    if exact_disease:
                        binder_sql += """
                            AND EXISTS (
                                SELECT 1
                                FROM binder_diseases bd
                                WHERE bd.binder_id = b.binder_id
                                  AND bd.disease_id = %s
                            )
                        """
                        binder_params.append(exact_disease["disease_id"])
                    binder_sql += """
                        ORDER BY
                            CASE
                                WHEN b.clinical_status ILIKE 'approved%%' THEN 1
                                WHEN b.clinical_status ILIKE 'phase 4%%' THEN 2
                                WHEN b.clinical_status ILIKE 'phase 3%%' THEN 3
                                WHEN b.clinical_status ILIKE 'phase 2%%' THEN 4
                                WHEN b.clinical_status ILIKE 'phase 1%%' THEN 5
                                ELSE 6
                            END,
                            b.binder_name
                        LIMIT %s;
                    """
                    binder_params.append(limits["binders"] + 1)

                elif exact_disease:
                    binder_sql = """
                        SELECT
                            b.binder_id,
                            b.binder_name,
                            b.binder_type,
                            b.sequence,
                            b.clinical_status,
                            b.mechanism_of_action,
                            b.approval_status,
                            bm.modality_name,
                            1 AS rank_score
                        FROM binder_diseases bd
                        JOIN binders b ON bd.binder_id = b.binder_id
                        LEFT JOIN binder_modalities bm ON b.modality_id = bm.modality_id
                        WHERE bd.disease_id = %s
                    """
                    binder_params = [exact_disease["disease_id"]]
                    if binder_type:
                        binder_sql += " AND b.binder_type = %s"
                        binder_params.append(binder_type)
                    if clinical_status:
                        binder_sql += " AND b.clinical_status = %s"
                        binder_params.append(clinical_status)
                    binder_sql += """
                        ORDER BY
                            CASE
                                WHEN b.clinical_status ILIKE 'approved%%' THEN 1
                                WHEN b.clinical_status ILIKE 'phase 4%%' THEN 2
                                WHEN b.clinical_status ILIKE 'phase 3%%' THEN 3
                                WHEN b.clinical_status ILIKE 'phase 2%%' THEN 4
                                WHEN b.clinical_status ILIKE 'phase 1%%' THEN 5
                                ELSE 6
                            END,
                            b.binder_name
                        LIMIT %s;
                    """
                    binder_params.append(limits["binders"] + 1)

                else:
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
                        like_term, like_term, like_term, like_term, like_term, like_term, like_term, like_term,
                    ]
                    if binder_type:
                        binder_sql += " AND b.binder_type = %s"
                        binder_params.append(binder_type)
                    if clinical_status:
                        binder_sql += " AND b.clinical_status = %s"
                        binder_params.append(clinical_status)
                    binder_sql += " ORDER BY rank_score, b.binder_name LIMIT %s;"
                    binder_params.append(limits["binders"] + 1)

                cur.execute(binder_sql, tuple(binder_params))
                binder_rows = cur.fetchall() or []
                has_more_binders = len(binder_rows) > limits["binders"]
                results["binders"] = binder_rows[:limits["binders"]]
                results["binders_pagination"] = {
                    "showing": len(results["binders"]),
                    "has_next": has_more_binders,
                    "next_limit": limits["binders"] + PAGE_SIZE,
                }

            # -----------------------------
            # Protein results
            # -----------------------------
            if "proteins" in sections_to_run:
                if exact_disease:
                    protein_sql = """
                        SELECT DISTINCT
                            p.protein_id,
                            p.protein_name,
                            p.uniprot_accession,
                            p.organism_name,
                            p.functional_description,
                            g.gene_symbol,
                            g.gene_name,
                            1 AS rank_score
                        FROM protein_diseases pd
                        JOIN proteins p ON pd.protein_id = p.protein_id
                        LEFT JOIN genes g ON p.gene_id = g.gene_id
                        WHERE pd.disease_id = %s
                        ORDER BY p.protein_name
                        LIMIT %s;
                    """
                    protein_params = [exact_disease["disease_id"], limits["proteins"] + 1]
                else:
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
                        WHERE (
                            p.protein_name ILIKE %s
                            OR p.uniprot_accession ILIKE %s
                            OR p.functional_description ILIKE %s
                            OR g.gene_symbol ILIKE %s
                            OR g.gene_name ILIKE %s
                        )
                        ORDER BY rank_score, p.protein_name
                        LIMIT %s;
                    """
                    protein_params = [
                        search_term, search_term,
                        like_term, like_term, like_term, like_term, like_term,
                        limits["proteins"] + 1,
                    ]
                cur.execute(protein_sql, tuple(protein_params))
                results["proteins"] = cur.fetchall() or []

            # -----------------------------
            # Trial results
            # -----------------------------
            if "trials" in sections_to_run:
                if exact_disease:
                    trial_sql = """
                        SELECT DISTINCT
                            ct.trial_id,
                            ct.nct_id,
                            ct.trial_title,
                            ct.condition_name,
                            ct.phase,
                            ct.recruitment_status,
                            1 AS rank_score
                        FROM trial_diseases td
                        JOIN clinical_trials ct ON td.trial_id = ct.trial_id
                        WHERE td.disease_id = %s
                        ORDER BY ct.trial_title
                        LIMIT %s;
                    """
                    trial_params = [exact_disease["disease_id"], limits["trials"] + 1]
                else:
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
                        WHERE (
                            ct.trial_title ILIKE %s
                            OR ct.condition_name ILIKE %s
                            OR ct.nct_id ILIKE %s
                            OR ct.brief_summary ILIKE %s
                        )
                        ORDER BY rank_score, ct.trial_title
                        LIMIT %s;
                    """
                    trial_params = [search_term, search_term, like_term, like_term, like_term, like_term, limits["trials"] + 1]
                cur.execute(trial_sql, tuple(trial_params))
                results["trials"] = cur.fetchall() or []

            # -----------------------------
            # Disease results
            # -----------------------------
            if "diseases" in sections_to_run:
                if exact_disease:
                    disease_sql = """
                        SELECT
                            d.disease_id,
                            d.disease_name,
                            d.disease_category,
                            d.description,
                            1 AS rank_score
                        FROM diseases d
                        WHERE d.disease_id = %s
                        LIMIT %s;
                    """
                    disease_params = [exact_disease["disease_id"], limits["diseases"] + 1]
                else:
                    disease_sql = """
                        SELECT
                            d.disease_id,
                            d.disease_name,
                            d.disease_category,
                            d.description,
                            CASE WHEN LOWER(COALESCE(d.disease_name, '')) = LOWER(%s) THEN 1 ELSE 2 END AS rank_score
                        FROM diseases d
                        WHERE d.disease_name ILIKE %s
                           OR d.description ILIKE %s
                        ORDER BY rank_score, d.disease_name
                        LIMIT %s;
                    """
                    disease_params = [search_term, like_term, like_term, limits["diseases"] + 1]
                cur.execute(disease_sql, tuple(disease_params))
                results["diseases"] = cur.fetchall() or []

    # Trim non-binder sections to their requested limits so the existing template
    # pagination still works normally.
    for section in ["proteins", "trials", "diseases"]:
        limit = limits[section]
        if len(results.get(section, [])) > limit:
            results[section] = results[section][:limit]

    results["binders"] = decorate_binder_records(results.get("binders", []))
    results = reorder_results_by_intent(results, intent)
    results = enrich_search_results(results)
    results["binder_class_breakdown"] = summarize_binder_classes(results.get("binders", []))
    return results

def get_page_param(param_name: str, default: int = 1) -> int:
    try:
        value = int(request.args.get(param_name, default))
        return max(value, 1)
    except (TypeError, ValueError):
        return default

def paginate_records(records, page: int, per_page: int = PAGE_SIZE):
    records = list(records or [])
    page = max(int(page or 1), 1)
    start = (page - 1) * per_page
    end = start + per_page
    return records[start:end], {
        "page": page,
        "per_page": per_page,
        "total": len(records),
        "showing": min(end, len(records)),
        "has_prev": page > 1,
        "has_next": end < len(records),
        "prev_page": max(page - 1, 1),
        "next_page": page + 1,
    }



def build_db_pagination(page: int, total: int, per_page: int = PAGE_SIZE):
    page = max(int(page or 1), 1)
    total = int(total or 0)
    showing = min(page * per_page, total)
    return {
        "page": page,
        "per_page": per_page,
        "total": total,
        "showing": showing,
        "has_prev": page > 1,
        "has_next": showing < total,
        "prev_page": max(page - 1, 1),
        "next_page": page + 1,
    }

def fetch_db_page(cur, data_sql: str, params=(), page: int = 1, count_sql: str | None = None, count_params=None, per_page: int = PAGE_SIZE):
    """Fetch one page without running expensive COUNT(*) queries.

    The previous version executed a COUNT query for each section before loading
    the 5 visible rows. On pages with large joins, those COUNT(*) queries can be
    slower than the actual page query. This version fetches one extra row
    instead. If the extra row exists, the template shows a Show More link.
    """
    page = max(int(page or 1), 1)
    per_page = max(int(per_page or PAGE_SIZE), 1)
    offset = (page - 1) * per_page
    data_sql = data_sql.strip().rstrip(';')

    cur.execute(data_sql + "\nLIMIT %s OFFSET %s;", tuple(params) + (per_page + 1, offset))
    fetched_rows = cur.fetchall() or []
    has_next = len(fetched_rows) > per_page
    rows = fetched_rows[:per_page]
    showing = offset + len(rows)

    pagination = {
        "page": page,
        "per_page": per_page,
        "total": f"{showing}+" if has_next else showing,
        "showing": showing,
        "has_prev": page > 1,
        "has_next": has_next,
        "prev_page": max(page - 1, 1),
        "next_page": page + 1,
    }
    return rows, pagination


def build_page_url(param_name: str, page_value: int):
    args = request.args.to_dict(flat=True)
    if page_value <= 1:
        args.pop(param_name, None)
    else:
        args[param_name] = page_value
    view_args = dict(request.view_args or {})
    return url_for(request.endpoint, **view_args, **args)

@app.context_processor
def inject_pagination_helpers():
    return {"pagination_url": build_page_url, "page_size": PAGE_SIZE}

def fetch_binder_browser(query=None, binder_type=None, clinical_status=None, disease_name=None, sort='clinical_status', page=1, per_page=PAGE_SIZE):
    """Return a paginated binder-first browse table with sponsor filters."""
    query = (query or '').strip()
    like_query = f"%{query}%"
    page = max(int(page or 1), 1)
    per_page = max(int(per_page or PAGE_SIZE), 1)
    offset = (page - 1) * per_page
    related_join_needed = bool(query or disease_name)

    sql = """
        WITH base_binders AS (
            SELECT DISTINCT
                b.binder_id,
                b.binder_name,
                b.binder_type,
                b.sequence,
                b.clinical_status,
                b.approval_status,
                b.mechanism_of_action,
                b.binder_description,
                bm.modality_name
            FROM binders b
            LEFT JOIN binder_modalities bm ON b.modality_id = bm.modality_id
    """
    if related_join_needed:
        sql += """
            LEFT JOIN protein_binders pb ON b.binder_id = pb.binder_id
            LEFT JOIN proteins p ON pb.protein_id = p.protein_id
            LEFT JOIN genes g ON p.gene_id = g.gene_id
            LEFT JOIN binder_diseases bd ON b.binder_id = bd.binder_id
            LEFT JOIN diseases d ON bd.disease_id = d.disease_id
        """
    sql += " WHERE 1=1 "

    params = []
    if query:
        sql += """
            AND (
                b.binder_name ILIKE %s
                OR b.binder_type ILIKE %s
                OR b.clinical_status ILIKE %s
                OR b.mechanism_of_action ILIKE %s
                OR b.binder_description ILIKE %s
                OR COALESCE(g.gene_symbol, '') ILIKE %s
                OR COALESCE(p.protein_name, '') ILIKE %s
                OR COALESCE(d.disease_name, '') ILIKE %s
            )
        """
        params.extend([like_query] * 8)
    if binder_type:
        sql += " AND b.binder_type = %s"
        params.append(binder_type)
    if clinical_status:
        sql += " AND b.clinical_status = %s"
        params.append(clinical_status)
    if disease_name:
        sql += " AND d.disease_name = %s"
        params.append(disease_name)

    if sort == 'name_desc':
        sql += " ORDER BY b.binder_name DESC"
    elif sort == 'binder_type':
        sql += " ORDER BY b.binder_type ASC NULLS LAST, b.binder_name ASC"
    else:
        sql += " ORDER BY b.binder_name ASC"

    sql += """
            LIMIT %s OFFSET %s
        )
        SELECT
            b.binder_id, b.binder_name, b.binder_type, b.sequence,
            b.clinical_status, b.approval_status, b.mechanism_of_action,
            b.binder_description, b.modality_name,
            (SELECT COUNT(DISTINCT pb.protein_id) FROM protein_binders pb WHERE pb.binder_id = b.binder_id) AS target_count,
            (SELECT COUNT(DISTINCT bd.disease_id) FROM binder_diseases bd WHERE bd.binder_id = b.binder_id) AS disease_count,
            (SELECT COUNT(DISTINCT bt.trial_id) FROM binder_trials bt WHERE bt.binder_id = b.binder_id) AS trial_count,
            (SELECT COUNT(DISTINCT bs.structure_id) FROM binder_structures bs WHERE bs.binder_id = b.binder_id) AS direct_structure_count,
            (
                SELECT STRING_AGG(DISTINCT COALESCE(g.gene_symbol, p.protein_name), ', ')
                FROM protein_binders pb
                JOIN proteins p ON pb.protein_id = p.protein_id
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                WHERE pb.binder_id = b.binder_id
            ) AS target_preview,
            (
                SELECT STRING_AGG(DISTINCT d.disease_name, ', ')
                FROM binder_diseases bd
                JOIN diseases d ON bd.disease_id = d.disease_id
                WHERE bd.binder_id = b.binder_id
            ) AS disease_preview
        FROM base_binders b
    """
    params.extend([per_page + 1, offset])

    if sort == 'name_desc':
        sql += " ORDER BY b.binder_name DESC"
    elif sort == 'binder_type':
        sql += " ORDER BY b.binder_type ASC NULLS LAST, b.binder_name ASC"
    elif sort == 'target_count':
        sql += " ORDER BY target_count DESC, b.binder_name ASC"
    elif sort == 'trial_count':
        sql += " ORDER BY trial_count DESC, b.binder_name ASC"
    else:
        sql += """
            ORDER BY
                CASE
                    WHEN b.clinical_status ILIKE 'approved%%' THEN 1
                    WHEN b.clinical_status ILIKE 'max phase 4%%' THEN 2
                    WHEN b.clinical_status ILIKE 'phase 4%%' THEN 2
                    WHEN b.clinical_status ILIKE 'phase 3%%' THEN 3
                    WHEN b.clinical_status ILIKE 'max phase 3%%' THEN 3
                    WHEN b.clinical_status ILIKE 'phase 2%%' THEN 4
                    WHEN b.clinical_status ILIKE 'max phase 2%%' THEN 4
                    WHEN b.clinical_status ILIKE 'phase 1%%' THEN 5
                    WHEN b.clinical_status ILIKE 'max phase 1%%' THEN 5
                    ELSE 6
                END, b.binder_name ASC
        """
    sql += ";"

    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute(sql, tuple(params))
            rows = cur.fetchall()

    has_next = len(rows) > per_page
    rows = rows[:per_page]
    rows = decorate_binder_records(rows)
    for row in rows:
        row['has_sequence'] = bool(safe_string(row.get('sequence')))
        row['has_direct_structure'] = (row.get('direct_structure_count') or 0) > 0
        row['target_preview_list'] = [x.strip() for x in safe_string(row.get('target_preview')).split(',') if x.strip()][:4]
        row['disease_preview_list'] = [x.strip() for x in safe_string(row.get('disease_preview')).split(',') if x.strip()][:4]

    pagination = {
        "page": page, "per_page": per_page, "has_prev": page > 1,
        "has_next": has_next, "prev_page": max(page - 1, 1),
        "next_page": page + 1, "showing": offset + len(rows), "total": None,
    }
    return rows, pagination


@app.route("/", methods=["GET", "POST"])
def index():
    query = request.values.get("q", "").strip()
    binder_type = request.values.get("binder_type", "").strip() or None
    clinical_status = request.values.get("clinical_status", "").strip() or None
    disease_name = request.values.get("disease_name", "").strip() or None
    binder_sort = request.values.get("binder_sort", "relevance").strip() or "relevance"
    protein_sort = request.values.get("protein_sort", "relevance").strip() or "relevance"
    trial_sort = request.values.get("trial_sort", "relevance").strip() or "relevance"
    disease_sort = request.values.get("disease_sort", "relevance").strip() or "relevance"
    try:
        binders_limit = max(int(request.values.get("binders_limit", 5)), 5)
    except (TypeError, ValueError):
        binders_limit = 5
    results = None; error = None; summary = None; route_hint = None; effective_query = query
    filter_options = get_filter_options(); structure_results = None; structure_error = None
    if request.method == "POST":
        uploaded_file = request.files.get("structure_file")
        if not uploaded_file or not uploaded_file.filename:
            structure_error = "Please choose a PDB or mmCIF file before running structure similarity search."
        elif not allowed_structure_filename(uploaded_file.filename):
            structure_error = "Unsupported structure file type. Please upload a .pdb, .cif, or .mmcif file."
        else:
            try: structure_results = run_structure_similarity_search(uploaded_file.filename, uploaded_file.read())
            except Exception as e: structure_error = str(e)
    elif query:
        try:
            route_decision = decide_search_route(query, binder_type=binder_type, clinical_status=clinical_status, disease_name=disease_name)
            if route_decision["mode"] == "redirect": return redirect(url_for(route_decision["endpoint"], **route_decision["values"]))
            effective_query = route_decision["effective_query"]; disease_name = route_decision["auto_disease_name"]; route_hint = route_decision["route_hint"]
            if route_decision.get("intent") == "sequence":
                results = fetch_sequence_similarity_results(effective_query, binder_type=binder_type, clinical_status=clinical_status, disease_name=disease_name)
            else:
                # Only ask the database for the records needed for the current page.
                # This prevents broad searches like "lung cancer" from loading every
                # related binder/trial/protein before the template renders.
                section_limits = {
                    "proteins": get_page_param("proteins_page") * PAGE_SIZE,
                    "trials": get_page_param("trials_page") * PAGE_SIZE,
                    "diseases": get_page_param("diseases_page") * PAGE_SIZE,
                }
                results = search_database(
                    effective_query,
                    binder_type=binder_type,
                    clinical_status=clinical_status,
                    disease_name=disease_name,
                    binder_limit=binders_limit,
                    section_limits=section_limits,
                    intent_override=route_decision.get("intent"),
                )
            results = apply_result_sorting(results, binder_sort, protein_sort, trial_sort, disease_sort)
            for section_name, page_param in [
                ("proteins", "proteins_page"),
                ("trials", "trials_page"),
                ("diseases", "diseases_page"),
            ]:
                section_page = get_page_param(page_param)
                paged_items, pagination = paginate_records(results.get(section_name, []), section_page)
                results[section_name] = paged_items
                results[f"{section_name}_pagination"] = pagination
                results[f"{section_name}_page_param"] = page_param
            results["binders_page_param"] = "binders_limit"
            summary = build_free_summary(query, results, binder_type=binder_type, clinical_status=clinical_status, disease_name=disease_name)
        except Exception as e:
            error = str(e)
    return render_template("index.html", query=query, effective_query=effective_query, results=results, error=error, summary=summary, route_hint=route_hint, filter_options=filter_options, selected_binder_type=binder_type, selected_clinical_status=clinical_status, selected_disease_name=disease_name, selected_binder_sort=binder_sort, selected_protein_sort=protein_sort, selected_trial_sort=trial_sort, selected_disease_sort=disease_sort, structure_results=structure_results, structure_error=structure_error)

@app.route("/api/search", methods=["GET"])

def api_search():
    query = request.args.get("q", "").strip()
    binder_type = request.args.get("binder_type", "").strip() or None
    clinical_status = request.args.get("clinical_status", "").strip() or None
    disease_name = request.args.get("disease_name", "").strip() or None

    if not query:
        return jsonify({"error": "Missing search query"}), 400

    try:
        search_payload = execute_search(
            query,
            binder_type=binder_type,
            clinical_status=clinical_status,
            disease_name=disease_name,
        )
        route_decision = search_payload["route_decision"]

        if route_decision["mode"] == "redirect":
            return jsonify({
                "redirect": True,
                "endpoint": route_decision["endpoint"],
                "values": route_decision["values"]
            })

        return jsonify({
            "redirect": False,
            "route_hint": route_decision["route_hint"],
            "effective_query": search_payload["effective_query"],
            "summary": search_payload["summary"],
            "results": search_payload["results"]
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500



@app.route("/binders", methods=["GET"])
def binders_browser():
    query = request.args.get("q", "").strip()
    binder_type = request.args.get("binder_type", "").strip() or None
    clinical_status = request.args.get("clinical_status", "").strip() or None
    disease_name = request.args.get("disease_name", "").strip() or None
    sort = request.args.get("sort", "clinical_status").strip() or "clinical_status"
    page = get_page_param("page")

    filter_options = get_filter_options_safe()
    binders = []
    binder_pagination = {"page": page, "per_page": PAGE_SIZE, "has_prev": False, "has_next": False, "prev_page": 1, "next_page": 2, "showing": 0, "total": None}
    error = None
    try:
        binders, binder_pagination = fetch_binder_browser(
            query=query,
            binder_type=binder_type,
            clinical_status=clinical_status,
            disease_name=disease_name,
            sort=sort,
            page=page,
            per_page=PAGE_SIZE,
        )
    except Exception as exc:
        error = str(exc)

    binder_type_breakdown = summarize_binder_classes(binders)
    status_breakdown = summarize_counter(Counter((b.get("clinical_status") or b.get("approval_status") or "Unknown") for b in binders))

    return render_template(
        "binders.html",
        binders=binders,
        binders_pagination=binder_pagination,
        query=query,
        error=error,
        filter_options=filter_options,
        selected_binder_type=binder_type,
        selected_clinical_status=clinical_status,
        selected_disease_name=disease_name,
        selected_sort=sort,
        binder_type_breakdown=binder_type_breakdown,
        status_breakdown=status_breakdown,
    )


@app.route("/protein/<int:protein_id>")
def protein_detail(protein_id):
    binders_page = get_page_param("binders_page")
    trials_page = get_page_param("trials_page")
    structures_page = get_page_param("structures_page")
    diseases_page = get_page_param("diseases_page")

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

            diseases, diseases_pagination = fetch_db_page(
                cur,
                """
                SELECT d.disease_id,d.disease_name,d.disease_category,pd.tag_reason
                FROM protein_diseases pd
                JOIN diseases d ON pd.disease_id = d.disease_id
                WHERE pd.protein_id = %s
                ORDER BY d.disease_name
                """,
                (protein_id,),
                diseases_page,
                """
                SELECT COUNT(*) AS total
                FROM protein_diseases pd
                WHERE pd.protein_id = %s
                """,
            )

            binders, binders_pagination = fetch_db_page(
                cur,
                """
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
                    b.binder_name
                """,
                (protein_id,),
                binders_page,
                """
                SELECT COUNT(*) AS total
                FROM protein_binders pb
                WHERE pb.protein_id = %s
                """,
            )

            trials, trials_pagination = fetch_db_page(
                cur,
                """
                SELECT ct.trial_id,ct.nct_id,ct.trial_title,ct.phase,ct.recruitment_status,ct.condition_name
                FROM protein_trials pt
                JOIN clinical_trials ct ON pt.trial_id = ct.trial_id
                WHERE pt.protein_id = %s
                ORDER BY ct.trial_title
                """,
                (protein_id,),
                trials_page,
                """
                SELECT COUNT(*) AS total
                FROM protein_trials pt
                WHERE pt.protein_id = %s
                """,
            )

            structures, structures_pagination = fetch_db_page(
                cur,
                """
                SELECT s.structure_id,s.pdb_id,s.structure_title,s.experimental_method,s.resolution
                FROM protein_structures ps
                JOIN structures s ON ps.structure_id = s.structure_id
                WHERE ps.protein_id = %s
                ORDER BY s.pdb_id
                """,
                (protein_id,),
                structures_page,
                """
                SELECT COUNT(*) AS total
                FROM protein_structures ps
                WHERE ps.protein_id = %s
                """,
            )

            # Important performance fix:
            # Only build binding-site annotations for the binders visible on the current page.
            # Previously this used every linked binder, which made EGFR and disease-related pages slow.
            binding_site_annotations, binding_site_mode = load_binding_annotations_for_protein(cur, protein, binders)

            # Lightweight aggregate summaries for all linked binders without fetching every binder row.
            cur.execute("""
                SELECT
                    COALESCE(b.binder_type, 'Unknown') AS binder_type,
                    COUNT(*) AS count
                FROM protein_binders pb
                JOIN binders b ON pb.binder_id = b.binder_id
                WHERE pb.protein_id = %s
                GROUP BY COALESCE(b.binder_type, 'Unknown')
                ORDER BY count DESC, binder_type;
            """, (protein_id,))
            binder_type_breakdown = [
                {"label": row["binder_type"], "count": row["count"]}
                for row in cur.fetchall()
            ]

            cur.execute("""
                SELECT
                    COALESCE(NULLIF(b.clinical_status, ''), NULLIF(b.approval_status, ''), 'Unknown') AS status,
                    COUNT(*) AS count
                FROM protein_binders pb
                JOIN binders b ON pb.binder_id = b.binder_id
                WHERE pb.protein_id = %s
                GROUP BY COALESCE(NULLIF(b.clinical_status, ''), NULLIF(b.approval_status, ''), 'Unknown')
                ORDER BY count DESC, status;
            """, (protein_id,))
            binder_status_breakdown = [
                {"label": row["status"], "count": row["count"]}
                for row in cur.fetchall()
            ]

    binders = decorate_binder_records(binders)
    structure_candidates = build_structure_candidates(
        structures,
        fallback_uniprot=protein.get("uniprot_accession"),
        fallback_title=protein.get("protein_name"),
    )
    default_structure = structure_candidates[0] if structure_candidates else None
    binding_region_summary = build_binding_region_summary(binding_site_annotations, binding_site_mode)
    protein_completeness = build_record_completeness(
        protein,
        [
            ("protein_name", "Protein name"),
            ("uniprot_accession", "UniProt accession"),
            ("organism_name", "Organism"),
            ("sequence_length", "Sequence length"),
            ("subcellular_location", "Subcellular location"),
            ("functional_description", "Functional description"),
        ],
    )

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
        default_structure=default_structure,
        binding_site_annotations=binding_site_annotations,
        binding_site_mode=binding_site_mode,
        binding_region_summary=binding_region_summary,
        protein_completeness=protein_completeness,
        display_value=display_value,
        binders_pagination=binders_pagination,
        trials_pagination=trials_pagination,
        structures_pagination=structures_pagination,
        diseases_pagination=diseases_pagination,
    )



@app.route("/binder/<int:binder_id>")
def binder_detail(binder_id):
    with get_connection() as conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""
                SELECT
                    b.binder_id,b.binder_name,b.binder_type,b.sequence,b.clinical_status,b.binder_description,
                    b.mechanism_of_action,b.approval_status,b.developer_company,bm.modality_name
                FROM binders b
                LEFT JOIN binder_modalities bm ON b.modality_id = bm.modality_id
                WHERE b.binder_id = %s;
            """, (binder_id,))
            binder = cur.fetchone()
            if not binder:
                return "Binder not found", 404
            binder = decorate_binder_record(binder)

            cur.execute("""
                SELECT
                    p.protein_id,p.protein_name,p.uniprot_accession,p.sequence_length,
                    g.gene_symbol,g.gene_name,p.functional_description,pb.interaction_type,
                    COUNT(DISTINCT pt.trial_id) AS linked_trial_count
                FROM protein_binders pb
                JOIN proteins p ON pb.protein_id = p.protein_id
                LEFT JOIN genes g ON p.gene_id = g.gene_id
                LEFT JOIN protein_trials pt ON pt.protein_id = p.protein_id
                WHERE pb.binder_id = %s
                GROUP BY p.protein_id,p.protein_name,p.uniprot_accession,p.sequence_length,g.gene_symbol,g.gene_name,p.functional_description,pb.interaction_type
                ORDER BY p.protein_name;
            """, (binder_id,))
            proteins = cur.fetchall()

            cur.execute("""
                SELECT ct.trial_id,ct.nct_id,ct.trial_title,ct.phase,ct.recruitment_status,ct.condition_name,ct.sponsor_name
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
                SELECT s.structure_id,s.pdb_id,s.structure_title,s.experimental_method,s.resolution
                FROM binder_structures bs
                JOIN structures s ON bs.structure_id = s.structure_id
                WHERE bs.binder_id = %s
                ORDER BY md5(%s || '-' || COALESCE(s.pdb_id, s.structure_id::text));
            """, (binder_id, str(binder_id)))
            structures = cur.fetchall()

            cur.execute("""
                SELECT d.disease_id,d.disease_name,d.disease_category,bd.tag_reason
                FROM binder_diseases bd
                JOIN diseases d ON bd.disease_id = d.disease_id
                WHERE bd.binder_id = %s
                ORDER BY d.disease_name;
            """, (binder_id,))
            diseases = cur.fetchall()

            cur.execute("""
                SELECT DISTINCT b2.binder_id,b2.binder_name,b2.binder_type,b2.clinical_status,bm2.modality_name
                FROM binders b2
                LEFT JOIN binder_modalities bm2 ON b2.modality_id = bm2.modality_id
                WHERE b2.binder_id <> %s
                  AND (
                    EXISTS (
                        SELECT 1
                        FROM protein_binders pb1
                        JOIN protein_binders pb2 ON pb1.protein_id = pb2.protein_id
                        WHERE pb1.binder_id = %s AND pb2.binder_id = b2.binder_id
                    )
                    OR EXISTS (
                        SELECT 1
                        FROM binder_diseases bd1
                        JOIN binder_diseases bd2 ON bd1.disease_id = bd2.disease_id
                        WHERE bd1.binder_id = %s AND bd2.binder_id = b2.binder_id
                    )
                  )
                ORDER BY b2.binder_name
                LIMIT 8;
            """, (binder_id, binder_id, binder_id))
            related_binders = cur.fetchall()

            cur.execute("""
                SELECT *
                FROM (
                    SELECT DISTINCT
                        p.protein_id,
                        p.protein_name,
                        p.uniprot_accession,
                        g.gene_symbol,
                        s.structure_id,
                        s.pdb_id,
                        s.structure_title,
                        s.experimental_method,
                        s.resolution,
                        md5(%s || '-' || COALESCE(s.pdb_id, s.structure_id::text)) AS structure_sort_key
                    FROM protein_binders pb
                    JOIN proteins p ON pb.protein_id = p.protein_id
                    LEFT JOIN genes g ON p.gene_id = g.gene_id
                    JOIN protein_structures ps ON p.protein_id = ps.protein_id
                    JOIN structures s ON ps.structure_id = s.structure_id
                    WHERE pb.binder_id = %s
                ) ranked_structures
                ORDER BY structure_sort_key, protein_name;
            """, (str(binder_id), binder_id))
            target_structures = cur.fetchall()

            binder_binding_maps, binder_binding_mode = build_binder_binding_maps(cur, binder, proteins)

    target_count=len(proteins); trial_count=len(trials); disease_count=len(diseases); structure_count=len(structures)
    sequence_length=len((binder.get("sequence") or "").replace("\n", "").replace(" ", "")) if binder.get("sequence") else 0
    related_binders = decorate_binder_records(related_binders)
    trial_phase_breakdown = summarize_counter(Counter((t.get("phase") or "Unknown") for t in trials))
    disease_category_breakdown = summarize_counter(Counter((d.get("disease_category") or "Uncategorized") for d in diseases))
    proteins, proteins_pagination = paginate_records(proteins, get_page_param("proteins_page"))
    trials, trials_pagination = paginate_records(trials, get_page_param("trials_page"))
    structures, structures_pagination = paginate_records(structures, get_page_param("structures_page"))
    diseases, diseases_pagination = paginate_records(diseases, get_page_param("diseases_page"))
    related_binders, related_binders_pagination = paginate_records(related_binders, get_page_param("related_binders_page"))

    binder_story_parts = [
        f'{binder.get("binder_name", "This binder")} is currently classified as {binder.get("binder_type") or "an unclassified binder"}.',
        f'It is linked to {target_count} target{"s" if target_count != 1 else ""}, {trial_count} clinical trial{"s" if trial_count != 1 else ""}, {disease_count} disease tag{"s" if disease_count != 1 else ""}, and {structure_count} directly linked structure{"s" if structure_count != 1 else ""}.',
    ]
    if binder.get("clinical_status") or binder.get("approval_status"):
        binder_story_parts.append(f'Clinical status: {binder.get("clinical_status") or binder.get("approval_status")}.')
    if proteins:
        protein_names = ", ".join([(p.get("gene_symbol") or p.get("protein_name") or "Unknown target") for p in proteins[:3]])
        binder_story_parts.append(f'Primary linked targets include {protein_names}.')
    binder_story = " ".join(binder_story_parts)

    target_structure_candidates=[]
    for row in target_structures:
        pdb_id=safe_string(row.get("pdb_id"))
        if not pdb_id:
            continue
        gene_or_name=safe_string(row.get("gene_symbol")) or safe_string(row.get("protein_name")) or "Target"
        method=safe_string(row.get("experimental_method")) or "Experimental structure"
        resolution=row.get("resolution")
        resolution_text=f" · Resolution: {resolution}" if resolution is not None else ""
        target_structure_candidates.append({
            "source_type":"pdb","source_id":pdb_id.upper(),"label":f"{gene_or_name} · PDB {pdb_id.upper()}",
            "title":safe_string(row.get("structure_title")) or f"{gene_or_name} structure",
            "subtitle":f"{method}{resolution_text}","viewer_url":build_molstar_pdb_viewer_url(pdb_id),
            "external_url":f"https://www.rcsb.org/structure/{pdb_id.upper()}","origin":f"Target-linked structure ({gene_or_name})"
        })
    fallback_uniprot=None; fallback_title=None
    if proteins:
        first_target=proteins[0]; fallback_uniprot=first_target.get("uniprot_accession"); fallback_title=first_target.get("protein_name")
    direct_candidates=build_structure_candidates(structures)
    structure_candidates=dedupe_structure_candidates(direct_candidates + target_structure_candidates)
    if not structure_candidates and fallback_uniprot:
        structure_candidates = build_structure_candidates([], fallback_uniprot=fallback_uniprot, fallback_title=fallback_title)
    default_structure=structure_candidates[0] if structure_candidates else None
    binder_3d_annotations = [
    ann
    for mapping in binder_binding_maps
    for ann in mapping.get("annotations", [])
    ]

    binding_region_summary = build_binding_region_summary(
        binder_3d_annotations,
        binder_binding_mode,
    )
    binder_classification = build_binder_classification(
        binder,
        proteins=proteins,
        trials=trials,
        diseases=diseases,
        structures=structures,
    )
    binder_visualization = build_binder_visualization(
        binder,
        proteins=proteins,
        diseases=diseases,
        trials=trials,
        related_binders=related_binders,
    )
    # Sponsor-facing completeness model.
    # This treats a binder as complete only when the page can show its core identity,
    # sequence/status metadata, target links, disease-tag context, clinical evidence,
    # and a direct or fallback structure status.
    binder_completeness_input = dict(binder)
    binder_completeness_input.update({
        "sequence": binder.get("sequence"),
        "linked_targets": "Available" if target_count > 0 else None,
        "disease_tags": "Available" if disease_count > 0 else None,
        "clinical_trials": "Available" if trial_count > 0 else None,
        "structure_availability": "Available" if structure_candidates else None,
    })

    binder_completeness = build_record_completeness(
        binder_completeness_input,
        [
            ("binder_name", "Binder name"),
            ("binder_type", "Binder type"),
            ("sequence", "Sequence / FASTA"),
            ("clinical_status", "Clinical status"),
            ("modality_name", "Modality"),
            ("linked_targets", "Linked target(s)"),
            ("disease_tags", "Disease tag(s)"),
            ("clinical_trials", "Clinical trial evidence"),
            ("structure_availability", "Structure availability"),
            ("binder_description", "Description"),
        ],
    )

    structure_status = {
        "direct_count": structure_count,
        "viewer_count": len(structure_candidates or []),
        "has_direct_structure": structure_count > 0,
        "has_fallback_structure": bool(structure_candidates) and structure_count == 0,
        "message": (
            "Direct binder structure available." if structure_count > 0
            else "No direct binder structure is linked; using target-linked or predicted fallback structure when available." if structure_candidates
            else "No direct, target-linked, or predicted structure is currently available."
        ),
    }

    sequence_note = binder_sequence_message(binder) if not binder.get("sequence") else None

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
        default_structure=default_structure,
        binder_binding_maps=binder_binding_maps,
        binder_binding_mode=binder_binding_mode,
        binding_region_summary=binding_region_summary,
        binder_classification=binder_classification,
        binder_visualization=binder_visualization,
        binder_3d_annotations=binder_3d_annotations,
        binder_completeness=binder_completeness,
        structure_status=structure_status,
        sequence_note=sequence_note,
        display_value=display_value,
        proteins_pagination=proteins_pagination,
        trials_pagination=trials_pagination,
        structures_pagination=structures_pagination,
        diseases_pagination=diseases_pagination,
        related_binders_pagination=related_binders_pagination,
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

    binders = decorate_binder_records(binders)
    binder_type_breakdown = summarize_binder_classes(binders)
    trial_context_counts = {
        "binders": len(binders),
        "proteins": len(proteins),
        "diseases": len(diseases),
        "related_trials": len(related_trials)
    }
    binders, binders_pagination = paginate_records(binders, get_page_param("binders_page"))
    proteins, proteins_pagination = paginate_records(proteins, get_page_param("proteins_page"))
    diseases, diseases_pagination = paginate_records(diseases, get_page_param("diseases_page"))
    related_trials, related_trials_pagination = paginate_records(related_trials, get_page_param("related_trials_page"))

    return render_template(
        "trial_detail.html",
        trial=trial,
        diseases=diseases,
        binders=binders,
        proteins=proteins,
        related_trials=related_trials,
        binder_type_breakdown=binder_type_breakdown,
        trial_context_counts=trial_context_counts,
        binders_pagination=binders_pagination,
        proteins_pagination=proteins_pagination,
        diseases_pagination=diseases_pagination,
        related_trials_pagination=related_trials_pagination
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

    binders = decorate_binder_records(binders)
    binder_type_breakdown = summarize_binder_classes(binders)
    protein_gene_breakdown = summarize_counter(
        Counter((p.get("gene_symbol") or p.get("protein_name") or "Unknown target") for p in proteins[:12])
    )
    binders, binders_pagination = paginate_records(binders, get_page_param("binders_page"))
    proteins, proteins_pagination = paginate_records(proteins, get_page_param("proteins_page"))
    trials, trials_pagination = paginate_records(trials, get_page_param("trials_page"))
    related_diseases, related_diseases_pagination = paginate_records(related_diseases, get_page_param("related_diseases_page"))

    return render_template(
        "disease_detail.html",
        disease=disease,
        proteins=proteins,
        trials=trials,
        binders=binders,
        related_diseases=related_diseases,
        binder_type_breakdown=binder_type_breakdown,
        protein_gene_breakdown=protein_gene_breakdown,
        binders_pagination=binders_pagination,
        proteins_pagination=proteins_pagination,
        trials_pagination=trials_pagination,
        related_diseases_pagination=related_diseases_pagination
    )


if __name__ == "__main__":
    debug_mode = os.getenv("FLASK_DEBUG", "False").lower() == "true"
    app.run(debug=debug_mode)