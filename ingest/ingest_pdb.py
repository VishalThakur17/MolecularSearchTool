import requests
from ingest.db import get_connection

# RCSB Search API
RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
# RCSB Data API
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry"


def search_pdb_by_gene(gene_symbol: str, max_results: int = 10):
    """
    Search PDB entries using an unstructured full-text query.
    This is the simplest stable alpha approach for gene symbols like EGFR, KRAS, ALK.
    """
    payload = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": gene_symbol
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": max_results
            }
        }
    }

    response = requests.post(RCSB_SEARCH_URL, json=payload, timeout=45)
    response.raise_for_status()
    data = response.json()
    result_set = data.get("result_set", [])
    return [r.get("identifier") for r in result_set if r.get("identifier")]


def fetch_pdb_entry(pdb_id: str):
    url = f"{RCSB_ENTRY_URL}/{pdb_id}"
    response = requests.get(url, timeout=45)
    response.raise_for_status()
    return response.json()


def upsert_source(cur, source_name: str, source_url: str):
    cur.execute(
        """
        INSERT INTO sources (source_name, source_url)
        VALUES (%s, %s)
        ON CONFLICT (source_name)
        DO UPDATE SET source_url = EXCLUDED.source_url
        RETURNING source_id;
        """,
        (source_name, source_url),
    )
    return cur.fetchone()[0]


def get_protein_id_by_gene(cur, gene_symbol: str):
    cur.execute(
        """
        SELECT p.protein_id
        FROM proteins p
        JOIN genes g ON p.gene_id = g.gene_id
        WHERE g.gene_symbol = %s
        LIMIT 1;
        """,
        (gene_symbol,),
    )
    row = cur.fetchone()
    return row[0] if row else None


def safe_get(dct, *keys, default=None):
    current = dct
    for key in keys:
        if not isinstance(current, dict):
            return default
        current = current.get(key)
        if current is None:
            return default
    return current


def normalize_pdb_date(date_str):
    """
    Keep only YYYY-MM-DD values for PostgreSQL DATE compatibility.
    """
    if not date_str:
        return None

    date_str = str(date_str).strip()

    if len(date_str) >= 10:
        return date_str[:10]

    return None


def upsert_structures_for_gene(gene_symbol: str, max_results: int = 10):
    """
    Search RCSB PDB by gene symbol and map structures into:
    - structures
    - protein_structures
    """
    pdb_ids = search_pdb_by_gene(gene_symbol, max_results=max_results)

    if not pdb_ids:
        print(f"[PDB] No structures found for {gene_symbol}")
        return

    with get_connection() as conn:
        with conn.cursor() as cur:
            source_id = upsert_source(cur, "RCSB PDB", "https://www.rcsb.org/")
            protein_id = get_protein_id_by_gene(cur, gene_symbol)

            if not protein_id:
                print(f"[PDB] No protein found in DB for gene {gene_symbol}. Skipping structure linking.")
                return

            inserted = 0

            for pdb_id in pdb_ids:
                try:
                    entry = fetch_pdb_entry(pdb_id)
                except Exception as e:
                    print(f"[PDB] Failed entry fetch for {pdb_id}: {e}")
                    continue

                title = safe_get(entry, "struct", "title", default=None)

                experimental_method = None
                exptl = entry.get("exptl")
                if isinstance(exptl, list) and exptl:
                    experimental_method = exptl[0].get("method")

                resolution = None
                rcsb_entry_info = entry.get("rcsb_entry_info") or {}
                resolution_list = rcsb_entry_info.get("resolution_combined")
                if isinstance(resolution_list, list) and resolution_list:
                    try:
                        resolution = float(resolution_list[0])
                    except Exception:
                        resolution = None

                raw_deposition_date = safe_get(entry, "rcsb_accession_info", "deposit_date", default=None)
                deposition_date = normalize_pdb_date(raw_deposition_date)

                structure_file_url = f"https://files.rcsb.org/download/{pdb_id}.cif"

                cur.execute(
                    """
                    INSERT INTO structures (
                        pdb_id,
                        structure_title,
                        experimental_method,
                        resolution,
                        deposition_date,
                        structure_file_url,
                        source_id
                    )
                    VALUES (%s, %s, %s, %s, %s, %s, %s)
                    ON CONFLICT (pdb_id)
                    DO UPDATE SET
                        structure_title = EXCLUDED.structure_title,
                        experimental_method = EXCLUDED.experimental_method,
                        resolution = EXCLUDED.resolution,
                        deposition_date = EXCLUDED.deposition_date,
                        structure_file_url = EXCLUDED.structure_file_url,
                        source_id = EXCLUDED.source_id
                    RETURNING structure_id;
                    """,
                    (
                        pdb_id,
                        title,
                        experimental_method,
                        resolution,
                        deposition_date,
                        structure_file_url,
                        source_id,
                    ),
                )
                structure_id = cur.fetchone()[0]

                cur.execute(
                    """
                    INSERT INTO protein_structures (
                        protein_id,
                        structure_id,
                        evidence_summary
                    )
                    VALUES (%s, %s, %s)
                    ON CONFLICT (protein_id, structure_id)
                    DO UPDATE SET
                        evidence_summary = EXCLUDED.evidence_summary;
                    """,
                    (
                        protein_id,
                        structure_id,
                        f"Imported from RCSB PDB search for {gene_symbol}",
                    ),
                )

                inserted += 1

    print(f"[PDB] Upserted {inserted} structures for {gene_symbol}")