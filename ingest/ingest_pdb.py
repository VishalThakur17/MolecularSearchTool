import time
import requests
from ingest.db import get_connection

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry"


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
    if not date_str:
        return None
    date_str = str(date_str).strip()
    return date_str[:10] if len(date_str) >= 10 else None


def search_pdb_by_gene(gene_symbol: str, max_results: int = 50):
    payload = {
        "query": {
            "type": "group",
            "logical_operator": "or",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                        "operator": "exact_match",
                        "value": gene_symbol,
                    },
                },
                {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": gene_symbol
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": max_results,
            },
            "results_content_type": ["experimental"],
        },
    }

    response = requests.post(RCSB_SEARCH_URL, json=payload, timeout=45)
    response.raise_for_status()

    data = response.json()
    result_set = data.get("result_set", [])

    pdb_ids = []
    seen = set()

    for row in result_set:
        pdb_id = row.get("identifier")
        if pdb_id and pdb_id.upper() not in seen:
            seen.add(pdb_id.upper())
            pdb_ids.append(pdb_id.upper())

    return pdb_ids


def fetch_pdb_entry(pdb_id: str):
    response = requests.get(f"{RCSB_ENTRY_URL}/{pdb_id}", timeout=45)
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
        WHERE UPPER(g.gene_symbol) = UPPER(%s)
        LIMIT 1;
        """,
        (gene_symbol,),
    )
    row = cur.fetchone()
    return row[0] if row else None


def upsert_structures_for_gene(gene_symbol: str, max_results: int = 50):
    try:
        pdb_ids = search_pdb_by_gene(gene_symbol, max_results=max_results)
    except Exception as e:
        print(f"[PDB] Search failed for {gene_symbol}: {e}")
        return 0

    if not pdb_ids:
        print(f"[PDB] No structures found for {gene_symbol}")
        return 0

    inserted = 0

    with get_connection() as conn:
        with conn.cursor() as cur:
            source_id = upsert_source(cur, "RCSB PDB", "https://www.rcsb.org/")
            protein_id = get_protein_id_by_gene(cur, gene_symbol)

            if not protein_id:
                print(f"[PDB] No protein found in DB for gene {gene_symbol}. Skipping.")
                return 0

            for pdb_id in pdb_ids:
                try:
                    entry = fetch_pdb_entry(pdb_id)
                    time.sleep(0.05)
                except Exception as e:
                    print(f"[PDB] Failed entry fetch for {pdb_id}: {e}")
                    continue

                title = safe_get(entry, "struct", "title")

                experimental_method = None
                exptl = entry.get("exptl")
                if isinstance(exptl, list) and exptl:
                    experimental_method = exptl[0].get("method")

                resolution = None
                resolution_list = (entry.get("rcsb_entry_info") or {}).get("resolution_combined")
                if isinstance(resolution_list, list) and resolution_list:
                    try:
                        resolution = float(resolution_list[0])
                    except Exception:
                        resolution = None

                deposition_date = normalize_pdb_date(
                    safe_get(entry, "rcsb_accession_info", "deposit_date")
                )

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
                    DO UPDATE SET evidence_summary = EXCLUDED.evidence_summary;
                    """,
                    (
                        protein_id,
                        structure_id,
                        f"Imported from expanded RCSB PDB search for {gene_symbol}",
                    ),
                )

                inserted += 1

    print(f"[PDB] Upserted {inserted} structures for {gene_symbol}")
    return inserted