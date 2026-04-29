import requests
from ingest.db import get_connection

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


def fetch_uniprot_by_gene(gene_symbol: str):
    """Fetch one reviewed human UniProtKB record for a gene symbol."""
    params = {
        "query": f"gene_exact:{gene_symbol} AND organism_id:9606 AND reviewed:true",
        "format": "json",
        "size": 1,
    }

    response = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=30)
    response.raise_for_status()
    data = response.json()
    results = data.get("results", [])
    return results[0] if results else None


def extract_comment_text(record: dict, comment_type: str) -> str | None:
    """Extract readable UniProt comment text, such as FUNCTION or SUBCELLULAR LOCATION."""
    values: list[str] = []

    for comment in record.get("comments", []) or []:
        if comment.get("commentType") != comment_type:
            continue

        for text_obj in comment.get("texts", []) or []:
            value = text_obj.get("value")
            if value:
                values.append(value.strip())

        # Subcellular location uses a nested structure, not always the generic texts field.
        for location_obj in comment.get("subcellularLocations", []) or []:
            location = location_obj.get("location", {}) or {}
            topology = location_obj.get("topology", {}) or {}
            orientation = location_obj.get("orientation", {}) or {}

            location_parts = []
            for obj in (location, topology, orientation):
                value = obj.get("value")
                if value:
                    location_parts.append(value.strip())

            if location_parts:
                values.append("; ".join(location_parts))

    if not values:
        return None

    # Preserve order while removing duplicates.
    seen = set()
    unique_values = []
    for value in values:
        key = value.lower()
        if key not in seen:
            seen.add(key)
            unique_values.append(value)

    return " ".join(unique_values)


def extract_protein_name(record: dict, fallback: str) -> str:
    protein_desc = record.get("proteinDescription", {}) or {}
    recommended = protein_desc.get("recommendedName", {}) or {}
    full_name = recommended.get("fullName", {}) or {}

    if isinstance(full_name, dict) and full_name.get("value"):
        return full_name["value"]

    submission_names = protein_desc.get("submissionNames", []) or []
    if submission_names:
        full_name = submission_names[0].get("fullName", {}) or {}
        if isinstance(full_name, dict) and full_name.get("value"):
            return full_name["value"]

    return fallback


def extract_gene_name(record: dict, fallback: str) -> str:
    genes = record.get("genes", []) or []
    if not genes:
        return fallback

    gene_obj = genes[0].get("geneName") or {}
    if isinstance(gene_obj, dict) and gene_obj.get("value"):
        return gene_obj["value"]

    return fallback


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


def upsert_gene_and_protein(gene_symbol: str):
    record = fetch_uniprot_by_gene(gene_symbol)
    if not record:
        print(f"[UniProt] No result found for {gene_symbol}")
        return None

    accession = record.get("primaryAccession")
    organism = (record.get("organism", {}) or {}).get("scientificName", "Homo sapiens")

    protein_name = extract_protein_name(record, gene_symbol)
    gene_name = extract_gene_name(record, gene_symbol)

    sequence_obj = record.get("sequence", {}) or {}
    sequence = sequence_obj.get("value")
    sequence_length = sequence_obj.get("length")

    functional_description = extract_comment_text(record, "FUNCTION")
    subcellular_location = extract_comment_text(record, "SUBCELLULAR LOCATION")

    with get_connection() as conn:
        with conn.cursor() as cur:
            source_id = upsert_source(cur, "UniProt", "https://www.uniprot.org")

            cur.execute(
                """
                INSERT INTO genes (gene_symbol, gene_name)
                VALUES (%s, %s)
                ON CONFLICT (gene_symbol)
                DO UPDATE SET gene_name = EXCLUDED.gene_name
                RETURNING gene_id;
                """,
                (gene_symbol, gene_name),
            )
            gene_id = cur.fetchone()[0]

            cur.execute(
                """
                INSERT INTO proteins (
                    gene_id,
                    uniprot_accession,
                    protein_name,
                    organism_name,
                    sequence,
                    sequence_length,
                    functional_description,
                    subcellular_location,
                    source_id
                )
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (uniprot_accession)
                DO UPDATE SET
                    gene_id = EXCLUDED.gene_id,
                    protein_name = EXCLUDED.protein_name,
                    organism_name = EXCLUDED.organism_name,
                    sequence = COALESCE(EXCLUDED.sequence, proteins.sequence),
                    sequence_length = COALESCE(EXCLUDED.sequence_length, proteins.sequence_length),
                    functional_description = COALESCE(EXCLUDED.functional_description, proteins.functional_description),
                    subcellular_location = COALESCE(EXCLUDED.subcellular_location, proteins.subcellular_location),
                    source_id = EXCLUDED.source_id
                RETURNING protein_id;
                """,
                (
                    gene_id,
                    accession,
                    protein_name,
                    organism,
                    sequence,
                    sequence_length,
                    functional_description,
                    subcellular_location,
                    source_id,
                ),
            )
            protein_id = cur.fetchone()[0]

    print(
        f"[UniProt] Upserted {gene_symbol} -> protein_id={protein_id} | "
        f"length={sequence_length or 'missing'} | "
        f"subcellular={'yes' if subcellular_location else 'missing'} | "
        f"function={'yes' if functional_description else 'missing'}"
    )
    return protein_id
