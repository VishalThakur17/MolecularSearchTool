import requests
from ingest.db import get_connection

UNIPROT_SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"


def fetch_uniprot_by_gene(gene_symbol: str):
    params = {
        "query": f"gene_exact:{gene_symbol} AND organism_id:9606",
        "format": "json",
        "size": 1,
    }

    response = requests.get(UNIPROT_SEARCH_URL, params=params, timeout=30)
    response.raise_for_status()
    data = response.json()
    results = data.get("results", [])
    return results[0] if results else None


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
    organism = record.get("organism", {}).get("scientificName", "Homo sapiens")

    protein_name = gene_symbol
    protein_desc = record.get("proteinDescription", {})
    recommended = protein_desc.get("recommendedName", {})
    full_name = recommended.get("fullName", {})
    if isinstance(full_name, dict):
        protein_name = full_name.get("value", gene_symbol)

    sequence_obj = record.get("sequence", {}) or {}
    sequence = sequence_obj.get("value")
    sequence_length = sequence_obj.get("length")

    gene_name = gene_symbol
    genes = record.get("genes", [])
    if genes:
        gene_obj = genes[0].get("geneName")
        if isinstance(gene_obj, dict):
            gene_name = gene_obj.get("value", gene_symbol)

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
                    source_id
                )
                VALUES (%s, %s, %s, %s, %s, %s, %s, %s)
                ON CONFLICT (uniprot_accession)
                DO UPDATE SET
                    gene_id = EXCLUDED.gene_id,
                    protein_name = EXCLUDED.protein_name,
                    organism_name = EXCLUDED.organism_name,
                    sequence = EXCLUDED.sequence,
                    sequence_length = EXCLUDED.sequence_length,
                    functional_description = EXCLUDED.functional_description,
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
                    f"Imported from UniProt for {gene_symbol}",
                    source_id,
                ),
            )
            protein_id = cur.fetchone()[0]

    print(f"[UniProt] Upserted {gene_symbol} -> protein_id={protein_id}")
    return protein_id