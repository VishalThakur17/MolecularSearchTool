import requests
from ingest.db import get_connection

CTGOV_V2_URL = "https://clinicaltrials.gov/api/v2/studies"


def fetch_trials(term: str, page_size: int = 10):
    params = {
        "query.term": term,
        "pageSize": page_size,
        "format": "json",
    }

    response = requests.get(CTGOV_V2_URL, params=params, timeout=45)
    response.raise_for_status()
    data = response.json()
    return data.get("studies", [])


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


def get_or_create_disease(cur, disease_name: str):
    cur.execute(
        """
        INSERT INTO diseases (disease_name, disease_category, description)
        VALUES (%s, %s, %s)
        ON CONFLICT (disease_name)
        DO UPDATE SET disease_name = EXCLUDED.disease_name
        RETURNING disease_id;
        """,
        (disease_name, "Oncology", f"Imported/linked disease record for {disease_name}"),
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


def first_list_value(value, default=None):
    if isinstance(value, list) and value:
        return value[0]
    return value if value is not None else default


def normalize_ctgov_date(date_str):
    """
    Convert ClinicalTrials.gov partial dates into PostgreSQL-safe YYYY-MM-DD dates.

    Examples:
    - 2036-11    -> 2036-11-01
    - 2027       -> 2027-01-01
    - 2025-08-03 -> 2025-08-03
    """
    if not date_str:
        return None

    date_str = str(date_str).strip()

    if len(date_str) == 10:
        # Already YYYY-MM-DD
        return date_str

    if len(date_str) == 7:
        # YYYY-MM
        return f"{date_str}-01"

    if len(date_str) == 4:
        # YYYY
        return f"{date_str}-01-01"

    return None


def extract_trial_fields(study: dict):
    protocol = study.get("protocolSection", {})

    identification = protocol.get("identificationModule", {})
    status = protocol.get("statusModule", {})
    conditions = protocol.get("conditionsModule", {})
    design = protocol.get("designModule", {})
    description = protocol.get("descriptionModule", {})
    sponsor = protocol.get("sponsorCollaboratorsModule", {})

    nct_id = identification.get("nctId")
    title = identification.get("briefTitle")
    condition_name = first_list_value(conditions.get("conditions"), default=None)
    phase = first_list_value(design.get("phases"), default=None)
    recruitment_status = status.get("overallStatus")
    study_type = design.get("studyType")
    brief_summary = description.get("briefSummary")
    sponsor_name = safe_get(sponsor, "leadSponsor", "name", default=None)

    start_date_struct = status.get("startDateStruct", {})
    completion_date_struct = status.get("completionDateStruct", {})

    raw_start_date = start_date_struct.get("date")
    raw_completion_date = completion_date_struct.get("date")

    start_date = normalize_ctgov_date(raw_start_date)
    completion_date = normalize_ctgov_date(raw_completion_date)

    return {
        "nct_id": nct_id,
        "trial_title": title,
        "condition_name": condition_name,
        "phase": phase,
        "recruitment_status": recruitment_status,
        "study_type": study_type,
        "brief_summary": brief_summary,
        "sponsor_name": sponsor_name,
        "start_date": start_date,
        "completion_date": completion_date,
        "trial_url": f"https://clinicaltrials.gov/study/{nct_id}" if nct_id else None,
    }


def upsert_trials_for_target(gene_symbol: str, disease_name: str, max_trials: int = 10):
    search_term = f'{gene_symbol} AND "{disease_name}"'
    studies = fetch_trials(search_term, page_size=max_trials)

    with get_connection() as conn:
        with conn.cursor() as cur:
            source_id = upsert_source(
                cur,
                "ClinicalTrials.gov",
                "https://clinicaltrials.gov",
            )

            disease_id = get_or_create_disease(cur, disease_name)
            protein_id = get_protein_id_by_gene(cur, gene_symbol)

            inserted = 0

            for study in studies:
                fields = extract_trial_fields(study)

                if not fields["nct_id"] or not fields["trial_title"]:
                    continue

                cur.execute(
                    """
                    INSERT INTO clinical_trials (
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
                        completion_date,
                        source_id
                    )
                    VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
                    ON CONFLICT (nct_id)
                    DO UPDATE SET
                        trial_title = EXCLUDED.trial_title,
                        condition_name = EXCLUDED.condition_name,
                        phase = EXCLUDED.phase,
                        recruitment_status = EXCLUDED.recruitment_status,
                        study_type = EXCLUDED.study_type,
                        sponsor_name = EXCLUDED.sponsor_name,
                        brief_summary = EXCLUDED.brief_summary,
                        trial_url = EXCLUDED.trial_url,
                        start_date = EXCLUDED.start_date,
                        completion_date = EXCLUDED.completion_date,
                        source_id = EXCLUDED.source_id
                    RETURNING trial_id;
                    """,
                    (
                        fields["nct_id"],
                        fields["trial_title"],
                        fields["condition_name"],
                        fields["phase"],
                        fields["recruitment_status"],
                        fields["study_type"],
                        fields["sponsor_name"],
                        fields["brief_summary"],
                        fields["trial_url"],
                        fields["start_date"],
                        fields["completion_date"],
                        source_id,
                    ),
                )
                trial_id = cur.fetchone()[0]

                cur.execute(
                    """
                    INSERT INTO disease_trials (disease_id, trial_id)
                    VALUES (%s, %s)
                    ON CONFLICT (disease_id, trial_id) DO NOTHING;
                    """,
                    (disease_id, trial_id),
                )

                if protein_id:
                    cur.execute(
                        """
                        INSERT INTO protein_trials (protein_id, trial_id, evidence_summary)
                        VALUES (%s, %s, %s)
                        ON CONFLICT (protein_id, trial_id)
                        DO NOTHING;
                        """,
                        (
                            protein_id,
                            trial_id,
                            f"Linked by ingestion query for {gene_symbol} and {disease_name}",
                        ),
                    )

                inserted += 1

    print(f"[Trials] Upserted {inserted} trials for {gene_symbol} / {disease_name}")