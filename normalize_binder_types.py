import os
import psycopg2
import psycopg2.extras
from dotenv import load_dotenv

# Load database credentials from .env
load_dotenv()


FULL_TAXONOMY = [
    "Small Molecule",
    "IgG",
    "Bispecific Antibody",
    "ADC",
    "VHH",
    "Fc Fusion",
    "Peptide",
    "Other",
]


def normalize_binder_type(raw_type, name="", description="", mechanism=""):
    text = " ".join([
        raw_type or "",
        name or "",
        description or "",
        mechanism or ""
    ]).lower()

    # ADC must come before Small Molecule because ADCs contain "drug"
    if any(keyword in text for keyword in [
        "adc",
        "antibody-drug conjugate",
        "antibody drug conjugate",
        "drug conjugate",
        "vedotin",
        "deruxtecan",
        "emtansine",
        "mafodotin",
        "duocarmycin",
        "duocarmazine",
        "tesirine",
        "ozogamicin",
        "maytansinoid",
        "payload"
    ]):
        return "ADC"

    if "bispecific" in text or "bi-specific" in text:
        return "Bispecific Antibody"

    if "fc fusion" in text or "fusion protein" in text:
        return "Fc Fusion"

    if "nanobody" in text or "vhh" in text or "single-domain antibody" in text:
        return "VHH"

    if "peptide" in text:
        return "Peptide"

    if any(k in text for k in [
        "monoclonal antibody",
        "antibody",
        "mab",
        "imab",
        "zumab",
        "umab",
        "ximab",
        "omab"
    ]):
        return "IgG"

    # Small Molecule should come after biologic categories
    if any(keyword in text for keyword in [
        "small molecule",
        "inhibitor",
        "kinase inhibitor",
        "tyrosine kinase",
        "oral drug",
        "chemical",
        "compound",
        "small-molecule",
        "drug"
    ]):
        return "Small Molecule"

    return "Other"


def get_connection():
    return psycopg2.connect(
        host=os.getenv("DB_HOST"),
        user=os.getenv("DB_USER"),
        password=os.getenv("DB_PASSWORD"),
        dbname=os.getenv("DB_NAME"),
        port=os.getenv("DB_PORT", 5432)
    )


def main():
    conn = get_connection()

    with conn:
        with conn.cursor(cursor_factory=psycopg2.extras.RealDictCursor) as cur:
            cur.execute("""
                SELECT binder_id, binder_name, binder_type,
                       description, mechanism_of_action
                FROM binders;
            """)
            binders = cur.fetchall()

            updated_count = 0

            for b in binders:
                new_type = normalize_binder_type(
                    b["binder_type"],
                    b["binder_name"],
                    b.get("description"),
                    b.get("mechanism_of_action")
                )

                cur.execute("""
                    UPDATE binders
                    SET binder_type = %s
                    WHERE binder_id = %s;
                """, (new_type, b["binder_id"]))

                updated_count += 1

    conn.close()

    print("Binder types normalized with FULL taxonomy + Small Molecule separation.")
    print(f"Updated binders: {updated_count}")
    print("Allowed taxonomy:")
    for binder_type in FULL_TAXONOMY:
        print(f"   - {binder_type}")


if __name__ == "__main__":
    main()