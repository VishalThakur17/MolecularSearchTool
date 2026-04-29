from ingest.db import get_connection


def ensure_binder_structures_table(cur):
    cur.execute(
        """
        CREATE TABLE IF NOT EXISTS binder_structures (
            binder_id INT NOT NULL REFERENCES binders(binder_id) ON DELETE CASCADE,
            structure_id INT NOT NULL REFERENCES structures(structure_id) ON DELETE CASCADE,
            evidence_summary TEXT,
            PRIMARY KEY (binder_id, structure_id)
        );
        """
    )

    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS idx_binder_structures_binder_id
        ON binder_structures(binder_id);
        """
    )

    cur.execute(
        """
        CREATE INDEX IF NOT EXISTS idx_binder_structures_structure_id
        ON binder_structures(structure_id);
        """
    )


def rebuild_binder_structure_links(cur, max_structures_per_binder=5):
    ensure_binder_structures_table(cur)

    cur.execute("DELETE FROM binder_structures;")

    cur.execute(
        """
        WITH unique_links AS (
            SELECT DISTINCT
                pb.binder_id,
                ps.structure_id,
                s.pdb_id,
                s.resolution
            FROM protein_binders pb
            JOIN protein_structures ps
                ON pb.protein_id = ps.protein_id
            JOIN structures s
                ON ps.structure_id = s.structure_id
        ),
        ranked_structures AS (
            SELECT
                binder_id,
                structure_id,
                ROW_NUMBER() OVER (
                    PARTITION BY binder_id
                    ORDER BY
                        CASE WHEN resolution IS NULL THEN 1 ELSE 0 END,
                        resolution ASC,
                        md5(binder_id::text || '-' || COALESCE(pdb_id, structure_id::text))
                ) AS rn
            FROM unique_links
        )
        INSERT INTO binder_structures (
            binder_id,
            structure_id,
            evidence_summary
        )
        SELECT
            binder_id,
            structure_id,
            'Linked through shared target protein with diversified structure selection'
        FROM ranked_structures
        WHERE rn <= %s
        ON CONFLICT (binder_id, structure_id)
        DO NOTHING;
        """,
        (max_structures_per_binder,),
    )

    print(
        f"[Finalize] Rebuilt binder_structures with up to "
        f"{max_structures_per_binder} structures per binder."
    )


def print_summary(cur):
    cur.execute("SELECT COUNT(*) FROM binder_structures;")
    binder_structure_count = cur.fetchone()[0]

    cur.execute("SELECT COUNT(*) FROM structures;")
    structure_count = cur.fetchone()[0]

    cur.execute("SELECT COUNT(*) FROM protein_structures;")
    protein_structure_count = cur.fetchone()[0]

    print("\n=== RELATIONSHIP SUMMARY ===")
    print(f"Structures: {structure_count}")
    print(f"Protein-Structure links: {protein_structure_count}")
    print(f"Binder-Structure links: {binder_structure_count}")


def main():
    with get_connection() as conn:
        with conn.cursor() as cur:
            rebuild_binder_structure_links(cur, max_structures_per_binder=5)
            print_summary(cur)
        conn.commit()

    print("\n[Finalize] Done.")


if __name__ == "__main__":
    main()