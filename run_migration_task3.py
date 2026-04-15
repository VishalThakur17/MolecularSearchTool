from ingest.db import get_connection

sql = """
CREATE TABLE IF NOT EXISTS protein_diseases (
    protein_id INT NOT NULL REFERENCES proteins(protein_id) ON DELETE CASCADE,
    disease_id INT NOT NULL REFERENCES diseases(disease_id) ON DELETE CASCADE,
    source_id INT NULL REFERENCES sources(source_id),
    tag_reason TEXT,
    PRIMARY KEY (protein_id, disease_id)
);

CREATE TABLE IF NOT EXISTS binder_diseases (
    binder_id INT NOT NULL REFERENCES binders(binder_id) ON DELETE CASCADE,
    disease_id INT NOT NULL REFERENCES diseases(disease_id) ON DELETE CASCADE,
    source_id INT NULL REFERENCES sources(source_id),
    tag_reason TEXT,
    PRIMARY KEY (binder_id, disease_id)
);

CREATE TABLE IF NOT EXISTS trial_diseases (
    trial_id INT NOT NULL REFERENCES clinical_trials(trial_id) ON DELETE CASCADE,
    disease_id INT NOT NULL REFERENCES diseases(disease_id) ON DELETE CASCADE,
    source_id INT NULL REFERENCES sources(source_id),
    tag_reason TEXT,
    PRIMARY KEY (trial_id, disease_id)
);

CREATE INDEX IF NOT EXISTS idx_protein_diseases_protein_id ON protein_diseases(protein_id);
CREATE INDEX IF NOT EXISTS idx_protein_diseases_disease_id ON protein_diseases(disease_id);

CREATE INDEX IF NOT EXISTS idx_binder_diseases_binder_id ON binder_diseases(binder_id);
CREATE INDEX IF NOT EXISTS idx_binder_diseases_disease_id ON binder_diseases(disease_id);

CREATE INDEX IF NOT EXISTS idx_trial_diseases_trial_id ON trial_diseases(trial_id);
CREATE INDEX IF NOT EXISTS idx_trial_diseases_disease_id ON trial_diseases(disease_id);
"""

with get_connection() as conn:
    with conn.cursor() as cur:
        cur.execute(sql)
    conn.commit()

print("Priority 3 migration completed successfully")