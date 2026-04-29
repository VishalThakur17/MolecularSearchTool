-- 004_full_binder_model_sync.sql
-- Safe schema/data sync for the sponsor binder model.
-- This file is designed to be rerunnable on an existing development database.

-- ============================================================
-- 1) Core binder columns
-- ============================================================
ALTER TABLE binders
ADD COLUMN IF NOT EXISTS binder_type TEXT,
ADD COLUMN IF NOT EXISTS sequence TEXT,
ADD COLUMN IF NOT EXISTS clinical_status TEXT,
ADD COLUMN IF NOT EXISTS sequence_source TEXT,
ADD COLUMN IF NOT EXISTS structure_notes TEXT,
ADD COLUMN IF NOT EXISTS data_completeness_score INTEGER;

-- ============================================================
-- 2) Relationship tables and defensive column syncing
-- ============================================================
CREATE TABLE IF NOT EXISTS binder_diseases (
    binder_id INT NOT NULL REFERENCES binders(binder_id) ON DELETE CASCADE,
    disease_id INT NOT NULL REFERENCES diseases(disease_id) ON DELETE CASCADE,
    PRIMARY KEY (binder_id, disease_id)
);

ALTER TABLE binder_diseases
ADD COLUMN IF NOT EXISTS source_id INT NULL REFERENCES sources(source_id),
ADD COLUMN IF NOT EXISTS tag_reason TEXT;

CREATE INDEX IF NOT EXISTS idx_binder_diseases_binder_id ON binder_diseases(binder_id);
CREATE INDEX IF NOT EXISTS idx_binder_diseases_disease_id ON binder_diseases(disease_id);

CREATE TABLE IF NOT EXISTS binder_trials (
    binder_id INT NOT NULL REFERENCES binders(binder_id) ON DELETE CASCADE,
    trial_id INT NOT NULL REFERENCES clinical_trials(trial_id) ON DELETE CASCADE,
    PRIMARY KEY (binder_id, trial_id)
);

ALTER TABLE binder_trials
ADD COLUMN IF NOT EXISTS relationship_type TEXT DEFAULT 'associated evidence',
ADD COLUMN IF NOT EXISTS source_id INT NULL REFERENCES sources(source_id),
ADD COLUMN IF NOT EXISTS notes TEXT;

CREATE INDEX IF NOT EXISTS idx_binder_trials_binder_id ON binder_trials(binder_id);
CREATE INDEX IF NOT EXISTS idx_binder_trials_trial_id ON binder_trials(trial_id);

CREATE TABLE IF NOT EXISTS binder_structures (
    binder_id INT NOT NULL REFERENCES binders(binder_id) ON DELETE CASCADE,
    structure_id INT NOT NULL REFERENCES structures(structure_id) ON DELETE CASCADE,
    PRIMARY KEY (binder_id, structure_id)
);

ALTER TABLE binder_structures
ADD COLUMN IF NOT EXISTS relationship_type TEXT DEFAULT 'direct binder structure',
ADD COLUMN IF NOT EXISTS source_id INT NULL REFERENCES sources(source_id),
ADD COLUMN IF NOT EXISTS notes TEXT;

CREATE INDEX IF NOT EXISTS idx_binder_structures_binder_id ON binder_structures(binder_id);
CREATE INDEX IF NOT EXISTS idx_binder_structures_structure_id ON binder_structures(structure_id);

-- ============================================================
-- 3) Normalize clinical status
-- ============================================================
UPDATE binders
SET clinical_status = CASE
    WHEN clinical_status IS NOT NULL AND TRIM(clinical_status) <> '' THEN clinical_status
    WHEN approval_status ILIKE '%approved%' THEN 'Approved'
    WHEN approval_status ILIKE '%max phase 4%' THEN 'Approved'
    WHEN approval_status ILIKE '%max phase 3%' THEN 'Phase 3'
    WHEN approval_status ILIKE '%max phase 2%' THEN 'Phase 2'
    WHEN approval_status ILIKE '%max phase 1%' THEN 'Phase 1'
    WHEN approval_status ILIKE '%max phase 0%' THEN 'Preclinical'
    ELSE COALESCE(NULLIF(approval_status, ''), 'Unknown')
END;

-- ============================================================
-- 4) Normalize binder type from modality/name where possible
-- ============================================================
UPDATE binders b
SET binder_type = CASE
    WHEN b.binder_type IS NOT NULL AND TRIM(b.binder_type) <> '' THEN b.binder_type
    WHEN bm.modality_name ILIKE '%adc%' OR b.binder_name ILIKE '%vedotin%' OR b.binder_name ILIKE '%deruxtecan%' OR b.binder_name ILIKE '%emtansine%' THEN 'ADC'
    WHEN bm.modality_name ILIKE '%bispecific%' OR b.binder_name ILIKE '%bispecific%' THEN 'Bispecific Antibody'
    WHEN bm.modality_name ILIKE '%vhh%' OR bm.modality_name ILIKE '%nanobody%' OR b.binder_name ILIKE '%nanobody%' THEN 'VHH'
    WHEN bm.modality_name ILIKE '%peptide%' OR b.binder_name ILIKE '%peptide%' THEN 'Peptide'
    WHEN bm.modality_name ILIKE '%antibody%' OR bm.modality_name ILIKE '%igg%' OR LOWER(b.binder_name) LIKE '%mab' THEN 'IgG'
    WHEN bm.modality_name ILIKE '%small molecule%' THEN 'Small Molecule'
    ELSE COALESCE(NULLIF(bm.modality_name, ''), 'Other')
END
FROM binder_modalities bm
WHERE b.modality_id = bm.modality_id;

UPDATE binders
SET binder_type = CASE
    WHEN binder_type IS NOT NULL AND TRIM(binder_type) <> '' THEN binder_type
    WHEN binder_name ILIKE '%vedotin%' OR binder_name ILIKE '%deruxtecan%' OR binder_name ILIKE '%emtansine%' THEN 'ADC'
    WHEN binder_name ILIKE '%bispecific%' THEN 'Bispecific Antibody'
    WHEN binder_name ILIKE '%nanobody%' THEN 'VHH'
    WHEN binder_name ILIKE '%peptide%' THEN 'Peptide'
    WHEN LOWER(binder_name) LIKE '%mab' THEN 'IgG'
    ELSE 'Other'
END
WHERE binder_type IS NULL OR TRIM(binder_type) = '';

-- ============================================================
-- 5) Backfill disease tags from target disease associations
-- ============================================================
INSERT INTO binder_diseases (binder_id, disease_id, source_id, tag_reason)
SELECT DISTINCT
    pb.binder_id,
    pd.disease_id,
    pd.source_id,
    'Backfilled from linked target disease association'
FROM protein_binders pb
JOIN protein_diseases pd ON pb.protein_id = pd.protein_id
ON CONFLICT (binder_id, disease_id) DO UPDATE SET
    source_id = COALESCE(EXCLUDED.source_id, binder_diseases.source_id),
    tag_reason = COALESCE(binder_diseases.tag_reason, EXCLUDED.tag_reason);

-- ============================================================
-- 6) Backfill clinical evidence from target clinical trial links
-- ============================================================
INSERT INTO binder_trials (binder_id, trial_id, relationship_type, source_id, notes)
SELECT DISTINCT
    pb.binder_id,
    pt.trial_id,
    'inferred from linked target clinical trial',
    pt.source_id,
    'Backfilled from protein_trials through protein_binders'
FROM protein_binders pb
JOIN protein_trials pt ON pb.protein_id = pt.protein_id
ON CONFLICT (binder_id, trial_id) DO UPDATE SET
    relationship_type = COALESCE(binder_trials.relationship_type, EXCLUDED.relationship_type),
    source_id = COALESCE(EXCLUDED.source_id, binder_trials.source_id),
    notes = COALESCE(binder_trials.notes, EXCLUDED.notes);

-- ============================================================
-- 7) Add a structure status note. Direct structures are optional, but
-- the UI can still show target-linked structures as fallback.
-- ============================================================
UPDATE binders b
SET structure_notes = CASE
    WHEN EXISTS (SELECT 1 FROM binder_structures bs WHERE bs.binder_id = b.binder_id)
        THEN 'Direct binder structure linked.'
    WHEN EXISTS (
        SELECT 1
        FROM protein_binders pb
        JOIN protein_structures ps ON ps.protein_id = pb.protein_id
        WHERE pb.binder_id = b.binder_id
    )
        THEN 'No direct binder structure linked; target-linked structure fallback is available.'
    ELSE 'No direct or target-linked structure currently available.'
END;

UPDATE binders
SET sequence_source = COALESCE(sequence_source, CASE WHEN sequence IS NOT NULL AND TRIM(sequence) <> '' THEN 'Ingested source sequence' ELSE NULL END);

-- ============================================================
-- 8) Compute data completeness score
-- ============================================================
WITH binder_counts AS (
    SELECT
        b.binder_id,
        CASE WHEN COALESCE(NULLIF(b.binder_name, ''), '') <> '' THEN 1 ELSE 0 END AS has_name,
        CASE WHEN COALESCE(NULLIF(b.binder_type, ''), '') <> '' THEN 1 ELSE 0 END AS has_type,
        CASE WHEN COALESCE(NULLIF(b.sequence, ''), '') <> '' THEN 1 ELSE 0 END AS has_sequence,
        CASE WHEN COALESCE(NULLIF(b.clinical_status, ''), NULLIF(b.approval_status, ''), '') <> '' THEN 1 ELSE 0 END AS has_status,
        CASE WHEN EXISTS (SELECT 1 FROM protein_binders pb WHERE pb.binder_id = b.binder_id) THEN 1 ELSE 0 END AS has_target,
        CASE WHEN EXISTS (SELECT 1 FROM binder_diseases bd WHERE bd.binder_id = b.binder_id) THEN 1 ELSE 0 END AS has_disease,
        CASE WHEN EXISTS (SELECT 1 FROM binder_trials bt WHERE bt.binder_id = b.binder_id) THEN 1 ELSE 0 END AS has_trial,
        CASE WHEN EXISTS (SELECT 1 FROM binder_structures bs WHERE bs.binder_id = b.binder_id)
               OR EXISTS (
                    SELECT 1
                    FROM protein_binders pb
                    JOIN protein_structures ps ON ps.protein_id = pb.protein_id
                    WHERE pb.binder_id = b.binder_id
               )
             THEN 1 ELSE 0 END AS has_structure_or_fallback
    FROM binders b
)
UPDATE binders b
SET data_completeness_score = ROUND((
    bc.has_name + bc.has_type + bc.has_sequence + bc.has_status +
    bc.has_target + bc.has_disease + bc.has_trial + bc.has_structure_or_fallback
) * 100.0 / 8.0)
FROM binder_counts bc
WHERE b.binder_id = bc.binder_id;

-- ============================================================
-- 9) Audit view for quick inspection in pgAdmin
-- ============================================================
CREATE OR REPLACE VIEW binder_completeness_audit AS
SELECT
    b.binder_id,
    b.binder_name,
    b.binder_type,
    b.clinical_status,
    CASE WHEN b.sequence IS NOT NULL AND TRIM(b.sequence) <> '' THEN 'Yes' ELSE 'No' END AS has_sequence,
    COUNT(DISTINCT pb.protein_id) AS target_count,
    COUNT(DISTINCT bd.disease_id) AS disease_tag_count,
    COUNT(DISTINCT bt.trial_id) AS clinical_trial_count,
    COUNT(DISTINCT bs.structure_id) AS direct_structure_count,
    CASE WHEN EXISTS (
        SELECT 1
        FROM protein_binders pb2
        JOIN protein_structures ps2 ON ps2.protein_id = pb2.protein_id
        WHERE pb2.binder_id = b.binder_id
    ) THEN 'Yes' ELSE 'No' END AS has_target_structure_fallback,
    b.data_completeness_score
FROM binders b
LEFT JOIN protein_binders pb ON pb.binder_id = b.binder_id
LEFT JOIN binder_diseases bd ON bd.binder_id = b.binder_id
LEFT JOIN binder_trials bt ON bt.binder_id = b.binder_id
LEFT JOIN binder_structures bs ON bs.binder_id = b.binder_id
GROUP BY b.binder_id, b.binder_name, b.binder_type, b.clinical_status, b.sequence, b.data_completeness_score;
