ALTER TABLE binders
ADD COLUMN IF NOT EXISTS binder_type TEXT,
ADD COLUMN IF NOT EXISTS sequence TEXT,
ADD COLUMN IF NOT EXISTS clinical_status TEXT;

-- Backfill binder_type from your existing modality table
UPDATE binders b
SET binder_type = bm.modality_name
FROM binder_modalities bm
WHERE b.modality_id = bm.modality_id
  AND (b.binder_type IS NULL OR b.binder_type = '');

-- Backfill clinical_status from the existing approval_status field
UPDATE binders
SET clinical_status = approval_status
WHERE (clinical_status IS NULL OR clinical_status = '')
  AND approval_status IS NOT NULL;