-- =========================
-- COMPLETE BINDER DATA MODEL
-- =========================

-- Add missing columns to binders
ALTER TABLE binders
ADD COLUMN IF NOT EXISTS binder_type TEXT,
ADD COLUMN IF NOT EXISTS sequence TEXT,
ADD COLUMN IF NOT EXISTS clinical_status TEXT,
ADD COLUMN IF NOT EXISTS description TEXT;

-- =========================
-- Binder Structures
-- =========================
CREATE TABLE IF NOT EXISTS binder_structures (
    binder_structure_id SERIAL PRIMARY KEY,
    binder_id INTEGER REFERENCES binders(binder_id) ON DELETE CASCADE,
    pdb_id TEXT,
    structure_url TEXT,
    structure_source TEXT,
    notes TEXT
);

-- =========================
-- Binder Diseases
-- =========================
CREATE TABLE IF NOT EXISTS binder_diseases (
    binder_id INTEGER REFERENCES binders(binder_id) ON DELETE CASCADE,
    disease_id INTEGER REFERENCES diseases(disease_id) ON DELETE CASCADE,
    tag_reason TEXT,
    PRIMARY KEY (binder_id, disease_id)
);

-- =========================
-- Fix binder_trials schema
-- =========================
ALTER TABLE binder_trials
ADD COLUMN IF NOT EXISTS relationship_type TEXT,
ADD COLUMN IF NOT EXISTS source_id INTEGER,
ADD COLUMN IF NOT EXISTS notes TEXT;