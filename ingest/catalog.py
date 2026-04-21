# ingest/catalog.py

GENES_TO_INGEST = [
    "EGFR",
    "ERBB2",
    "KRAS",
    "ALK",
    "PIK3CA",
    "MET",
    "ROS1",
    "RET",
    "ESR1",
    "BRCA1",
    "BRCA2",
]

# Keep this mostly aligned to your sponsor focus on lung and breast cancer.
TRIAL_TARGETS = [
    {"gene": "EGFR", "disease": "Lung Cancer"},
    {"gene": "EGFR", "disease": "Non-Small Cell Lung Cancer"},
    {"gene": "ERBB2", "disease": "Breast Cancer"},
    {"gene": "ERBB2", "disease": "HER2-Positive Breast Cancer"},
    {"gene": "KRAS", "disease": "Lung Cancer"},
    {"gene": "KRAS", "disease": "Non-Small Cell Lung Cancer"},
    {"gene": "ALK", "disease": "Lung Cancer"},
    {"gene": "ALK", "disease": "Non-Small Cell Lung Cancer"},
    {"gene": "PIK3CA", "disease": "Breast Cancer"},
    {"gene": "PIK3CA", "disease": "Metastatic Breast Cancer"},
    {"gene": "MET", "disease": "Lung Cancer"},
    {"gene": "MET", "disease": "Non-Small Cell Lung Cancer"},
    {"gene": "ROS1", "disease": "Lung Cancer"},
    {"gene": "ROS1", "disease": "Non-Small Cell Lung Cancer"},
    {"gene": "RET", "disease": "Lung Cancer"},
    {"gene": "RET", "disease": "Non-Small Cell Lung Cancer"},
    {"gene": "ESR1", "disease": "Breast Cancer"},
    {"gene": "ESR1", "disease": "Metastatic Breast Cancer"},
    {"gene": "BRCA1", "disease": "Breast Cancer"},
    {"gene": "BRCA1", "disease": "Triple Negative Breast Cancer"},
    {"gene": "BRCA2", "disease": "Breast Cancer"},
    {"gene": "BRCA2", "disease": "Triple Negative Breast Cancer"},
]