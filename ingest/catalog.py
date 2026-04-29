GENES_TO_INGEST = [
    "EGFR", "ERBB2", "KRAS", "ALK", "PIK3CA", "MET", "ROS1", "RET", "ESR1", "BRCA1", "BRCA2",

    "TP53", "BRAF", "NRAS", "PTEN", "AKT1", "MTOR", "CDK4", "CDK6", "MDM2",
    "MYC", "NTRK1", "NTRK2", "NTRK3", "FGFR1", "FGFR2", "FGFR3",
    "VEGFA", "PDGFRA", "KIT", "JAK2", "IDH1", "IDH2", "AR", "NBN",
]

DISEASES_TO_SEARCH = [
    "Lung Cancer",
    "Non-Small Cell Lung Cancer",
    "Breast Cancer",
    "Metastatic Breast Cancer",
    "Triple Negative Breast Cancer",
    "Colorectal Cancer",
    "Melanoma",
    "Prostate Cancer",
    "Ovarian Cancer",
    "Glioblastoma",
]

TRIAL_TARGETS = [
    {"gene": gene, "disease": disease}
    for gene in GENES_TO_INGEST
    for disease in DISEASES_TO_SEARCH
]
