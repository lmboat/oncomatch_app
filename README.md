# OncoMATCH - Optimizing Oncology Combination Therapy Prediction Through Genomic, Structural, and Network Analysis

# What's the problem?
### Despite advancements in precision medicine, identifying effective cancer treatments remains a challenge due to the complexity of linking genetic data to drug pathways. [Prior research](https://osf.io/preprints/biohackrxiv/c5wtr_v1) - has demonstrated a method for identifying potential colorectal cancer drugs by integrating genetic data with drug effectiveness scoring. However, this approach has not been widely applied to other cancers or optimized for clinical decision-making.

# How are we solving the problem?
### We aim to solve this problem by expanding the existing methodology to other cancers with similar classification systems, including bladder, ovarian, and small-cell lung cancer, thereby broadening its impact. Additionally, we plan to enhance diagnostics by leveraging Drug Central to map current treatments to genetic cancer data, and developing an evidence-based scoring system to assess the effectiveness of drug combinations in clinical settings. By advancing these efforts, we seek to improve precision medicine by making cancer treatment more personalized and data-driven.

# Workflow for OncoMatch
### OncoMatch is a web app that matches the cancer and gene mutation type to the most effective drug therapy based on the activity score of the drug against the specific gene target. We utilize the [COSMIC database](https://cancer.sanger.ac.uk/cosmic/browse/tissue?wgs=off&sn=ovary&ss=all&hn=all&sh=&in=t&src=tissue&all_data=n) - to identify the gene mutations for specific cancer types and the [Drug Central database](https://drugcentral.org/) - to identify the drug therapies and the specific genes they target.
<img width="807" alt="Screen Shot 2025-04-02 at 2 22 11 PM" src="https://github.com/user-attachments/assets/4a595176-0375-4500-9a64-d740717ec0cb" />

# Datasets for OncoMatch
### [Cosmic_CancerGeneCensus_v99_GRCh38.tsv](data/Cosmic/Cosmic_CancerGeneCensus_v99_GRCh38.tsv)
This data comes from a Gene Census located in the Cosmic database. Columns: GENE_SYMBOL, NAME, COSMIC_GENE_ID, CHROMOSOME, GENOME_START, GENOME_STOP, CHR_BAND, SOMATIC, GERMLINE, TUMOUR_TYPES_SOMATIC, TUMOUR_TYPES_GERMLINE, CANCER_SYNDROME, TISSUE_TYPE, MOLECULAR_GENETICS, ROLE_IN_CANCER, MUTATION_TYPES, TRANSLOCATION_PARTNER, OTHER_GERMLINE_MUT, OTHER_SYNDROME, TIER, SYNONYMS

### [EMA_Approved.csv](data/DrugCentral/raw_data/EMA_Approved.csv)
This data comes from Drug Central, containing a single column of drug names which have been approved by the (European Medicines Agency) EMA

### [FDA_Approved.csv](data/DrugCentral/raw_data/FDA_Approved.csv)
This data comes from Drug Central, containing a single column of drug names which have been approved by the (Food and Drug Administration) FDA

### [PDMA_Approved.csv](data/DrugCentral/raw_data/PMDA_Approved.csv)
This data comes from Drug Central, containing a single column of drug names which have been approved by the (Product Development and Management Association) PDMA

### [drug.target.interaction.tsv](data/DrugCentral/raw_data/drug.target.interaction.tsv)
This data comes from Drug Central, containing many columns including drug name, the gene it targets, and efficacy metrics/values
Columns: DRUG_NAME, STRUCT_ID, TARGET_NAME, TARGET_CLASS, ACCESSION, GENE, SWISSPROT, ACT_VALUE, ACT_UNIT, ACT_TYPE, ACT_COMMENT, ACT_SOURCE, RELATION, MOA, MOA_SOURCE, ACT_SOURCE_URL, MOA_SOURCE_URL, ACTION_TYPE, TRL, ORGANISM

### [structures.smiles.tsv](data/DrugCentral/raw_data/structures.smiles.tsv)
This data comes from Drug Central, containing many columns including SMILE data and ID
Columns: SMILES, InChl, InChlKey, ID, INN, CAS_RN

### [drug.structure.id.csv](data/DrugCentral/processed_data/drug.structure.id.csv)
This data comes from Drug Central, containing 2 columns mapping drug name to struct id
Columns: DRUG_NAME, STRUCT_ID

### [drug.target.interaction.fda.cosmic.cancer.type](data/DrugCentral/processed_data/drug.target.interaction.fda.cosmic.cancer.type.tsv)
This data comes from Drug Central, containing many columns including drug name, the gene it targets, efficacy metrics/values, and cancer type
Columns: DRUG_NAME, STRUCT_ID, TARGET_NAME, TARGET_CLASS, ACCESSION, GENE, SWISSPROT, ACT_VALUE, ACT_UNIT, ACT_TYPE, ACT_COMMENT, ACT_SOURCE, RELATION, MOA, MOA_SOURCE, ACT_SOURCE_URL, MOA_SOURCE_URL, ACTION_TYPE, TDL, ORGANISM, Driver_Gene, FDA_Approved, EMA_Apprived, PMDA_Approved, Bladder, Colon, NSCLC, Ovarian

### [drug.target.interaction.scores.tsv](data/DrugCentral/processed_data/drug.target.interaction.scores.tsv)
This data comes from Drug Central, containing many columns including drug name gene targeted, and precision score for each cancer type
Columns: DRUG_NAME, nGENE, BLADDER_CPS, BLADDER_CGPS, COLON_CPS, COLON_CGPS, NSCLC_CPS, NSCLC_CGPS, OVARIAN_CPS, OVARIAN_CGPS

### [gene.id.csv](data/DrugCentral/processed_data/gene.id.csv)
This data comes from Drug Central, containing 2 columns, mapping gene name to cancer type
Columns: GENE, CANCER

 # Sample Radar Plots for OncoMatch
![image](https://github.com/user-attachments/assets/7fb4a0ac-dce8-4675-9d52-820486dad6da)
