# OncoMatch - Optimizing Oncology Combination Therapy Prediction Through Genomic, Structural, and Network Analysis

# What's the problem?
### Despite advancements in precision medicine, identifying effective cancer treatments remains a challenge due to the complexity of linking genetic data to drug pathways. Prior research - https://osf.io/preprints/biohackrxiv/c5wtr_v1 - has demonstrated a method for identifying potential colorectal cancer drugs by integrating genetic data with drug effectiveness scoring. However, this approach has not been widely applied to other cancers or optimized for clinical decision-making.

# How are we solving the problem?
### We aim to solve this problem by expanding the existing methodology to other cancers with similar classification systems, including bladder, ovarian, and small-cell lung cancer, thereby broadening its impact. Additionally, we plan to enhance diagnostics by leveraging Drug Central to map current treatments to genetic cancer data, and developing an evidence-based scoring system to assess the effectiveness of drug combinations in clinical settings. By advancing these efforts, we seek to improve precision medicine by making cancer treatment more personalized and data-driven.

# Workflow for OncoMatch
### OncoMatch is a web app that matches the cancer and gene mutation type to the most effective drug therapy based on the activity score of the drug against the specific gene target. We utilize the COSMIC database - https://cancer.sanger.ac.uk/cosmic/browse/tissue?wgs=off&sn=ovary&ss=all&hn=all&sh=&in=t&src=tissue&all_data=n - to identify the gene mutations for specific cancer types and the Drug Central database - https://drugcentral.org/ - to identify the drug therapies and the specific genes they target.
<img width="813" alt="Screen Shot 2025-04-02 at 11 58 29 AM" src="https://github.com/user-attachments/assets/c6ca0c50-5a05-448a-9b79-5112c674cdb4" />

 
