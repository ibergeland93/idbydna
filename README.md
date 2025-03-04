# idbydna
This repo contains projects I worked on during my tenure at IDbyDNA (now Illumina).
All required data files for reproduction can be found in this G-Drive: https://drive.google.com/drive/folders/1ekrhX3TqLBGDJtwYqVPqI34E7dlnLH1W?usp=drive_link


## Project 1: Data Vault
The bioinformatics team was creating a "Data Vault" of approved clinical samples that could be used for product testing. The clinical samples were obtained from CLIA labs and were sequenced and analyzed in the IDbyDNA wet-lab or by other bioinformatics teams using the IDbyDNA analysis software ("Explify") which detected microorganisms (bacteria, virus, fungus, parasites) in a sample. Following a sequencing run, samples would be output as files for post-processing; however, only samples that met specific criteria could be included for testing.

### create_cp.tsv.py: 
is a python script that creates a resulting dataframe with the following columns: ['accession', 'file_path', 'taxids_expected',	'repids_expected', 'repids_included', 'repids_excluded','kingdom', 'prefix', 'split']), with only samples that had microorganisms detected at a certain level.

### map_repid_to_cp.tsv
reads two tab-separated (.tsv) datasets containing analytical data related to product releases. It then compares the accession values from both datasets, and if a match is found, it updates the corresponding row in df with a specific value (repid) from product_release_df. Finally, it saves the modified df to a new .tsv file. 

### check_duplicates_in_subfolder.py:
is a python script that compares files across two folders and checks for duplicates.

### map_log_to_blacklist.py
is a python script that reads from a "blacklist" of key words that indicate if a sample should be excluded from a batch.

### validate_rn_taxonomy.py
This script validates and reorganizes the Reporting Names Table for an Explify classification database release. It ensures the accuracy of taxonomy mappings, detects duplicate or incorrect taxonomic IDs, verifies consistency with the NCBI taxonomy database, and logs any errors found. Finally, it reorders taxonomic IDs based on relevance to improve data organization.

## Project 2: Respiratory Pathogen Panel (RPIP)
BaseSpace Sequence Hub is a cloud environment for the analysis, storage, and sharing of genomic data. During the COVID-19 pandemic, the bioinformatics team was developing an enrichment test that could detect hundreds of respiratory pathogens. See more info here: https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub/apps/explify-rpip-data-analysis.html

### map_basespace_path.py
As part of testing, I developed a script to compare results of an internally developed analysis pipeline against a tool developed by Illumina (Basespace). This script maps file locations for the same files. I would provide results to the Dir. of Bioinformatics to evaluate.

### map_off-profile.py
This script maps off-profile organisms by analyzing taxonomic lineage data from the NCBI taxonomy database. It reads a list of organisms and their corresponding taxonomic IDs, determines their genus-level classification, and applies exceptions for specific families and orders. It then identifies and records organisms that are not in the expected taxonomic profile while filtering out parent-child relationships to refine the final output.

## Project 3: Open-Wound Pathogen Panel (Pamplona)
The Bioinformatics team was developing a pathogen test for detecting pathogenic micro-organisms in open wounds. We conducted our production testing in AWS and our R&D testing locally via HPC.

### md5_check.py
This script compares MD5 checksums of gzip-compressed (.gz) files from two directories (pamp_asimov and pamp_aws) to verify data integrity. It extracts filenames, computes their MD5 hashes, and stores them in dictionaries for comparison. If a file exists in both directories, the script checks whether their checksums match and logs the results in pamp_lab_source__md5.txt. Any mismatches or missing files are flagged for further investigation.

## Project 4: Product Testing
Produce a concise yet comprehensive overview of product performance for strategic planning and stakeholder updates.

### summarize_LOD.py
This script automates Limit of Detection (LoD) testing by processing sample metadata, detection results, and expected outcomes to validate product accuracy. It identifies false positives/negatives, corrects deprecated data, and calculates detection performance metrics. By generating a structured LoD report, it ensures regulatory compliance, supports quality control, and streamlines product validation.

### summarize_profile.py 
This script parses and summarizes taxonomic data from a database release, specifically for the Explify classification database (Explify is what IDbyDNA called their pathogen detection software). It extracts organism profiles, classifies them by type (bacterial, viral, fungal, parasite), and generates summary reports. The script is designed to validate and organize taxonomy information for a specific database version used in microbial classification.

### file_size_check.py
This script compares file sizes between two datasets to ensure data integrity. It reads two tab-separated (.tsv) files containing initial and controlled file sizes, matches file paths, and checks for discrepancies. If a difference is found, it prints the mismatched file paths and their sizes.















