import pandas as pd


product_release_df = pd.read_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/scripts/pamp_lab_product_release_cp.tsv', sep='\t')
df = pd.read_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/all_cp_tsvs/pamp_lab.cp.tsv', sep='\t')

accession = product_release_df.loc[: , 'accession']
accession_df = df.loc[:, 'accession']

for i in range(len(accession_df)):
	for j in range(len(accession)):
		if accession_df[i] == accession[j]:
			repid = product_release_df.iloc[j, 3]
			df.iloc[i, 3] = repid

df.to_csv('/data/analysis_group2/data_vault/datasets/org_challenge_datasets/alloid_test_data/cd5_all_analytical/pamp_lab_mapped_cp.tsv', sep='\t')