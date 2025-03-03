import pandas as pd

accessions_file = pd.read_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/panels/rpip/rpip_prod2/temp.txt', sep='\t')
file_paths = pd.read_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/panels/rpip/rpip_prod2/all_paths.txt', sep='\t')
rp_prod_cp = pd.read_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/all_cp_tsvs/rpip_prod.cp.tsv', sep='\t')


#for i in range(len(file_paths)):

#	path = file_paths.iloc[i, 0]
#	file = file_paths.iloc[i, 0].split('/')[-1].split('.')[0]

#	file_paths.iloc[i,1] = file_paths.iloc[i, 0].split('/')[-1].split('.')[0]

#new_df = pd.merge(accessions_file, file_paths)
#print(new_df)


#new_df.to_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/scripts/map_basespace.tsv', sep='\t', index=False)

map_basespace = pd.read_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/scripts/map_basespace.tsv', sep='\t')

short_paths = []

for filepath in map_basespace['filepath']:

	short_path = '/'.join(filepath.split('/')[2:])
	short_paths.append(short_path)

map_basespace['short_paths'] = short_paths

short_paths = []
for filepath in rp_prod_cp['file_path']:

	short_path = filepath.split('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/panels/rpip/rpip_prod/')[1]
	short_paths.append(short_path)

rp_prod_cp['short_paths'] = short_paths

new_df = pd.merge(map_basespace, rp_prod_cp)
new_df = new_df.loc[:,['accession', 'filepath', 'file_path', 'short_paths']]

print(new_df)