import pandas as pd

logfile = pd.read_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/all_logfiles/logfile.other_lab_mixed.out', sep='\t')
blacklist = pd.read_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/all_blacklist/other_lab.blacklist.tsv', sep='\t')

retained_files = {}
for i in range(len(logfile)):
	
	flag = logfile.iloc[i, 0]
	filepath = logfile.iloc[i, 2]
	chosen_filepath = logfile.iloc[i, 3]

	if flag == 'file_size':
		retained_files[filepath] = flag
	elif flag == 'duplicate':
		retained_files[chosen_filepath] = flag
	elif flag == 'stop_word_in_path':
		retained_files[filepath] = flag

for key, value in retained_files.items():
	key_short = key.split('/')[-1]

	for i in range(len(blacklist)):
		accession = blacklist.iloc[i, 0]

		if accession == key_short:
			
			blacklist.iloc[i, 1] = value

blacklist.to_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/temp_blacklist/other_lab.blacklist.tsv', sep='\t', index=False)