import pandas as pd
import os
import glob
import gzip
from pathlib import Path

from pandas.io.parsers import read_csv


def get_readcount(filename):

	with gzip.open(filename, "r") as f:
		lines = f.readlines()  # read file as bytes
 
		readcount = 0
		for line in lines:
			readcount += 1

		readcount = readcount / 4

		return readcount

def compare_readcounts():

#	initial = pd.read_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/pamp_lab_initial_readcount.txt', sep='\t')
	initial = pd.read_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/scripts/output/pamp_lab_initial_file_sizes.txt', sep='\t')
	controlled = pd.read_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/scripts/output/pamp_lab_controlled_file_sizes.txt', sep='\t')
#	controlled = pd.read_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/pamp_lab_controlled_readcount.txt', sep='\t')

#	initial_dict = initial.set_index('file')['read_count'].to_dict()
#	controlled_dict = controlled.set_index('file')['read_count'].to_dict()

	initial_dict = initial.set_index('path')['file_size'].to_dict()
	controlled_dict = controlled.set_index('path')['file_size'].to_dict()

#	for file, read_count in initial_dict.items():
	for file, size in initial_dict.items():

#		for file_c, read_count_c in controlled_dict.items():
		for file_c, size_c in controlled_dict.items():

			if file == file_c:
				if initial_dict[file] == controlled_dict[file_c]:
				#	print(file, read_count)
					pass
				else:
						print(file, initial_dict[file], controlled_dict[file_c])

compare_readcounts()
#pamp = "/data/analysis_group2/data_vault/datasets/analytical/pamplona/enriched/aws/release_2_reprocess/*/postqual_fastqs/*gz"

#targets = []

#with open("pamp_lab_initial_readcount.txt", "a+") as o:
#	for target in glob.glob(pamp):

#		targets.append(target)

#	print('Targets finished.')
#	for target in targets:
#		readcount = get_readcount("{}".format(target))
#		o.write("{}\t{}\n".format(target, int(readcount)))

#	o.close()

