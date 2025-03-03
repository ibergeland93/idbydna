# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import glob
import argparse
import os

def parse_args():
	"""
	Calls arguments used at the command line.
	
	Args:
		required:
			- input_dir (str): Path to input directory.
			- prefix (str): type of data being soured
	Returns:
		- args (str): arg parser object.
	"""

	parser = argparse.ArgumentParser()
	parser.add_argument("-input_dir",
						type=str,
						help="input is either .txt or .xlsx",
						required=True)
	parser.add_argument("--compare_dir",
						type=str,
						help="input is either .txt or .xlsx",
						required=False)
	args = parser.parse_args()


	return args

def check_duplicates(input_dir):

	all_files_list = glob.glob("{}/*/*gz".format(input_dir))

	store_files_set = set()

	with open('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/to_remove/duplicates.logfile.tsv', 'a') as logfile:
		logfile.write("{}{}{}{}{}{}".format("older_file",'\t',"newer_file", '\t', "duplicate_in_compared_dir",'\n'))
		for file_long in all_files_list:

			file_short = file_long.split('/')[-1]

			if file_short in store_files_set:
				
				duplicates = glob.glob('{}/**/{}'.format(input_dir,file_short), recursive=True)
				f1 = duplicates[0]
				f2 = duplicates[1]

				logfile.write("{}{}{}{}{}{}".format(f1,'\t',f2, '\t', '','\n'))

			else:
				store_files_set.add(file_short)
		
		logfile.close()

	return store_files_set, logfile

def check_between_folders(compare_dir, store_files_set, logfile):

	all_files_list = glob.glob("{}/**/*.fastq.gz".format(compare_dir), recursive=True)

	with open('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/to_remove/duplicates.logfile.tsv', 'a') as logfile:

		for file_long in all_files_list:

			file_short = file_long.split('/')[-1]

			if file_short in store_files_set:
				pass
			else:
				
				logfile.write("{}{}{}{}{}{}".format(file_short,'\t','', '\t', file_long,'\n'))
		
		logfile.close()


def main():
	"""
	Main function -- gathers input, maps, produces json output.
	"""

	#gather required command line parameters
	args = parse_args()
	input_dir = args.input_dir
	compare_dir = args.compare_dir

	store_files_set, logfile = check_duplicates(input_dir)

	if compare_dir == compare_dir:
		check_between_folders(compare_dir, store_files_set, logfile)

if __name__ == "__main__":
    """
    main function that directs flow of code execution
    """

    main() 