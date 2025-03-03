# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os

import pandas as pd
import re
import datetime

from idbd_bio_utils import ReportingNames


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
	parser.add_argument("-prefix",
						type=str,
						help="input is type of data being sourced",
						required=True)
	parser.add_argument("--min_file_size",
						type=int,
						help="minimum threshold for file size. will issue a warning if less than given value.")
	args = parser.parse_args()


	return args

def get_accession(file):

	nuc = re.search('\-d\-[0-9]+|\-r\-[0-9]+',file)
	try:
		split_name = file.split(nuc[0])
		accession = split_name[0] + str(nuc[0])
	except:
		accession = file.split('.')[0]
	return accession

def find_latest_from_hash(current,compared):

	# Get time stamps
	current_timestamp = None

	try:
		nuc_seq_ind_current = re.search('[CAGT]+\-[CAGT]+?\(-|_)',current).span()[-1]
		current_timestamp = current[nuc_seq_ind_current:].split('_')[0] # in hex hash
	except:
		pass
	if current_timestamp == None:
		try:
			nuc_seq_ind_current = re.search(r'[0-9]{4}-[0-9]{2}-[0-9]{2}', current)
			current_timestamp = nuc_seq_ind_current[0]
		except:
			pass

	compared_timestamp = None
	try:
		nuc_seq_ind_compared = re.search('[CAGT]+\-[CAGT]+?\(-|_)',compared).span()[-1]	
		compared_timestamp = compared[nuc_seq_ind_compared:].split('_')[0] # in hex hash
	except:
		pass
	if compared_timestamp == None:
		try:
			nuc_seq_ind_compared = re.search(r'[0-9]{4}-[0-9]{2}-[0-9]{2}', compared)
			compared_timestamp = nuc_seq_ind_compared[0]
		except:
			pass

	if current_timestamp is not None and compared_timestamp is not None:
#		try: 
#			current_timestamp = datetime.fromtimestamp(int(current_timestamp,16))
#		except:
#			print(f"ERROR: current_timestamp not converted {current_timestamp}, from {nuc_seq_ind_current}, from {current}")
#		try: 
#			compared_timestamp = datetime.fromtimestamp(int(compared_timestamp,16))
#		except:
#			print(f"ERROR: compared_timestamp not converted {compared_timestamp}, from {nuc_seq_ind_compared}, from {compared}")

		if current_timestamp < compared_timestamp:
			return compared, current_timestamp, compared_timestamp
		else:
			return current, current_timestamp, compared_timestamp

	else:
		if 'post' in current.lower() and 'post' not in compared.lower():
			return current
		elif 'post' in compared.lower() and 'post' not in current.lower():
			return compared
		else:
			return current

def do_not_transfer(accession, acc_path_dict, filepath, logfile, min_file_size):
    """
    Specific rules to not add to cp.tsv
    """
    # check if it's a symbolic link (will break if this is not done first)
    if os.path.islink(filepath):
        if not os.path.exists(filepath):
            logfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('sym_link_broken', 'Symbolic link does not exist.', filepath, '.', '.', '.', '.'))
            return True

    # check if it's a small (or empty file)
    file_size = float(os.stat(filepath).st_size)
    if file_size<min_file_size:
        logfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('file_size', 'File size too small ({}) -- skip file.'.format(file_size), filepath, '.', '.', '.', '.'))
        return True
        
    # check if any of the stop words are in the file name    
    bad_fastqs = {'downsample','trim','blk','poscon','negcon', 'demo'}
    for i in bad_fastqs:
        if i.lower() in filepath.lower():
            logfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('stop_word_in_path', 'Stop word {} recognized -- skip file.'.format(i), filepath, '.', '.', '.', '.'))
            return True

def get_files(input_dir, prefix, repnames, logfile, min_file_size):

	list_of_files = list()
	for (dirpath, dirnames, filenames) in os.walk(input_dir):
		list_of_files += [os.path.join(dirpath, file) for file in filenames if 'fastq.gz' or '.fasta' or '.fa' or '.fq' in file]
	
	rows = []
	acc_path_dict = {}
	for filepath in list_of_files:
		file = filepath.split('/')[-1]
		accession = get_accession(file)

		if accession in acc_path_dict.keys():
			current = filepath
			compared = acc_path_dict[accession]

			try:
				new, current_timestamp, compared_timestamp = find_latest_from_hash(current,compared)
			except:
				new = find_latest_from_hash(current,compared)
				current_timestamp = 'NA'
				compared_timestamp = 'NA'
			if new == filepath:
				logfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('duplicate', 'Incoming filepath replaces existing', '.', new, acc_path_dict[accession], current_timestamp, compared_timestamp))
				acc_path_dict[accession] = new
			elif new == acc_path_dict[accession]:
				logfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('duplicate', 'Existing filepath retained', '.', new, filepath, compared_timestamp, current_timestamp,))
				acc_path_dict[accession] = new
		
		else:
			if do_not_transfer(accession, acc_path_dict, filepath, logfile, min_file_size):
				pass	
			else:
				try:
					acc_path_dict[accession] = os.readlink(filepath)
				except:
					acc_path_dict[accession] = filepath


	for accession, filepath in acc_path_dict.items():
		short_accession = '-'.join(accession.split('-')[-4:-2])
		row = [short_accession, filepath, '', '', '', '', '', prefix, 2]
		rows.append(row)
	df = pd.DataFrame(rows, columns = ['accession', 'file_path', 'taxids_expected',
										'repids_expected', 'repids_included', 'repids_excluded',
										'kingdom', 'prefix', 'split'])
	
	df.to_csv('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/all_cp_tsvs/{}.cp.tsv'.format(prefix), sep='\t', index=False)

	return(df)	


def main():
	"""
	Main function -- gathers input, maps, produces json output.
	"""

	#gather required command line parameters
	args = parse_args()
	input_dir = args.input_dir
	prefix = args.prefix
	min_file_size = args.min_file_size

	repnames = ReportingNames("/data/taxonomer2/ibergeland_work/cloned_repos/explify-config/reporting_names/explify_reporting_name_info_table.txt")
	print(repnames)
	with open('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/all_logfiles/logfile.{}.out'.format(prefix), 'a') as logfile:
		
		logfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('flag', 'notes', 'filepath', 'chosen_filepath', 'excluded_filepath', 'chosen_timestamp', 'excluded_timestamp'))
		get_files(input_dir, prefix, repnames, logfile, min_file_size)
	
	logfile.close()


if __name__ == "__main__":
    """
    main function that directs flow of code execution
    """

    main() 