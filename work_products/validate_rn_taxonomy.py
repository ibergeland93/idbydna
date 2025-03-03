import os
import sys
import argparse
import glob
import pandas as pd

from datetime import datetime

from idbd_bio_utils import NcbiTaxonomy

#bio-utils functions
sys.path.append('/data/taxonomer2/ibergeland_work/cloned_repos/idbd-bio-utils/idbd_bio_utils/functions/')

from bio_utils_functions import load_taxid_repname_dict


def import_NCBI_taxonomy():
	"""
	Imports NBCI taxonomy class from the idbd_bio_utils repo and stores the variables: 
	merges, nodes, names, & ncbi.

	"""
    
	ncbi_dir = "/data/analysis_group2/ncbi_tax/2021_12_23"
	merged = os.path.join(ncbi_dir, "merged.dmp")
	nodes = os.path.join(ncbi_dir, "nodes.dmp")
	names = os.path.join(ncbi_dir, "names.dmp")
	ntax = NcbiTaxonomy(merged = merged, nodes = nodes, names = names)

	return ntax

def parse_args():
	"""
	Calls arguments used at the command line.
	
	Args:
		required:
			- input (str): Path to input file.
		optional:
			- initialize_release (str): Use flag to initialize a new database release (ex. --initialize_release 5.8.1)

	Returns:
		- args (str): arg parser object.

	"""

	parser = argparse.ArgumentParser()
	parser.add_argument("--db_version",
						type=str,
						help="use this flag with db release version to initialize a new Explify Release into the DataVault")
	parser.add_argument("--rn_path",
						type=str,
						help="path to reporting names table")
	args = parser.parse_args()

	return args


def initialize_explify_release_properties(db_version):
	"""
	Args:
		- update_type (str): Specifies which DataVault property is being updated.
		- reporting_names_filename: The reporting names table filename (ie. explify_reporting_name_info_table.txt)
		- product_name (str): The corresponding profile / use-case (ie. UPIP, RPIP, PAMPLONA).
		- database_release_version (str): The database version that was utilized to process data (ex 5.8.1).
		- ncbi_dir (str): The ncbi taxonomy db release version used for this release (ex. 2021_05_21).
		- test_profile_filename (str): The name of text file storing profile information (ex. ExplifyRespEnrichmentPanel.txt)

	Returns:
		- explify_release_properties (dict): Dictionary object of database release folder resources.

	"""

	db_release_folder = '/data/analysis_group1/explify-classification-database/'
	release_dir_full_path = glob.glob(db_release_folder + '*/*/explify_classification_database_release-*{}'.format(db_version), recursive=True)[0]

	reporting_names_filename = 'explify_reporting_name_info_table.txt'
	product_name = release_dir_full_path.split('/')[4]
	ncbi_dir = glob.glob(release_dir_full_path + '/*/*.dmp', recursive=True)[0].split('/')[-2]
	test_profile_filename = glob.glob(release_dir_full_path + '/Explify*.txt', recursive=True)[0].split('/')[-1]
	
	explify_release_properties = {'reporting_names_table_path': release_dir_full_path + '/' + reporting_names_filename,
								  'product_name':product_name,
								  'database_release_version':db_version,
								  'ncbi_dir':ncbi_dir,
								  'test_profile_path': release_dir_full_path + '/' + test_profile_filename}

def load_repname_table(rn_table):
	"""
	Loads Reporting Names Table through provided path.

	Args:
		- rn_table (str): Path to Reporting Names Table.

	Returns:
		- repnames (class): Uses the idbio_utils ReportingNames method; 
				Class to fetch information from a reporting names info table for a given organism. 
				Especially useful for converting between reporting names, taxids, and reporting IDs.
		- repnmaes_df (dataframe): Loads Reporting Names Table as a DataFrame; Useful when NCBI_taxonomy class is used.

	"""

	repnames = load_taxid_repname_dict(rn_table, key='reporting_name')

	repnames_df = pd.read_csv(rn_table, sep="\t")

	return repnames, repnames_df

def check_unique_rn(repnames_df, logfile):
	"""
	Checks that a given reporting name is not duplicated in the Reporting Names Table.

	Args:
		- logfile (output file): Writes to reporting_name_error_log.txt if errors are found.
		- repnames_df (dataframe): Loads Reporting Names Table as a DataFrame; Useful when NCBI_taxonomy class is used.

	"""

	reporting_names = set()

	for repname in repnames_df['reporting_name']:
		if repname not in reporting_names:
			reporting_names.add(repname)
		else:
			error = 'Error: Duplicate reporting name {}'.format(repname)
			logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))

def check_unique_tx(repnames_df, logfile):
	"""
	Checks that a given taxid is not duplicated in the Reporting Names Table.

	Args:
		- logfile (output file): Writes to reporting_name_error_log.txt if errors are found.
		- repnmaes_df (dataframe): Loads Reporting Names Table as a DataFrame; Useful when NCBI_taxonomy class is used.

	"""

	taxids = set()

	for tx in repnames_df['taxids']:
		tx = tx.split(',')
		for i in range(len(tx)):
			if tx[i] not in taxids:
				taxids.add(tx[i])
			else:
				error = 'Error: Duplicate taxid {}'.format(tx[i])
				logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))

def check_unique_repid_and_int(repnames_df, logfile):
	"""
	Checks that a given reporting id is not duplicated in the Reporting Names Table.

	Args:
		- logfile (output file): Writes to reporting_name_error_log.txt if errors are found.
		- repnmaes_df (dataframe): Loads Reporting Names Table as a DataFrame; Useful when NCBI_taxonomy class is used.

	"""

	repids = set()

	for repid in repnames_df['reporting_id']:
		if repid not in repids:
			repids.add(repid)
		else:
			error = 'Error: Duplicate reporting id {}'.format(repid)
			logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))
	
	for repid in repnames_df['reporting_id']:
		try:
			repid = repid + 1
		except TypeError:
			error = 'TypeError: Repid ({}) is not type int.'.format(repid)
			logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))
		

def check_unique_compoundid(repnames_df, logfile):
	"""
	Checks that a given compoundID is not duplicated in the Reporting Names Table.

	Args:
		- logfile (output file): Writes to reporting_name_error_log.txt if errors are found.
		- repnmaes_df (dataframe): Loads Reporting Names Table as a DataFrame; Useful when NCBI_taxonomy class is used.

	"""

	compound_ids = set()

	for compound_id in repnames_df['compound_id']:
		if compound_id not in compound_ids:
			compound_ids.add(compound_id)
		else:
			error = 'Error: Duplicate compound id {}'.format(compound_id)
			logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))

def check_tx_in_compound_id(repnames_df, logfile, ntax):
	"""
	Checks that a compound IDs taxid component is found in the listed taxids for a given row in the Reporting Names Table.

	Args:
		- logfile (output file): Writes to reporting_name_error_log.txt if errors are found.
		- repnmaes_df (dataframe): Loads Reporting Names Table as a DataFrame; Useful when NCBI_taxonomy class is used.

	"""

	for i in range(len(repnames_df)):
		
		repname = repnames_df['reporting_name'][i]
		compound_id = repnames_df['compound_id'][i]

		repid = int(compound_id.split('_')[0])
		tx = int(compound_id.split('_')[1])
		updated_tx = int(ntax.get_updated_taxid(tx))

		if tx == updated_tx:
			pass
		else:
			updated_compound_id = '{}_{}'.format(repid, updated_tx)
			error = 'Error: {} CompoundID {} is outdated.. Updated taxid = {}. Updated CompoundID is: {}'.format(repname, compound_id, updated_tx, updated_compound_id)
			logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))

def verify_tx_maps_to_ncbi_taxonomy(repnames_df, logfile, ntax):
	"""
	Verifies that a given taxid maps to the most recent NCBI taxonomy.

	Args:
		- logfile (output file): Writes to reporting_name_error_log.txt if errors are found.
		- repnmaes_df (dataframe): Loads Reporting Names Table as a DataFrame; Useful when NCBI_taxonomy class is used.
		- ntax (class): NCBI taxonomy class from idbd_bio_utils

	"""

	for i in range(len(repnames_df)):
		
		taxids = repnames_df['taxids'][i].split(',')
		repname = repnames_df['reporting_name'][i]

		for tx in taxids:
			tx = int(tx)

			updated_tx = int(ntax.get_updated_taxid(tx))

			if tx != updated_tx:
				error = 'Error: {} -- taxonomy not up to date. tx {} != updated_tx {}'.format(repname, tx, updated_tx)
				logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))

def check_tabulation(repnames_df, logfile):
	"""
	Checks that there are no extra tabs in the Reporting Names Table for a given row.

	Args:
		- logfile (output file): Writes to reporting_name_error_log.txt if errors are found.
		- repnmaes_df (dataframe): Loads Reporting Names Table as a DataFrame; Useful when NCBI_taxonomy class is used.

	"""

	columns_length = len(repnames_df.columns)

	for i in range(len(repnames_df)):
		row_length = len(repnames_df.iloc[i, :])

		try:
			if row_length == columns_length:
				pass
		except pd.ParserError:
			error = 'Error: Row {} is not equal to column length'.format(i)
			logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error))
		
		for j in range(len(repnames_df.columns)):
			column = repnames_df.columns[j]
			reporting_name = repnames_df.iloc[i,0]
			if column == 'parent_reporting_name':
				value = repnames_df.iloc[i, j]
				if value != value:
					error = 'Error: Row {}; {} has semantic parent_reporting_names value {}. {}'.format(i, reporting_name, value)
					logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))
			if column == 'children_reporting_names':
				value = repnames_df.iloc[i, j]
				if value != value:
					error = 'Error: Row {}; {} has children_reporting_names value {}.{}'.format(i, reporting_name, value)
					logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))
#					repnames_df.iloc[i, j] = '.'
			if column == 'semantic_group':
				value = repnames_df.iloc[i, j]
				if value != value:
					error = 'Error: Row {}; {} has semantic group value {}.{}'.format(i, reporting_name, value)
					logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))
#					repnames_df.iloc[i, j] = '.'

def check_class_subclass_nucleic_acid(repnames_df, logfile):

	subclass_viral = ['plant_virus', 'protist_virus', 'protist_virus', 'fungal_virus', 'virophage', 'viral', 'protist_virus', 'fungal_virus', 'plant_virus', 'endogenous_virus', 'phage' , 'armored_rna']
	subclass_fungal = ['fungal']
	subclass_parasite = ['parasite']
	subclass_bacterial = ['bacterial']

	nucleic_acid_types = ['unclear', 'RNA', 'retro', 'DNA']

	for i in range(len(repnames_df)):

		class_type = repnames_df.iloc[i, repnames_df.columns.get_loc('class_type')]
		subclass = repnames_df.iloc[i, repnames_df.columns.get_loc('subclass')].split(',')[0]
		nucleic_acid = repnames_df.iloc[i, repnames_df.columns.get_loc('nucleic_acid')]
		
		if class_type == 'bacterial':
			if subclass not in subclass_bacterial:
				error = 'Error (row {}): Missing bacterial subclass; added to reporting names table.'.format(i)
				logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))
				repnames_df.iloc[i, repnames_df.columns.get_loc('subclass')] = 'bacterial'
		elif class_type == 'fungal':
			if subclass not in subclass_fungal:
				error = 'Error (row {}): Missing fungal subclass; added to reporting names table.'.format(i)
				logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))
				repnames_df.iloc[i, repnames_df.columns.get_loc('subclass')] = 'fungal'
		elif class_type == 'parasite':
			if subclass not in subclass_parasite:
				error = 'Error (row {}): Missing parasite subclass; added to reporting names table.'.format(i)
				logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))
				repnames_df.iloc[i, repnames_df.columns.get_loc('subclass')] = 'parasite'
		elif class_type == 'viral':
			if subclass not in subclass_viral:
				error = 'Error (row {}): Missing viral subclass; added to reporting names table.'.format(i)
				logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))
				repnames_df.iloc[i, repnames_df.columns.get_loc('subclass')] = subclass
			if nucleic_acid not in nucleic_acid_types:
				error = 'Error (row {}): Missing nucleic acid type.'.format(i)
				logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', '\t', error, '\n'))


def log_parent_child_semantic(repnames_df, logfile):
	
	for i in range(len(repnames_df)):

		repname = repnames_df.iloc[i, repnames_df.columns.get_loc('reporting_name')]
		parent_reporting_name = repnames_df.iloc[i, repnames_df.columns.get_loc('parent_reporting_name')]
		children_reporting_names = repnames_df.iloc[i, repnames_df.columns.get_loc('children_reporting_names')]
		semantic_group = repnames_df.iloc[i, repnames_df.columns.get_loc('semantic_group')]
		rname_sgroup_match_flag = '.'

		if parent_reporting_name != '.':
			logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format(repname, '\t', parent_reporting_name, '\t', children_reporting_names, '\t', semantic_group, '\t', rname_sgroup_match_flag, '\t', '\t', '\n'))
		elif children_reporting_names != '.':
			logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format(repname, '\t', parent_reporting_name, '\t', children_reporting_names, '\t', semantic_group, '\t', rname_sgroup_match_flag, '\t', '\t', '\n'))
		elif semantic_group != '.':
			if semantic_group == repname:
				rname_sgroup_match_flag = 'T'
				logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format(repname, '\t', parent_reporting_name, '\t', children_reporting_names, '\t', semantic_group, '\t', rname_sgroup_match_flag, '\t', '\t', '\n'))
			else:
				rname_sgroup_match_flag = 'F'
				logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format(repname, '\t', parent_reporting_name, '\t', children_reporting_names, '\t', semantic_group, '\t', rname_sgroup_match_flag, '\t', '\t', '\n'))


def sort_taxid_by_relevance(repnames_df, ntax):
	
	with open('/data/taxonomer2/ibergeland_work/cloned_repos/virtual-environments/idbd-bio-utils/idbd_bio_utils/scripts/validation_scripts/reporting_names_table/output_logs/reporting_name_tx_rank.txt', 'a') as logfile:
		
		logfile.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format('flag', 'reporting_name', 'class_type', 'flag_description', 'species_tx', 'non_species_tx', 'chosen_tx', 'chosen_rank', 'trim_first_20_taxids', '\n'))

		for i in range(len(repnames_df)):

			#gathers necessary data from reporting names table
			taxids = repnames_df.iloc[i, 2].split(',')
			reporting_name = repnames_df.iloc[i, 0]
			class_type = repnames_df.iloc[i, 4]

			#splits taxids into two lists -- species or no_species
			species = []
			no_species = []
			for tx in taxids:
				tx = int(tx)
				
				rank = ntax.get_rank(tx)
				
				if rank == 'species':
					species.append(tx)
				else:
					no_species.append(tx)
			
			#ranks and reorganizes taxids in the species list
			multi_species_flag = 0
			if len(species) != 0:
				multi_species_flag = 1

				if len(species) > 1:
					multi_species_flag = 2

					rank_sort_dict = {}
					for tx in species:
						lineage = ntax.get_lineage_lists(tx)
						len_lineage = len(lineage[0])
						
						rank_sort_dict[tx] = len_lineage

					sorted_ranks = sorted(rank_sort_dict.items(), key=lambda x: x[1], reverse=False)

					sorted_ranks_and_taxids = sorted(sorted_ranks, key=lambda item: (item[1], item[0]))

					new_species_taxid_order = []
					
					for j in range(len(sorted_ranks_and_taxids)):
						tx = sorted_ranks_and_taxids[j][0]
						new_species_taxid_order.append(str(tx))
				
				if len(new_species_taxid_order) > 20:
					readable_new_species_taxid_order = new_species_taxid_order[0:21]
				else:
					readable_new_species_taxid_order = new_species_taxid_order

			#ranks and reorganizes taxids in the no_species list
			multi_no_species_flag = 0
			if len(no_species) != 0:
				multi_no_species_flag = 1

				if len(no_species) > 1:
					multi_no_species_flag = 2

					rank_sort_dict = {}
					for tx in no_species:
						lineage = ntax.get_lineage_lists(tx)

						len_lineage = len(lineage[0])
						
						rank_sort_dict[tx] = len_lineage

					sorted_ranks = sorted(rank_sort_dict.items(), key=lambda x: x[1], reverse=False)

					sorted_ranks_and_taxids = sorted(sorted_ranks, key=lambda item: (item[1], item[0]))

					new_no_species_taxid_order = []
					
					for k in range(len(sorted_ranks_and_taxids)):
						tx = sorted_ranks_and_taxids[k][0]
						new_no_species_taxid_order.append(str(tx))
					
				if len(new_species_taxid_order) > 20:
					readable_new_species_taxid_order = new_species_taxid_order[0:21]
				else:
					readable_new_species_taxid_order = new_species_taxid_order


			#writes logfile output based on conditions
			#if multiple species taxids
			if multi_species_flag == 2:

				#and if multiple "non-species" taxids
				if multi_no_species_flag == 2:
					chosen_tx = new_ordered_taxids[0]
					chosen_rank = ntax.get_rank(int(chosen_tx))
					first_twenty = new_ordered_taxids[0:21]
					first_twenty_str = ','.join(first_twenty)

					error = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('MULTI-SPECIES + MULTI-NON-SPECIES', reporting_name, class_type,
															'mulitple species taxids & multiple non-species taxids exist',
															 readable_new_species_taxid_order, no_species, chosen_tx, chosen_rank, 
															 first_twenty_str)
					logfile.write(error)
				
				elif multi_no_species_flag == 1:
					new_ordered_taxids = new_species_taxid_order + no_species
					error = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('MULTI-SPECIES + NON-SPECIES', reporting_name, class_type,
															'mulitple species taxids & a singular non-species taxids exist',
															 readable_new_species_taxid_order, no_species, new_ordered_taxids)
					logfile.write(error)
					
				#if only species, replace with ordered list
				else:
					new_ordered_taxids = new_species_taxid_order
					error = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('MULTI-SPECIES', reporting_name, class_type,
															'multiple species taxids exist',
															 readable_new_species_taxid_order, '.', new_ordered_taxids)
					logfile.write(error)

			#if a single species taxid exits
			elif multi_species_flag == 1:
				
				#and if multiple "non-species" taxids exist, combine single species taxid with non-species ordered list
				if multi_no_species_flag == 2:
					new_ordered_taxids = species + new_no_species_taxid_order
					error = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('SPECIES + MULTI-NON-SPECIES', reporting_name, class_type,
															'a singular species taxid, and multiple non-species taxids exist',
															 readable_new_species_taxid_order, no_species, new_ordered_taxids)
					logfile.write(error)

				#and a single "non-species" taxid exists, combine single species taxid with single non-species
				elif multi_no_species_flag == 1:
					new_ordered_taxids = species + no_species
					error = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('SPECIES + NON-SPECIES', reporting_name, class_type,
															'a singular species taxid, and a singular non-species taxids exist',
															 readable_new_species_taxid_order, no_species, new_ordered_taxids)
					logfile.write(error)
				
				#if only single species
				else:
					new_ordered_taxids = species
					error = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('SPECIES', reporting_name, class_type,
															'a singular species taxid exists', readable_new_species_taxid_order, 
															'.', new_ordered_taxids)
			
			#if no species taxids, and only singular non-species taxid
			elif multi_no_species_flag == 1:
				new_ordered_taxids = no_species
				error = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('NO SPECIES', reporting_name, class_type,
															'a singular, non-species taxid exists', '.',
															 no_species, new_ordered_taxids)
			
			#if no species taxids, and multiple non-species taxids; replace with ordered list
			else:
				new_ordered_taxids = new_no_species_taxid_order
				error = '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format('MULTI-NON-SPECIES', reporting_name, class_type,
															'multiple non-species taxids exist', '.',
															 no_species, new_ordered_taxids)
				logfile.write(error)

			for m in range(len(new_ordered_taxids)):
				new_ordered_taxids[m] = str(new_ordered_taxids[m])

			new_taxids = ','.join(new_ordered_taxids)
			
			repnames_df.iloc[i, 2] = new_taxids

	logfile.close()

	repnames_df.to_csv('/data/taxonomer2/ibergeland_work/cloned_repos/virtual-environments/idbd-bio-utils/idbd_bio_utils/scripts/validation_scripts/reporting_names_table/output_logs/explify_reporting_name_info_table_taxids_reordered.txt', sep='\t', index=False)


def main():

	#import ncbi taxonomy
	ntax = import_NCBI_taxonomy()

	#gather required command line parameters
	args = parse_args()

	#gather optional command line parameters
	try:
		db_version = args.db_version
	except AttributeError:
		db_version = None

	try:	
		rn_path = args.rn_path
	except AttributeError:
		rn_path = None

	date = datetime.now().strftime("%Y_%m_%d")

	with open('/data/taxonomer2/ibergeland_work/cloned_repos/idbd-bio-utils/idbd_bio_utils/scripts/validation_scripts/reporting_names_table/output_logs/{}_reporting_name_error_log.txt'.format(date), 'a') as logfile:
#		logfile.write('{}{}{}{}{}{}{}{}{}{}{}{}'.format('reporting_name', '\t','parent_reporting_name', '\t', 'children_reporting_names', '\t', 'semantic_group', '\t', 'rname_sgroup_match_flag', '\t', 'flag_other', '\n'))
		if db_version is not None:
			explify_release_properties = initialize_explify_release_properties(db_version)


		repnames_dict, repnames_df = load_repname_table(rn_path)

		#perform Reporting Name Table checks
#		check_unique_rn(repnames_df, logfile)
#		check_unique_tx(repnames_df, logfile)
#		check_unique_repid_and_int(repnames_df, logfile)
#		check_unique_compoundid(repnames_df, logfile)
#		check_tx_in_compound_id(repnames_df, logfile, ntax)
#		verify_tx_maps_to_ncbi_taxonomy(repnames_df, logfile, ntax)
#		check_tabulation(repnames_df, logfile)
#		check_class_subclass_nucleic_acid(repnames_df, logfile)
#		log_parent_child_semantic(repnames_df, logfile)

#		logfile.close()

		#reorganize Reporting Name Table
		sort_taxid_by_relevance(repnames_df, ntax)


if __name__ == "__main__":
	"""
	main function that directs flow of code execution

	"""
    
	main()