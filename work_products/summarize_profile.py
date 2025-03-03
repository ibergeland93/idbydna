import os
import argparse
import glob
import json
import pandas as pd

from idbd_bio_utils import NcbiTaxonomy


def parse_args():
	"""
	Calls arguments used at the command line.
	
	Args:
		required:
			- db_version (str): Database release version (ex. 5.8.1) & will map all information.

	Returns:
		- args (str): arg parser object.

	"""

	parser = argparse.ArgumentParser()
	parser.add_argument("--db_version",
						type=str,
						help="use this flag with db release version to initialize a new Explify Release into the DataVault")
	args = parser.parse_args()

	return args
	

def import_NCBI_taxonomy():
	"""
	Imports NBCI taxonomy class from the idbd_bio_utils repo and stores the variables: 
	merges, nodes, names, & ncbi.

	"""
    
	ncbi_dir = "/data/analysis_group2/ncbi_tax/2021_12_23/"
	merged = os.path.join(ncbi_dir, "merged.dmp")
	nodes = os.path.join(ncbi_dir, "nodes.dmp")
	names = os.path.join(ncbi_dir, "names.dmp")
	ntax = NcbiTaxonomy(merged = merged, nodes = nodes, names = names)

	return ntax


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


	return explify_release_properties

def parse_organism_profile(org_profile, product_name, db_version, ntax):

	json_objs = []
	bacterial = []
	viral = []
	fungal = []
	parasite = []

	with open(org_profile, 'r', encoding='utf-16') as org_profile:
		with open('/data/taxonomer2/ibergeland_work/cloned_repos/virtual-environments/idbd-bio-utils/idbd_bio_utils/scripts/summarization_scripts/org_profile/profile_summary_output/{}_{}_organism_summary.tsv'.format(product_name, db_version), 'a') as outfile:
			outfile.write('{}\t{}\t{}\t{}\n'.format('reporting_name', 'genus_name', 'org_class', 'rep_id'))
			for json_obj in org_profile.readlines():
				json_objs.append(json_obj)

			for i in range(len(json_objs)):
				line = json.loads(json_objs[i])
				
				repname = line['reporting_name']
				rep_id = line['reporting_id']

				org_class = line['class_type']
				if org_class == 'bacterial':
					bacterial.append(repname)
				elif org_class == 'viral':
					viral.append(repname)
				elif org_class == 'fungal':
					fungal.append(repname)
				elif org_class == 'parasite':
					parasite.append(repname)

				taxids = line['taxids']
				for tx in taxids:
					genus_tx = ntax.get_genus_taxid(tx)
					genus_name = ntax.get_name(genus_tx)
			
				outfile.write('{}\t{}\t{}\t{}\n'.format(repname, genus_name, org_class, rep_id))
			
			with open('/data/taxonomer2/ibergeland_work/cloned_repos/virtual-environments/idbd-bio-utils/idbd_bio_utils/scripts/summarization_scripts/org_profile/profile_summary_output/{}_{}_class_counts.tsv'.format(product_name, db_version), 'a') as class_counts:
				class_counts.write('bacterial count: {} \n viral count: {} \n fungal count: {} \n parasite_count: {} \n'.format(
				len(bacterial), len(viral), len(fungal), len(parasite)))


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

	explify_release_properties = initialize_explify_release_properties(db_version)

	org_profile = explify_release_properties['test_profile_path']
	product_name = explify_release_properties['product_name']
	db_version = explify_release_properties['database_release_version']

	parse_organism_profile(org_profile, product_name, db_version, ntax)


if __name__ == "__main__":
	"""
	main function that directs flow of code execution

	"""
    
	main()