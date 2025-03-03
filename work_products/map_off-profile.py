## run at /data/taxonomer2/ibergeland_work/cloned_repos/virtual-environments/idbd-bio-utils/map_off-profile.py

import pandas as pd

from idbd_bio_utils import NcbiTaxonomy

tax_dir = "/data/analysis_group2/ncbi_tax/2020_08_25"
ntax = NcbiTaxonomy(merged=f"{tax_dir}/merged.dmp", nodes=f"{tax_dir}/nodes.dmp", names=f"{tax_dir}/names.dmp")

input = '/data/taxonomer2/ibergeland_work/cloned_repos/virtual-environments/explify-config/test_profiles/resources/organism/respiratory/rpp/mr_evidence/manual_review/explify_rpp_final_mr_reporting_names.txt'
rnt = '/data/taxonomer2/ibergeland_work/cloned_repos/virtual-environments/explify-config/reporting_names/explify_reporting_name_info_table.txt'

with open(input) as file:
	lines = [line.rstrip() for line in file]

reporting_names_dict = {}
with open(rnt) as rnt:
	for line in rnt:
		line = line.split('\t')
		reporting_name = line[0]
		taxids = line[2].split(',')
		reporting_names_dict[reporting_name] = {'taxids':taxids}

#family level for enterobacteriaceae and actinomycetaceae, order level for mucorales
family_exceptions = ['Enterobacteriaceae', 'Actinomycetaceae']
order_exceptions = ['Mucorales']

repname_genus_tx = {}

for org in lines:
	taxids = reporting_names_dict[org]['taxids']
	for tx in taxids:
		genus_tx = ntax.get_genus_taxid(int(tx))
		try:
			genus_name = ntax.get_name(int(genus_tx))
		except:
			genus_name = 'None'
		
		lineage_lists = ntax.get_lineage_lists(int(tx))
		
		txs = lineage_lists[0]
		names_in_lineage = lineage_lists[1]
		rank =lineage_lists[2]
		
		for val_f in family_exceptions:
			for val_o in order_exceptions:
				if val_f in names_in_lineage:
					genus_name = val_f
				elif val_o in names_in_lineage:
					genus_name = val_o
				else:
					pass
			
		if genus_name is not 'None':
			index = names_in_lineage.index(genus_name)
			rank = rank[index]
			t = txs[index]
			repname_genus_tx[org] = [t, rank]


repname_off_profile_orgs = {}

for prof_org, genus_tx_rank in repname_genus_tx.items():
	off_profile_orgs = []

	for reporting_name, taxids in reporting_names_dict.items():
		taxids = reporting_names_dict[reporting_name]['taxids']

		genus_tx = genus_tx_rank[0]
		genus_rank = genus_tx_rank[1]

		#if reporting name not on profile
		if reporting_name not in repname_genus_tx.keys():
			
			if 'Influenza A' in prof_org and 'Influeza A' in reporting_name:
				pass
				off_profile_orgs.append(reporting_name)
				repname_off_profile_orgs[prof_org] = off_profile_orgs
				print(prof_org, repname_off_profile_orgs[prof_org])
			elif 'Influenza B' in prof_org and 'Influeza B' in reporting_name:
				pass
				off_profile_orgs.append(reporting_name)
				repname_off_profile_orgs[prof_org] = off_profile_orgs
				print(prof_org, repname_off_profile_orgs[prof_org])
			elif 'Influenza C' in prof_org and 'Influeza C' in reporting_name:
				pass
				off_profile_orgs.append(reporting_name)
				repname_off_profile_orgs[prof_org] = off_profile_orgs
				print(prof_org, repname_off_profile_orgs[prof_org])
			
			#for each taxid in a given reporting name
			parent_names = set()
			children_names = set()
			for tx in taxids:
				try:
					off_prof_lineage_lists = ntax.get_lineage_lists(int(tx))	
				
					taxids_in_parent_lineage = off_prof_lineage_lists[0]
					child_taxids = ntax.get_children(int(tx))

					#create parent child mappings for viruses
					if '10239' in taxids_in_parent_lineage:
						for parent_tx in taxids_in_parent_lineage:
							parent_name = ntax.get_name(parent_tx)
							parent_names.add(parent_name)
						for child_tx in child_taxids:
							child_name = ntax.get_name(child_tx)
							children_names.add(child_name)

					if genus_tx in off_prof_lineage_lists[0]:
						index = off_prof_lineage_lists[0].index(genus_tx)
						rank = off_prof_lineage_lists[2][index]

						if genus_rank == rank:
							if reporting_name not in off_profile_orgs:
								off_profile_orgs.append(reporting_name)
							else:
								pass
							repname_off_profile_orgs[prof_org] = off_profile_orgs
				except:
					off_prof_lineage_lists = 'NA'
					pass
			
			#exclude halo orgs that are parent / children of repname
			for halo_org in off_profile_orgs:
				if halo_org in parent_names:
					off_profile_orgs.remove(halo_org)
					repname_off_profile_orgs[prof_org] = off_profile_orgs
				elif halo_org in children_names:
					off_profile_orgs.remove(halo_org)
					repname_off_profile_orgs[prof_org] = off_profile_orgs


with open('map-off-profile.txt', 'a') as f:
	for key, value in repname_off_profile_orgs.items():
		print('{}\t{}\n'.format(key, ','.join(value)))
		f.write('\n{}\t{}\n'.format(key, ','.join(value)))
