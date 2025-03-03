import argparse
import numpy as np
import pandas as pd


def _parse_args():
	"""
	Gathers arguments from command line input.

	"""

	parser = argparse.ArgumentParser(description="Create formatted LOD table from detections file")
	parser.add_argument(
		"--data_transfer",
		required=True,
		help="Path to Data Transfer file",
	)
	parser.add_argument(
		"--rnt",
		required=True,
		help="Path to reporting names file",
	)
	parser.add_argument(
		"--detections",
		required=True,
		help="Path to detections file",
	)
	parser.add_argument(
		"--genotype",
		required=False,
		default=False,
		action="store_true",
		help="Use flag if genotype reporting performance by genotype is desired",
	)
	parser.add_argument(
		"--outfile",
		required=False,
		default='LOD_table.xlsx',
		help="Path to output file",
	)

	return parser.parse_args()


def build_samples_meta(samples):
	"""
	Builds dictionary that stores metadata from the Samples dataframe.
	Adds accessions (keys) to dictionary only if they are relevant to LoD experiments.
	- Input (dataframe == samples): A spreadsheet from the Data Transfer Template that contains sample metadata.
	- Output (dictionary == metadata): A dictionary with accessions as keys and sample metadata as values.

	"""
	
	print(" Building a dictionary that stores metadata from samples dataframe ... ", flush=True)
	metadata = {}

	print(" >>> Adding accessions that meet LoD conditions ...")
	for idx, row in samples.iterrows():
		
		expt = row['Experiment Name']
		
		if ('lod' in expt.lower() or 'titration' in expt.lower()) and (not 'checks' in expt):
#		if 'rpip illumina contrived' in expt.lower():
			acc = row['Accession Enrichment']
			if acc in metadata:
				print(f'Accession {acc} is duplicated in Data Transfer sheet')
			metadata[acc] = {}
			try:
				for f in ['Accession Enrichment',
						'Batch ID Enrichment',
						'Run ID Enrichment',
						'Sample Name',
						'Replicate',
						'BFx Notes',
						'Sample Type',
						'Internal Control Used']:
					metadata[acc][f] = row[f]
				metadata[acc]['Experiment Name'] = expt
			except:
				for f in ['Accession Enrichment',
						'Batch ID Enrichment',
						'Run ID Enrichment',
						'Sample Name',
						'Replicate',
						'BFx Notes',
						'Sample Type',
						'Internal Control Used-1',
						'Internal Control Used-2']:
					metadata[acc][f] = row[f]
				metadata[acc]['Experiment Name'] = expt

	return metadata


def add_micro_meta_to_existing(micro, metadata):
	"""
	Adds microorganism metadata from existing samples metadata dictionary.
	If more than one org is expected in an accession, it will be gathered in a single string object per accession (sep= ';').
	- Input:
		- (dataframe == micro): A spreadsheet from the Data Transfer Template that contains microorganism metadata.
		- (dictionary == metadata): A dictionary with accessions as keys and sample metadata as values; Deprecated accessions replaced
	- Output:
		- (dictionary == metadata): Same dictionary as input; but with a mechanism for capturing multiple expected orgs in the same sample;
												In the case of multiple expecteds, a final list of each is appended as the dictionary value for
												both Quant Value and Microorganism Name.

	"""
	
	for idx, row in micro.iterrows():

		acc = row['Accession Enrichment']

		if acc in metadata:
			try:
				old_quant = metadata[acc]['Quant Value']
				new_quant = row['Quant Value']
				metadata[acc]['Quant Value'] = '{};{}'.format(old_quant, new_quant)
			except:
				metadata[acc]['Quant Value'] = row['Quant Value']

			try:
				old_micro = metadata[acc]['Microorganism Name']
				new_micro = row['Microorganism Name']
				metadata[acc]['Microorganism Name'] = '{};{}'.format(old_micro, new_micro)
			except:
				metadata[acc]['Microorganism Name'] = row['Microorganism Name']
	
	for accession, meta in metadata.items():

		try:
			quant_vals = meta['Quant Value'].split(';')
			exptd_orgs = meta['Microorganism Name'].split(';')
		except:

			try:
				quant_vals = meta['Quant Value']
				exptd_orgs = meta['Microorganism Name']
			except:
				acc = accession.split(' ')[0]

		metadata[accession]['Quant Value'] = quant_vals
		metadata[accession]['Microorganism Name'] = exptd_orgs
	
	return metadata


def replace_dep_accs(samples, metadata):
	"""
	Builds a dictionary that stores deprecated accessions and their replacement mapping.
	Builds another dictionary from the samples metadata dictionary and replaces outdated accessions.
	- Input:
		- (dataframe == samples): A spreadsheet from the Data Transfer Template that contains sample metadata.
		- (dictionary == metadata): A dictionary with accessions as keys and sample metadata as values.
	- Output: (dictionary == metadata_replaced): Same dictionary as input, but with deprecated accessions replaced

	"""
	
	print(" Building a dictionary that stores deprecated accessions and their replacement values from samples dataframe ... ", flush=True)
	deprecated = {}  # deprecated: replacement
	for idx, row in samples.iterrows():
		notes = row['BFx Notes']
		if isinstance(notes, str) and 'deprecated' in notes.lower():
			replacement_acc = 'IDBD-' + notes.split('IDBD-')[-1]
			deprecated_acc = row['Accession Enrichment']
			deprecated[deprecated_acc] = replacement_acc
		else:
			continue

	print(" Replacing deprecated accessions with new accessions in metadata dictionary ...", flush=True)
	metadata_replaced = {}
	for accession, data in metadata.items():
		if accession in deprecated.keys():
			metadata[accession]['Accession Enrichment'] = deprecated[accession] + ' -- (replaces {})'.format(accession)
			metadata_replaced[deprecated[accession] + ' -- (replaces {})'.format(accession)] = metadata[accession]
		else:
			metadata_replaced[accession] = metadata[accession]
	
	return metadata_replaced


def map_long_to_short_accession(detections, metadata_replaced):

	selected_accessions = set(metadata_replaced.keys())
	detections_accessions = list(set(detections['accession']))

	accession_converter = {}

	for acc in selected_accessions:
		try:
			index = [detections_accessions.index(i) for i in detections_accessions if acc in i]
			long_acc = detections_accessions[index[0]]
			accession_converter[acc] = long_acc
		except:
			metadata_replaced.pop(acc)
	
	return accession_converter


def add_detections_meta_to_existing(detections, rnt, metadata_replaced, genotype, accession_converter):
	"""
	Creates dictionary of reporting names to map parent / children relationships
	Adds detection data to existing accessions within the metadata dictionary.
	Appends cases where detected name is in list of expecteds.
	Appends cases where DT Microorganism Name & Class Type are NoneType ... (Indicates PC or NC)
	- Input:
		- (dataframe == detections): A spreadsheet containing detection / performance metadata for processed data
		- (text file == rnt): A text file containing reporting information
		- (dictionary == metadata_replaced): A dictionary with accessions as keys and sample metadata as values; 
				Deprecated accessions replaced; metadata for multiple expectations
		- (genotype == argument): command line argument indicating if genotype-level summarization is desired
	- Output: (dictionary == metadata_replaced): Same dictionary as input; but with each detection object
				added to a list of detection objects for a given accession.

	"""

	print(" Adding detection data to existing accessions within the metadata dictionary ...", flush=True)

	reporting_names_dict = {}

	with open(rnt) as rnt:
		for line in rnt:
			line = line.split('\t')
			reporting_name = line[0]
			parent = line[7]
			children = line[8]
			reporting_names_dict[reporting_name] = {'parent':parent, 'children':children}
	
	#create dictionary that maps accessions to batchIDs
	acc_batchID = {}
	for accession, obj in metadata_replaced.items():
		batch_id = obj['Batch ID Enrichment']

		if '--' in accession:
			acc = accession.split()[0]
		else:
			acc = accession
		try:
			acc = accession_converter[acc]
		except:
			pass

		acc_batchID[acc] = batch_id

	detections['Batch ID Enrichment'] = detections['accession'].map(acc_batchID)
	## add part mapping long form to short form

	#create dictionary that maps batch expectation by batchID
	batch_id_expt = {}
	for b in acc_batchID.values():
		detections_subset = detections.loc[detections['Batch ID Enrichment'] == b]
		expected_names = set(detections_subset.loc[:, 'expected_names'])
		
		expected_names_list = []
		for exp in expected_names:
			if exp == exp:
				if ';' in exp:
					exp = exp.split(';')
					for cntrl in exp:
						expected_names_list.append(cntrl)
				else:
					expected_names_list.append(exp)
				batch_id_expt[b] = expected_names_list
	
	store_sample_read_count = {}
	for accession, obj in metadata_replaced.items():
		batch_id = obj['Batch ID Enrichment']
		micro_name = obj['Microorganism Name']

		if '--' in accession:
			acc = accession.split()[0]
		else:
			acc = accession
		try:
			acc = accession_converter[acc]
		except:
			pass

		detections_subset = detections.loc[detections['accession'] == acc]

		if genotype:
			detections_subset = detections_subset.loc[:, ['expected_names', 'name', 'detected_on_profile', 'expected', 
									'predicted', 'sample_read_count', 'coverage', 'rpkm', 'absolute_quant_ratio',
									'tiered_logic_result', 'read_count', 'read_count_ratio', 'union_result', 'strict_result', 'simple_result']]
		elif 'read_count' in detections_subset.columns and 'read_count_ratio' in detections_subset.columns:
			detections_subset = detections_subset.loc[:, ['expected_names', 'name', 'detected_on_profile', 'expected', 
									'predicted', 'sample_read_count', 'coverage', 'rpkm', 'absolute_quant_ratio',
									'tiered_logic_result', 'read_count', 'read_count_ratio', 'simple_result']]
		else:
			detections_subset = detections_subset.loc[:, ['expected_names', 'name', 'detected_on_profile', 'expected', 
									'predicted', 'sample_read_count', 'coverage', 'rpkm', 'absolute_quant_ratio',
									'tiered_logic_result', 'simple_result']]
			detections_subset['union_result'] = detections_subset['simple_result']
			detections_subset['strict_result'] = detections_subset['simple_result']

		parent_rows = []
		genotype_rows = []
		other_rows = []
		acc_level_spl_read_count = set()

		#checking all detections for a given accession
		for i in range(len(detections_subset)):
			meta = detections_subset.iloc[i, :]
			expected_names = detections_subset.iloc[i, 0]
			detected_name = detections_subset.iloc[i, 1]
			sample_read_count = detections_subset.iloc[i, 5]

			#if genotype, collect simple_result, union_result, strict_result (only union and strict are used)
			if genotype:	
				simple_result = detections_subset.iloc[i, 14]
				union_result = detections_subset.iloc[i, 12]
				strict_result = detections_subset.iloc[i, 13]
			#else, use only simple result by setting union and strict equal to simple
			else:
				simple_result = detections_subset.iloc[i, 12]
				union_result = simple_result
				strict_result = simple_result

			#storing accession-level sample count (same for all detections) in case NA
			acc_level_spl_read_count.add(sample_read_count)

			if sample_read_count == 'None':
				acc_level_spl_read_count.remove('None')
				try:
					sample_read_count = list(acc_level_spl_read_count)[0]
				except:
					pass
				meta['sample_read_count'] = sample_read_count
			store_sample_read_count[acc] = sample_read_count

			#homogenizing detected names in the kmer detections file
			if '(unable to genotype)' in detected_name or '(unable to genotype) ' in detected_name:
				detected_name = detected_name.split('(unable to genotype)')[0].strip()
				meta['name'] = detected_name

			#if an accession has no expectation (ie water or negative control), filter to correct bin (genotype, parent, other)
			if expected_names == expected_names:
				expected_names = expected_names.split(';')

				if union_result == 'TP' and strict_result == 'TP':
					if reporting_names_dict[detected_name]['parent'] == '.' and reporting_names_dict[detected_name]['children'] == '.':
						other_rows.append(dict(meta))
					elif reporting_names_dict[detected_name]['parent'] != '.':
						genotype_rows.append(dict(meta))
					else:
						other_rows.append(dict(meta))

				elif union_result == 'TP' and (strict_result != strict_result):
					if reporting_names_dict[detected_name]['parent'] == '.' and reporting_names_dict[detected_name]['children'] == '.':
						other_rows.append(dict(meta))
					elif reporting_names_dict[detected_name]['parent'] == '.' and reporting_names_dict[detected_name]['children'] != '.':
						for e in expected_names:
							if e == detected_name:
								other_rows.append(dict(meta))
						parent_rows.append(dict(meta))
					
				elif union_result == 'FN' or strict_result == 'FN':
					if reporting_names_dict[detected_name]['parent'] == '.' and reporting_names_dict[detected_name]['children'] == '.':
						if strict_result == 'None':
							pass
						else:
							other_rows.append(dict(meta))
					elif reporting_names_dict[detected_name]['parent'] == '.' and reporting_names_dict[detected_name]['children'] != '.':
						parent_rows.append(dict(meta))
						if strict_result == 'None':
							pass
						else:
							other_rows.append(dict(meta))
					elif reporting_names_dict[detected_name]['parent'] != '.':
						genotype_rows.append(dict(meta))

			else:
				#check for batch-level contamination
				if union_result == 'FP' and strict_result == 'FP':
					for exp in batch_id_expt[batch_id]:
						if exp != 'UNKNOWN':
							if reporting_names_dict[exp]['parent'] == detected_name:
								parent_rows.append(dict(meta))
							elif detected_name in reporting_names_dict[exp]['children']:
								genotype_rows.append(dict(meta))
							elif detected_name == exp:
								other_rows.append(dict(meta))
		
		#if genotype flag is set; include genotype and non-parent rows 
		if genotype:
			combo = genotype_rows + other_rows 
		#else only include parent rows and non-parent rows
		else:
			combo = parent_rows + other_rows

		metadata_replaced[accession]['detections_meta'] = combo

	return metadata_replaced, reporting_names_dict, store_sample_read_count


def add_sample_meta_to_all_detections(metadata_replaced, genotype, add_sample_meta_to_all_detections):
	"""
	Gathers detetction objects (representing a detected org within an accession) and updates each with identical sample information.
	Appends each dictionary to a list of metadata objects to be retained.
	Converts expected quant values to scientific notation.
	- Input (dictionary == metadata_replaced): A dictionary with sample metadata and nested dictionary containing detections metadata.
	- Output (dictionary == all_metadata_df): dataframe with all detections for a given sample reported individually as rows

	"""

	all_metadata = {}
	for accession, metadata in metadata_replaced.items():
		
		retain = []

		#separate sample metadata and detections metadata 
		sample_meta = {key:value for (key,value) in metadata.items() if key != 'detections_meta'}
		detections_meta = {key:value for (key,value) in metadata.items() if key == 'detections_meta'}

		#add missing sample read count based on what is seen for other detections for a given accession
		for meta in detections_meta.values():
			if len(meta) != 0:
				for i in range(len(meta)):
					combine_meta = {}
					combine_meta.update(sample_meta)
					combine_meta.update(meta[i])
					retain.append(combine_meta)
					store_columns = meta[i].keys()
				
			#if an accession is missing a detection for an expected detection, add row of "NA"s
			elif genotype:
				try:
					store_sample_read_count = add_sample_meta_to_all_detections[accession]
				except:
					store_sample_read_count = 'NA'
				if 'read_count' in store_columns and 'read_count_ratio' in store_columns:
					combine_meta = {'expected_names': 'NA', 'name': 'NA',
									'detected_on_profile': 'NA', 'expected': 'NA',
									'predicted': 'False', 'sample_read_count': store_sample_read_count, 'read_count': 'NA',
									'coverage': 'NA', 'rpkm': 'NA','absolute_quant_ratio': 'NA', 'read_count_ratio': 'NA',
									'tiered_logic_result': 'NA', 'union_result': 'ND', 'strict_result': 'ND'}
				else:
					combine_meta = {'expected_names': 'NA', 'name': 'NA',
									'detected_on_profile': 'NA', 'expected': 'NA',
									'predicted': 'False', 'sample_read_count': store_sample_read_count,
									'coverage': 'NA', 'rpkm': 'NA','absolute_quant_ratio': 'NA',
									'tiered_logic_result': 'NA', 'union_result': 'ND', 'strict_result': 'ND', 'simple_result': 'ND'}
				combine_meta.update(sample_meta)
				retain.append(combine_meta)
			else:
				try:
					store_sample_read_count = add_sample_meta_to_all_detections[accession]
				except:
					store_sample_read_count = 'NA'
				combine_meta = {'expected_names': 'NA', 'name': 'NA',
									'detected_on_profile': 'NA', 'expected': 'NA',
									'predicted': 'False', 'sample_read_count': store_sample_read_count,
									'coverage': 'NA', 'rpkm': 'NA','absolute_quant_ratio': 'NA',
									'tiered_logic_result': 'NA', 'union_result': 'ND', 'strict_result': 'ND', 'simple_result': 'ND'}
				combine_meta.update(sample_meta)
				retain.append(combine_meta)

		all_metadata[accession] = retain

	#flatten dictionary so that sample metadata is appended to each detection object
	retain_all = []
	for accession, det_objs in all_metadata.items():
		for det_obj in det_objs:
			det_name = det_obj['name']
			quant_val = det_obj['Quant Value']
			microorganisms = det_obj['Microorganism Name']
			if isinstance(quant_val, list):
				indices = [i for i, x in enumerate(quant_val) if x == 'nan']
				if len(indices) == 0:
					if det_name in ['Saccharomyces cerevisiae', 'Lactobacillus fermentum']:
						pass
					else:
						retain_all.append(det_obj)
				for i in indices:
					exp_name = microorganisms[i]
					if exp_name == det_name:
						pass
					else:
						retain_all.append(det_obj)
			else:
				retain_all.append(det_obj)
	
	all_metadata_df = pd.DataFrame(retain_all)

	return all_metadata_df


def add_rep(row):
	"""
	Applies a singular replicate value, 'A', when a blank value is encountered in the Data Transfer Template.
	- Input (dataframe == all_metadata_df): a dataframe consisting of sample, microorganism, and detections metadata
	- Output (dataframe == all_metadata_df): same as input dataframe, but with updated 'Replicates' values

	"""
	
	rep = row.Replicate
	if pd.isnull(rep):
		return 'A'
	else:
		return rep


def check_IC(row):
	"""
	Cross-references if internal control T7 was spiked in the sample. If no, the new value for absolute ratio will be changed to 'No IC'.
	- Input (dataframe == all_metadata_df): a dataframe consisting of sample, microorganism, and detections metadata
	- Output (dataframe == all_metadata_df): same as input dataframe, but with updated 'absolute_quant_ratio' values

	"""

	abs_quant = row.absolute_quant_ratio
	try:
		ic_val = row['Internal Control Used']

		if 't7' not in ic_val.lower():
			ic_val = 'no IC'
			return ic_val
		else:
			return abs_quant
	
	except:
		ic_val_1 = row['Internal Control Used-1']
		ic_val_2 = row['Internal Control Used-2']

		try:
			if 't7' not in ic_val_1.lower() or 'ms2' not in ic_val_2.lower():
				ic_val = 'no IC'
				return ic_val
			else:
				return abs_quant
		except:
			return abs_quant


def mark_expected_detected(row):
	"""
	Applies new column to dataframe indicating 1 if the outcome is a True Positive (TP) or 0 otherwise.
	- Input (dataframe == all_metadata_df): a dataframe consisting of sample, microorganism, and detections metadata; 
											replicates added; duplicates removed
	- Output (dataframe == all_metadata_df): same as input dataframe, but with outcome labeled

	"""
	simple_result = row['simple_result']
	
	try:
		union_result = row['union_result']
	except:
		union_result = simple_result
	try:
		strict_result = row['strict_result']
	except:
		strict_result = simple_result
	
	genotype = row['genotype']

	#if genotype flag is set at the command line, return detection based on strict and union results
	if genotype:
		if strict_result == 'TP' or strict_result == 'FP':
			return 1
		elif (strict_result != strict_result) and union_result == 'TP':
			return 1
		else: 
			return 0
	#if parent, only use union_result in determining positivity
	else:
		if union_result == 'TP' or union_result == 'FP':
			return 1
		else:
			return 0


def modify_exp_quant(row):
	"""
	Converts spike in quant value to scientific notation.
	- Input: (all_metadata_df['quant_val'])
	"""

	quant_val = row['Quant Value']
	name = row['name']
	dtt_name = row['Microorganism Name']
	sample_type = row['Sample Type']
	sample_name = row['Sample Name']

	#homogenize quant val
	if isinstance(quant_val, list):
		for i in range(len(quant_val)):
			if name in dtt_name[i]:
				qv = quant_val[i]
			elif name == 'NA':
				qv = quant_val[i]
				if qv == 'Undetected':
					qv = 'NA'
				if '<' in qv:
					qv = int(qv.split('<')[-1])

	elif quant_val != quant_val:
		if 'plasma' in sample_type.lower() or 'plasma' in sample_name.lower():
			qv = 'Plasma'
		elif 'urine' in sample_type.lower() or 'urine' in sample_name.lower():
			qv = 'Urine'
		elif 'human' in sample_type.lower() or 'human' in sample_name.lower():
			qv = 'Human'
		else:
			qv = 'Water'
	
	else:
		qv = quant_val

	return qv


def modify_abs_quant(row):
	"""
	Homogenizes absolute quantification value to desired format given input.
	- Input: (all_metadata_df['absolute_quant_ratio'])

	"""

	absolute_quant_ratio = row['absolute_quant_ratio']

	if isinstance(absolute_quant_ratio, str):
		if absolute_quant_ratio == absolute_quant_ratio:
			if absolute_quant_ratio == 'None':
				qv = 'NA'
			elif absolute_quant_ratio == 'NA':
				qv = 'NA'
			elif absolute_quant_ratio == 'no IC':
				qv = 'no IC'
			else:
				qv = '{:.2e}'.format(float(absolute_quant_ratio))
			return qv
	else:
		try:
			qv = '{:.2e}'.format(float(absolute_quant_ratio))
		except:
			qv = 'NA'
		return qv

def split_replicate(row):
	
	sample_name = row['Sample Name'].split(' ')
	replicates = ['a', 'b', 'c']

	for element in sample_name:
		if element.lower() in replicates:
			sample_name.remove(element)
	
	sample_name = ' '.join(sample_name)
	return sample_name
	

def format_lod(lod_table):
	"""
	Formats the final summary table by batch ID, rearranges columns, and converts values to desired dtypes.
	All objects are appended to a list, so that data will be displayed separately by batch.
	- Input (dataframe == lod_table): same dataframe as above (all_metadata_df), except pivoted by Replicate, 
										and a new column 'pos_ct' added (values = 1 or 0 depending on prediction).
	- Output (dataframe == lod_table): same as input dataframe, but with reformatted

	"""

	dfs = []

	for batch in lod_table['Batch ID Enrichment'].unique():

		tmp = lod_table[lod_table['Batch ID Enrichment'] == batch].copy()

		for i in range(len(tmp)):
			quant_val = tmp.iloc[i, 1]
			if isinstance(quant_val, str):
				try:
					quant_val = float(quant_val)
					tmp.iloc[i, 1] = quant_val
				except:
					pass

		# Rename Quant Value column with experiment details
		experiment = tmp['Experiment Name'].value_counts().idxmax()
		samp_type = tmp['Sample Type'].value_counts().idxmax()
		new_qv = experiment + ' in ' + samp_type
		tmp = tmp.rename(columns={'Quant Value': new_qv,
									'Batch ID Enrichment':'Batch ID',
									'pos_ct': 'No. of positive replicates'})
		tmp.columns = [' '.join(c) if isinstance(c, tuple) else c for c in tmp.columns]

		tmp = tmp.sort_values(['Sample Name'], ascending = [False])
		controls = ['pos', 'neg', 'blank']
		tmp_num = []
		tmp_str = []
		for i in range(len(tmp)):
			quant_val = tmp.iloc[i, 1]
			sample_name = tmp.iloc[i, 4]
			name = tmp.iloc[i, 6]

			for val in controls:
				if val in sample_name.lower():
					try:
						quant_val = str('{} in '.format(name) + '{:.2e}'.format(float(quant_val)))
						tmp.iloc[i, 1] = quant_val
						tmp_str.append(tmp.iloc[i, :])
						break
					except:
						pass

			if isinstance(quant_val, str):
				try:
					quant_val = str('{} in '.format(name) + '{:.2e}'.format(float(quant_val)))
					tmp.iloc[i, 1] = quant_val
					if quant_val == quant_val:
						tmp_str.append(tmp.iloc[i, :])
				except:
					tmp_str.append(tmp.iloc[i, :])
				
			else:
				quant_val = str('{} in '.format(name) + '{:.2e}'.format(float(quant_val)))
				tmp.iloc[i, 1] = quant_val
				tmp_num.append(tmp.iloc[i, :])

		tmp_str = pd.DataFrame(tmp_str)
		tmp_num = pd.DataFrame(tmp_num)

		try:
			tmp_str = tmp_str.sort_values(by=['accession A'], ascending = [True])
			tmp_num = tmp_num.sort_values(by=['Sample Name'], ascending = [False])
			tmp_final = pd.concat([tmp_num, tmp_str])
		except:
			tmp_final = tmp

		tmp_final = pd.concat([tmp_num, tmp_str])

		tmp_final = tmp_final.replace([None, '', 'None'], 'NA', regex=True)

		dfs.append(tmp_final)

	return dfs


def write_formatted_lod(dfs, outfile):
	"""
	Writes formatted lod table to Excel spreadsheet.
	- Input:
		- list of dataframes (dfs)
		- outfile path (outfile)

	"""

	pd_writer = pd.ExcelWriter(outfile, engine='xlsxwriter')

	if len(dfs) > 0:
		startrow = 0
		for df in dfs:
			df.to_excel(pd_writer, engine="xlsxwriter", startrow=startrow, sheet_name='formatted_lod', index=False)
			startrow += (df.shape[0] + 2)
			workbook = pd_writer.book
			worksheet = pd_writer.sheets['formatted_lod']

		merge_format = workbook.add_format({'align': 'center', 'valign': 'vcenter'})

		df_start = 1
		for df in dfs:
			df = df.reset_index(drop=True)
			df_end = df_start + df.shape[0] - 1
			worksheet.merge_range(df_start, 0, df_end, 0, df.loc[0, 'Batch ID'], merge_format)
			df_start = df_end + 3
	else:
		print('NO LOD TABLES PRESENT!')

	pd_writer.save()


def main():
	"""
	Main executables for producing the LoD performance summary output spreadsheet

	"""
	args = _parse_args()
	dt_file = args.data_transfer
	rnt = args.rnt
	detections_file = args.detections
	outfile = args.outfile
	genotype = args.genotype

	pd.set_option('display.max_rows', None)

	print(" Reading in samples as dataframe from {}".format(dt_file), flush=True)
	samples = pd.read_excel(open(dt_file, 'rb'), sheet_name='Samples', header=1)
	metadata = build_samples_meta(samples)

	print(" Reading in microorganisms as dataframe from {}".format(dt_file), flush=True)
	micro = pd.read_excel(open(dt_file, 'rb'), sheet_name='Microorganisms', header=1)
	metadata_replaced = add_micro_meta_to_existing(micro, metadata)

	print(" Replacing deprecated accessions within the metadata dictionary ... ", flush=True)
	metadata_replaced = replace_dep_accs(samples, metadata)

	print(" Switching attention to the k-mer detections file ... reading in detections as dataframe ... ***", flush=True)
	detections = pd.read_csv(detections_file, sep='\t', low_memory=False)
	
	accession_converter = map_long_to_short_accession(detections, metadata_replaced)

	metadata_replaced, reporting_names_dict, store_sample_read_count = add_detections_meta_to_existing(detections, rnt, metadata_replaced, genotype, accession_converter)
	all_metadata_df = add_sample_meta_to_all_detections(metadata_replaced, genotype, store_sample_read_count)
	
	print(" Adding in blank replicates with A ... (indicates single replicate)", flush=True)
	all_metadata_df['Replicate'] = all_metadata_df.apply(add_rep, axis=1)
	all_metadata_df['absolute_quant_ratio'] = all_metadata_df.apply(check_IC, axis=1)

	#set flag if genotype reporting is desired or not
	if genotype:
		all_metadata_df['genotype'] = 1
	else:
		all_metadata_df['genotype'] = 0

	print(" Adding a new column w value == 1 if TP or SUPERSEDED...", flush=True)
	all_metadata_df['pos_ct'] = all_metadata_df.apply(mark_expected_detected, axis=1)
	print(" Homogenizing the values in quant_val ... ", flush=True)
	all_metadata_df['quant_val'] = all_metadata_df.apply(modify_exp_quant, axis=1)
	all_metadata_df['absolute_quant_ratio'] = all_metadata_df.apply(modify_abs_quant, axis=1)
	all_metadata_df['Sample Name'] = all_metadata_df.apply(split_replicate, axis=1)

	all_metadata_df = all_metadata_df.rename(columns={'Accession Enrichment': 'accession'})

	all_metadata_df = all_metadata_df[all_metadata_df['quant_val'] != 'nan']
	all_metadata_df = all_metadata_df.drop(columns=['Quant Value', 'Microorganism Name', 'expected_names',
																'detected_on_profile', 'expected', 'predicted', 
																'Run ID Enrichment', 'BFx Notes', 'tiered_logic_result'], axis=1).reset_index(drop=True)
	all_metadata_df = all_metadata_df.drop_duplicates()

	if 'read_count' in all_metadata_df.columns and 'read_count_ratio' in all_metadata_df.columns:
		if genotype:
			wide_table = pd.pivot(all_metadata_df, columns=['Replicate'],
						values=['rpkm', 'coverage', 'union_result', 'strict_result', 'absolute_quant_ratio', 'sample_read_count', 'accession', 'name', 'read_count', 'read_count_ratio'],
                       	index=['Batch ID Enrichment', 'Experiment Name', 'Sample Name','Sample Type', 'quant_val', 'accession', 'name'])
			col_order = np.array(['Batch ID Enrichment', 'quant_val', 'Experiment Name', 'Sample Type', 'Sample Name', 'pos_ct', ('name', 'A'), ('name', 'B'), ('name', 'C'), ('union_result', 'A'), ('union_result', 'B'), ('union_result', 'C'), ('strict_result', 'A'), ('strict_result', 'B'), ('strict_result', 'C'),	('coverage', 'A'),('coverage', 'B'), ('coverage', 'C'), ('absolute_quant_ratio', 'A'), ('absolute_quant_ratio', 'B'), ('absolute_quant_ratio', 'C'), ('read_count_ratio', 'A'), ('read_count_ratio', 'B'), ('read_count_ratio', 'C'), ('sample_read_count', 'A'), ('sample_read_count', 'B'), ('sample_read_count', 'C'), ('read_count', 'A'), ('read_count', 'B'), ('read_count', 'C'), ('rpkm', 'A'), ('rpkm', 'B'), ('rpkm', 'C'), ('accession', 'A'), ('accession', 'B'), ('accession', 'C')], dtype=object)
		else:
			wide_table = pd.pivot(all_metadata_df, columns=['Replicate'],
						values=['rpkm', 'coverage', 'simple_result', 'absolute_quant_ratio', 'sample_read_count', 'accession', 'name', 'read_count', 'read_count_ratio'],
                    	index=['Batch ID Enrichment', 'Experiment Name', 'Sample Name','Sample Type', 'quant_val', 'accession', 'name'])
			col_order = np.array(['Batch ID Enrichment', 'quant_val', 'Experiment Name', 'Sample Type', 'Sample Name', 'pos_ct', ('name', 'A'), ('name', 'B'), ('name', 'C'), ('simple_result', 'A'), ('simple_result', 'B'), ('simple_result', 'C'), ('coverage', 'A'),('coverage', 'B'), ('coverage', 'C'),('absolute_quant_ratio', 'A'), ('absolute_quant_ratio', 'B'), ('absolute_quant_ratio', 'C'), ('read_count_ratio', 'A'), ('read_count_ratio', 'B'), ('read_count_ratio', 'C'), ('sample_read_count', 'A'), ('sample_read_count', 'B'), ('sample_read_count', 'C'), ('read_count', 'A'), ('read_count', 'B'), ('read_count', 'C'), ('rpkm', 'A'), ('rpkm', 'B'), ('rpkm', 'C'), ('accession', 'A'), ('accession', 'B'), ('accession', 'C')], dtype=object)
	else:
		wide_table = pd.pivot(all_metadata_df, columns=['Replicate'],
						values=['rpkm', 'coverage', 'simple_result', 'absolute_quant_ratio', 'sample_read_count', 'accession', 'name'],
                       	index=['Batch ID Enrichment', 'Experiment Name', 'Sample Name','Sample Type', 'quant_val', 'accession', 'name'])
		col_order = np.array(['Batch ID Enrichment', 'quant_val', 'Experiment Name',
							'Sample Type', 'Sample Name', 'pos_ct', ('name', 'A'), ('name', 'B'), ('name', 'C'), ('simple_result', 'A'),
							('simple_result', 'B'), ('simple_result', 'C'), ('coverage', 'A'), ('coverage', 'B'), ('coverage', 'C'),
							('absolute_quant_ratio', 'A'), ('absolute_quant_ratio', 'B'), ('absolute_quant_ratio', 'C'), ('sample_read_count', 'A'),
							('sample_read_count', 'B'), ('sample_read_count', 'C'), ('rpkm', 'A'), ('rpkm', 'B'), ('rpkm', 'C'),
							('accession', 'A'), ('accession', 'B'), ('accession', 'C')], dtype=object)
	wide_table.to_excel('/data/analysis_group2/data_vault/datasets/org_challenge_datasets/rpip/troubleshoot/wide_table.xlsx')
	#separate clinical from contrived to handle replicate collapse
	searchfor = ['clinical', 'Clinical']
	clin_metadata_df = wide_table[wide_table.index.get_level_values('Sample Type').str.contains('|'.join(searchfor))]
	if not clin_metadata_df.empty:
		summary_clin = pd.DataFrame(all_metadata_df.groupby(['Batch ID Enrichment', 'name', 'quant_val', 'Sample Name'])['pos_ct'].sum())
		clin_metadata_df = clin_metadata_df.join(summary_clin).reset_index()

		contrived_metadata_df = wide_table[~wide_table.index.get_level_values('Sample Type').str.contains('|'.join(searchfor))]
		summary_contrived = pd.DataFrame(all_metadata_df.groupby(['Batch ID Enrichment', 'name', 'quant_val'])['pos_ct'].sum())
		contrived_metadata_df = contrived_metadata_df.join(summary_contrived).reset_index()
		contrived_metadata_df = contrived_metadata_df.groupby(['Batch ID Enrichment', 'name', 'quant_val']).first().reset_index()

		# join back
		contrived_metadata_df = contrived_metadata_df.groupby(['Batch ID Enrichment', 'name', 'quant_val']).first().reset_index()
		wide_table = contrived_metadata_df.append(clin_metadata_df, ignore_index=True)
	else:
		summary = pd.DataFrame(all_metadata_df.groupby(['Batch ID Enrichment', 'name', 'quant_val', 'Sample Name'])['pos_ct'].sum())
		wide_table = wide_table.join(summary).reset_index()
		wide_table = wide_table.groupby(['Batch ID Enrichment', 'name', 'quant_val', 'Sample Name']).first().reset_index()
	
	
	print(" Rearranging columns order in Pivot Table ...", flush=True)
	wide_table = wide_table[col_order]

	dfs = format_lod(wide_table)

	write_formatted_lod(dfs, outfile)

if __name__ == "__main__":

    main()
