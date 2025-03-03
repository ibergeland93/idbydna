import hashlib
import os
import glob
import gzip

import pandas as pd

def get_md5(filename):
    with open(filename, "rb") as f:
        bytes = f.read()  # read file as bytes
        return hashlib.md5(bytes).hexdigest()


pamp_asimov = "/data/analysis_group2/data_vault/datasets/org_challenge_datasets/alloid_test_data/cd5_all_analytical/dataset/*210729B0*gz"
pamp_aws = "/data/analysis_group2/data_vault/datasets/org_challenge_datasets/alloid_test_data/cd5_all_analytical/aws_data_processing/0.4.3/aws_output/210731-2-1/postqual_fastqs/*gz"

with open("pamp_lab_source__md5.txt", "w+") as o:

    directories = [pamp_asimov, pamp_aws]
    pamp_asimov_dict = {}
    pamp_aws_dict = {}

    for path in directories:

        pamp_asimov = '/'.join(pamp_asimov.split('/')[0:9])
        pamp_aws = '/'.join(pamp_aws.split('/')[0:11])

        for target in glob.glob(path):
            with gzip.open(target, 'rb') as f_in:
                target = target.split('/')[-1].split('.')[0]
#                unzipped_file = f_in.read()

#                with open('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/scripts/md5_output/{}.postQual.fastq'.format(target), 'wb') as f_out:
#                    f_out.write(unzipped_file)
#                    f_out.close()
            

                md5 = get_md5('/data/analysis_group2/data_vault/datasets/analytical/all_sourced_analytical/scripts/md5_output/{}.postQual.fastq'.format(target))

                if pamp_asimov in path:
                    print('asimov', target, md5)
                    pamp_asimov_dict[target] = md5
                else:
                    print('aws', target, md5)
                    pamp_aws_dict[target] = md5


        for target, md5 in pamp_asimov_dict.items():

            if target in pamp_aws_dict.keys():

                if pamp_asimov_dict[target] != pamp_aws_dict[target]:
                    print('Checksums differ: {}\t{}\t{}\n'.format(target,md5,pamp_aws_dict[target]))
                elif pamp_asimov_dict[target] == pamp_aws_dict[target]:
                    o.write("{}\t{}\t{}\t{}\n".format(target, md5, target, pamp_aws_dict[target]))
                else:
                    print('bad', target)
    
    o.close()
