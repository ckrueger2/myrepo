import os
import argparse
import numpy as np
import pandas as pd
import argparse
import sys
import hail as hl
import subprocess

#arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="import hail table from phenotype accession number")
    parser.add_argument('--phecode', help='path to input file', required='True')
    parser.add_argument('--pop', help='path to input file', required='True')
    return parser.parse_args(args)

args = check_arg(sys.argv[1:])

phenotype_id = args.phecode
pop = args.pop

#check if phecode file exists
path = f"gs://fc-aou-datasets-controlled/AllxAll/v1/ht/ACAF/{pop}/phenotype_{phenotype_id}_ACAF_results.ht"

result = os.system(f"gsutil -u $GOOGLE_PROJECT ls {path} > /dev/null 2>&1")

if result == 0:
    print(f"Phenotype {phenotype_id} is in the All of Us database")
else:
    print(f"Phenotype {phenotype_id} is not in the All of Us database; enter valid phenotype ID")
    sys.exit(1)

#PULL TABLE FROM ALL OF US DATABASE   
#initialize hail
hl.init()

#define bucket to save to
bucket = os.getenv('WORKSPACE_BUCKET')
bucket # gs://fc-secure-bb61452f-d5e2-4d26-9227-6a9444241af8/

#not necessary, shows files in bucket searched
os.system(f"gsutil -u $GOOGLE_PROJECT ls gs://fc-aou-datasets-controlled/AllxAll/v1/ht/ACAF/{pop}/phenotype_{phenotype_id}_ACAF_results.ht")

#find hail table and save to variable
ht = hl.read_table(f"gs://fc-aou-datasets-controlled/AllxAll/v1/ht/ACAF/{pop}/phenotype_{phenotype_id}_ACAF_results.ht")

#save full table to bucket for S-PrediXcan input
ht_path = f'{bucket}/data/{pop}_full_{phenotype_id}.tsv'
ht.export(ht_path)

#show first few lines of hail table
ht.show(20)

#summary of table
globals_dict = ht.globals.collect()[0]

#extract specific fields
n_cases = globals_dict['n_cases']
n_controls = globals_dict['n_controls']
heritability = globals_dict['heritability']
meta_pop = globals_dict['meta_pop']
phenoname = globals_dict['phenoname']

#print the fields
print(f"Number of cases: {n_cases}")
print(f"Number of controls: {n_controls}")
print(f"Heritability: {heritability}")
print(f"Population(s): {meta_pop}")
print(f"Phenotype name: {phenoname}")

#table dimentions
rows, cols = ht.count(), len(ht.row)
print(f"Table dimensions: {rows} rows x {cols} columns")

#FILTER BY PVALUE
#set thresholds
pval = 0.05
max_snps = 2000000
min_pval = 5e-8

#intial thinning of SNPs
significant_snps = ht.filter(ht.Pvalue < pval)
num_snps = significant_snps.count()

#continue to thin until less than max_snps threshold
while num_snps > max_snps and pval > min_pval:
    #make threshold 5x more stringent
    pval /= 5
    if pval < min_pval:
        pval = min_pval
    #thin again
    significant_snps = ht.filter(ht.Pvalue < pval)
    num_snps = significant_snps.count()
    if pval == min_pval:
        break

#print number of SNPs and pvalue in filtered table
print(f"Number of SNPs: {num_snps} at pvalue < {pval}")

#sort SNPs
sig_snps_sorted = significant_snps.order_by(significant_snps.Pvalue)

#print top SNPs
sig_snps_sorted.show(10)

#save filtered file
filtered_path = f'{bucket}/data/{pop}_filtered_{phenotype_id}.tsv'
sig_snps_sorted.export(filtered_path)

#CHECK IF FILES ARE SAVED TO BUCKET
#full file
try:
    check_filtered = subprocess.check_output(
        f"gsutil ls {bucket}/data/ | grep {ht_path}", 
        shell=True, 
        stderr=subprocess.DEVNULL
    )
    #if command succeeded 
    print("Full file successfully saved to bucket.\n")
except subprocess.CalledProcessError:
    #if command failed
    sys.exit(f"ERROR: File '{ht_path}' was not found in {bucket}/data/.\n")
    
#pvalue file
try:
    check_filtered = subprocess.check_output(
        f"gsutil ls {bucket}/data/ | grep {filtered_path}", 
        shell=True, 
        stderr=subprocess.DEVNULL
    )
    #if command succeeded 
    print("Pvalue filtered file successfully saved to bucket.\n")
except subprocess.CalledProcessError:
    #if command failed
    sys.exit(f"ERROR: File '{filtered_path}' was not found in {bucket}/data/.\n")
