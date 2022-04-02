#!/usr/bin/env python

import argparse
import re
import csv
from Bio.Seq import Seq
import os

######################### ARGPARSE #########################
parser = argparse.ArgumentParser(description='Identifiy potential ORFs within extended nucleotide sequences.')
parser.add_argument('-c', '--csv_file', help='csv file output from LitiContigPull.py (required)')
parser.add_argument('-o', '--output_file', help='output file')
args = parser.parse_args()

######################### FUNCTIONS #########################

def find_longest_orf(nt_seq: str) -> str:
    '''Find longest orf seq for input nt seq (input must be 5'->3').'''
    # for each dict entry (strain), find longest orf
    orf_re_pattern = "(?=(ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)))"
    all_orfs = re.findall(orf_re_pattern, nt_seq.upper())
    longest_orf = max(all_orfs, key=len)
    return longest_orf

def store_gene_dict(filename: str, gene_name: str) -> list:
    '''Store dictionary of strain info for specified gene in file.'''
    gene_dict = []
    with open(filename, "r") as gene_h:
        next(gene_h) # ignore header
        for line in gene_h: # store dict entry for each line, if gene matches
            line = line.strip("\n")
            strain,acc,gene,seq=line.split(",")
            if gene == gene_name: # ensure we're analyzing only the desired gene
                print('Finding', gene,'ORFs in:', strain)
                # find longest orf
                long_orf_nt = find_longest_orf(seq)
                # translate longest orf 
                long_orf_aa = Seq(long_orf_nt).translate(to_stop=True)
                # store in list of dictionaries
                gene_dict.append({"File":strain,"Contig":acc,"Gene":gene,"Seq":long_orf_nt,"Prot":long_orf_aa,"Length":len(long_orf_aa)})
    return gene_dict

######################### MAIN #########################

## read input file and store sequence dictionary for each gene
khr_dict = store_gene_dict(args.csv_file, "KHR")
khs_dict = store_gene_dict(args.csv_file, "KHS")
                        
## write output file
csv_columns = ['File','Contig','Gene','Seq','Prot','Length']
with open(args.output_file, 'w') as output_handle:
    writer = csv.DictWriter(output_handle, fieldnames=csv_columns, lineterminator = '\n')
    writer.writeheader()
    for data in khr_dict:
        writer.writerow(data)
    for data in khs_dict:
        writer.writerow(data)
