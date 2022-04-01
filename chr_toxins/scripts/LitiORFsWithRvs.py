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
def store_gene_dict(filename: str, gene_name: str) -> dict:
    '''Store dictionary of strain info for specified gene in file.'''
    gene_dict = {}
    with open(filename, "r") as gene_h:
        gene_seq = ''
        gene_header = ''
        for line in gene_h:
            line = line.strip("\n")
            strain,acc,gene,seq=line.split(",")
            if gene == gene_name:
                gene_dict[strain] = seq
    return gene_dict

######################### MAIN #########################
## read input file
input_handle = open(args.csv_file,'r')
contigdata = csv.reader(input_handle)
next(contigdata) # ignores header

## set global variables
tempdict={}
finaldict=[]

## find orfs
for i in contigdata: 
    strain = i[0]
    contig = i[1]
    gene= i[2]
    seq = Seq(i[3])
    StrSeq = str(seq)

    ## find ORFs 
    print('Getting ORFs of:', strain)
    for startMatch in re.finditer('ATG', StrSeq): # look for start codons
        # store nt seq from start codon to end of seq
        remaining = seq[startMatch.start():] 
        # look for downstream stop codon
        for stopMatch in re.finditer('TAA|TGA|TAG',str(remaining)): 
            # save seq from first start to first stop codon
            substring = remaining[:stopMatch.end()] 
            # proceed only if stop is in-frame and > 20 aa
            if len(substring) % 3 == 0 and len(substring)>20: 
                # translate nt seq to aa
                subprot=substring.translate(to_stop=True)
                # identify KHR match
                if "MGHLAI" in subprot or "DSIYLK" in subprot or "DSIYLK" in subprot: 
                    tempdict = {"File":strain,"Contig":contig,"Gene":gene,"Seq":substring,"Prot":subprot}
                    finaldict.append(tempdict)  
                    print('{} aa: {}...{}'.format(len(subprot),subprot[0:10],subprot[len(subprot)-11:len(subprot)-1]))
                    print("stored KHR match") 
                # identify KHS match
                elif "YSLFLST" in subprot or "VIVAALV" in subprot or "QEFATSW" in subprot: 
                    tempdict = {"File":strain,"Contig":contig,"Gene":gene,"Seq":substring,"Prot":subprot}
                    finaldict.append(tempdict)    
                    print('{} aa: {}...{}'.format(len(subprot),subprot[0:10],subprot[len(subprot)-11:len(subprot)-1]))
                    print("stored KHS match") 
                break
                        
## write output file
csv_columns = ['File','Contig','Gene','Seq','Prot']
with open(args.output_file, 'w') as output_handle:
    writer = csv.DictWriter(output_handle, fieldnames=csv_columns, lineterminator = '\n')
    writer.writeheader()
    for data in finaldict:
        writer.writerow(data)

## close files used    
input_handle.close()
output_handle.close()