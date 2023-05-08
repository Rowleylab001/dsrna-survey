#!/usr/bin/env python

import argparse
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Align
import re
import csv

######################### ARGPARSE ########################
parser = argparse.ArgumentParser(description='Ensures BLAST hits capture extended proteins and categorizes translated proteins based on reference sequences.')
parser.add_argument('-f', '--fasta_dir', help='folder containing original fasta files (required)')
parser.add_argument('-b', '--blast_db', help='fasta file containing reference sequences used in BLAST')
parser.add_argument('-c', '--csv_file', help='csv file output from 01_chrtox_blastn.sh (required)')
parser.add_argument('-o', '--output_file', help='output filename')
args = parser.parse_args()

######################## FUNCTIONS ########################
def revc (seq: str) -> str:
	'''converts DNA seq to reverse complement'''
	seq=seq.upper()
	nucd = {"A":  "T", "T": "A", "C": "G", "G": "C", "N": "N"}
	rseq=seq.translate(str.maketrans(nucd))[::-1] 
	return rseq

def make_blast_db_dict(blast_db_file: str) -> dict:
    '''Create a dictionary containing blast database reference sequences'''
    blast_dict={}
    for gene in SeqIO.parse(blast_db_file, "fasta"):
        acc = gene.id.split(".")[0]
        prot = str(gene.seq.translate(to_stop=True))
        blast_dict[acc] = prot
    return blast_dict

def find_longest_orf(nt_seq: str) -> str:
    '''Find longest orf seq for input nt seq (input must be 5'->3').'''
    # find longest orf
    longest_orf = "" # returns empty string if no orfs found
    orf_re_pattern = "(?=(ATG(?:(?!TAA|TAG|TGA)...)*(?:TAA|TAG|TGA)))"
    all_orfs = re.findall(orf_re_pattern, nt_seq.upper())
    if all_orfs != []: # if not empty list (i.e. no orfs found)
        longest_orf = max(all_orfs, key=len)
    return longest_orf

def calc_pid_shortseq(seq1: str, seq2: str) -> float:
    '''Calculate percent identity after aligning two protein sequences'''
    # pid = 100 * (identical positions) / (length shorter sequence)
    # will return 100% if number of aligned matches = length of shortest seq
    aligner = Align.PairwiseAligner()
    shortest_len = len(seq1) if len(seq1) < len(seq2) else len(seq2)
    pid = round(
        100*aligner.align(seq1, seq2).score/shortest_len,
        1)
    return pid

def calc_pid_refseq(qseq: str, refseq: str) -> float:
    '''Calculate percent identity after aligning two protein sequences'''
    # pid = 100 * (identical positions) / (length reference sequence)
    # will return 100% if number of aligned matches = length of shortest seq
    aligner = Align.PairwiseAligner()
    pid = round(
        100*aligner.align(qseq, refseq).score/len(refseq),
        1)
    return pid

def pid_series(q_col: pd.Series, ref_col: pd.Series) -> pd.Series:
    return pd.Series([
        calc_pid_refseq(qseq, ref_seq)
        for (qseq, ref_seq) in zip(q_col, ref_col)
        ])

def categorize(qseq: str, refseq: str) -> str:
    '''Categorize a mutation based on reference sequence.'''
    if qseq==refseq: 
        return 'canonical'
    elif len(qseq)==len(refseq):
        return 'mutant'
    elif len(qseq)<len(refseq):
        return 'truncated'
    else:
        return 'extended'

def categorize_series(df: pd.DataFrame) -> pd.Series:
    return pd.Series([
        categorize(qseq, refseq)
        for (qseq, refseq) in zip(df['ext_aaseq'], df['ref_aaseq'])
        ])

def print_ext_fasta_seq(fa_file: str, contig_acc: str, start_pos: int, end_pos: int)->str:
    f = open(fa_file, 'r')
    extension_num = 50
    for record in SeqIO.parse(f, "fasta"):
        print(record.id)
        print(contig_acc)
        if record.id == contig_acc:
            # determine extended start position
            start = start_pos - extension_num
            if start_pos < extension_num: # if extension exceeds start of contig 
                start = 0
            # determine extended stop position
            stop = end_pos + extension_num    
            if end_pos > len(record.seq) - extension_num:
                stop = len(record.seq)
                
            #A second nested if statement needs to go above for STOP, where it compares to length of seq.
            subseq = record.seq[start:stop]
            # if gene is in reverse orientation, reverse sequence 
            if gene_start > gene_stop:
                subseq = subseq.reverse_complement()
            subseq.id = record.id

def extend_BLASTn_hits(basepath: str, blast_filename: str) -> pd.Series:
    # import blast csv file as pandas dataframe
    blast_df = pd.read_csv(blast_filename)
    blast_df = blast_df.reset_index()  # make sure indexes pair with number of rows
    ext_seq_nt_list=[] # initialize list for storing extended sequences
    for index, row in blast_df.iterrows():
        ext_seq_nt="" # initialize extended sequence 
        # convert nt positions to python positions 
        hit_seq_start = int(row['qstart'])-1
        hit_seq_stop = int(row['qend'])-1
        gene_start = int(row['sstart'])-1
        gene_stop = int(row['send'])-1 
        # open contig file and extract extended nucleotide sequence
        print("Opening: " + row['filename'])
        openme = basepath + "/" + row['filename']
        handle = open(openme)    #Opens an assembly file as "handle" depending on the strain variable.
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == row['contig']:
                start=hit_seq_start # initialize start position
                stop=hit_seq_stop # initialize stop position
                #Parses the fasta file, dictates sequence to pull, and places into a dictionary.
                if int(hit_seq_start) < 50:
                    start = 0
                else:
                    start = int(hit_seq_start) - 50
                if int(hit_seq_stop) > int(len(record.seq))-50:
                    stop = len(record.seq)
                else:
                    stop = int(hit_seq_stop) + 50
                #A second nested if statement needs to go above for STOP, where it compares to length of seq.
                ext_seq_nt = record.seq[start:stop]
                # if gene is in reverse orientation, reverse sequence 
                if gene_start > gene_stop:
                    #print(record.id+" is reversed!")
                    ext_seq_nt = ext_seq_nt.reverse_complement()
                break # stops looping after it stores info from contig file
        ext_seq_nt_list.append(str(ext_seq_nt))   
    # save extended sequences to new column
    blast_df['ext_seq_nt'] = ext_seq_nt_list
    blast_df=blast_df.drop(columns=['index'])
    return(blast_df)

########################## MAIN ###########################

# store pandas dataframe containing blast hit info and extended hits 
# (hits extended 50bp to either side of hit start and stop locations)
blast_pd = extend_BLASTn_hits(basepath=args.fasta_dir, blast_filename=args.csv_file)

# store blast database as dictionary
blast_dict = make_blast_db_dict(args.blast_db)

# create new column for revcomp seqs
blast_pd['ntseq'] = blast_pd['qseq_nt'].str.replace('-','') # delete gaps
blast_pd['ntseq'] = np.where( # return reverse complement if gene is reversed
    blast_pd['sstart'] > blast_pd['send'], 
    np.array([Seq(x).reverse_complement() for x in blast_pd['ntseq']], dtype=object), 
    blast_pd['ntseq'])

# create new column to print orf nt seq
blast_pd['orfseq'] = blast_pd['ntseq'].apply(lambda x: 
    find_longest_orf(str(x)))

# create new column to designate if blast hit contained complete orf
blast_pd['hit_v_orf'] = np.where( 
    len(blast_pd['ntseq']) >= len(blast_pd['orfseq']), 
    'good', 
    'incomplete_hit')

# create new column for original translated protein seqs
blast_pd['og_aaseq'] = np.where( 
    blast_pd['ntseq'] != '', 
    blast_pd['ntseq'].apply(lambda x: 
    ''.join(Seq(x).translate(to_stop=True))), 
    '')

# create new column for extended translated protein seqs
blast_pd['ext_aaseq'] = np.where( 
    blast_pd['orfseq'] != '', 
    blast_pd['orfseq'].apply(lambda x: 
    ''.join(Seq(x).translate(to_stop=True))), 
    '')

# create new column for best aaseq
blast_pd['best_aaseq'] = np.where( 
    len(blast_pd['ntseq']) >= len(blast_pd['orfseq']), 
    blast_pd['og_aaseq'], 
    blast_pd['ext_aaseq'])

# create new column for translated protein seq lengths
blast_pd['best_aalen'] = blast_pd['best_aaseq'].str.len()

# create column for reference seqs
blast_pd['ref_aaseq'] = blast_pd['match_acc'].map(blast_dict)

# create new column for best_aaseq/ref percent identity 
blast_pd['best_aapid_ref'] = pid_series(blast_pd['best_aaseq'], blast_pd['ref_aaseq'])

# create new column for mutation type
blast_pd['mutation_type'] = categorize_series(blast_pd)

# Write DataFrame to CSV
blast_pd.to_csv(args.output_file, index=False)