#!/usr/bin/env python

from Bio import SeqIO
import csv
import argparse

######################### ARGPARSE #########################
parser = argparse.ArgumentParser(description='BLAST sanity check to ensure results were not fragmented.')
parser.add_argument('-f', '--fasta_dir', help='folder containing original fasta files (required)')
parser.add_argument('-c', '--csv_file', help='csv file output from liti_chr_aa_convert.R (required)')
parser.add_argument('-o', '--output_file', help='output file')
args = parser.parse_args()

######################### MAIN #########################
contigs=[]                        
seqdict={}
totaldict=[]
#Open contig list file as your data file.
data = open(args.csv_file, 'r')
contigfile = csv.reader(data) #Redefines the data file after reading it.
#Headers are dumb. We don't need them where we're going.
header = next(contigfile) 
if header != None:              
    for i in contigfile:
        #Defining variables to use from the data file for our assembly files.
        strain = i[0]
        hit_contig = i[1]
        hit_seq_start = i[7] 
        hit_seq_stop = i[8]
        gene = i[14]
        #Doublechecking strain is correct /w a print.
        print("Opening:" + strain)
        basepath = args.fasta_dir
        filename = strain    
        openme = basepath + "/" + filename
        handle = open(openme)    #Opens an assembly file as "handle" depending on the strain variable.
        for record in SeqIO.parse(handle, "fasta"):
            if record.id == hit_contig:
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
                subseq = record.seq[start:stop]
                subseq.id = record.id
                seqdict = {
                        "strain":strain,
                        "acc":subseq.id,
                        "gene":gene,
                        "seq":subseq
                        }
                totaldict.append(seqdict)
        handle.close()  #Close before exporting data.
data.close()

#Exporting of the total dictionary list to a csv for access.
with open(args.output_file, 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=['strain','acc','gene','seq'])
    writer.writeheader()
    for data in totaldict:
        writer.writerow(data)