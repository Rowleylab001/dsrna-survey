#!/bin/bash

## prevent code from continuing to run if errors are encountered
set -e 

## BLASTn pipeline
## Authors: Angela Crabtree, Mason Shipley

########## OPTIONS ############
while getopts "f:b:o:th" opt; do
	case ${opt} in
		f) fasta_dir=$OPTARG ;;
		b) blast_db=$OPTARG ;;
		o) usroutdir=$OPTARG ;; 
		t) mode="test" ;;
		h)
			printf "\n\n-------------------------------------------------------\n"
			printf "options:\n"
			printf "\n"
			printf "   -f [arg]	directory containing fasta files (required)\n"
			printf "   -b [arg]	BLASTn reference fasta file (required)\n"
			printf "   -o [arg]	destination of output folder\n"
			printf "   -t		test (ensures required CL apps are working)\n"
			printf "   -h		help\n"
			printf "-------------------------------------------------------\n\n\n"
			exit 0
			;;
	esac
done

########## ESTABLISH BACKGROUND LOGISTICS ############

## Set filename variables
	scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"   # location of this script file
	if [ -n "$mode" ]; then # check if the $mode variable is an empty string
		printf "\n\tOutput files will be saved to your current directory.\n"
		usroutdir=$(pwd)
		fasta_dir=${scriptdir}/../data/test
		blast_db=${scriptdir}/../data/reference_seqs/chromTox_nt_db.fasta
	fi	
    blast_output=${usroutdir}/01_blastn_output.csv
    chr_categories=${usroutdir}/02_chr_categories.csv
    chr_quality_check=${usroutdir}/03_chr_compare.csv

###################### BLASTn ######################
## Run BLASTn to find KHR and KHS in contigs
printf "\n\n*** Performing BLASTn using local BLAST Database ***\n\n"
chmod a+x 1_chrtox_blastn.sh
./1_chrtox_blastn.sh $fasta_dir $blast_db $blast_output

###################### EXPAND HITS & CATEGORIZE ######################
# Obtain nucleotide sequences 50 bp to either side of the BLASTn hit and 
# identify KHR/KHS-like ORFs from within the extended sequence. 
# This ensures BLASTn did not exclude nucleotide sequences on either end 
# due to sequence dissimilarity.
printf "\n\n*** Expanding BLASTn hits and categorizing mutations ***\n\n"
chmod a+x 2_blast_qc.py
./2_blast_qc.py\
 -f $fasta_dir\
 -b $blast_db\
 -c $blast_output\
 -o $chr_categories

###################### CHECK RESULTS ######################
# Check if there were any longer KHR or KHS ORFs missed by BLASTn. 
printf "\n\n*** Checking the BLAST hits match largest expanded ORFs ***\n\n"
chmod a+x 3_compare_orfs.R
./3_compare_orfs.R $chr_categories $chr_quality_check

###################### MERGE RESULTS ######################
# Merge results with killer assay data



printf "\nAll done!\n\n" 