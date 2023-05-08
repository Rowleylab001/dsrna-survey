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
    chr_extended=${usroutdir}/03_chr_extended.csv
    chr_revised=${usroutdir}/04_chr_revised.csv
    chr_quality_check=${usroutdir}/05_chr_compare.csv

###################### BLASTn ######################
## Run BLASTn to find KHR and KHS in contigs
printf "\n\n*** Performing BLASTn using local BLAST Database ***\n\n"
chmod a+x 01_chrtox_blastn.sh
./01_chrtox_blastn.sh $fasta_dir $blast_db $blast_output

###################### TRANSLATE & CATEGORIZE ######################
# Translate BLASTn hits to protein sequences and give designations depending 
# on if they’re truncated, and whether the full-length proteins are canonical 
# or mutated. 
printf "\n\n*** Translating and categorizing BLASTn output ***\n\n"
chmod a+x 02_chr_aa_convert.R
./02_chr_aa_convert.R $blast_output $blast_db $chr_categories

###################### EXPAND HITS ######################
# Obtain nucleotide sequences 50 bp to either side of the BLASTn hit and 
# identify KHR/KHS-like ORFs from within the extended sequence. 
# This ensures BLASTn did not exclude nucleotide sequences on either end 
# due to sequence dissimilarity.
printf "\n\n*** Expanding BLAST hit areas ***\n\n"
chmod a+x 03_extend_hits.py
./03_extend_hits.py -f $fasta_dir -c $chr_categories -o $chr_extended

###################### IDENTIFY NEW ORFs ######################
# Identifiy longest ORFs within the extended nucleotide sequences
printf "\n\n*** Finding longest ORFs in expanded hits ***\n\n"
chmod a+x 04_identify_orfs.py
./04_identify_orfs.py -c $chr_extended -o $chr_revised

###################### CHECK RESULTS ######################
# Merge new list of extended sequences with original BLASTn results to compare 
# protein sequences. We can then see if there was any KHR or KHS genes erroneously 
# incomplete by the BLASTn. 
printf "\n\n*** Checking the BLAST hit ORFs match largest expanded ORFs ***\n\n"
chmod a+x 05_compare_orfs.R
./05_compare_orfs.R $chr_categories $chr_revised $chr_quality_check

printf "\nAll done!\n\n" 