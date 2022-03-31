#!/bin/bash

## Liti BLASTn pipeline
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
			printf "\noptions:\n"
			printf "\n"
			printf "   -f [arg]	directory containing fasta files (required)\n"
			printf "   -b [arg]	BLASTn reference fasta file (required)\n"
			printf "   -o [arg]	destination of output folder\n"
			printf "   -t		test (ensures required CL apps are working)\n"
			printf "   -h		help\n"
			printf "\n-------------------------------------------------------\n\n\n"
			exit 0
			;;
	esac
done

########## ESTABLISH BACKGROUND LOGISTICS ############

## SET UP DIRECTORIES
	scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"   # location of this script file
	if [$mode == "test"]; then
		printf "\n\tOutput files will be saved to your current directory.\n"
		usroutdir=$(pwd)
		fasta_dir=${scriptdir}/../data/Liti_contigs/_local/GENOMES_ASSEMBLED
		blast_db=${scriptdir}/../data/reference_seqs/chromTox_nt_db.fasta
	fi	
    blast_output="${usroutdir}"/liti_blastn_output.csv

###################### BLASTn ######################

# make necessary database (only needs to be done once, but this way no one has to check)
	makeblastdb \
        -in $blast_db \
        -dbtype nucl \
        -parse_seqids

# perform blast on each file in directory
printf "\n\n*** Performing BLASTn using local BLAST Database ***\n\n"
for filename in "$fasta_dir"/*.fa; do
    base_filename="$(basename ${filename})"
	printf "\nblasting ${base_filename}\n\n"
    temp_blast_output="${usroutdir}/${base_filename}.csv"

    ## run BLASTn to find KHR and KHS in Liti collection contigs
	blastn -db $blast_db \
		-query $filename \
		-max_target_seqs 1 \
		-outfmt "10 qacc sacc stitle pident length qstart qend sstart send" \
		-out "${temp_blast_output}"

	# adds a column for the sample name
    awk -v FOO=${base_filename} \
        '{print FOO","$0}' "${temp_blast_output}" \
        > "${temp_blast_output}".tmp
    # append to main results file
    cat "${temp_blast_output}".tmp >> "${blast_output}".tmp
    # remove intermediary files to prevent file bloat
    rm "${temp_blast_output}".tmp
    rm "${temp_blast_output}"
done

# adds column headers (may not work on non-Mac OS)
sed "1i\\
filename,contig,match_acc,match_description,pident,length,qstart,qend,sstart,send
" "${blast_output}".tmp > "${blast_output}"
# remove temp file
rm "${blast_output}".tmp