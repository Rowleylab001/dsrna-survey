#!/bin/bash

## FGSC BLASTx pipeline
## Authors: Angela Crabtree

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
			printf "   -b [arg]	BLASTx reference fasta file (required)\n"
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
		fasta_dir=${scriptdir}/../data/FGSC_chrom/_local/FGSC_chr_files
		blast_db=${scriptdir}/../data/reference_seqs/chromTox_aa_db.fasta
	fi	
    blast_output="${usroutdir}"/FGSC_blastx_output.csv

###################### BLASTn ######################

# make necessary database (only needs to be done once, but this way no one has to check)
	makeblastdb \
        -in $blast_db \
        -dbtype prot \
        -parse_seqids

# perform blast on each file in directory
printf "\n\n*** Performing BLASTx using local BLAST Database ***\n\n"
for filename in "$fasta_dir"/*.fa; do
    base_filename="$(basename ${filename})"
	printf "\nblasting ${base_filename}\n\n"
    temp_blast_output="${usroutdir}/${base_filename}.csv"

    ## run BLASTx to find KHR and KHS in FGSC collection chromosomes
	blastx -db $blast_db \
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
