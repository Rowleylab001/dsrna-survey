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

## Set filename variables
	scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"   # location of this script file
	if [ -n $mode ]; then
		printf "\n\tOutput files will be saved to your current directory.\n"
		usroutdir=$(pwd)
		fasta_dir=${scriptdir}/../data/Liti_contigs/_local/GENOMES_ASSEMBLED
		blast_db=${scriptdir}/../data/reference_seqs/chromTox_nt_db.fasta
	fi	
    blast_output=${usroutdir}/liti_blastn_output.csv
    chr_categories=${usroutdir}/liti_chr_categories.csv
    chr_extended=${usroutdir}/liti_chr_extended.csv
    chr_revised=${usroutdir}/liti_chr_revised.csv
    chr_quality_check=${usroutdir}/liti_chr_compare.csv

###################### BLASTn ######################

# Make necessary database (only needs to be done once, but this way no one has to check)
	makeblastdb \
        -in $blast_db \
        -dbtype nucl \
        -parse_seqids

# Perform blast on each file in specified directory
printf "\n\n*** Performing BLASTn using local BLAST Database ***\n\n"
for filename in ${fasta_dir}/*.fa; do
    base_filename="$(basename ${filename})"
    sample="$(echo ${base_filename} | sed -e 's/_*[[:digit:]]*.re.fa$//')"
	printf "\nblasting ${base_filename}\n\n"
    temp_blast_output="${usroutdir}/${base_filename}.csv"

    ## run BLASTn to find KHR and KHS in Liti collection contigs
    ## options here: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.options_common_to_all_blast/
	blastn \
        -query $filename \
        -db $blast_db \
        -evalue 1e-1 \
		-max_target_seqs 1 \
		-outfmt "10 qacc sacc stitle evalue pident length qstart qend sstart send qseq" \
		-out $temp_blast_output

	# adds a column for the sample name
    awk -v FOO=${base_filename} \
        '{print FOO","$0}' $temp_blast_output \
        > ${temp_blast_output}.tmp
    # append to main results file
    cat ${temp_blast_output}.tmp >> ${blast_output}.tmp
    # remove intermediary files to prevent file bloat
    rm ${temp_blast_output}.tmp
    rm $temp_blast_output
done

# adds column headers (may not work on non-Mac OS)
sed "1i\\
filename,contig,match_acc,match_description,eval,pid,ntlen,qstart,qend,sstart,send,qseq_nt
" ${blast_output}.tmp > ${blast_output}
# remove temp file
rm ${blast_output}.tmp

###################### TRANSLATE & CATEGORIZE ######################
# Translate BLASTn hits to protein sequences and give designations depending 
# on if they’re truncated, and whether the full-length proteins are canonical 
# or mutated. 
printf "\n\n*** Translating and categorizing BLASTn output ***\n\n"
RScript liti_chr_aa_convert.R $blast_output $blast_db $chr_categories

###################### EXPAND HITS ######################
# Obtain nucleotide sequences 50 bp to either side of the BLASTn hit and 
# identify KHR/KHS-like ORFs from within the extended sequence. 
# This ensures BLASTn did not exclude nucleotide sequences on either end 
# due to sequence dissimilarity.
printf "\n\n*** Expanding BLAST hit areas ***\n\n"
python3 LitiContigPull.py -f $fasta_dir -c $chr_categories -o $chr_extended

###################### REVISE ORFs ######################
# Identifiy potential ORFs within the extended nucleotide sequences
printf "\n\n*** Checking for ORFs in expanded hits ***\n\n"
python3 LitiORFsWithRvs.py -c $chr_extended -o $chr_revised

###################### CHECK RESULTS ######################
# Merge new list of extended sequences with original BLASTn results to compare 
# protein sequences. We can then see if there was any KHR or KHS genes erroneously 
# incomplete by the BLASTn. 
printf "\n\n*** Checking the BLAST hit ORFs match largest expanded ORFs ***\n\n"
RScript RLitiORFComparison.R $chr_categories $chr_revised $chr_quality_check

printf "\nAll done!\n\n" 