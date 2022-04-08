#!/bin/bash

########## SET FILENAME VARIABLES ############

	scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"   # location of this script file
	fasta_dir=$1
	blast_db=$2
	blast_output=$3
	usroutdir="$(dirname ${blast_output})"
	tmp_blast_db=${usroutdir}/blast_db_tmp.fasta

###################### BLASTn ######################

# Make necessary database (only needs to be done once, but this way no one has to check)
	# first, reformat headers if they contain commas (prevents downstream errors)
	tr ',' ';' < $blast_db | tr ':' ' ' > $tmp_blast_db
	# generate blast database index files
	makeblastdb \
        -in $tmp_blast_db \
        -dbtype nucl \
        -parse_seqids

# Perform blast on each file in specified directory
for filename in ${fasta_dir}/*.fa*; do
    base_filename="$(basename ${filename})"
    sample="$(echo ${base_filename} | sed -e 's/_*[[:digit:]]*.re.fa$//')"
	printf "\nblasting ${base_filename}\n\n"
    temp_blast_output=${blast_output}.txt

    ## run BLASTn to find KHR and KHS in contigs
    ## options here: https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.options_common_to_all_blast/
	blastn \
        -query $filename \
        -db $tmp_blast_db \
        -evalue 1e-1 \
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
" ${blast_output}.tmp > $blast_output

# remove temp files
rm ${blast_output}.tmp
rm $tmp_blast_db
rm ${tmp_blast_db}.*