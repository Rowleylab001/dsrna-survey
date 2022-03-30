#!/bin/env sh

## de novo assembly pipeline
## Author: Angela Crabtree

########## OPTIONS ############
while getopts "s:1:2:r:o:th" OPTION; do
	case $OPTION in
		s) SAMPLE=$OPTARG ;;
		1) FREAD=$OPTARG ;;
		2) RREAD=$OPTARG ;;
		o) OUTDIR=$OPTARG ;; 
		t) MODE="test" ;;
		h)
			printf "\n\n-------------------------------------------------------\n"
			printf "\nOptions:\n"
			printf "\n"
			printf "   -s [arg]	sample name (required)\n"
			printf "   -1 [arg]	forward read file (required)\n"
			printf "   -2 [arg]	reverse read file (required)\n"
			printf "   -o [arg]	destination of output folder (required)\n"
			printf "   -t		test (ensures required CL apps are working)\n"
			printf "   -h		help\n"
			printf "\n-------------------------------------------------------\n\n\n"
			exit 0
			;;
	esac
done

########## ESTABLISH BACKGROUND LOGISTICS ############

## LOAD MODULES SINCE WE'RE USING THIS ON A SERVER
	module load java
	module load python
	module load spades
	module load R
	module load ncbi-blast

## SET UP DIRECTORIES
	S_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"   # location of this script file
	if [[ ${MODE} = "test" ]]; then
		printf "\n\tOutput files will be saved to your current directory.\n"
		OUTDIR=$(pwd)
		SAMPLE="TEST_STRAIN"
		FREAD="${S_DIR}"/test/NCYC190_S3_L001_R1_001.fastq.gz
		RREAD="${S_DIR}"/test/NCYC190_S3_L001_R2_001.fastq.gz
	fi	
	mkdir -p "${OUTDIR}"/output/"${SAMPLE}"/denovo/stats   # store stats
	WORKDIR="${OUTDIR}"/output/"${SAMPLE}"/denovo   # store sample output
	APPDIR=~/bin # this is the folder where fastp and bbmap are stored

############################## READ CLEANING ###################################

## GENERATE CONTIGS AFTER CLEANING READS WITH FASTP
	printf "\n\n*** Cleaning reads with Fastp ***\n\n"
	# run fastp, pretty much just for the .html report
	${APPDIR}/fastp/fastp \
		-i ${FREAD} \
		-I ${RREAD} \
		-o "${WORKDIR}"/${SAMPLE}_R1_filt.fastq.gz \
		-O "${WORKDIR}"/${SAMPLE}_R2_filt.fastq.gz \
		-h "${WORKDIR}"/stats/${SAMPLE}_fastp.html -c -f 6 -F 6 -w 16
	rm fastp.json

############################## READ FILTERING ###################################

## SPECIFY REFERENCE GENOMES TO FILTER OUT
	if [[ ${MODE} = "test" ]]; then
		REF="${S_DIR}"/refseq/Scer.fasta.gz # this is the genome to filter out in the test condition
	else
		SCER_REF="${S_DIR}"/refseq/Scer.fasta.gz
		SPARA_REF="${S_DIR}"/refseq/Spara.fasta.gz
		HSAP_REF="${S_DIR}"/refseq/Hsap.fasta.gz
		REF=${SCER_REF},${SPARA_REF},${HSAP_REF} # these are the genomes to filter out
	fi

## REMOVE READS THAT MAP TO REFERENCE GENOME(S)
	printf "\n\n*** Filtering S. cerevisiae, S. paradoxus, and H. sapiens reads with BBsplit ***\n\n"
	cd "${WORKDIR}" # you MUST be in the same working directory as the reads or bbsplit won't work
	sh ${APPDIR}/bbmap/bbsplit.sh -da \
		in1=${SAMPLE}_R1_filt.fastq.gz \
		in2=${SAMPLE}_R2_filt.fastq.gz \
		ref=${REF} \
		basename=dirty_reads_%.fastq.gz \
		outu1=clean1.fastq.gz \
		outu2=clean2.fastq.gz
	rm -r "${WORKDIR}"/ref
	
############################## ASSEMBLE CONTIGS ###################################

## ASSEMBLE USING SPADES DEFAULT PARAMETERS, + AUTO COVERAGE CUTOFF
	printf "\n\n*** Assembling contigs with SPAdes ***\n\n"
 	spades.py \
 		--cov-cutoff 5 \
 		-1 "${WORKDIR}"/clean1.fastq.gz \
 		-2 "${WORKDIR}"/clean2.fastq.gz \
 		-o "${WORKDIR}"/${SAMPLE}_spades_output	
 	cp "${WORKDIR}"/${SAMPLE}_spades_output/contigs.fasta \
 		"${WORKDIR}"/${SAMPLE}_contigs.fasta

############################## STATS & GRAPHS ###################################

	printf "\n\n*** Helpful Stats ***\n\n"
	READWC1=$(echo "$(gzip -dc ${FREAD} | wc -l | cut -d" " -f1)/4" | bc)
	READWC2=$(echo "$(gzip -dc ${RREAD} | wc -l | cut -d" " -f1)/4" | bc)
	T_READ_I=$(echo "$READWC1 + $READWC2"| bc)
	printf "\n    Number of original reads: $T_READ_I \n"
	
	READWC3=$(echo "$(gzip -dc ${WORKDIR}/clean1.fastq.gz | wc -l | cut -d" " -f1)/4" | bc)
	READWC4=$(echo "$(gzip -dc ${WORKDIR}/clean2.fastq.gz | wc -l | cut -d" " -f1)/4" | bc)
	T_READ_F=$(echo "$READWC3 + $READWC4"| bc)
	printf "\n    Number of polished reads: $T_READ_F \n"
	
	T_READ_HS=$(echo "$(gzip -dc ${WORKDIR}/dirty_reads_Hsap.fastq.gz | wc -l | cut -d" " -f1)/4" | bc)
	T_READ_SC=$(echo "$(gzip -dc ${WORKDIR}/dirty_reads_Scer.fastq.gz | wc -l | cut -d" " -f1)/4" | bc)
	T_READ_SP=$(echo "$(gzip -dc ${WORKDIR}/dirty_reads_Spara.fastq.gz | wc -l | cut -d" " -f1)/4" | bc)
	READWC5=$(printf %.2f $(echo "($T_READ_HS/$T_READ_I)*100"| bc -l))
	printf "\n    Percent of reads from human contamination: $READWC5 "
	READWC6=$(printf %.2f $(echo "($T_READ_SC/$T_READ_I)*100"| bc -l))
	printf "\n    Percent of reads from S. cerevisiae contamination: $READWC6 "
	READWC7=$(printf %.2f $(echo "($T_READ_SP/$T_READ_I)*100"| bc -l))
	printf "\n    Percent of reads from S. paradoxus contamination: $READWC7 \n\n"

## GRAPH READ DEPTH 
	printf "\n*** RScript - Contig Graph ***\n"
	cd "${S_DIR}"
	chmod a+x CL_contig_graph.R
	./CL_contig_graph.R \
		"${WORKDIR}"/${SAMPLE}_contigs.fasta

## LOCAL BLASTN
	echo "\n\n*** Performing blastn Search of Narnavirus Database ***\n\n"
	# make necessary database (only needs to be done once, but this way no one has to check)
	makeblastdb -in "${S_DIR}"/refseq/dsrna.fasta -dbtype nucl -parse_seqids
	blastn -db "${S_DIR}"/refseq/dsrna.fasta \
		-query "${WORKDIR}"/${SAMPLE}_contigs.fasta \
		-max_target_seqs 1 \
		-outfmt "10 qacc sacc stitle pident length qstart qend sstart send qseq" \
		-out "${WORKDIR}"/${SAMPLE}_dsrna_blast.csv
	# Format the output file
	echo "\n\n*** Formatting Blast Output ***\n\n"
	# adds a column for the sample name
	sed -i'.og' -e "s/^/${SAMPLE}, /" "${WORKDIR}"/${SAMPLE}_dsrna_blast.csv
	#rm "${WORKDIR}"/*.og
	# adds column headers
	gsed -i'.og2' -e '1s/^/sample,qacc,sacc,stitle,pident,length,qstart,qend,sstart,send,qseq\n/' \
		"${WORKDIR}"/${SAMPLE}_dsrna_blast.csv
	#rm "${WORKDIR}"/*.og2
	
	# also output in html format
	blastn -db "${S_DIR}"/refseq/dsrna.fasta \
		-query "${WORKDIR}"/${SAMPLE}_contigs.fasta \
		-max_target_seqs 1 \
		-html \
		-out "${WORKDIR}"/${SAMPLE}_dsrna_blast.html
