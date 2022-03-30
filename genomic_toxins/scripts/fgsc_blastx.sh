#!/bin/bash

## BLASTx pipeline
## Author: Angela Crabtree

## load conda environment
conda activate dsrna-survey

## assign variables
d=/projects/bgmp/acrabtre/bioinfo/Bi623/QAA/star
f_read=/projects/bgmp/acrabtre/bioinfo/Bi623/QAA/trim/Undetermined.1.trim.fastq.gz
r_read=/projects/bgmp/acrabtre/bioinfo/Bi623/QAA/trim/Undetermined.2.trim.fastq.gz
prfx="Undetermined_Mmus_"
samfile=$d/${prfx}Aligned.out.sam

cd $d

## make BLAST database (generates indexes for BLAST to use)
/usr/bin/time -v STAR \
    --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir $d/mmus \
    --genomeFastaFiles $d/mmus/Mus_musculus.GRCm39.dna.primary_assembly.fa \
    --sjdbGTFfile $d/mmus/Mus_musculus.GRCm39.104.gtf

## run BLASTx to find KHR and KHS in FGSC collection
/usr/bin/time -v STAR --runThreadN 8 --runMode alignReads \
    --outFilterMultimapNmax 3 \
    --outSAMunmapped Within KeepPairs \
    --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --readFilesCommand zcat \
    --readFilesIn $f_read $r_read \
    --genomeDir $d/mmus \
    --outFileNamePrefix $d/$prfx

## run python script to get counts of mapped and unmapped reads
python /projects/bgmp/acrabtre/bioinfo/Bi623/QAA/scripts/inspector_sam.py $samfile > $d/${prfx}map_counts.txt

## run HTseq to count reads that map to features 
htseq-count --stranded=no $samfile $d/mmus/Mus_musculus.GRCm39.104.gtf > $d/${prfx}hts.genecounts
htseq-count --stranded=yes $samfile $d/mmus/Mus_musculus.GRCm39.104.gtf > $d/${prfx}hts_str.genecounts

######################## RUN STAR ON SAMPLE 2 ########################

## assign variables
d=/projects/bgmp/acrabtre/bioinfo/Bi623/QAA/star
f_read=/projects/bgmp/acrabtre/bioinfo/Bi623/QAA/trim/2_2B_control.1.trim.fastq.gz
r_read=/projects/bgmp/acrabtre/bioinfo/Bi623/QAA/trim/2_2B_control.2.trim.fastq.gz
prfx="2_2B_control_Mmus_"
samfile=$d/${prfx}Aligned.out.sam

## run STAR to assemble reads to reference genome
/usr/bin/time -v STAR --runThreadN 8 --runMode alignReads \
    --outFilterMultimapNmax 3 \
    --outSAMunmapped Within KeepPairs \
    --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --readFilesCommand zcat \
    --readFilesIn $f_read $r_read \
    --genomeDir $d/mmus \
    --outFileNamePrefix $d/$prfx

## run python script to get counts of mapped and unmapped reads
python /projects/bgmp/acrabtre/bioinfo/Bi623/QAA/scripts/inspector_sam.py $samfile > $d/${prfx}map_counts.txt

## run HTseq to count reads that map to features 
htseq-count --stranded=no $samfile $d/mmus/Mus_musculus.GRCm39.104.gtf > $d/${prfx}hts.genecounts
htseq-count --stranded=yes $samfile $d/mmus/Mus_musculus.GRCm39.104.gtf > $d/${prfx}hts_str.genecounts






########## ESTABLISH BACKGROUND LOGISTICS ############

## LOAD MODULES (SERVER USE ONLY)
local=$(module load 2>&1 >/dev/null | grep 'command not found' | wc -l)
if [ $local != 1 ]; then
	module load ncbi-blast
fi

## SET UP DIRECTORIES
	scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"   # location of this script file
	if [ -n "$mode" ]; then
		printf "\n\tOutput files will be saved to your current directory.\n"
		usroutdir=$(pwd)
		sample="TEST_STRAIN"
		fread=${scriptdir}/test/NCYC190_S3_L001_R1_001.fastq.gz
		rread=${scriptdir}/test/NCYC190_S3_L001_R2_001.fastq.gz
	fi	

	refseqdir="${scriptdir}"/../big_data
	outdir="${usroutdir}"/"${sample}"/denovo
	mkdir -p ${outdir}/stats ${outdir}/img ${outdir}/polished-reads ${outdir}/contigs
	blast_db=${scriptdir}/refseq/dsrna.fasta

## INPUT FILES & DIRECTORIES
seq_dir=

## ESTABLISH OUTPUT FILE NAMES & DIRECTORIES
	fastp_html=${outdir}/stats/${sample}_fastp.html
	contig_file=${outdir}/${sample}_contigs.fasta
	contig_img=${outdir}/img/${sample}_contigs.jpeg
	blast_output=${outdir}/${sample}_dsrna_blast.csv
	report_file=${outdir}/${sample}_denovo_report.md
	report_html=${outdir}/${sample}_denovo_report.html
	bbsplitR1=${outdir}/polished-reads/${sample}_R1_bbsplit.fastq.gz
	bbsplitR2=${outdir}/polished-reads/${sample}_R2_bbsplit.fastq.gz

############################## BLASTx ###################################

## Local blastx
	printf "\n\n*** Performing BLASTx ***\n\n"
	# make necessary database (only needs to be done once, but this way no one has to check)
	makeblastdb -in $blast_db -dbtype prot -parse_seqids
	blastn -db $blast_db \
		-query $contig_file \
		-max_target_seqs 1 \
		-outfmt "10 qacc sacc stitle pident length qstart qend sstart send" \
		-out $blast_output
	# Format the output file
	printf "\n\n*** Formatting Blast Output ***\n\n"
	# adds a column for the sample name
	sed -i'.og' -e "s/^/${sample}, /" $blast_output
	# adds column headers (doesn't work on linux for some reason)
	sed -i '1 i\sample,contig,match_acc,match_description,pident,length,qstart,qend,sstart,send' "$blast_output"