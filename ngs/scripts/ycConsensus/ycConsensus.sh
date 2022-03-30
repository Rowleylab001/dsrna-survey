#!/bin/bash

## alignment pipeline
## Authors: Angela Crabtree, Chris Large, Caiti S. Heil

########## OPTIONS ############
while getopts "s:1:2:r:o:th" OPTION; do
	case $OPTION in
		s) SAMPLE=$OPTARG ;;
		1) FREAD=$OPTARG ;;
		2) RREAD=$OPTARG ;;
		r) REF=$OPTARG ;;
		o) OUTDIR=$OPTARG ;; 
		t) MODE="test" ;;
		h)
			printf "\n\n-------------------------------------------------------\n"
			printf "\nOptions:\n"
			printf "\n"
			printf "   -s [arg]	sample name (required)\n"
			printf "   -1 [arg]	forward read file (required)\n"
			printf "   -2 [arg]	reverse read file (required)\n"
			printf "   -r [arg]	reference genome file (required)\n"
			printf "   -o [arg]	destination of output folder (required)\n"
			printf "   -t		test (ensures required CL apps are working)\n"
			printf "   -h		help\n"
			printf "\n-------------------------------------------------------\n\n\n"
			exit 0
			;;
	esac
done

## (for server only) loading modules
#module load muscle
#module load bwa
#module load samtools
#module load bcftools
#module load java
#module load python/3.6.7
#module load R

## set up directories
S_DIR=$( dirname `readlink -e -- "$0"` ) # stores the path to your script directory
if [ "${MODE}" = "test" ]; then
	printf "\n\tOutput files will be saved to your current directory.\n"
	OUTDIR=$(pwd)
	REF="${S_DIR}"/test/M1.fasta
	SAMPLE="TEST_STRAIN"
	FREAD="${S_DIR}"/test/NCYC190_S3_L001_R1_001.fastq.gz
	RREAD="${S_DIR}"/test/NCYC190_S3_L001_R2_001.fastq.gz
fi	

mkdir -p "${OUTDIR}"/output/"${SAMPLE}"/troubleshooting   # make directories
WORKDIR="${OUTDIR}"/output/"${SAMPLE}"   # store sample output
REFNAME=$(basename "${REF}" | cut -d"." -f1)   # store name of reference file, w/o extension

############################## ALIGNMENT ###################################

### Align reads with bwa
printf "\n*** BWA - mem -R ***\n"
bwa mem \
	-R '@RG\tID:'${SAMPLE}'\tSM:'${SAMPLE}'\tLB:1' \
	"${REF}" \
	"${FREAD}" \
	"${RREAD}" \
	> "${WORKDIR}"/${SAMPLE}_${REFNAME}.sam \
	2> "${WORKDIR}"/troubleshooting/bwa_stderr.txt
	# if index files are missing, automatically create them
	if grep -q "fail to locate the index files" "${WORKDIR}"/troubleshooting/bwa_stderr.txt; 
	then
		## create index files
		printf "\n*** BWA - creating index files ***\n"
		bwa index "${REF}"
		samtools faidx "${REF}"
		java -Dpicard.useLegacyParser=false -jar /mnt/c/Users/Admin/bin/picard.jar CreateSequenceDictionary \
			-REFERENCE "${REF}" 

		### Align reads with BWA
		printf "\n*** BWA - mem -R ***\n"
		bwa mem -R '@RG\tID:'${SAMPLE}'\tSM:'${SAMPLE}'\tLB:1' \
		"${REF}" \
		"${FREAD}" \
		"${RREAD}" \
		> "${WORKDIR}"/${SAMPLE}_${REFNAME}.sam 
		2> "${WORKDIR}"/troubleshooting/bwa_stderr.txt
	fi

############################## REFORMATTING ###################################

## converts .sam to .bam
printf "\n*** Samtools - View ***\n"
samtools view \
	-bS "${WORKDIR}"/${SAMPLE}_${REFNAME}.sam \
	-o "${WORKDIR}"/${SAMPLE}_${REFNAME}.bam \

## produces a sorted .bam file that can then be used by other programs (like bcftools)
printf "\n*** Samtools - Sort ***\n"
samtools sort \
	"${WORKDIR}"/${SAMPLE}_${REFNAME}.bam \
	-o "${WORKDIR}"/${SAMPLE}_${REFNAME}_sorted.bam 

## produces a BAM index file (.bai) that could be useful for variant calling, etc.
printf "\n*** Samtools - Index ***\n"
samtools index \
	"${WORKDIR}"/${SAMPLE}_${REFNAME}_sorted.bam 

## Print stats on how well the alignment worked
printf "\n*** Samtools - Flagstat ***\n"
samtools flagstat \
	"${WORKDIR}"/${SAMPLE}_${REFNAME}.bam \
	> "${WORKDIR}"/troubleshooting/samtools_flagstat.txt

## Print read depth information
printf "\n*** Samtools - Read Depth ***\n"
samtools depth \
	"${WORKDIR}"/${SAMPLE}_${REFNAME}_sorted.bam \
	> "${WORKDIR}"/${SAMPLE}_${REFNAME}_read_depth.txt

## Graph read depth 
printf "\n*** RScript - Read Depth Graph ***\n"
cd "${S_DIR}"
chmod a+x read_depth_graphs.R
Rscript read_depth_graphs.R \
	"${WORKDIR}"/${SAMPLE}_${REFNAME}_read_depth.txt

############################## MAKE CONSENSUS SEQUENCE ################################

mkdir -p "${WORKDIR}/bcftools"
## call variants with bcftools
printf "\n*** BCFtools - mpileup ***\n"
bcftools mpileup \
	-Ou -f "${REF}" \
	"${WORKDIR}"/${SAMPLE}_${REFNAME}_sorted.bam \
	| bcftools call \
	-mv -Oz -o "${WORKDIR}"/bcftools/${SAMPLE}_${REFNAME}_calls.vcf.gz 

## normalize indels
printf "\n*** BCFtools - norm ***\n"
bcftools norm \
	-f "${REF}" \
	"${WORKDIR}"/bcftools/${SAMPLE}_${REFNAME}_calls.vcf.gz -Ob -o \
	"${WORKDIR}"/bcftools/${SAMPLE}_${REFNAME}_calls.norm.bcf

## filter adjacent indels within 5bp
printf "\n*** BCFtools - filter ***\n"
bcftools filter \
	--IndelGap 5 "${WORKDIR}"/bcftools/${SAMPLE}_${REFNAME}_calls.norm.bcf \
	-Ob -o "${WORKDIR}"/bcftools/${SAMPLE}_${REFNAME}_calls.norm.flt-indels.bcf

## index normalized, filtered variants for use in making consensus
printf "\n*** BCFtools - index ***\n"
bcftools index \
	"${WORKDIR}"/bcftools/${SAMPLE}_${REFNAME}_calls.norm.flt-indels.bcf

## generate consensus file using reference file and variant call file (vcf/bcf)
printf "\n*** BCFtools - consensus ***\n"
cat "${REF}" | bcftools consensus \
	"${WORKDIR}"/bcftools/${SAMPLE}_${REFNAME}_calls.norm.flt-indels.bcf \
	> "${WORKDIR}"/${SAMPLE}_${REFNAME}_consensus.fasta

#"\n*** Calculating Nucleotide Alignment Identity ***\n"	
## generate the truncated version of this file, based on the first and last 
## positions listed in the read depth file (I made a python program to do this)
printf "\n*** Python - trimming consensus ***\n"
cd "${S_DIR}"
chmod a+x assemblyRange.py
./assemblyRange.py \
	"${WORKDIR}"/${SAMPLE}_${REFNAME}_read_depth.txt \
	"${WORKDIR}"/${SAMPLE}_${REFNAME}_consensus.fasta \
	"${WORKDIR}"/${SAMPLE}_${REFNAME}_consensus_trunc.fasta
rm "${WORKDIR}"/${SAMPLE}_${REFNAME}_consensus.fasta # I don't want this file available, just in case

## generate a reference sequence matching the lengths of the truncated file (for id calc)
./assemblyRange.py \
	"${WORKDIR}"/${SAMPLE}_${REFNAME}_read_depth.txt \
	"${REF}" \
	"${WORKDIR}"/${SAMPLE}_${REFNAME}_ref_trunc.fasta

## produce a pairwise alignment file using muscle
printf "\n*** Creating Pairwise Alignment File ***\n"
awk '{print}' "${REF}" "${WORKDIR}"/${SAMPLE}_${REFNAME}_consensus_trunc.fasta \
	> "${WORKDIR}"/${SAMPLE}_${REFNAME}_muscle_input.fasta
muscle \
	-in "${WORKDIR}"/${SAMPLE}_${REFNAME}_muscle_input.fasta \
	-out "${WORKDIR}"/${SAMPLE}_${REFNAME}_aligned.fasta \
	&> "${WORKDIR}"/troubleshooting/muscle_output.txt
	
## print percent identity of alignment
printf "\n*** Alignment Info ***\n"
printf "\n\tSample: ${SAMPLE}\n"
printf "\n\tReference: ${REFNAME}\n"
printf "\n\tPercent of reads mapped:   "
sed -n 5p "${WORKDIR}"/troubleshooting/samtools_flagstat.txt
printf "\n"
cd "${S_DIR}"
chmod a+x percent_id.py
NTNUM=$(wc -l "${WORKDIR}"/${SAMPLE}_${REFNAME}_read_depth.txt | awk '{ print $1 }')
python3 percent_id.py \
	"${WORKDIR}"/${SAMPLE}_${REFNAME}_consensus_trunc.fasta \
	"${WORKDIR}"/${SAMPLE}_${REFNAME}_ref_trunc.fasta \
	"${REF}" ${NTNUM}

############################## FILE CLEANUP ###################################

# Remove files you don't care about for now unless testing
cd "${WORKDIR}"
rm "${WORKDIR}"/${SAMPLE}_${REFNAME}_ref_trunc.fasta 
rm "${WORKDIR}"/${SAMPLE}_${REFNAME}_muscle_input.fasta 
rm -r "${WORKDIR}"/bcftools
rm "${WORKDIR}"/*.bam
rm "${WORKDIR}"/*.bam.bai
rm "${WORKDIR}"/*.sam

printf "\n\tAll done!\n\n"
printf "*********************************************************************\n"
