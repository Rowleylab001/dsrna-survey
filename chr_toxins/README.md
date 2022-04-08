
Analysis of KHR and KHS Toxins in Liti & FGSC Collections (WGS)
===============================================================


## Data

*FGSC Collection*

- `chr_toxins/data/FGSC_chrom/Strope_Table_S19_AccNo.xlsx` - Contains accession numbers for all FGSC strains assembled by Strope et al. These accession numbers are available in the supplementary file "Table _S19". Strope assembled and submitted complete chromosomes under these accessions. 

- `chr_toxins/scripts/fgsc_fasta_download.py` - Run script to automatically download Chromosomes 5 and 9 from FGSC collection for all strains in the collection. These are the chromosomes carrying KHR and KHS. There should be a total of 186 FASTA files after completion. You will prompted to enter your own email (so the Entrez database server can identify you) and you must have Biopython installed. The files will be deposited in your current working directory. To run script, use the following command: 

		$ python3 fgsc_fasta_download.py

If the above script produces an error, it is very likely that Biopython needs to be installed. 

*Liti Collection*

- `chr_toxins/data/Liti_contigs/Liti1000 std_name to sra accession.xlsx` - Matches strain names to SRA accession numbers (study id# ERP014555) for raw reads of all Liti strains sequenced by Peter et al. The de novo assembled contigs used in this analysis were downloaded from: http://1002genomes.u-strasbg.fr/files/ under the "1011Assemblies.tar.gz" file. Caution: the file is 3.7Gb zipped, >12Gb unzipped), which is VERY large. 

*KHR & KHS Reference Sequences*

The KHS toxin sequence published by Goto et al. on NCBI (GenBank: S77712.1) was refuted by a later paper by Frank and Wolf who published an updated sequence for KHS which should fix the problems in the Goto sequence. Hence, the correct [KHS sequence](https://www.ncbi.nlm.nih.gov/nuccore/AAFW02000048.1?report=fasta&log$=seqview&from=551404&to=552456) encodes protein SCY_1690 from strain YJM789. 

> *"We believe that there are errors in the original report of the sequence of KHS1, from *S. cerevisiae* strain no. 115 (17), because it differs from all other sequences by many frameshifts and by a 1.4-kb inversion that ends at a restriction site used in cloning."* (Frank & Wolf 2009)

The KHR sequence, published by Goto et al., is the ORF within [D00905.1](https://www.ncbi.nlm.nih.gov/nuccore/D00905.1?report=fasta&log$=seqview&from=139&to=1029) as the published sequence included flanking sequences. 

- `chr_toxins/data/reference_seqs/chromTox_nt_db.fasta` - Contains KHR and KHS reference sequences downloaded from NCBI:

| Gene       | Seq Type   | NCBI Accession               | Link                                        |
| :--------- | :--------- | :--------------------------- | :------------------------------------------ |
| KHS1       | Nucleotide | AAFW02000048.1:551404-552456 | https://www.ncbi.nlm.nih.gov/nuccore/AAFW02000048.1?report=fasta&log$=seqview&from=551404&to=552456 |
| KHR1       | Nucleotide | D00905.1:139-1029            | https://www.ncbi.nlm.nih.gov/nuccore/D00905.1?report=fasta&log$=seqview&from=139&to=1029                     | 



## Workflow

Josephine Boyer performed exploratory analysis of KHS & KHR toxins in FGSC collection WGS data using local BLASTx and manual filtering. The data was later analyzed in addition to the Liti collection using the `chrtox_pipeline.sh` script. The scripts used therein were initially written by Mason Shipley and further developed by Angela Crabtree. 

	$ chmod a+x chrtox_pipeline.sh
	$ ./chrtox_pipeline.sh -h

	-------------------------------------------------------
	options:

	   -f [arg]	directory containing fasta files (required)
	   -b [arg]	BLASTn reference fasta file (required)
	   -o [arg]	destination of output folder
	   -t		test (ensures required CL apps are working)
	   -h		help
	-------------------------------------------------------

	$ ./chrtox_pipeline.sh -f [fasta folder] -b [BLASTn database fasta file] -o [output folder]


The `chrtox_pipeline.sh` script uses the following scripts and workflow:

1. [**01_chrtox_blastn.sh**](scripts/01_chrtox_blastn.sh)

	a. Script performs BLASTn using local BLAST database. 

	b. A local BLASTn was performed on files using a [nucleotide database of KHR & KHS](data/reference_seqs/chromTox_nt_db.fasta). Note that the word "KHR" or "KHS" need to be included in the header for the downstream R scripts to work properly. 

2. [**02_chr_aa_convert.R**](scripts/02_chr_aa_convert.R)

	a. Script translates and categorizes BLASTn output. 

	b. Resultant BLASTn hits were translated to protein sequences and given designations depending on if the proteins were truncated, and whether the full-length proteins were canonical or mutated. These designations were listed in a column called "type" as "truncated", "canonical", or "mutant". The "mutant" designation was given if the protein was the same length as the canonical protein and the protein sequence differed from canonical. KHS and KHR canonical proteins are translated from the nucleotide sequences in the blastn nucleotide database file `chr_toxins/data/reference_seqs/chromTox_nt_db.fasta`. 

3. [**03_extend_hits.py**](scripts/03_extend_hits.py)

	a. Script extends BLAST hit nucleotide sequences. 

	b. Nucleotide sequences 50 bp to either side of the BLASTn hit were extracted from original sequence files. The purpose of this was to ensure that the BLASTn did not exclude nucleotide sequences on either end due to sequence dissimilarity. 

4. [**04_identify_orfs.py**](scripts/04_identify_orfs.py)

	a. Script identifies the longest ORFs within extended nucleotide sequences. 

	b. The longest ORF within each extended hit sequence was found. 

5. [**05_compare_orfs.R**](scripts/05_compare_orfs.R)

	a. Script compares the blast hit protein sequences with the extended hit protein sequences. 

	b. The protein sequences from the original blast and the extended orfs were compared and designated as "Match" if the two agreed and "Mismatch" if a different protein was found in the extended search. This determines whether there were any KHR or KHS genes erroneously incomplete by the BLASTn. 


## References

Frank, A Carolin, and Kenneth H Wolfe. “Evolutionary capture of viral and plasmid DNA by yeast nuclear chromosomes.” *Eukaryotic cell* vol. 8,10 (2009): 1521-31. doi:10.1128/EC.00110-09

Peter J, De Chiara M, Friedrich A, et al. Genome evolution across 1,011 Saccharomyces cerevisiae isolates. *Nature*. 2018;556(7701):339-344. doi:10.1038/s41586-018-0030-5

Strope PK, Skelly DA, Kozmin SG, et al. The 100-genomes strains, an S. cerevisiae resource that illuminates its natural phenotypic and genotypic variation and emergence as an opportunistic pathogen. *Genome Res*. 2015;25(5):762-774. doi:10.1101/gr.185538.114



