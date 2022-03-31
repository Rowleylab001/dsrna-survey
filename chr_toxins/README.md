
Analysis of KHR and KHS Toxins in Liti & FGSC Collections (WGS)
===============================================================


## Data

*FGSC Collection*

- `genomic_toxins/data/FGSC_chrom/Strope_Table_S19_AccNo.xlsx` - Contains accession numbers for all FGSC strains assembled by Strope et al. These accession numbers are available in the supplementary file "Table _S19". Strope assembled and submitted complete chromosomes under these accessions. 

- `genomic_toxins/scripts/fgsc_fasta_download.py` - Run script to automatically download Chromosomes 5 and 9 from FGSC collection for all strains in the collection. These are the chromosomes carrying KHR and KHS. There should be a total of 186 FASTA files after completion. You will prompted to enter your own email (so the Entrez database server can identify you) and you must have Biopython installed. The files will be deposited in your current working directory. To run script, use the following command: 

		$ python3 fgsc_fasta_download.py

*Liti Collection*

- `genomic_toxins/data/Liti_contigs/Liti1000 std_name to sra accession.xlsx` - Matches strain names to SRA accession numbers (study id# ERP014555) for raw reads of all Liti strains sequenced by Peter et al. The de novo assembled contigs used in this analysis were downloaded from: http://1002genomes.u-strasbg.fr/files/ under the "1011Assemblies.tar.gz" file. Caution: the file is 3.7Gb zipped, >12Gb unzipped), which is VERY large. 

*KHR & KHS Reference Sequences*

- `genomic_toxins/data/reference_seqs/` - Contains KHR and KHS reference sequences downloaded from NCBI (

| Gene       | Seq Type   | NCBI Accession           | Link                                            |
| :--------- | :--------- | :----------------------- | :---------------------------------------------- |
| KHS1       | Protein    | EDN63163.1               | https://www.ncbi.nlm.nih.gov/protein/EDN63163.1 |
| KHR1       | Protein    | BAA00751.1 or EDN61464.1 | https://www.ncbi.nlm.nih.gov/protein/EDN61464.1 |
| KHS1       | Nucleotide | ???                      |                                                 |
| KHR1       | Nucleotide | D00905.1                 | https://www.ncbi.nlm.nih.gov/nuccore/D00905.1   |



## Workflow

#### Part 1: FGSC BLASTx

Josephine Boyer performed preliminary analysis of KHS & KHR toxins in FGSC collection WGS data. 

1. Download the FGSC chromosomes 5 and 9 using `fgsc_fasta_download.py` detailed above. 

2. Perform a local BLASTx using the protein sequences as the database. BLASTx results were manually filtered based on sequence similarity to canonical sequences. 

#### Part 2: Liti BLASTn

Analysis on the Liti collection was performed mainly by Mason Shipley and further developed/analyzed by Angela Crabtree. 

1. BLASTn on Liti contigs using a nucleotide database of KHR & KHS. KHS was identified from some pre-existing publication and the KHR sequence was the ORF within [D00905.1](https://www.ncbi.nlm.nih.gov/nuccore/D00905.1)

> BLAST nucleotide database: <button onclick="window.open('seqs/chromTox_nt_db.fasta')">View</button> 

> `liti_blastn.sh` script to perform BLASTn on single contig file: <button onclick="window.open('scripts/Mason/liti_blastn.sh')">View</button>

> `Liti_runall.sh` script to run liti_blastn.sh script on all liti files: <button onclick="window.open('scripts/Mason/Liti_runall.sh')">View</button>

2. Translate BLASTn hits to protein sequences and give designations depending on if they're truncated, and whether the full-length proteins are canonical or mutated. These designations are listed in a column called "type" as "truncated", "canonical", or "mutant". Note that the "mutant" designation is given only if the protein is the same length as the canonical protein. KHS and KHR canonical proteins are from [EDN63163.1](https://www.ncbi.nlm.nih.gov/protein/EDN63163.1) and [BAA00751.1](https://www.ncbi.nlm.nih.gov/protein/BAA00751.1) from NCBI, respectively. 

> `KHSKHR.R` script to translate and compare sequences from BLASTn output: <button onclick="window.open('scripts/Mason/KHSKHR.R')">View</button>

3. Obtain nucleotide sequences 50 bp to either side of the BLASTn hit and identify KHR/KHS-like ORFs from within the extended sequence. The purpose of this is to ensure that the BLASTn did not exclude nucleotide sequences on either end due to sequence dissimilarity. 

> `LitiContigPull.py` script identifies potential ORFs within given BLASTn output. Note that the BLASTn output has been modified to only include the columns for strain name, contig name, blast hit start position, and blast hit end position. This script uses the location information of blast hit start & stop to find what positions to grab the extended nucleotide sequence from. This includes an exception for when the exted sequence positions have exceeded the perimeter of the contig sequence. The output csv contains 3 columns: strain name, contig name, and extended nucleotide sequence. 

> <button onclick="window.open('scripts/Mason/LitiContigPull.py')">View</button>

> `LitiORFsWithRvs.py` script identifies potential ORFs within the extended nucleotide sequences from `LitiContigPull.py` output. When an ATG codon is found within the extended sequence, the script scans downstream for a stop codon. If it finds one, it checks to see if it is in frame and, if so, it will translate it into a protein. It repeats this for the reverse complement of the extended sequence and returns a csv with 4 columns: strain name, contig name, ORF nucleotide sequence, and translated protein sequence. 

> <button onclick="window.open('scripts/Mason/LitiORFsWithRvs.py')">View</button>

4. Merge new list of extended sequences with original BLASTn results in order to compare protein sequences. We can then see if there was any KHR or KHS genes erroneously incomplete by the BLASTn. The R script merges the two lists and makes a new column to compare whether the two protein sequences match. Mismatches were filtered and delted using Excel. 

> `RLitiORFComparison.R` script merges the two lists compares whether the protein sequences match between BLASTn and extended search: <button onclick="window.open('scripts/Mason/RLitiORFComparison.R')">View</button>



## References

Peter J, De Chiara M, Friedrich A, et al. Genome evolution across 1,011 Saccharomyces cerevisiae isolates. *Nature*. 2018;556(7701):339-344. doi:10.1038/s41586-018-0030-5

Strope PK, Skelly DA, Kozmin SG, et al. The 100-genomes strains, an S. cerevisiae resource that illuminates its natural phenotypic and genotypic variation and emergence as an opportunistic pathogen. *Genome Res*. 2015;25(5):762-774. doi:10.1101/gr.185538.114

