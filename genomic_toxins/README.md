
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

| Gene       | Seq Type   | NCBI Accession |
| :--------- | :--------- | :------------- |
| KHS1       | Protein    | EDN63163.1     |
| KHR1       | Protein    | BAA00751.1     |
| KHS1       | Nucleotide | ???     |
| KHR1       | Nucleotide | ???     |


## References

Peter J, De Chiara M, Friedrich A, et al. Genome evolution across 1,011 Saccharomyces cerevisiae isolates. *Nature*. 2018;556(7701):339-344. doi:10.1038/s41586-018-0030-5

Strope PK, Skelly DA, Kozmin SG, et al. The 100-genomes strains, an S. cerevisiae resource that illuminates its natural phenotypic and genotypic variation and emergence as an opportunistic pathogen. *Genome Res*. 2015;25(5):762-774. doi:10.1101/gr.185538.114

