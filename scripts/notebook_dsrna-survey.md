### Project: The Paper Formerly Known as "MDPI"
Analysis of KHR and KHS Toxins in Liti & FGSC Collections (WGS)
======================================================

## Objectives

1. Identify the chromosomal toxins, KHS and KHR, within the Liti and FGSC collections

2. Determine whether there is a statistical relationship between functional KHS/KHR and the killing of hypersensitive lawns (Y-2046 & ATCC2001) 

## Methods

#### *Sequence Assembly of dsRNA Totiviruses & Satellites*

Raw Illumina reads were cleaned and filtered with fastp (v. ___), genomic contaminants (NCBI accessions: ___ for S. cerevisiae, ___ for S. paradoxus, ___ for H. sapiens) were removed using bbsplit (samtools v. ____), and de novo contigs were assembled with SPAdes assembler (v. ____). The custom script is available on GitHub (insert github repo link). Resultant contigs were screened for totivirus and satellite sequences with local BLASTn (v. ____) using a database of published totiviruses and satellites (Table S___) downloaded from NCBI's nucleotide database (accessed on XX Month 202X). 

#### *Comparative Analysis of dsRNA Totiviruses & Satellites among Sequenced Strains*



#### *Comparative Analysis of KHR1 and KHS1 among Screened Collections*

FGSC genome assemblies were downloaded from NCBI's Entrez database (accessed on XX Month 202X) using accession numbers listed in the Strope et al. 2015 Supplemental [Table S19](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4417123/bin/supp_gr.185538.114_Table_S19_AccNo.xlsx)(Strope et al. 2015). 

Liti Collection genome assemblies were downloaded from http://1002genomes.u-strasbg.fr/files/1011Assemblies.tar.gz (accessed on XX Month 202X)(Peter et al. 2018). 

#### *Analysis of Killer Phenotype Data*

Heatmaps were constructed using custom R scripts. 







## References (AMA style)

Peter J, De Chiara M, Friedrich A, et al. Genome evolution across 1,011 Saccharomyces cerevisiae isolates. *Nature.* 2018;556(7701):339-344. doi:10.1038/s41586-018-0030-5

Strope PK, Skelly DA, Kozmin SG, et al. The 100-genomes strains, an S. cerevisiae resource that illuminates its natural phenotypic and genotypic variation and emergence as an opportunistic pathogen. *Genome Res.* 2015;25(5):762-774. doi:10.1101/gr.185538.114



-----------------------------------------------------------------------------------------

**Josie's Workflow** 

> She performed preliminary analysis of KHS & KHR toxins in FGSC collection WGS data. She accomplished this by downloading the FGSC genomes and performing a local BLASTx using the protein sequences as the database, and later manually filtered the samples based on sequence similarity to canonical sequences. 

**Mason's workflow:** 

1. BLASTn on Liti contigs using a nucleotide database of KHR & KHS. KHS was identified by Paul from some pre-existing publication? and the KHR sequence was the ORF within [D00905.1](https://www.ncbi.nlm.nih.gov/nuccore/D00905.1)

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

* * * *

## Data

*Online*

* [Liti Raw Reads on SRA](https://trace.ncbi.nlm.nih.gov/Traces/sra/?study=ERP014555) - study id# ERP014555

* [1011Assemblies.tar.gz](http://1002genomes.u-strasbg.fr/files/1011Assemblies.tar.gz) - de novo assemblies of the Liti1011 collection (caution: 3.7Gb zipped, >12Gb unzipped)

*Located in* `data/`

* `KHR_ref.fasta` & `KHS_ref.fasta` - KHR & KHS canonical protein seqs

* `KHR_liti_raw.fasta` & `KHS_liti_raw.fasta` - KHR & KHS protein seqs from Mason after his data cleaning process

* `KHR_liti_all_full.fasta` & `KHS_liti_all_full.fasta` - contains fastas of only full-length KHR & KHS proteins (for mutational analysis)

**Metadata**

        data/meta/
        
* `FGSC_accession.xlsx` - Lists accession numbers for FGSC collection genomes (available on NCBI)

* `FGSC_acc_to_strain_name.csv` - Pairs FGSC accession numbers to strain name and chromosome (for Chr5 and Chr9 only, since that's where KHR & KHS are located)

* `Josie_FGSC_KHR_KHS_summary_edited.xlsx` - Josie's data output from her analysis of KHR/KHS presence in FGSC collection (may not be updated)

* `Liti1000 std_name to sra accession.xlsx` - Pairs Liti accession numbers to strain name and some other info

**Analysis**

        data/analysis/

* `mutation_analysis_KHR_full.csv` & `mutation_analysis_KHR_full.csv` - output of msa_varcall.py script, calling all amino acid mutations between each Liti fasta and the canonical reference sequence

* * * *

## Scripts

        scripts/

* `khs_khr_filt.R` - Used for filtering csv-converted fasta sequences based on length of seqs; should probably add on a section for formatting straight from fasta and converting back to fasta after filtering (it was annoying to regex it in Notepad++)

* `lollipop_prot.R` - Used for generating lollipop figures based on output from `msa_varcall.py` script

* `merge_khs_khr.R` - Used for merging KHS & KHR data to existing killer data (needs to be performed on updated KHR/KHS data)

* `msa_varcall.py` - Used for calling mutations between a reference sequence and a file contaning 1+ fasta-formatted (potentially) mutated sequences; this script takes 1 fasta-formatted reference file and 1 fasta-formatted mutant file as input and outputs a csv containing a row for each mutation with columns for id/accession of mutant, reference nt/aa, position of mutation, and mutated nt/aa
        
## Notes

- Lots of mutations in KHR & KHS. Not many truncations that look functional (i.e. few truncations longer than 50 or 100 aa for KHR and KHS, respectively. 

| Metric                               | KHR    | KHS    |
| :----------------------------------- | :----- | :----- |
| canonical length                     | 296 aa | 350 aa |
| total #sequences curated from liti\* | 801    | 1142   |
| sequences > 50 aa                    | 590    | 1131   |
| sequences > 100 aa                   | n/a    | 635    |
| sequences = 296 aa                   | 582    | 619    |

> *\*Liti sequences were curated by Mason using the workflow detailed in the "Methods" section, above.*

![](img/KHR_mutations.jpeg){ width=80% }

![](img/KHS_mutations.jpeg){ width=80% }
