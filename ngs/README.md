
Assembly and Identification of Killer Toxins from NGS
======================================================

## NGS Data

Raw FASTQ sequences for this study can be found in the SRA under the following accessions: 

| Strain         | SRA Accession* | 
| :------------- | :------------ | 
| YJM1077        | SRXxxxxxxx    | 
| YJM1287        | SRXxxxxxxx    | 
| YJM1307        | SRXxxxxxxx    | 
| YJM1621        | SRXxxxxxxx    | 
| YO1482         | SRXxxxxxxx    | 
| YO1490         | SRXxxxxxxx    | 
| YO1619         | SRXxxxxxxx    | 
| YO1622         | SRXxxxxxxx    | 

*"SRXxxxxxxx" denotes that this data has yet to be submitted. 



## Workflow

1. Use `ngs/scripts/ycDeNovo/ycDeNovo.sh` script to assemble de novo contigs from FASTQ files (paired end sequencing required). For usage, refer to `ngs/scripts/ycDeNovo/README.docx`. 

2. Perform BLASTn using top contigs to determine potential toxin candidates to use as reference sequences for the next step. 

3. Use `ngs/scripts/ycConsensus/ycConsensus.sh` script with suspected toxin reference sequences to assess whether the toxin is likely to be encoded by the sequenced satellite dsRNA. For script usage, refer to `ngs/scripts/ycConsensus/README.docx`. If the reads mapped with relatively high coverage across the majority of the reference sequence, it was considered to be a dsRNA element of the yeast. 




