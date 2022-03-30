
Figures for dsRNA Survey Paper
===============================

## Data

Main spreadsheet includes merged data for killing phenotype, alternative strain names, plate locations, toxin types (if available, from previous publications or our own NGS data), collection (FGSC, Liti, etc.), results of genomic toxin BLAST searches (KHR & KHS), and agarose gel data manually examined for presence of dsRNAs. This spreadsheet is located in: 

	`killer_assays/data/killer_assays.csv`

## Workflow

Workflow for figure generation is arranged from top to bottom in R Markdown script located in:

	`figures/scripts/Figures.Rmd`

A rendered HTML report named `Figures.html` is located in the same directory. 

SVG images are located in here:

	`figures/img`

Note that the heatmaps need to be made manually within R Studio because of a niche formatting
error which prevents their rendering during knitting unless "lhei" argument is removed from plot code. 

