
dsRNA Survey Project
====================

- Current Manuscript: Google Doc [Link](https://docs.google.com/document/d/12VYn31xDQXrQYVJ30M-PAtjWb1uwu5kt/edit?usp=sharing&ouid=106641077409198019548&rtpof=true&sd=true)

- GitHub [Link]()

- Local Folder:

	$ pwd
	/Users/angela/rowley/dsrna-survey

- Conda Environment: dsrna-survey.yaml

	$ conda activate dsrna-survey
	
---------------------------------------------------------------------------------------------------------

**************** 2022-02-23 ****************

- Setup new folder and file system

	$ mkdir /Users/angela/rowley/dsrna-survey
	$ cd /Users/angela/rowley/dsrna-survey
	$ mkdir data test_data scripts output img

- examined existing conda environments

	$ conda info --env
	
- Setting up conda environment

	$ conda create --name dsrna-survey python=3.8
	$ conda activate dsrna-survey

- Exported current environment to YAML file

	$ conda env export --name dsrna-survey > dsrna-survey.yaml
	
- updated conda

	$ conda update -n base -c defaults conda
	
- setup folders for scripts

	$ cd scripts/
	$ mkdir dsrna-seq chrom-seq figures supp private
	
**************** 2022-03-28 ****************

Decided to take a lower-effort approach to tying up loose ends and just put 
all the R scripts to generate tables and graphs together in one Rmd

Discovered the numbers were off when recreating Table 1 (the summary of dsRNA/killer prevalence).
This is likely due to the fact that these exact permutations were not requested and instead 
were summed post hoc. The data set has also been altered since the original numbers were calculated. 

**************** 2022-03-29 ****************

#### Bubble Plots

Worked on Fig. 3ab (bubble plots). They're looking awesome now!
However, the figures in the manuscript need to be changed because there is a 357-strain point
when Y-2046 is included. This is presumably strains that only kill Y-2046, since 
643 - 357 = 286 strains. 

* All strains:

	> max(pcaDF.agg.combo$num_strains)	[1] 357

	> sum(pcaDF.agg.combo$num_strains)	[1] 643

* No Y-2046:

	> max(pcaDF.agg.combo$num_strains)	[1] 58

	> sum(pcaDF.agg.combo$num_strains)	[1] 286

#### Genome-encoded toxins

??? in Table 3 including strains that have either KHR or KHS data (should this be BOTH?)

- 1 strain off original on Table 3 -KHS1/-KHR1 killers (not too bad!)

??? Do we need to update our raw data sheet for public use? (if so, I need to redo Rmd)

>> Next time: work on supplementary figures



