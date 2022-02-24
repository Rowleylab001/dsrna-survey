dsRNA Survey Project
====================

> Current Manuscript: Google Doc [Link](https://docs.google.com/document/d/12VYn31xDQXrQYVJ30M-PAtjWb1uwu5kt/edit?usp=sharing&ouid=106641077409198019548&rtpof=true&sd=true)

> GitHub [Link]()

> Local Folder:

	$ pwd
	/Users/angela/rowley/dsrna-survey

> Conda Environment: dsrna-survey.yaml

	$ conda activate dsrna-survey
	
-----------------------------------------------------------------

**************** 2022-02-23 ****************

## Setup new folder and file system

	$ mkdir /Users/angela/rowley/dsrna-survey
	$ cd /Users/angela/rowley/dsrna-survey
	$ mkdir data test_data scripts output img

## examined existing conda environments

	$ conda info --env
	
## Setting up conda environment

	$ conda create --name dsrna-survey python=3.8
	$ conda activate dsrna-survey

## Exported current environment to YAML file

	$ conda env export --name dsrna-survey > dsrna-survey.yaml
	
## updated conda

	$ conda update -n base -c defaults conda
	
## setup folders for scripts

	$ cd scripts/
	$ mkdir dsrna-seq chrom-seq figures supp private
	
## 




