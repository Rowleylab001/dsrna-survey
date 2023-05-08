#!/usr/bin/env Rscript

main <- function() {
  ###################### LIBRARIES ###################### 
  ## install bioconductor installer, if needed
  if (!requireNamespace("BiocManager")) install.packages("BiocManager")
  ## install and load required R packages
  if(!require(Biostrings)) BiocManager::install("Biostrings", dependencies = TRUE)
  library(Biostrings)
  if(!require(msa)) BiocManager::install("msa", dependencies = TRUE)
  library(msa)
  if(!require(readr)) install.packages("readr", dependencies = TRUE)
  library(readr)
  if(!require(seqinr)) install.packages("seqinr", dependencies = TRUE)
  library(seqinr)
  if(!require(ape)) BiocManager::install("ape", dependencies = TRUE)
  library(ape)
  
  ###################### ARGPARSE ###################### 
  
  ## parse arguments
  args <- commandArgs(trailingOnly = TRUE)
  ## any fasta file containing more than one sequence
  fasta_file <- args[1]  
  data_dir <- dirname(fasta_file)
  setwd(data_dir)
  # option for output filename 
  if (is.na(args[2])==FALSE){
    output_file_nwk = args[2]
  } else {
    output_file_nwk <- paste0(data_dir, "/", 
      unlist(strsplit(basename(fasta_file), "[.]"))[1], "_tree.nwk") 
  }
  # option for gene name
  if (is.na(args[3])==FALSE){
    output_str <- args[3]
  } else {
    output_str <- ""
  }
  
  ###################### FUNCTIONS ###################### 
  
  # determine if fasta file contains DNA sequences (if FALSE, then it's protein)
  is_fasta_dna = function(filename){
    ## import blast output file first sequence
    fseq=toupper(as.character(readr::read_lines(fasta_file, n_max = 2)[2]))
    is_dna = grepl("^[ATGCN]*$", fseq)
    return(is_dna)
  }
  
  # save pdf of MSA visualization
  save_msa_viz = function(msa_object, outfile){
    msaPrettyPrint(toxinAln, output="pdf", 
                   file=msa_viz_file, 
                   showNames="left", 
                   shadingMode="similar", shadingColors="blues", 
                   showLogo="none", askForOverwrite=FALSE)
    dev.off()  # close the file
    print("Wrote MSA visualization to file: \n   ")
    print(outfile)
  }
  
  save_tree_data = function(msa_object, outfile_nwk, tree_img_file=""){
    # create phylogenetic tree based on neighbor-joining
    msa_object2 <- msaConvert(msa_object, type="seqinr::alignment") # -> seqinr format
    d <- dist.alignment(msa_object2, "identity")
    myTree <- nj(d) # construct tree 
    # write tree info to Newick file
    ape::write.tree(myTree, file=outfile_nwk)
    print("Wrote Tree info to Newick file: \n   ")
    print(outfile_nwk)
    
    # save dendrogram image
    if (tree_img_file!=""){ # if tree image file is specified, write file
      pdf(tree_img_file)
      plot(myTree, main=paste("Phylogenetic Tree (Neighbor Joining)"))
      dev.off()  # close the file
      print("Wrote Dendrogram to file: \n   ")
      print(tree_img_file)
    }
  }

  ###################### MAIN ###################### 
  
  if (is_fasta_dna(fasta_file)==FALSE){ # if the sequence is a protein sequences
    
    # perform MSA
    mySequences <- readAAStringSet(fasta_file)
    toxinAln=msa(mySequences)
    
    # save MSA visualization
    msa_viz_file=paste0(data_dir,"/",output_str, "msa_aa.pdf")
    save_msa_viz(toxinAln, msa_viz_file)
    
    # create phylogenetic tree and write tree info and image to files
    # (tree is made by neighbor-joining method)
    output_dendrogram = paste0(dirname(outfile_nwk), "tree_aa.pdf") 
    save_tree_data(toxinAln, output_file_nwk, output_dendrogram)
    
  } else { # the file contains DNA sequences
    
    # perform MSA
    mySequences <- readDNAStringSet(fasta_file)
    toxinAln=msa(mySequences)
    
    # save MSA visualization
    msa_viz_file=paste0(data_dir,"/",output_str, "msa_nt.pdf")
    save_msa_viz(toxinAln, msa_viz_file)
    
    # create phylogenetic tree and write tree info and image to files
    # (tree is made by neighbor-joining method)
    output_dendrogram = paste0(dirname(outfile_nwk), "tree_nt.pdf") 
    save_tree_data(toxinAln, output_file_nwk, output_dendrogram)
  }
}
main()
