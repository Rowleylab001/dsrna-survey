#!/usr/bin/env Rscript

main <- function() {
  ###################### LIBRARIES ###################### 
  ## install bioconductor installer, if needed
  if (!requireNamespace("BiocManager")) install.packages("BiocManager")
  ## install and load required R packages
  if(!require(Biostrings)) BiocManager::install("Biostrings", dependencies = TRUE)
  library(Biostrings)
  if(!require(reshape2)) install.packages("reshape2", dependencies = TRUE)
  library(reshape2)
  if(!require(ggplot2)) install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
  if(!require(stringr)) install.packages("stringr", dependencies = TRUE)
  library(stringr)
  if(!require(stringi)) install.packages("stringi", dependencies = TRUE)
  library(stringi)

  ###################### ARGPARSE ###################### 
  
  ## parse arguments
  args <- commandArgs(trailingOnly = TRUE)
  # blast input file
  blast_file <- args[1]
  data_dir <- gsub("(\\w+)/\\w+\\.csv", "\\1", blast_file)[1]
  setwd(data_dir)
  # blast db file
  database_file <- args[2]
  # option for output filename 
  if (is.na(args[3])==FALSE){
    output_file <- args[3]
  } else {
    output_file <- paste0(data_dir, '/liti_chr_categories.csv')
  }
  
  ###################### MAIN ###################### 
  ## Read data 
  hits.df <- read.csv(blast_file, header=T)
  # delete any spaces added into nt seq by BLAST
  hits.df$qseq_nt <- lapply(hits.df$qseq_nt, gsub, pattern = "-", replacement = "", fixed = TRUE)
  # create new column to store converted aa seq, reverse if necessary
  hits.df$aaseq = as.character(hits.df$qseq_nt)
  # translate from nt to aa (and also make sure all nt and aa seqs are 5'->3')
  for (i in 1:as.numeric(nrow(hits.df))){
    if (hits.df$send[i] < hits.df$sstart[i]){
      hits.df$qseq_nt[i] = as.character(
        Biostrings::reverseComplement(DNAString(as.character(hits.df$qseq_nt[i]))))
    }
    hits.df$aaseq[i] = as.character(
      Biostrings::translate(DNAString(as.character(hits.df$qseq_nt[i]))))
  }
  # make a column to list all full-length protein sequences (multiple hits)
  # it determines this by there being a "*" at end of seq, representing stop codon
  hits.df$aalst = stri_extract_all(hits.df$aaseq, regex="^[A-Z]++\\*$")
  # truncate aa seq at internal stop codon
  hits.df$aaseq = str_split_fixed(hits.df$aaseq, '\\*', 2)[,1] 
  # creates gene column that identifies if a hit is KHS/KHR from description
  hits.df$gene = gsub(".*(KHR*S*).*", "\\1", hits.df$match_description)
  hits.df$gene = as.factor(hits.df$gene)
  # add column for aa length
  hits.df$aalen=nchar(hits.df$aaseq)
  # add column to display whether aa is full-length
  for (i in 1:as.numeric(nrow(hits.df))){
    if (hits.df$gene[i] == "KHS"){
      if (hits.df$aalen[i] == 350){
        hits.df$full[i] = "yes"
      } else {
        hits.df$full[i] = "no"
      }
    }
    if (hits.df$gene[i] == "KHR"){
      if (hits.df$aalen[i] == 296){
        hits.df$full[i] = "yes"
      } else {
        hits.df$full[i] = "no"
      }
    }
  }

  # function to convert blast db file into dataframe
  fasta_to_df = function(fasta_file){ 
    fasta_biostring <- readDNAStringSet(fasta_file)
    faheader=names(fasta_biostring)
    seq = paste(fasta_biostring)
    df <- data.frame(faheader, seq)
    df$acc = gsub("(\\S+)[.].*", "\\1", df$faheader)
    df$gene = gsub(".*(KHR*S*).*", "\\1", df$faheader)
    for (i in 1:as.numeric(nrow(df))){
      df$prot[i] = as.character(
        Biostrings::translate(DNAString(as.character(df$seq[i]))))
    }
    df$prot=str_split_fixed(df$prot, '\\*', 2)[,1]
    return(df)
  }
  # define KHR & KHS
  fa_bio.df = fasta_to_df(database_file)
  KHR=as.character(na.omit(ifelse(fa_bio.df$gene=="KHR",fa_bio.df$prot,NA)))
  KHS=as.character(na.omit(ifelse(fa_bio.df$gene=="KHS",fa_bio.df$prot,NA)))

  ## define what criteria defines categories
  hits.df$type[hits.df$aaseq == KHR | hits.df$aaseq == KHS] = "canonical"
  hits.df$type[hits.df$aaseq != KHR & hits.df$aaseq != KHS] = "mutant"
  hits.df$type[hits.df$full == "no"] = "truncated"
  # prepare data to collapse all identical proteins for MSA data & plot
  hits.df$num_hits = 1
  hits.df.collapsed = aggregate(num_hits ~ gene+type, data=hits.df, FUN=sum)
  hits.df.msa = aggregate(num_hits ~ gene+type+aaseq, data=hits.df, FUN=sum)
  hits.df.msa$type=as.character(hits.df.msa$type)
  hits.df.collapsed.khs.mut = hits.df.msa[which(hits.df.msa$type=="mutant" & hits.df.msa$gene=="KHS"),]
  hits.df.collapsed.khr.mut = hits.df.msa[which(hits.df.msa$type=="mutant" & hits.df.msa$gene=="KHR"),]
  
  ## plot hits 
  myPalette <- c("#377EB8", "#4DAF4A", "#FA3232")
  p = ggplot(data=hits.df.collapsed, aes(x=gene, y=num_hits, fill=type)) +
    geom_bar(stat="identity", position="stack", color="grey")  +
    theme_classic(base_size = 28) +
    labs(y="Strain Count", x="") +
    theme(axis.text.x = element_text(angle=315, hjust=0.2, vjust=0.5),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major.y = element_line(
            size = 0.5, linetype = 'solid', colour = "grey"),
          plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm"),
          legend.position = "right", legend.title=element_text(size=22),
          legend.text=element_text(size=22),
          legend.background = element_rect(
            fill="white", size=0.5, linetype="solid", colour ="black")) +
    scale_fill_manual("Protein Type", values = myPalette) +
    ylim(0, length(unique(hits.df$filename)))
  # save plot to SVG file
  img_file = paste0(data_dir, "/Liti_blastn_chr_hits.svg") # filename
  ggsave(
    img_file,
    plot = p,
    width = 5, height = 4)
  
  ## write data to output files
  hits.df$full = as.character(hits.df$full)
  hits.df$aalst = as.character(hits.df$aalst)
  hits.df$qseq_nt = as.character(hits.df$qseq_nt)
  write.csv(hits.df, file=output_file, row.names = FALSE)
  
  ## write files for MSA to determine KHR/KHS protein functionality
  hits.df.msa.khs = hits.df.msa[which(hits.df.msa$gene=="KHS"),]
  write.csv(hits.df.msa.khs, file=paste0(data_dir,"/Liti-KHS-msa.csv"), 
            row.names = FALSE)
  hits.df.msa.khr = hits.df.msa[which(hits.df.msa$gene=="KHR"),]
  write.csv(hits.df.msa.khs, file=paste0(data_dir,"/Liti-KHR-msa.csv"), 
            row.names = FALSE)
  }
main()
