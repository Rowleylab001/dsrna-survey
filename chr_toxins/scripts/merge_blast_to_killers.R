#!/usr/bin/env Rscript

main <- function() {
  ###################### LIBRARIES ###################### 
  ## install bioconductor installer, if needed
  if (!requireNamespace("BiocManager")) install.packages("BiocManager")
  ## install and load required R packages
  if(!require(Biostrings)) BiocManager::install("Biostrings", dependencies = TRUE)
  library(Biostrings)
  if(!require(ggplot2)) install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
  if(!require(tidyr)) install.packages("tidyr", dependencies = TRUE)
  library(tidyr)
  library(rstudioapi)
  
  ###################### FILES ######################
  
  script.dir = dirname(rstudioapi::getSourceEditorContext()$path)
  fgsc_names_file = paste0(script.dir,"/../data/reference_info/FGSC_accession.csv")
  fgsc_hits_file = "/Users/angela/rowley/dsrna-survey/chr_toxins/output/FGSC/02_chr_categories.csv"
  liti_hits_file = "/Users/angela/rowley/dsrna-survey/chr_toxins/output/Liti/02_chr_categories.csv"
  killer_file = "/Users/angela/rowley/dsrna-survey/killer_assays/data/killer_assays.csv"
  
  output_file = paste0(dirname(killer_file),'/killer_assays_new.csv')
  
  ###################### MAIN ###################### 
  
  # import and reformat the fgsc strain name reference sheet
  fgsc_names = read.csv(fgsc_names_file, header=FALSE)
  colnames(fgsc_names) = as.character(fgsc_names[2,])
  fgsc_names = fgsc_names[-c(1:2),]
  fgsc_names = fgsc_names %>% 
    tidyr::pivot_longer(cols=Chr01:`2micron`, names_to="Chr", values_to="accession")
  # get unique_name for fgsc collection
  processed_hits_fgsc = read.csv(fgsc_hits_file)
  processed_hits_fgsc$filestrain=gsub(
    "([a-zA-Z0-9]+).*\\..+","\\1", processed_hits_fgsc$filename)
  processed_hits_fgsc$unique_name <- fgsc_names$Strain[match(
    processed_hits_fgsc$filestrain, fgsc_names$accession)]
  processed_hits_fgsc = subset (processed_hits_fgsc, select = -filestrain)
  
  # get unique_name for liti collection
  processed_hits_liti = read.csv(liti_hits_file)
  processed_hits_liti$unique_name=gsub(
    "([a-zA-Z0-9]+).*\\..+","\\1", processed_hits_liti$filename)
  
  # combine liti and fgsc information
  processed_hits = rbind(processed_hits_liti, processed_hits_fgsc)
  # sort lists based on percent aa identity to reference seqs, then take highest score
  # (should represent most canonical genomic toxin for that strain)
  processed_hits$gene = gsub(
    ".+ - (\\w+) - .+", "\\1", processed_hits$match_description)
  data_ordered <- processed_hits[order(processed_hits$best_aapid_ref, decreasing = TRUE), ] 
  data_ordered_khs <- data_ordered[which(data_ordered$gene=='KHS'), ] 
  best_khs_by_strain <- data_ordered_khs[!duplicated(data_ordered_khs$unique_name), ] 
  data_ordered_khr <- data_ordered[which(data_ordered$gene=='KHR'), ] 
  best_khr_by_strain <- data_ordered_khr[!duplicated(data_ordered_khr$unique_name), ]
  best_toxin_by_strain = rbind(best_khs_by_strain, best_khr_by_strain)
  # only retain subsetted data (khr_type, khr_prot, khs_type, khs_prot)
  best_toxin_by_strain = subset (
    best_toxin_by_strain, 
    select = c(unique_name, gene, mutation_type, best_aaseq))
  best_toxin_by_strain = best_toxin_by_strain %>% 
    tidyr::pivot_wider(names_from = gene, values_from = c(mutation_type, best_aaseq))
  colnames(best_toxin_by_strain)=c("unique_name","khs_type","khr_type","khs_prot","khr_prot")
  best_toxin_by_strain <- best_toxin_by_strain[, c("unique_name","khr_type","khr_prot","khs_type","khs_prot")]
  
  # import killer assay data and make row for unique_name
  killer_df = read.csv(killer_file)
  unique_name_v = as.character()
  for (row in 1:nrow(killer_df)) {
    if (is.na(killer_df$std_name[row])){ # if there's no std_name (not Liti)
      unique_name_v = append(unique_name_v, killer_df$strain_name[row])
    } else {
      unique_name_v = append(unique_name_v, killer_df$std_name[row])
    }
  }
  killer_df$unique_name = unique_name_v
  killer_df = subset (killer_df, select = -c(khr_type,khr_prot,khs_type,khs_prot))
  
  # merge with killer data
  full_killer_blast_df = merge(killer_df, best_toxin_by_strain, 
                               all.x = TRUE, by="unique_name")
  full_killer_blast_df = subset (full_killer_blast_df, select = -unique_name)
  # write new killer data to file
  write.csv(full_killer_blast_df, file=output_file, row.names = FALSE)
}
