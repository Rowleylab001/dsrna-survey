## download required packages, if necessary
if (!requireNamespace("RColorBrewer")) install.packages("RColorBrewer")
if (!requireNamespace("gplots")) install.packages("gplots")

## import libraries
library(RColorBrewer)
library(gplots) #functions: heatmap.2
library(ggplot2) #functions: ggsave

######################## DATA TRANSFORM FUNCTIONS ######################## 

KillerList = function(mx){
  
}

########################### PLOTTING FUNCTIONS ########################### 
canonical_heatmap = function(df, remove_lawns){
  #' Plots canonical heatmap
  #' 
  #' Creates a heatmap of strains containing canonical toxins
  #' 
  #' @param df Dataframe with killer strains information
  #' @param remove Vector of names of lawns to remove
  #' @usage canonical_heatmap(df)
  
  # specify which columns contain lawn killing data
  lawn_columns = 9:(21-length(remove_lawns))
  # retain only rows with canonical toxins
  canonical.df <- df[which(df$strain_name=="BJH001" |    
                             df$strain_name=="NCYC1001" |    
                             df$strain_name=="MS300c" |    
                             df$strain_name=="DSM 70459" ),]
  ## create numeric matrix of killer data only
  killerMatrix <- as.matrix(canonical.df[, lawn_columns])
  killerMatrix <- matrix(as.numeric(killerMatrix), # Convert to numeric matrix
                         ncol = ncol(killerMatrix))
  ## rename rows and columns
  killer_toxin_type <- canonical.df$toxin
  lawn_strains <- colnames(canonical.df[, lawn_columns])
  rownames(killerMatrix) <- killer_toxin_type
  colnames(killerMatrix) <- lawn_strains
  myColors = c('oldlace','steelblue')
  heatmap.2(killerMatrix, dendrogram = 'none', Rowv = TRUE, Colv = FALSE,
            margins = c(8, 6), lwid = c(0.2,5), lhei = c(0.2,5), 
            trace = 'none', key = FALSE, col = myColors,
            xlab='Lawns', ylab="Killers", colsep = 0:ncol(killerMatrix), 
            rowsep=0:nrow(killerMatrix),
            sepcolor='white', sepwidth=c(0.02, 0.02), 
            main="Toxin Controls")
}

gg_canonical_heatmap = function(df){
  #' Plots canonical heatmap
  #' 
  #' Creates a heatmap of strains containing canonical toxins
  #' 
  #' @param df Dataframe with killer strains information
  #' @usage gg_canonical_heatmap(df)
  
  # retain only rows with canonical toxins
  canonical.df <- df[which(df$strain_name=="BJH001" |    
                                    df$strain_name=="NCYC1001" |    
                                    df$strain_name=="MS300c" |    
                                    df$strain_name=="DSM 70459" ),]
  # killer toxin in last column and lawn data in other columns
  toxin = canonical.df$toxin
  canonical.df = canonical.df[, lawn_columns]
  canonical.df$toxin = toxin
  # put all killer phenotype values in one column
  canonical.df.long = canonical.df %>% pivot_longer(cols = -toxin)
  # formatting for plot
  color_palette = c('oldlace','steelblue')
  label_text = c('no effect', 'killing phenotype')
  theme_change <- theme(
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 90, hjust=1)
  )
  # plot
  ggplot(canonical.df.long, aes(name, toxin, fill=factor(value))) + 
    geom_tile(color = "black") +
    scale_fill_manual(values = color_palette, name = "", labels = label_text) +
    labs(main="Canonical Toxins", x="Lawn", y="Canonical Toxins") +
    theme_change
}

collapsed_heatmap(killer.df) = function(df){
  #' Plots collapsed heatmap
  #' 
  #' Creates a heatmap of killers screened, collapsing those with the same pattern
  #' 
  #' @param df Dataframe with killer strains information
  #' @usage collapsed_heatmap(df)
  
  killer.df = subset(killer.df, select = -c(Y2046) )
  killer.df = subset(killer.df, select = -c(ATCC2001) )
  killer.df <- killer.df[which(killer.df$Y1088!=0 |  
                                 killer.df$DSM70459!=0 |
                                 killer.df$CYC1058!=0 |  
                                 killer.df$Y27106!=0 | 
                                 killer.df$NCYC1006!=0 |
                                 killer.df$NCYC777!=0 | 
                                 killer.df$NBRC1802!=0 | 
                                 killer.df$NBRC1815!=0 |
                                 killer.df$Y5509!=0 | 
                                 killer.df$MS300c!=0 |
                                 killer.df$BY4741!=0),]
  Toxin <- killer.df$toxin
  Strain <- killer.df$strain_name
  killer.df$num_strains = 1
  killer.df.collapsed = aggregate(num_strains ~ Y1088+DSM70459+CYC1058+Y27106+NCYC1006+NCYC777+NBRC1802+NBRC1815+Y5509+MS300c+BY4741, 
                                  data=killer.df, FUN=sum)
  killer.df.agg.strain = aggregate(strain_name ~ Y1088+DSM70459+CYC1058+Y27106+NCYC1006+NCYC777+NBRC1802+NBRC1815+Y5509+MS300c+BY4741, 
                                   data=killer.df, paste, collapse=",")
  killer.df.collapsed$strain_name = killer.df.agg.strain$strain_name
  killer.df.collapsed$group = seq(1:nrow(killer.df.collapsed))
  K1 = as.numeric(grep('BJH001', killer.df.collapsed$strain_name))
  K2 = as.numeric(grep('NCYC1001', killer.df.collapsed$strain_name))
  K28 = as.numeric(grep('MS300c', killer.df.collapsed$strain_name))
  Klus = as.numeric(grep('DSM 70459', killer.df.collapsed$strain_name))
  K21 = as.numeric(grep('T21.4 OS40', killer.df.collapsed$strain_name))
  K45 = as.numeric(grep('N-45 OS78', killer.df.collapsed$strain_name))
  K62 = as.numeric(grep('Q62.5 OS169', killer.df.collapsed$strain_name))
  K74a = as.numeric(grep('Q74.4 OS294', killer.df.collapsed$strain_name))
  K74b = as.numeric(grep('Y8.5 OS299', killer.df.collapsed$strain_name))
  
  killer.df.collapsed$group[K1]=paste(K1, "(K1)")
  killer.df.collapsed$group[K2]=paste(K2, "(K2)")
  killer.df.collapsed$group[K28]=paste(K28, "(K28)")
  killer.df.collapsed$group[Klus]=paste(Klus, "(Klus)")
  killer.df.collapsed$group[K21]=paste(K21, "(K21)")
  killer.df.collapsed$group[K45]=paste(K45, "(K45)")
  killer.df.collapsed$group[K62]=paste(K62, "(K62)")
  killer.df.collapsed$group[K74a]=paste(K74a, "(K74)")
  killer.df.collapsed$group[K74b]=paste(K74b, "(K74)")
  killer.df.collapsed$group = paste("Group ", killer.df.collapsed$group, 
                                    ", ", killer.df.collapsed$num_strains, " strains", sep="")
  
  nlawns = ncol(killer.df.collapsed)-3
  killer.df.collapsed[,1:nlawns] <- lapply(killer.df.collapsed[,1:nlawns], as.integer)
  killerMatrix <- as.matrix(killer.df.collapsed[,1:nlawns])
  rownames(killerMatrix) <- killer.df.collapsed$group
  myColors = c('oldlace','steelblue')
  dev.off()
  heatmap.2(killerMatrix, dendrogram = 'row', Rowv = TRUE, Colv = FALSE,
            trace = 'none', margins = c(8, 10), lwid = c(0.2,5),
            lhei = c(0.2,5), key = FALSE, col = myColors,
            xlab='Lawns', ylab="Killer Groups", colsep = 0:ncol(killerMatrix), 
            rowsep=0:nrow(killerMatrix),
            sepcolor='white', sepwidth=c(0.01, 0.01), 
            main="Collapsed Heatmap")
}