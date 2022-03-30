#!/usr/bin/env Rscript

main <- function() {
  dir.create(Sys.getenv("R_LIBS_USER"), recursive = TRUE)  # create personal library
  .libPaths(Sys.getenv("R_LIBS_USER"))  # add to the path
  ## open libraries
  packages = c('ggplot2','stringr')
  # check.packages function: install and load multiple R packages.
  # Check to see if packages are installed. Install them if they are not, then load them into the R session.
  check.packages <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
      install.packages(new.pkg, dependencies = TRUE, repo="http://cran.us.r-project.org")
    sapply(pkg, require, character.only = TRUE)
  }
  msg.trap <- capture.output( suppressMessages( check.packages(packages)) ) # silent mode
  
  ## parse arguments
  args <- commandArgs(trailingOnly = TRUE)
  ## parse arguments
  file <- args[1]  ## needed for importing data
  # output option (argument #2)
  if (is.na(args[2])==FALSE){
    gextension = args[2]
  } else {
    output_file=paste(str_split(file, "[.]")[[1]][1],'.jpeg',sep="")
  }
  # graph height option (argument #3)
  if (is.na(args[3])==FALSE){
    grheight = args[3]
  } else {
    grheight = 3
  }
  # graph width option (argument #4)
  if (is.na(args[4])==FALSE){
    grwidth = args[4]
  } else {
    grwidth = 5
  }
  
  
  #****************** MAIN SCRIPT *******************

  ## Read data from pre-formatted dataset.
  depth.df <- read.table(file, sep="\t", header = TRUE)
  names(depth.df) = c('sample','pos','depth')
  depth.df$depth[depth.df$depth == 0] <- NA
  depth.df = na.omit(depth.df)
  filename = sapply(str_split(file, "[/]"), tail, 1)
  filename_v = as.vector(unlist(str_split(filename, "_")))
  referencename = filename_v[length(filename_v)-2]
  samplename = paste(filename_v[1:(length(filename_v)-3)], collapse=" ")
  
  ## make sure graph will show the whole toxin genome
  ref = as.character(depth.df[1,1])
  if (ref == 'NC_001782.1') { # M1
    tox_len = 1801 
  } else if (ref == 'MF957266.1') { # M2
    tox_len = 1723
  } else if (ref == 'KJ796682.1') { # M28
    tox_len = 1038
  } else if (ref == 'GU723494.2') { # Mlus
    tox_len = 2033
  } else if (ref == 'AF039063.1') { # 20s
    tox_len = 2514
  } else if (ref == 'U90136.1') { # 23s
    tox_len = 2891
  } else if (ref == 'MN267768.1') { # I329
    tox_len = 2509
  }
  
  ############## AREA PLOT #################
  if (exists("tox_len")){
    ggplot(depth.df, aes(x=pos, y=depth)) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      geom_area(fill = "lightblue") +
      theme_classic() +
      labs(x="Position (nt)", y="Coverage", 
           title=paste(samplename,"Reads Mapped to",referencename)) +
      scale_x_continuous(breaks = seq(0, tox_len, by=200), limits=c(0,tox_len))
    ggsave(filename=output_file, height=as.numeric(grheight),
           width=as.numeric(grwidth))
  } else {
    ggplot(depth.df, aes(x=pos, y=depth)) +
      scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      geom_area(fill = "lightblue") +
      theme_classic() +
      labs(x="Position (nt)", y="Coverage", 
           title=paste(samplename,"Reads Mapped to",referencename))
    ggsave(filename=output_file, height=as.numeric(grheight),
           width=as.numeric(grwidth))
  }

}
main()
