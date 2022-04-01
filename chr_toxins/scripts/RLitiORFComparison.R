
##Import the original Liti R script file from khr/khs R script as well as the generated ORF file
table1 <- read.csv('/Users/lab/Downloads/LitiKillList.csv', header = TRUE)
table2 <- read.csv('/Users/lab/Downloads/AlignedORFsForFischer.csv', header = TRUE)
mergedTbl = merge(table1,table2,by='Strain',all=TRUE)
for(i in 1:nrow(mergedTbl)){
  if ((mergedTbl$aaseq[i]) != (mergedTbl$Prot[i])){ #This matches the new ORF aa with the original aa
    mergedTbl$ORFMatch[i] <- "Mismatch"
    }else {
      mergedTbl$ORFMatch[i] <- "Match"
}}
write.csv(mergedTbl, file="/Users/lab/Downloads/AlignedORFsNoMin.csv", row.names = FALSE)
###print out a new merged file with A Match/Mismatch file. Most will be matches. Filter out in excel 
### and can move to a new csv if desired.####