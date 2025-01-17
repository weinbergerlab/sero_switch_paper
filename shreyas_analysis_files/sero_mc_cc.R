rm(list=ls())

library(tidyr)
library(dplyr)
library(reshape2)
#library(svMisc) #for progress

#GENERATE sero_cc_unique.csv from serotype_mlst_mat.csv
#This dataset will contain unique serotype-MLST combinations that were observed in data from PubMLST
#Load serotype*mlst matrix
sero_cc <- read.csv("./shreyas_analysis_files/serotype_cluster_mat.csv")
#Remove X from column names
colnames(sero_cc) <- substring(colnames(sero_cc),2)
#Add column name to first column which will be assigned as id for the melt command
colnames(sero_cc)[1] <- "Serotype"
#Remove NT
#sero_cc <- sero_cc[-99,]

#Convert the data frame to binary
#Convert the matrix to binary
#The first column which contains the names of the serotypes should not be changed while converting the data frame to binary
#Copy the sero_cc data frame from second column onwards into a new one
sero_tmp <- sero_cc[,2:ncol(sero_cc)]
sero_tmp[sero_tmp != 0] <- 1 
sero_cc[,2:ncol(sero_cc)] <- sero_tmp

#Make a list of shared mlst between a pair of serotypes by converting mlst_NbyN to wide format
sero_cc_long <- melt(sero_cc, id="Serotype")
colnames(sero_cc_long) <- c("stA","CC","frequency")

#Remove rows where frequency is 0
sero_cc_unique <- sero_cc_long[-(which(sero_cc_long$frequency == 0)),]
#Convert serotype values to character datatype
sero_cc_unique$stA <- as.character(sero_cc_unique$stA)
#Save the data in new csv file
#write.csv(sero_cc_unique,"sero_cc_unique.csv")

#Make a dataset with serotype pairs and shared MLST between them
#Assign serotype names as row names
rownames(sero_cc) <- sero_cc$Serotype
#Remove column with rownames
sero_cc <- sero_cc[,-1]
#Make a table with all possible combinations of serotypes for pairwise comparisons
sero_pairs_cc <- combn(rownames(sero_cc),2)
#Transform the table
sero_pairs_cc <- data.frame(t(sero_pairs_cc))
colnames(sero_pairs_cc) <- c("Sero1","Sero2")
sero_pairs_cc$Sero1 <- as.character(sero_pairs_cc$Sero1)
sero_pairs_cc$Sero2 <- as.character(sero_pairs_cc$Sero2)
sero_pairs_cc$shared_cc <- 0

#Calculate the number of shared cc
#To calculate the number of shared cc, number of cc that are associated with both serotypes in a pair will be calculated
#There are several MLST (almost 372 out of 641)in the data that are associated with only 1 serotype
length(which(colSums(sero_cc) == 1)) #372

#To make computation faster, include only those CC that are associated with more than one serotype
sero_cc <- sero_cc[,which(colSums(sero_cc) > 1)]

for(i in 1:nrow(sero_pairs_cc)){
  pb   <- txtProgressBar(1, nrow(sero_pairs_cc), style=3)
  Sys.sleep(0.01)
  sero_pairs_cc$shared_cc[i] <- length(which(as.numeric(sero_cc[sero_pairs_cc[i,1],])+as.numeric(sero_cc[sero_pairs_cc[i,2],]) == 2))
  setTxtProgressBar(pb, i)
}

#Save results in csv
#write.csv(sero_pairs_cc,"monte_carlo/CC/sero_pairs_cc.csv")
#write.csv(sero_pairs_cc,"index/sero_pairs_cc.csv")
#Sampling without replacement for 1000 iterations
samples <- list()
#There are a total of 12812 unique serotype-mlst combinations.
#The following code randomly samples one of the numbers between 
#1 and 12812 with each value occuring only once (that's why called sampling without replacement)
#This is done for 1000 iterations and the output is a list of results from those iterations
for(i in 1:1000){
  samples[[i]] <- sample(1:nrow(sero_cc_unique), size = nrow(sero_cc_unique), replace = FALSE)
}

# Monte Carlo (MC)
#Through the values obtained from resampling from the steps above, perform pairwise operation similar to the previous step,
#and store the values in 1000 columns of a data frame
mc_output <- sero_pairs_cc
colnames(mc_output) <- c("stA","stB","shared_cc")
for (j in 1:1000){
  pb_mc   <- txtProgressBar(1, 1000, style=3)
  Sys.sleep(0.01)
  tmp_df1 <- sero_cc_unique[,1:2]
  for(k in 1:nrow(sero_cc_unique)){
    tmp_df1$stA[k] <- as.character(sero_cc_unique$stA[samples[[j]][k]])
  }
  tmp_df2 <- tmp_df1
  colnames(tmp_df2) <- c("stB","CC")
  df_merged <- merge(x = tmp_df1, y = tmp_df2, by = "CC")
  #Delete columns where bpth serotypes in a pair are same
  df_merged <- df_merged[which(df_merged$stA != df_merged$stB),]
  merged_count <- df_merged %>% group_by(stA,stB) %>% summarize(n())
  mc_output <- left_join(mc_output, merged_count, by = c("stA" = "stA", "stB" = "stB"))
  colnames(mc_output)[j+3] <- paste("Run",j,sep = "_")
  setTxtProgressBar(pb_mc, j)
}

#Change NA values to 0
mc_output[is.na(mc_output) == T] <- 0
#Save results
write.csv(mc_output, "./shreyas_analysis_files/monte_carlo/CC/1000_run.csv")
#Add a column with values for upper 97.5 percentile (percent quantile) for each row
#mc_output[,1004] <- apply(mc_output[,4:ncol(mc_output)], 1, quantile, probs = 0.025)

quants <- t(apply(mc_output[,4:ncol(mc_output)], 1, quantile, probs = c(0.025, 0.975)))



#Make a data frame containing serotype pairs, frequency of that pairing, and 2.5%ile value
sero_quant <- cbind(mc_output[,c(1:3)],quants)
#Remove serotype pairs with no shared MLST
sero_quant <- sero_quant[-which(sero_quant$shared_cc == 0),]

#Make a final data frame for serotype pairs where the observed frequency is greater than 97.5 %ile value
#sero_fin <- sero_quant[(sero_quant$V1004 < sero_quant$shared_mlst),]
sero_fin <- sero_quant
#Assign column names and store the results in a separate file
colnames(sero_fin) <- c("Serotype1","Serotype2","Freq","025tile", "975tile")
#write.csv(sero_fin, "monte_carlo/CC/975tile_cc.csv")
write.csv(sero_fin, "./shreyas_analysis_files/monte_carlo/CC/975tile_cc.csv")

#Filter sero_fin with threshold of 3 on the basis of number of shared cc
sero_fin_wt <- sero_fin[sero_fin$Freq > 2,]
#Save results
#write.csv(sero_fin_wt,"monte_carlo/CC/975tile_cc_wt.csv")
write.csv(sero_fin_wt,"./shreyas_analysis_files/monte_carlo/CC/975tile_cc_wt.csv")
