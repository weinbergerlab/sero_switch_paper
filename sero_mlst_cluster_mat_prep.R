rm(list = ls())
#This code is being written to identify the cause for differences in serotype associations from MLST and clusters
#in monte carlo and market basket analysis
#Load required libraries
library(dplyr)
library(reshape2)

#The following code is to fix the serotype names for the MLST definitions obtained from PubMLST on 03/22/2019
mlst_records <- read.csv("data_prep/sero_mlst_032219.csv")
nrow(mlst_records) #39725
#Make new dataframe with serotype and mlst info
sero_st <- mlst_records[,c(8,32:39)]

#FIX SEROTYPE NAMES
#Make serotype column to character
sero_st[,1] <- as.character(sero_st[,1])

#Copy sero_st into a new data frame
#Remove the serotype definitions that have more than 4 characters in their name
#since these are the ones most likely non typeable or unclear serotype names
validsero_st <- sero_st[which(nchar(sero_st[,1])<5),]
colnames(validsero_st)[1] <- "Var1"
validsero_st$Var1 <- as.character(validsero_st$Var1)

#Fix the serotype names
validsero_st[which(validsero_st$Var1 == "15 C"),1] <- "15B/C"
validsero_st[which(validsero_st$Var1 == "15B"),1] <- "15B/C"
validsero_st[which(validsero_st$Var1 == "15C"),1] <- "15B/C"
validsero_st[which(validsero_st$Var1 == "15BC"),1] <- "15B/C"
validsero_st[which(validsero_st$Var1 == "15c"),1] <- "15B/C"
validsero_st[which(validsero_st$Var1 == "17 F"),1] <- "17F"
validsero_st[which(validsero_st$Var1 == "18c"),1] <- "18C"
validsero_st[which(validsero_st$Var1 == "19 F"),1] <- "19F"
validsero_st[which(validsero_st$Var1 == "19a"),1] <- "19A"
validsero_st[which(validsero_st$Var1 == "23 F"),1] <- "23F"
validsero_st[which(validsero_st$Var1 == "23a"),1] <- "23A"
validsero_st[which(validsero_st$Var1 == "33 C"),1] <- "33C"
validsero_st[which(validsero_st$Var1 == "35f"),1] <- "35F"
validsero_st[which(validsero_st$Var1 == "6 B"),1] <- "6B"
validsero_st[which(validsero_st$Var1 == "7f"),1] <- "7F"
validsero_st[which(validsero_st$Var1 == "nt"),1] <- "NT"
validsero_st[which(validsero_st$Var1 == "06C"),1] <- "6C/D"
validsero_st[which(validsero_st$Var1 == "06A"),1] <- "6A"
validsero_st[which(validsero_st$Var1 == "06N"),1] <- "6N"
validsero_st[which(validsero_st$Var1 == "019F"),1] <- "19F"
validsero_st[which(validsero_st$Var1 == "015A"),1] <- "15A"
validsero_st[which(validsero_st$Var1 == "017F"),1] <- "17F"
validsero_st[which(validsero_st$Var1 == "011A"),1] <- "11A"
validsero_st[which(validsero_st$Var1 == "016F"),1] <- "16F"
validsero_st[which(validsero_st$Var1 == "006B"),1] <- "6B"
validsero_st[which(validsero_st$Var1 == "035B"),1] <- "35B"
validsero_st[which(validsero_st$Var1 == "022A"),1] <- "22A"
validsero_st[which(validsero_st$Var1 == "015B"),1] <- "15B/C"
validsero_st[which(validsero_st$Var1 == "023F"),1] <- "23F"
validsero_st[which(validsero_st$Var1 == "006E"),1] <- "6E"
validsero_st[which(validsero_st$Var1 == "015C"),1] <- "15B/C"
validsero_st[which(validsero_st$Var1 == "019A"),1] <- "19A"
validsero_st[which(validsero_st$Var1 == "019B"),1] <- "19B"
validsero_st[which(validsero_st$Var1 == "006A"),1] <- "6A"
validsero_st[which(validsero_st$Var1 == "007C"),1] <- "7C"
validsero_st[which(validsero_st$Var1 == "009N"),1] <- "9N"
validsero_st[which(validsero_st$Var1 == "09N"),1] <- "9N"
validsero_st[which(validsero_st$Var1 == "011B"),1] <- "11B"
validsero_st[which(validsero_st$Var1 == "012F"),1] <- "12F"
validsero_st[which(validsero_st$Var1 == "028A"),1] <- "28A"
validsero_st[which(validsero_st$Var1 == "018A"),1] <- "18A"
validsero_st[which(validsero_st$Var1 == "023B"),1] <- "23B"
validsero_st[which(validsero_st$Var1 == "007F"),1] <- "7F"
validsero_st[which(validsero_st$Var1 == "06B"),1] <- "6B"
validsero_st[which(validsero_st$Var1 == "06D"),1] <- "6C/D"
validsero_st[which(validsero_st$Var1 == "6D"),1] <- "6C/D"
validsero_st[which(validsero_st$Var1 == "6C"),1] <- "6C/D"
validsero_st[which(validsero_st$Var1 == "7B"),1] <- "7B/C"
validsero_st[which(validsero_st$Var1 == "7C"),1] <- "7B/C"
validsero_st[which(validsero_st$Var1 == "19A"),1] <- "19A"

#Make a vector with the serotypes that will not be included in the analysis because of ambiguous serotype calls
#The following serotypes don't exist so they will also be excluded
#10,11,12,15,16,17, 18,19,20,22,23,24,25,28,33,35,41,47,6,7,9
invalid_Sero <- c("15AF","19_","1A5","6A/B","6A/C","6AB","6ABC","6Bii","9L/N","9N/L","9NL","9V/A","N/A","N17F","ND","R","NT")
not_exist <- c("10","11","12","15","16","17","18","19","20","22","23","24","24","25","28","33","35","41","47","6","7","9","5F","53","62","")
validsero_st <- validsero_st[!(validsero_st$Var1 %in% invalid_Sero),]
validsero_st <- validsero_st[!(validsero_st$Var1 %in% not_exist),]

#There is a problem with serotype 19A since it occurs twice among the unique serotypes
# > sort(unique(validsero_st$Var1))
# [1] "1"     "10A"   "10B"   "10C"   "10F"   "11A"   "11B"   "11C"   "11D"   "11E"   "11F"   "12A"   "12B"   "12F"   "13"    "14"    "15A"  
# [18] "15B/C" "15F"   "16A"   "16F"   "17A"   "17F"   "18A"   "18B"   "18C"   "18F"   "19A"   "19B"   "19C"   "19F"   "19А"   "2"     "20A"  
# [35] "20B"   "21"    "22A"   "22F"   "23A"   "23B"   "23F"   "24A"   "24B"   "24F"   "25A"   "25F"   "27"    "28A"   "28F"   "29"    "3"    
# [52] "31"    "32A"   "32F"   "33A"   "33B"   "33C"   "33D"   "33F"   "34"    "35A"   "35B"   "35C"   "35F"   "36"    "37"    "38"    "39"   
# [69] "4"     "40"    "41A"   "41F"   "42"    "43"    "44"    "45"    "46"    "47A"   "47F"   "48"    "5"     "6A"    "6B"    "6C/D"  "6E"   
# [86] "6F"    "6G"    "7A"    "7B/C"  "7F"    "8"     "9A"    "9L"    "9N"    "9V"   

#Workaround
#To take care of this problem, make a new temporary data frame without the serotype 19A
tmp_without_19a <- validsero_st[-which(validsero_st$Var1 == "19A"),]
nrow(tmp_without_19a) #31346
#check for unique serotypes
# > sort(unique(tmp_without_19a$Var1))
# [1] "1"     "10A"   "10B"   "10C"   "10F"   "11A"   "11B"   "11C"   "11D"   "11E"   "11F"   "12A"   "12B"   "12F"   "13"    "14"    "15A"  
# [18] "15B/C" "15F"   "16A"   "16F"   "17A"   "17F"   "18A"   "18B"   "18C"   "18F"   "19B"   "19C"   "19F"   "19А"   "2"     "20A"   "20B"  
# [35] "21"    "22A"   "22F"   "23A"   "23B"   "23F"   "24A"   "24B"   "24F"   "25A"   "25F"   "27"    "28A"   "28F"   "29"    "3"     "31"   
# [52] "32A"   "32F"   "33A"   "33B"   "33C"   "33D"   "33F"   "34"    "35A"   "35B"   "35C"   "35F"   "36"    "37"    "38"    "39"    "4"    
# [69] "40"    "41A"   "41F"   "42"    "43"    "44"    "45"    "46"    "47A"   "47F"   "48"    "5"     "6A"    "6B"    "6C/D"  "6E"    "6F"   
# [86] "6G"    "7A"    "7B/C"  "7F"    "8"     "9A"    "9L"    "9N"    "9V"  

#There is still a 19A and this is causing the problem down the with 19A-19A pairings
#test1 <- tmp_without_19a[which(tmp_without_19a$Var1 == "19А"),] #serotype copied from above output
#test2 <- tmp_without_19a[which(tmp_without_19a$Var1 == "19A"),] #serotype manually typed
#rm(test1,test2)

#Hard to figure out the difference. Fix the serotype name in the data.
#Within the brackets inside the which function, make sure to copy the serotype name from above and not manually enter it.
tmp_without_19a$Var1[which(tmp_without_19a$Var1 == "19А")] <- "19A"

#Make another temporary dataframe with 19A from validsero_st
tmp_with_19a <- validsero_st[which(validsero_st$Var1 == "19A"),]

#Combine both dataframes
validsero_st <- rbind(tmp_with_19a, tmp_without_19a)

###################eBURST Input####################
#The following code prepares date that cen be used used as input in eBURST for assignment of clonal clusters
#Assign validsero_st to sero_st
sero_st <- validsero_st
#The alleles for ddl are factor. Convert them to int
sero_st$ddl <- as.integer(sero_st$ddl)
#Some of the MLST profiles have NA because of unidentified alleles
length(which(is.na(sero_st$ST..MLST.) == TRUE))
#Drop records where MLSTs are NA 
sero_st <- sero_st[-which(is.na(sero_st$ST..MLST.) == TRUE),]
#Make separate data frame for MLST profiles and save it in a csv file
mlst_profiles <- sero_st[,c(9,2:8)]
colnames(mlst_profiles)[1] <- "MLST"
#Save mlst_profiles in a csv file
write.csv(mlst_profiles,"data_prep/mlst_4eBURST.csv")

#eBURST accepts tab delimited text file without row and column names
write.table(mlst_profiles,"data_prep/eBURST_input.txt",sep = "\t",col.names = F,row.names = F)

#Inside Java jnlp program for eBURST, click "Open Ref" to load the eburst input file and click "Compute" 
#button in Analysis panel
#Save the output as "eBURST_output.txt"
#Run the R script "eBURST_output_data_scraping.R" with "eBURST_output.txt"
#The resulting "eBURST_MLST_clusters.csv" file will be used to create serotype*cluster matrix to be used in
#market basket analysis
##################################################

#Make a dataframe containing the details of serotypes and number of rows that are there for a single serotype
#There are several serotypes with multiple rows
sero_cols_mlst <- as.data.frame(table(sero_st$Var1))
write.csv(sero_cols_mlst,"data_prep/serotype_composition_mlst.csv")

######################SEROTYPE*MLST MATRIX PREP#########################
#Make serotype*MLST matrix from the sero_st data frame
sero_mlst <- as.data.frame(table(sero_st$Var1, sero_st$ST..MLST.))
#Convert validsero_st to wide format
sero_mlst <- dcast(sero_mlst, Var1~Var2, value.var = "Freq")
colnames(sero_mlst)[1] <- "Serotype"

#Assign serotype names as row names
rownames(sero_mlst) <- sero_mlst$Serotype
#Remove the first column
sero_mlst <- sero_mlst[,-1]

#Store the results in a new file that can be used for SDI and other association analysis
#write.csv(sero_mlst, "data_prep/serotype_mlst_mat.csv")

######################SEROTYPE*CLUSTER MATRIX PREP######################
#Load eBURST cluster information for mlst from "eBURST_MLST_clusters.csv"
eburst_output <- read.csv("data_prep/eBURST_MLST_clusters.csv")
#remove first column with row numbers
eburst_output <- eburst_output[,-1]
#Make new data frmae cluster_profiles with data from mlst_profiles data frame
cluster_profiles <- sero_st
#Add column for cluster to sero_clus
cluster_profiles$Cluster <- eburst_output$Cluster[match(cluster_profiles$ST..MLST., eburst_output$ST)]
#Identify MLST of isolates for which the clsuters have not been assigned
length(which(is.na(cluster_profiles$Cluster)==TRUE)) #3478
#Drop the isolates/records that do not have cluster assignment
cluster_profiles <- cluster_profiles[-which(is.na(cluster_profiles$Cluster)==TRUE),]
#Make a dataframe containing the details of serotypes and number of rows that are there for a single serotype
#There are several serotypes with multiple rows
sero_cols_clus <- as.data.frame(table(cluster_profiles$Var1))

#Make serotype*cluster matrix from the sero_st data frame
sero_clus <- as.data.frame(table(cluster_profiles$Var1, cluster_profiles$Cluster))
#Convert validsero_st to wide format
sero_clus <- dcast(sero_clus, Var1~Var2, value.var = "Freq")
colnames(sero_clus)[1] <- "Serotype"

#Assign serotype names as row names
rownames(sero_clus) <- sero_clus$Serotype
#Remove the first column
sero_clus <- sero_clus[,-1]

#Store the results in a new file that can be used for SDI and other association analysis
write.csv(sero_clus, "data_prep/serotype_cluster_mat.csv")
