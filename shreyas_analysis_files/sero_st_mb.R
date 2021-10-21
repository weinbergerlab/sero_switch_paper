rm(list = ls())

#The following code is to conduct Market-Basket analysis on serotype*mlst matrix
#the first step is to get the data in the correct format of the transaction datatype
#In this format, MLSTs are like customers who are shopping for the optimum polysacchride composition
#in the serotypes
library(arules)
library(arulesViz)

#Load file containing serotype*mlst matrix
sero_st <- read.csv("./shreyas_analysis_files/serotype_mlst_mat.csv")
#Assign serotypes as row names
rownames(sero_st) <- sero_st$X
#Remove first column containing serotypes
sero_st <- sero_st[,-1]

#Remove X from column names
colnames(sero_st) <- substring(colnames(sero_st),2)
#Convert the matrix to binary
sero_st[sero_st != 0] <- 1 

#Transpose data
sero_st <- data.frame(t(sero_st))
#Remove X from column names
colnames(sero_st) <- substring(colnames(sero_st),2)

#Change the column name for 15B.C to 15B/C
colnames(sero_st)[which(colnames(sero_st) == "15B.C")] <- "15B/C"
colnames(sero_st)[which(colnames(sero_st) == "7B.C")] <- "7B/C"
colnames(sero_st)[which(colnames(sero_st) == "6C.D")] <- "6C/D"
sero_st <- sero_st[which(rowSums(sero_st) != 0),]
sero_st <- sero_st[-which(rowSums(sero_st) == 1),]

#Convert dataset to transcation data type for market basket analysis
ps_mat <- as.matrix(sero_st)
ps_tran <- as(ps_mat,"transactions")
summary(ps_tran)

#The plot shows the frequency of serotypes
itemFrequencyPlot(ps_tran, topN = 25, type='absolute')
itemFrequency(ps_tran, type='absolute')

#Since we are interested in the pairwise interactions between serotypes,
#rules with length of two will be generated to identify patterns of serotype co-occurence
rules_2 <- apriori(ps_tran,parameter = list(supp = 0.0001, conf = 0.0001,minlen=2, maxlen=2))
summary(rules_2)
#Data frame for rules with 2 serotypes sorted by support
rules_2sero <- inspect(sort(rules_2, by ="support"))

#Default plot of rules_2 from arules
plot(rules_2, method = "graph")

#The following code adds column for p-values from Fishers Exact test
#Make a dataframe with quality values for rules
r2_qual <- inspect(rules_2)
#Include interest measures
r2_qual$ST_chisq <- interestMeasure(rules_2, measure = c("chiSquared"), significance=TRUE, transactions = ps_tran)
r2_qual$ST_fisher <- interestMeasure(rules_2, measure = c("fishersExactTest"), transactions = ps_tran)

#r2_qual_sig <- r2_qual[which(r2_qual$count > 9),]

#Drop duplicates from final results in r2_qual
#Function to remove duplicated combinations occuring on even number of rows
#rem_dup <- function(comb_df){
#  i <- 1:nrow(comb_df)
#  j <- which(i %% 2 == 0)
#  comb_df <- comb_df[j,]
#}

#Remove duplicate rules
rules_all_2sero <- r2_qual
#Make separate data frame for serotype pairs associated with 3 or more MLST
rules_sig_2sero <- rules_all_2sero[rules_all_2sero$count > 2,]

#Save results in csv
write.csv(rules_all_2sero, "./shreyas_analysis_files/MB/rules_all_2serotypes_st.csv")
write.csv(rules_sig_2sero, "./shreyas_analysis_files/MB/rules_sig_2serotypes_st.csv")

#Save results in csv
#write.csv(rules_all_2sero, "index/rules_all_2serotypes_st.csv")
#write.csv(rules_sig_2sero, "index/rules_sig_2serotypes_st.csv")





