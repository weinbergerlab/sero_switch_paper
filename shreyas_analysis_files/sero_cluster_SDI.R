rm(list=ls())
#Load the prepped data file into the environment
sero_st <- read.csv("data_prep/serotype_cluster_mat.csv")
#Assign serotype names as row names
rownames(sero_st) <- sero_st$X
#Remove column containing the serotypes
sero_st <- sero_st[,-1]
#Remove X from column names
colnames(sero_st) <- substring(colnames(sero_st),2)
#SDI values for each serotype (each row) will be calculated
#Calculate big_N by adding values for the whole row
sero_st$big_N <- rowSums(sero_st)
#Remove serotypes with number of isolates less than 10
sero_st <- sero_st[sero_st$big_N > 10,]
#SDI values will be calculated for two different formulae
#Make separate variables for small_n for both approaches
sero_st$small_n1 <- 0
sero_st$small_n2 <- 0
sero_st$sdi1 <- 0
sero_st$sdi2 <- 0

#Compute SDI values for each serotype with respect to MLST
for (i in 1:nrow(sero_st)){
  for (j in 1:581){
    if(sero_st[i,j] > 0){
      sero_st$small_n1[i] <- sero_st$small_n1[i] + (sero_st[i,j]/sero_st$big_N[i])^2
      sero_st$small_n2[i] <- sero_st$small_n2[i] + (sero_st[i,j]*(sero_st[i,j]-1))
    }
    else{
      sero_st$small_n1[i] <- sero_st$small_n1[i]
      sero_st$small_n2[i] <- sero_st$small_n2[i]
    }
  }
  sero_st$sdi1[i] <- 1-sero_st$small_n1[i]
  sero_st$sdi2[i] <- 1-(sero_st$small_n2[i]/(sero_st$big_N[i]*(sero_st$big_N[i] - 1)))
}
sero_sdi <- data.frame(rownames(sero_st),sero_st$sdi1,sero_st$sdi2,sero_st$big_N)
sero_sdi$sero_st.sdi1 <- sero_st$sdi1*(sero_st$big_N/(sero_st$big_N - 1))
write.csv(sero_sdi, "sdi_cc.csv")


