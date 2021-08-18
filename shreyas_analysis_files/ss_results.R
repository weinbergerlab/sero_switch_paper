
mc_st <- read.csv("./shreyas_analysis_files/monte_carlo/MLST/975tile_st.csv")
mc_cc <- read.csv("./shreyas_analysis_files/monte_carlo/CC/975tile_cc.csv")
mb_st <- read.csv("./shreyas_analysis_files/MB/rules_all_2serotypes_st.csv")
mb_cc <- read.csv("./shreyas_analysis_files/MB/rules_all_2serotypes_cc.csv")

#Fix Serotype names to make them readable
mb_st$lhs <- gsub('\\{','',mb_st$lhs)
mb_st$lhs <- gsub('\\}','',mb_st$lhs)
mb_st$rhs <- gsub('\\{','',mb_st$rhs)
mb_st$rhs <- gsub('\\}','',mb_st$rhs)

mb_cc$lhs <- gsub('\\{','',mb_cc$lhs)
mb_cc$lhs <- gsub('\\}','',mb_cc$lhs)
mb_cc$rhs <- gsub('\\{','',mb_cc$rhs)
mb_cc$rhs <- gsub('\\}','',mb_cc$rhs)

#Make serotype pair names same format in both MC and MC
mc_st$pair <- paste(mc_st$Serotype1, "-",mc_st$Serotype2)
mc_cc$pair <- paste(mc_cc$Serotype1, "-",mc_cc$Serotype2)

mb_st$pair <- paste(mb_st$lhs,"-",mb_st$rhs)
mb_cc$pair <- paste(mb_cc$lhs,"-",mb_cc$rhs)

#Identify unique pairs in all output files
pair <- list(mc_st$pair, mc_cc$pair, mb_st$pair, mb_cc$pair)
pair <- Reduce(intersect,pair)

#Combine results from MC output for both ST and CC
mc_merge <- merge(mc_st,mc_cc,by = "pair")
mc_merge <- mc_merge[,c(1,3:6,10,11)]
colnames(mc_merge) <- c("sero_pair","Serotype1","Serotype2","ST_Obs","ST_Exp","CC_Obs","CC_Exp")
#Prepare index
#Weight = 1 if the number of isolates for a pair are 3 or more
#Weight = 1 if the frequency of observed isolates is higher than the values obtained from monte carlo (MC)
mc_merge$mc_st_index_freq <- 0
mc_merge$mc_st_index_freq[which(mc_merge$ST_Obs > 2)] <- 1
mc_merge$mc_st_index_sig <- 0
mc_merge$mc_st_index_sig[which(mc_merge$ST_Obs>mc_merge$ST_Exp)] <- 1
mc_merge$mc_cc_index_freq <- 0
mc_merge$mc_cc_index_freq[which(mc_merge$CC_Obs > 2)] <- 1
mc_merge$mc_cc_index_sig <- 0
mc_merge$mc_cc_index_sig[which(mc_merge$CC_Obs>mc_merge$CC_Exp)] <- 1

#Combine results from MB output for both ST and CC
mb_merge <- merge(mb_st,mb_cc,by = "pair")
#mb_merge <- mb_merge[,c(1,6:11,16:21)]
mb_merge <- mb_merge[,c('pair','support.x','confidence.x','lift.x','count.x',"ST_chisq",   "ST_fisher",
                        'support.y','confidence.y','lift.y','count.y',"CC_chisq",   "CC_fisher"  )]

colnames(mb_merge) <- c("sero_pair","ST_support","ST_confidence","ST_lift","ST_count","ST_chisq","ST_fisher","CC_support","CC_confidence","CC_lift","CC_count","CC_chisq","CC_fisher")

#Prepare index for ST
#Weight = 0.25 if confidence measure values are in the top quartile. These can be adjusted by adding "probs=" to the quantile function
mb_merge$mb_st_index_con <- 0
confidence_threshold <- quantile(mb_merge$ST_confidence)
mb_merge$mb_st_index_con[which(mb_merge$ST_confidence > confidence_threshold[4])] <- 0.25
#Weight = 0.25 if chi square values are significant
mb_merge$mb_st_index_chi <- 0
mb_merge$mb_st_index_chi[which(mb_merge$ST_chisq < 0.05)] <- 0.25
#Weight = 0.25 if lift > 1
mb_merge$mb_st_index_lift <- 0
mb_merge$mb_st_index_lift[which(mb_merge$ST_lift > 1)] <- 0.25
#Weight = 0.25 if fisher's exact test is significant
mb_merge$mb_st_index_fisher <- 0
mb_merge$mb_st_index_fisher[which(mb_merge$ST_fisher < 0.05)] <- 0.25

#Prepare index for CC
#Weight = 0.25 if confidence measure values are in the top quartile. These can be adjusted by adding "probs=" to the quantile function
mb_merge$mb_cc_index_con <- 0
confidence_threshold <- quantile(mb_merge$CC_confidence)
mb_merge$mb_cc_index_con[which(mb_merge$CC_confidence > confidence_threshold[4])] <- 0.25
#Weight = 0.25 if chi square values are significant


mb_merge$mb_cc_index_chi <- 0
mb_merge$mb_cc_index_chi[which(mb_merge$CC_chisq < 0.05)] <- 0.25

#Weight = 0.25 if lift > 1
mb_merge$mb_cc_index_lift <- 0
mb_merge$mb_cc_index_lift[which(mb_merge$CC_lift > 1)] <- 0.25
#Weight = 0.25 if fisher's exact test is significant
mb_merge$mb_cc_index_fisher <- 0
mb_merge$mb_cc_index_fisher[which(mb_merge$CC_fisher < 0.05)] <- 0.25

#SUMMARY OF INDEX
#Essentially we have assigned a weight of 1 to the count/frequency of isolates for a serotype pair
#Weight of 1 to both ST and CC for significant outputs from Monte Carlo
#In Market Basket weights assigned from outputs for ST and CC are still 1 
#Thus the score in the index ranges from 0 to 5

#Merge results for MC and MB
ss_results <- merge(mc_merge,mb_merge,by = "sero_pair")

#Make a separate table for index
ss_index <- ss_results[,c(1:3,8:11,24:31)]
ss_index$total <- rowSums(ss_index[,4:15])

write.csv(ss_results,"./shreyas_analysis_files/ss_results.csv")
write.csv(ss_index, "./shreyas_analysis_files/ss_index.csv")

# length(ss_index$total)
# length(which(ss_index$total >= 1))
# length(which(ss_index$total >= 2))
# length(which(ss_index$total >= 3))
# length(which(ss_index$total >= 4))
# length(which(ss_index$total >= 5))
# 
# length(which(ss_index$total > 5.75))
# test <- ss_index[which(ss_index$mc_st_index_freq == 1),]
# test <- ss_index[order(-ss_index$total),]
# test <- test[,c(1,14)]

