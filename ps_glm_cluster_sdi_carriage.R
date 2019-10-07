rm(list = ls())

library("ggplot2")
library("dplyr")
library("stargazer")

#set working directory
setwd("D:\\ss_paper")

#load file with carriage data
carriage <- read.csv("sdi/pneumo_master_carriage.csv",stringsAsFactors = FALSE)
#The serotype 15 in carriage dataset has 2 rwos for 15B and 15C. Combine those rows.
carriage_15BC <- carriage[c(which(carriage$ST == "15B"),which(carriage$ST == "15C")),]
carriage_15BC <- data.frame(t(c("15B/C",colSums(carriage_15BC[2:ncol(carriage_15BC)]))))
colnames(carriage_15BC) <- colnames(carriage)
carriage <- carriage[-c(which(carriage$ST == "15B"),which(carriage$ST == "15C")),]
carriage <- rbind(carriage,carriage_15BC)
#Change serotype 6A/C to 6A
carriage$ST[which(carriage$ST == "6A/C")] <- "6A" 
#Unlist the columns with carriage values and convert the data type to integer
carriage[,2:ncol(carriage)] <- as.integer(unlist(carriage[,2:ncol(carriage)]))
#Convert Na to 0
carriage[is.na(carriage)] <- 0

#Load the data containing sdi values for clusters and obtained from eBURST into the environment
sdi_st <- read.csv("sdi/sdi_cc.csv")
#Select columns with serotypes, SDI, and big_N
sdi_st <- sdi_st[,c(2,3,5)]
colnames(sdi_st) <- c("Serotype","SDI","Frequency")
#Sort SDI values
sdi_st <- sdi_st[order(-sdi_st$SDI),]
#Make barplot
sdi_barplot <- sdi_st
#colnames(sdi_barplot) <- c("Serotype","SDI","MLST Frequency")
p<-ggplot(data=sdi_barplot, aes(x=reorder(Serotype,-SDI), y=SDI)) +
  geom_point(alpha = 0.5, aes(size = Frequency)) + labs(x = "Serotype")
p+coord_flip() + theme(panel.background = element_blank(), axis.line = element_line(size = 0.5), axis.text = element_text(face = "bold")) 
#rownames(sdi_barplot) <- sdi_barplot$Serotype


#Load polysaccharide composition data
ps_composition <- read.csv("sdi/ps_composition.csv",stringsAsFactors = FALSE)
#Fix the serotype names
ps_composition$Serotype <- as.character(ps_composition$Serotype)
ps_composition[which(ps_composition$Serotype == "15B"),1] <- "15B/C"
ps_composition[which(ps_composition$Serotype == "15C"),1] <- "15B/C"
ps_composition$Serotype <- as.factor(ps_composition$Serotype)
#Combine multiple rows of single serotypes by adding the values
ps_composition <- aggregate(ps_composition[,2:ncol(ps_composition)], by=list(Serotype=ps_composition$Serotype), FUN=sum)
ps_tmp <- ps_composition[which(ps_composition$Serotype == "15B/C"),2:ncol(ps_composition)]
ps_tmp[ps_tmp == 2] <- 1
ps_composition[which(ps_composition$Serotype == "15B/C"),2:ncol(ps_composition)] <- ps_tmp
#Add a column for presence or absence of N-Acetylated sugars
nacs <- ps_composition %>% select(Serotype ,FucNAc, GalNAc, GlcNAc, ManNAc, PneNAc)
nacs$NAc <- rowSums(nacs[,2:6])
nacs$NAc[nacs$NAc > 0] <- 1 
ps_composition$NAc <- nacs$NAc

#Add a column with SDI values for clusters in the polysaccharide composition data frame by matching for serotypes
ps_composition$sdi <- sdi_st$SDI[match(ps_composition$Serotype, sdi_st$Serotype)]
ps_composition$freq <- sdi_st$freq[match(ps_composition$Serotype, sdi_st$Serotype)]
#Which rows have missing values for SDI?
sdi_missing <- ps_composition$Serotype[which(is.na(ps_composition$sdi)==TRUE)]
#Remove the serotypes (from above) for which sdi values re not available
ps_composition <- ps_composition[!(ps_composition$Serotype %in% sdi_missing),]
#Remove serotypes with SDI value 0
#ps_composition <- ps_composition[!(ps_composition$sdi == 0),] #dont need this because the isolate number threshold was 10 in mat_prep
#Perform logit transformation of sdi values
ps_composition$tr_sdi <- log(ps_composition$sdi/(1 - ps_composition$sdi))
#ps_composition_sdi <- ps_composition_sdi[ps_composition_sdi$freq > 25,]

#Process the carriage data by taking percentage for each value in a column, followed by 
#averaging across rows and final square root transformation
#for(i in 2:ncol(carriage)){
#  carriage[,i] <- sapply(carriage[,i], function(x) (x/sum(carriage[,i])*100))
#}
#carriage$Ave <- rowSums(carriage[,2:ncol(carriage)])/(ncol(carriage)-1)
#carriage$Ave <- sqrt(carriage$Ave)
               
#Add carriage$Ave to the ps_composition dataset
#ps_composition$Carriage <- carriage$Ave[match(ps_composition$Serotype,carriage$ST)]

#Add Seeman carriage values from England to ps_composition dataset
ps_composition$Carriage <- carriage$SLEEMANCARRN[match(ps_composition$Serotype,carriage$ST)]

#Remove rows for which carriage values are NA
ps_composition <- na.omit(ps_composition)
#ps_composition[is.na(ps_composition)] <- 0
#ps_composition <- ps_composition[-which(ps_composition$tr_sdi == Inf),]
#plot(ps_composition$sdi, ps_composition$Carriage)

#include only those sugars for which have more than 5 serotypes
ps_composition1 <- ps_composition[,2:26]
ps_composition1[(nrow(ps_composition1)+1),] <- colSums(ps_composition1) 
less_than_3 <- colnames(ps_composition1[,which((ps_composition1[nrow(ps_composition1),]) < 3)])
ps_composition <- dplyr::select(ps_composition, -one_of(less_than_3))

#Remove serotypes with carriage values 0
#ps_composition <- ps_composition[!ps_composition$Carriage == 0,]

#Regression analysis
fdat2 <- ps_composition[,2:(ncol(ps_composition)-3)]
reg.ace.nb1 <- vector("list", ncol(fdat2)) 
aic.ace.nb1 <- 1: ncol(fdat2)
pred.ace.nb1<-matrix(nrow = nrow(fdat2), ncol = ncol(fdat2))
reg.coef.nb1 <- 1:ncol(fdat2)
for (i in 1:ncol(fdat2)){
  if(nlevels(as.factor(fdat2[,i])) > 1){
    reg.ace.nb1[[i]]<-glm(ps_composition$tr_sdi ~ ps_composition$Carriage + as.factor(fdat2[,i]),family = gaussian)
    aic.ace.nb1[[i]]<-reg.ace.nb1[[i]]$aic
    pred.ace.nb1[,i]<-predict(reg.ace.nb1[[i]])
    reg.coef.nb1[[i]] <- reg.ace.nb1[i][[1]]$coefficients[3][[1]]
  }
  else{
    #AIC values cannot be calculated since there is only one level
    reg.ace.nb1[[i]]<-"1 level"
    #Randomly assign value of 1000 since minimum vaue of aic will be used for computing weights (refer wgt)
    aic.ace.nb1[[i]]<-1000
    #Predicted value will be same as existing value
    pred.ace.nb1[,i]<-fdat2[1,i]
    #Regression coefficient will be the same
    reg.coef.nb1[[i]] <- reg.ace.nb1[i][[1]]$coefficients[3][[1]]
  }
}
#Compute weights for each position
#wgt_sugar <- exp(-0.5*(aic.ace.nb1-min(aic.ace.nb1)))/sum(exp(-0.5*(aic.ace.nb1-min(aic.ace.nb1)))) #Akaike Weights

#Compute weight using model likelihood
wgt_sugar <- exp(-0.5*(aic.ace.nb1-min(aic.ace.nb1))) #Model Likelihood

#The columns which had one level had constant values for weight.
#Replace the constant values with 0 (for better visualization of higher weights in a plot)
wgt1_sugar <- wgt_sugar
#wgt1_sugar[wgt1_sugar == wgt1_sugar[1]] <- 0
pos1 <- colnames(fdat2)
#Make a dataframe containing positions and weights
ps_table <- data.frame(pos1,aic.ace.nb1,reg.coef.nb1) 
names(ps_table) <- c("Component","AIC","Regression.Coefficient")
#Sort ps_table on the basis of aic values
ps_table <- ps_table[order(ps_table$AIC),]
tabular((Component + 1)~(n=1), AIC, Regression.Coefficient, ps_table)
write.csv(ps_table,"glm_AIC_reg_coef_cc.csv")

null_model <- glm(ps_composition$tr_sdi ~ ps_composition$Carriage,family = gaussian)
null_model$aic #166.8775 #169.571
#wgt_null_model <- exp(-0.5*(null_model$aic-min(aic.ace.nb1)))/sum(exp(-0.5*(aic.ace.nb1-min(aic.ace.nb1))))
#wgt_null_model #0.1238 #0.3467
#MAKE SCATTER PLOT OF SDI AND CARRIAGE
#plot(ps_composition$tr_sdi, ps_composition$Carriage)

# In the model, instead of using the null model threshold, use the value of exp(-0.5*(2)) as threshold.
#The reasoning behind it is that we are using the delta of 2 for AIC directly in computing the threshold.
z <- ggplot(data=ps_table, aes(x=ps_table$Component, y=ps_table$Weight, group=1)) +
  geom_line() + theme(axis.text.x = element_text(angle=75, hjust=1, size = 10),panel.background = element_blank()) +
  geom_hline(yintercept = exp(-0.5*(2)), linetype = 2, color = "blue")+
  theme(axis.line.x = element_line(color="black", size = 0.5), axis.line.y = element_line(color="black", size = 0.5))+
  labs(x = "Polysaccharide", y = "Akaike Weights")
z

#Make a plot to understand as to how the serotype diversity relates to the polysaccharide composition
ps_sdi <- ps_composition[,1:11]
ps_freq_sdi <- data.frame(colnames(ps_composition)[2:15], colSums(ps_composition[,2:15]))
rownames(ps_freq_sdi) <- 1:nrow(ps_freq_sdi)
colnames(ps_freq_sdi) <- c("CPS","Frequency")
ps_freq_sdi$sdi_mean <- NA
ps_freq_sdi$carr_mean <- NA

for (i in 2:15){
  tmp_df <- data.frame(ps_composition[,i], ps_composition$sdi, ps_composition$Carriage)
  tmp_df <- tmp_df[which(tmp_df[,1]==1),]
  ps_freq_sdi[(i-1),] <- c(colnames(ps_composition)[i],sum(tmp_df[,1]),round(mean(tmp_df[,2]),3), round(mean(tmp_df[,3]),3))
}
plot(ps_freq_sdi$sdi_mean, ps_freq_sdi$carr_mean)
