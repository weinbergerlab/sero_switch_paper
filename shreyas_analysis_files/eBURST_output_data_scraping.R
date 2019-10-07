#The goal of the code is to scrape the results from eBURST clustering/grouping tool for MLSTs
rm(list=ls())

#The "eBURST_output.txt" contains mlst clusters made by the eburst program. The goal of this code is to extract the information about the clusters
#from the text file. Reading the whole text file into the memory can slow down the program. To read the text file line by line, below is the code
#The current file is not so big to worry about the memory, but the following code block with "readLines" function should be used for similar tasks in the future.
#Create and open connection to the file
filecon <- file("data_prep/eBURST_output.txt", "r")
#Read all the lines from file and store all lines into a vector
grp_lines <- readLines(filecon)
#Close connection
close(filecon)
#length(grp_lines) 

#Extract information about all the groups/clusters from the file
#Setting value to true provides verbose output
grp_info <- grep("Group ",grp_lines, value = TRUE)
#Same command with value=FALSE (DEFAULT) provides ROW INDEXES 
grprow_info <- grep("Group ",grp_lines)

#Create vectors to store cluster information from grp_info
grp_no <- vector()
isolates <- vector()
seq_types <- vector()
founder <- vector()

#Loop through the grp_info vector to extract values
for (i in 1:length(grp_info)){
  #Store values in respective vectors after splitting the line (element within vector)
  temp_vec <- strsplit(grp_info[i]," | ")[[1]]
  grp_no[i] <- temp_vec[2]
  isolates[i] <- as.numeric(temp_vec[7])
  seq_types[i] <- as.numeric(temp_vec[14])
  founder[i] <- as.character(temp_vec[21])
}
#Remove ":" from group names
grp_no <- gsub(":","",grp_no)
#Add all extracted values for each cluster to a data frame
grp_df <- cbind.data.frame(grp_no,isolates,seq_types,founder)
#From the file text, extract row indices for the lines containing cluster information
strow_info <- as.numeric(grep("ST   ",grp_lines))

#len_vec <- vector() #Stores size of each cluster (refer loop below)
#Separate data frames will be made for each cluster, and all clusters will be saved in a list
clus_df <- list() 
#loop through the cluster data to make a final data frame containing values with respect to sequence types
for (i in 1:length(strow_info)){
    #Make a temporary vector to store all the lines for a single cluster
    temp_vec <- grp_lines[strow_info[i]:(strow_info[i]+grp_df$seq_types[i])]
    #Store the size of the cluster (vector created above)
    #len_vec[i] <- length(temp_vec)
    #Declare individual vectors that will store individual information from each line
    ST <- vector()
    Freq <- vector()
    SLV <- vector()
    DLV <- vector()
    TLV <- vector()
    SAT <- vector()
    Distance <- vector()
    Group <- vector()
    Subgrp <- vector()
    Cluster <- vector()
    #Add a nested loop that scans each line within the vector for a cluster to extract 
    #the information for a sequence type within a cluster
    #The first line contains the headings so the loop was started at 2 (see above)
    for (j in 2:length(temp_vec)){
      #Since the loop was started at 2, decrease the index by 1 so as to adjust for it
      j <- j-1
      #temp_line <- vector()
      #Store values to respective vectors after splitting the by \t
      temp_line <- strsplit(temp_vec[j], "\t")[[1]]
      ST[j] <- gsub(" ","", temp_line[1])
      Freq[j] <- temp_line[3]
      SLV[j] <- temp_line[4]
      DLV[j] <- temp_line[5]
      TLV[j] <- temp_line[6]
      SAT[j] <- temp_line[7]
      Distance[j] <- temp_line[8]
      Group[j] <- temp_line[9]
      Subgrp[j] <- temp_line[10]
      Cluster[j] <- i
    }
    #Make a data frame for all the values for a single cluster
    tmp_df <- cbind.data.frame(ST,Freq,SLV,DLV,TLV,SAT,Distance,Group,Subgrp,Cluster)
    #Remove rows with missing values
    #tmp_df <- na.omit(tmp_df)
    #Create a cluster name based on the cluster number
    #df_name <- paste("Cluster",i,sep = "")
    #Assign cluster name (from the line above) to the data frame
    #assign(df_name, tmp_df)
    #Add the cluster data frame to the list (created before the loop)
    clus_df[[i]] <- tmp_df
}
#Combine all cluster data frames from the list into a single data frame
clus_df1 <- do.call(rbind, clus_df)
#The final data frame still contains the headings with variable names
#Remove those lines and save it into a csv file
clus_df1 <- clus_df1[-which(clus_df1$ST == "ST"),]
write.csv(clus_df1, "data_prep/eBURST_MLST_clusters.csv")
