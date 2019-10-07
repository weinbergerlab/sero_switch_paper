rm(list = ls())
library(igraph)

#set working directory
setwd("D:\\ss_paper")
#Load network data file
network_data <- read.csv("index/ss_index.csv")
network_data <- network_data[which(network_data$total >= 5),]

#Reorder the sequence of columns in network data with lhs, rhs and total
#the index total values will be assigned as weights in the network
network_data <- network_data[,c(3,4,17)]
colnames(network_data) <- c("from","to","weight")
network_data$from <- as.character(network_data$from)
network_data$to <- as.character(network_data$to)

#check if there are any pairs where serotypes are the same
network_data[which(network_data$from == network_data$to),]
#Check if the interactions are within or outside serogroup
network_data$in_sero <- "inside"
network_data$from_serogroup <- gsub(c("[A-Z]"),"",network_data$from)
network_data$from_serogroup <- gsub("/","",network_data$from_serogroup)
network_data$to_serogroup <- gsub(c("[A-Z]"),"",network_data$to)
network_data$to_serogroup <- gsub("/","",network_data$to_serogroup)

for(i in 1:nrow(network_data)){
  if (network_data$from_serogroup[i] != network_data$to_serogroup[i]){
    network_data$in_sero[i] <- "outside"
  }
}

#Make a new data frame with details of serotypes
uni_from <- unique(network_data$from)
uni_to <- unique(network_data$to)
serogroup <- data.frame(unique(c(uni_from,uni_to)))
colnames(serogroup)[1] <- "serotype"
#Remove alphabets from serotype names to assign serogroups
serogroup$serogroup <- gsub("[A-Z]","",serogroup$serotype)
serogroup$serogroup <- gsub("/","",serogroup$serogroup)

#Step by step
nodes <- serogroup
links <- network_data
nrow(nodes)
length(unique(nodes$serotype))
nrow(links)
nrow(unique(links[,c("from","to")]))

net <- graph.data.frame(links, nodes, directed = F)
net <- simplify(net, remove.multiple = F, remove.loops = T)

E(net)
V(net)$serogroup

#Set node size on the basis of number of degrees(links) for a serotype
deg <- degree(net, mode = "all")
V(net)$size <- 15
#Set wdge width based on weight
E(net)$length <- 30
E(net)$width <- E(net)$weight

l <- layout.auto(net)
plot(net, edge.color=c("red","green")[(E(net)$in_sero=="inside")+1],layout = l)

write_graph(net, "plots/index_5_or_more.graphml", format = "graphml")
