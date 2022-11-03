library(cluster)
library(dplyr)
library(reshape2)
library(readr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("./Functions.R")

center=read_delim("Centers.txt", delim = "\t", 
                  escape_double = FALSE, trim_ws = TRUE) %>% data.frame()

#unflag following lines to use personalized clusters

#predefined=read_delim("Centers.txt", delim = "\t",
#                      escape_double = FALSE, trim_ws = TRUE) %>% data.frame()
#center=GetCenter(predefined)
#write.csv(center,file="center.csv",row.names=F)

seeds=read_delim("Seeds.txt", delim = "\t", 
                 escape_double = FALSE, trim_ws = TRUE) %>% data.frame()

test=read_delim("Matrix.txt", delim = "\t", 
                escape_double = FALSE, trim_ws = TRUE) %>% data.frame()

cluster1=cluster2=GetCluster(test,center,seeds)
colnames(cluster1)=c("Sample","Cluster")
write.csv(cluster1,file="cluster.csv",row.names=F)


#unflag folloWing lines to generate clusters including TP53
#TP=read_delim("TP53.txt", delim = "\t", 
#              escape_double = FALSE, trim_ws = TRUE) %>% data.frame()
#cluster2[which(TP$TP53==1),2]="TP53"
#colnames(cluster2)=c("Sample","Cluster")
#write.csv(cluster2,file="cluster_TP53.csv",row.names=F)