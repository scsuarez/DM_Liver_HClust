# DrugMatrix-Hierarchical cluster/Dyntreecut processing script (DoD BHSAI)
# Modified by SCS

library(annotate)
library(affy)
library("Biobase")
library("rat2302.db")
library("dynamicTreeCut")
library("moduleColor")
library(cluster)

###########
#HC-Pearson
###########
# load table
FC9K<-as.matrix(read.delim("liver_fc_pc25_nsFilter_rma_qc.txt", row.names=1, header=TRUE, sep="\t"))
# transpose table. Since "cor" function operates on columns
tFC9K<-t(FC9K)
# calculate correlation matrix
# changed complete.pairwise.obs to everything
# had to use complete.pairwise.obs due to presence of one column of NaN
## unknown at the time of initial code running, now fixed in DM_liver_preprocess
CORFC9K_P <- cor(tFC9K, y=NULL, use= "everything", method = "pearson")
# CORFC9K_P<-cor(tFC9K, method="pearson")
# convert correlation into distance using "1-cor" and then get the distance matrix
DISTFC9K_P<-as.dist(1-CORFC9K_P)
# finally do hierarchical clustering using distance matrix and "average" linkage method
HCFC9K_PA<-hclust(DISTFC9K_P, method="average")
# plot dendrogram
plot(HCFC9K_PA,labels=FALSE)
##################
# DYNAMIC TREE CUT
##################
DISTFC9K_Pm<-as.matrix(DISTFC9K_P)
DHCUT_FC9KPA2<-cutreeDynamic(dendro= HCFC9K_PA, distM= DISTFC9K_Pm, method="hybrid", cutHeight=0.95, minClusterSize=16,deepSplit=3, pamStage=TRUE,  respectSmallClusters=TRUE)
table(DHCUT_FC9KPA2)
#******OUT***********
FC9KDHCUT_PA2<-data.frame(FC9K,DHCUT_FC9KPA2)#FCDATAFRAME_WITHCLUSTINFO
write.table(FC9KDHCUT_PA2,file="FC9KDHCUT_PA2out.txt",sep="\t",row.names=TRUE)
#************************
###############################

#initial Cluster analysis for SOT
Clusters <- FC9KDHCUT_PA2[ ,647]
names(Clusters) <- rownames(FC9KDHCUT_PA2)
ClusterSort <- sort(Clusters)

#convert Cluster Probe IDs to gene information
PROBES<- as.character(names(ClusterSort))
# Use select and chip.db to extract ids (instead of get/mget).Note: one to many mapping occurs
AOUT <- select(rat2302.db, PROBES, c("SYMBOL","ENTREZID", "GENENAME"))
DUP_AOUT <- AOUT[duplicated(AOUT$PROBEID)|duplicated(AOUT$PROBEID,fromLast=TRUE),]
DUP_AOUT_PROBE <- DUP_AOUT$PROBEID
AOUT_UNIQ <- AOUT[!(AOUT$PROBEID %in% DUP_AOUT_PROBE),]
# in duplicated, remove LOC probes + duplicate copy of probes
# Note: LOC probes are probes that map to gene symb starting with LOC... (duplicate).Same probe also map to another gene symb
DUP_AOUT1<-DUP_AOUT[!grepl("LOC",DUP_AOUT$SYMBOL)&!duplicated(DUP_AOUT$PROBEID),]
#rbind DUP_AOUT1 to AOUT_UNIQ
PROB_RID_OUT <- rbind(AOUT_UNIQ,DUP_AOUT1)
# bind our cluster info to gene information
ClusterGeneSort <- cbind(PROB_RID_OUT, ClusterSort)
write.table(ClusterGeneSort,file="ClusterGeneSort.txt",sep="\t",row.names=FALSE)





