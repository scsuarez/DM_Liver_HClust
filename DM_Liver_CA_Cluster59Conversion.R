# load packages needed for probe conversion to genes
library("rat2302.db")
library(annotate)
# file for conversion of Rat entrez gene id into human entrez gene id
# utilizes Bioconductor package biomaRt
library("biomaRt")
# utilizes complicated call structure, but this is the order you do it in
# generic - can be used for any of the list functions of biomaRt
# ensembl_us_east <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="useast.ensembl.org")
# if you want to use any of the functions of biomaRt you need to know, 
## the mart (database) you want to query
### obtained by listMarts(host="useast.ensembl.org")
## the dataset 
### available using listDatasets()
## the host
### default host ="www.ensembl.org" is no longer functional
### use a mirror like host="useast.ensembl.org"
# for joining and df maniulation
library(plyr)
# for human gene annotation
library("org.Hs.eg.db")
# making database simplification easier
library(sqldf)

# use the ensembl database for human genes for one part of our query
human <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org")
# use the ensembl database for rat genes for the first part of our query
rat<-useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="rnorvegicus_gene_ensembl", host="useast.ensembl.org")

# select info from fold change matrix for Cluster 59, only cluster with absolute value of average Z score for conditions of interest above 1
# means average of gene z scores in the cluster is greater than 1 SD from the gene mean
sd_Clust <- subset(FC9K_CSF, FC9K_CSF$DHCUT_FC9KPA2 = "59")

# convert probes to Rat Genes
PROBES<- as.character(row.names(sd_Clust))
# Use select and chip.db to extract ids (instead of get/mget).Note: one to many mapping occurs if we had not filtered for entrezid duplicates
OUT <- select(rat2302.db, PROBES, "ENTREZID")
colnames(OUT) <- c("PROBES", "ENTREZ_GENE_ID_RAT")
# get information on 2 linked datasets - homology mapping
# attributes available by listAttributes()
# filters available by listFilters()
RHIDs <- getLDS(attributes=c("entrezgene"), filters="entrezgene", values= OUT$ENTREZID, mart=rat, attributesL=c("entrezgene"), martL=human, verbose = TRUE, uniqueRows = FALSE)
# set names of input and output columns
colnames(RHIDs) <- c("ENTREZ_GENE_ID_RAT", "ENTREZ_GENE_ID_HUMAN")
# join the 2 data frames including all values in count.UPSIG.entrez and matching values in RHIDs_UPSIG
RH <- join(OUT, RHIDs, by = "ENTREZ_GENE_ID_RAT", type = "left")
RH[ ,3] <- as.character(RH[ ,3])
# get gene annotation information
h.entrezid <- as.character(RH$ENTREZ_GENE_ID_HUMAN[!is.na(RH$ENTREZ_GENE_ID_HUMAN)])
h.entrezsymbol  <- as.character(unlist(mget(h.entrezid, envir=org.Hs.egSYMBOL, ifnotfound=NA)))
h.entrezgene  <- as.character(unlist(mget(h.entrezid, envir=org.Hs.egGENENAME, ifnotfound=NA)))

# bind entrez id, symbol, and gene information into 1 df
h.entrezinfo <- as.data.frame(cbind(h.entrezid, h.entrezsymbol, h.entrezgene))
# set colnames of df so it can be joined
colnames(h.entrezinfo) <- c("ENTREZ_GENE_ID_HUMAN", "SYMBOL", "GENENAME")
# set column classes to character for join
for(i in 1:dim(h.entrezinfo)[2]){
  h.entrezinfo[ ,i] <- as.character(h.entrezinfo[ ,i])
}
# join df with human entrez information to our table with probes, frequency, rat entrezid
RHcomplete <- join(RH, h.entrezinfo, by = "ENTREZ_GENE_ID_HUMAN", type = "left")
# use SQL and sqldf to subset for distinct cases (no repeated information, which happens for some reason with join)
RHcomplete <- sqldf('SELECT DISTINCT * FROM [RHcomplete]')

# write significant up genes with counts to delimited text file
write.table(RHcomplete ,file = "Analysis/DM_HClust59_entrez.txt",sep="\t",row.names=FALSE)

# write probe and human entrez id to file for DAVID web analysis
write(RHcomplete$PROBEID, file = "Analysis/DM_HClust59_probes.txt")
write(h.entrezid, file = "Analysis/DM_HClust59_hentrezid.txt")