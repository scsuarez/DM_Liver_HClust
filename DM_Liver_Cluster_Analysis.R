# Code for Z-score creation and Activation Score Analysis
# should be run after DM_Liver_HC_SCS.R

# read in files if not already in environment
# either need liver_fc_pc25_nsFilter_rma_qc.txt
## original log matrix out
# or FC9KDHCUT_PA2out.txt
## log matrix text file with cluster information added

# for ease use FC9KDHCUT_PA2out.txt as FC9KDHCUT_PA2

# sort by cluster
FC9K_ClusterSort <- FC9KDHCUT_PA2[order(DHCUT_FC9KPA2), ]
# create a frequency table 
ClustCount <- table(FC9K_ClusterSort[ ,646])
# make the frequency table into a data frame
CC <- as.data.frame(ClustCount)
# rename the columns of that DF
colnames(CC) <- c("Cluster", "Freq")

# function to create a new data frame based on FC9K_ClusterSort cluster identifier that replicates the frequency of cluster in a new data frame
expand_cc <- function(x){
  LS <- sapply(1:nrow(x), function(i) x[rep(i, each = x$Freq[i]), ], simplify = FALSE) 
  DF <- do.call("rbind", LS)
  rownames(DF) <- NULL
return (DF)
}

# use above function to create a data frame that has equal number of rows to our FC matrix and has the cluster identifier 
ClustFreqMatch <- expand_cc(CC)
rownames(ClustFreqMatch) <- rownames(FC9K_ClusterSort)







