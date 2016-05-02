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
FC9K_CSF <- cbind(FC9K_ClusterSort, ClustFreqMatch$Freq)
FC9KZ <- data.matrix(FC9K_CSF[ ,1:645])

# Zscore matrix creation
# rowSds is dependent on genefilter package
library(genefilter)
# normalize log ratio values and calculate z score
z_score <- (FC9KZ-rowMeans(FC9KZ))/rowSds(FC9KZ)

# Cluster Activation Score Calculation 
# select columns that are Pathology Severity Score >= 1
pathcond <- c(colnames(z_score)[137], colnames(z_score)[142], colnames(z_score)[196], colnames(z_score)[203], colnames(z_score)[268], colnames(z_score)[353], colnames(z_score)[528])
# for all 9 conditions of interest
## 7 are present in our FC matrix
# CARBON.TETRACHLORIDE.1175.mg.kg.3.d.CORN.OIL.100...ORAL.GAVAGE
# CARBON.TETRACHLORIDE.400.mg.kg.7.d.CORN.OIL.100...ORAL.GAVAGE
# CLOTRIMAZOLE.178.mg.kg.5.d.CORN.OIL.100...ORAL.GAVAGE
# CLOTRIMAZOLE.89.mg.kg.3.d.CORN.OIL.100...ORAL.GAVAGE
# ERLOTINIB.58.mg.kg.1.d.WATER.100...ORAL.GAVAGE
# HYDRAZINE.45.mg.kg.5.d.WATER.100...ORAL.GAVAGE
# PHENOBARBITAL.54.mg.kg.1.d.WATER.100...ORAL.GAVAGE"

# 2 are absent
# CLOTRIMAZOLE.178.mg.kg.3.d.CORN.OIL.100...ORAL GAVAGE
## all arrays for this chemical exposure were judged outliers by ArrayQualityMetrics
# MICONAZOL.920.mg.kg.5.d.CORN.OIL.100...ORAL.GAVAGE
## all arrays for this chemical exposure were judged outliers by ArrayQualityMetrics

# subset the z_score matrix for our path conditions of interest
AC_Raw <- z_score[ ,pathcond]
# sum the rows of AC_Raw and add our cluster identifier column and frequency to our subset z_score matrix
AC_Clust <- cbind(rowSums(AC_Raw), ClustFreqMatch)
AC_Calc <- aggregate(AC_Clust[,1], list(Cluster = AC_Clust[,2]), sum)
AC_Calc <- cbind(AC_Calc, CC[ ,2, drop = FALSE])
AC_FINAL <- (1/(AC_Calc$Freq*dim(AC_Raw)[2]))*AC_Calc$x
AC_FINAL <- cbind(CC, AC_FINAL)
