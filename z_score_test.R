# z_score
# load in fc_matrix if not already loaded
# fc_final <- as.matrix(read.delim("Workspaces_And_Objects/liver_fc_pc25_nsFilter_rma_qc.txt", row.names=1, header=TRUE, sep="\t"))
# same, but for other workspaces
fc_final <- as.matrix(read.delim("liver_fc_pc25_nsFilter_rma_qc.txt", row.names=1, header=TRUE, sep="\t"))
# if already in DM_Liver_Preprocess and environment is loaded from DM_Liver_Preprocess_18Apr16.RData
fc_z <- fc_final
# skip if importing from file
# fc_z <- fc_z[ ,2:647]
# identifying which cells contain the value NaN
# NaN prevents any numerical calculations
# fc_num <- apply(fc_z, 2, is.nan)
# fc_nan <- which(fc_num, arr.ind = TRUE)
# removing column which contains all NaN
# due to fact that there is no corresponding drug matrix control
# TREATMENT_CHEMICAL == "ROFLUMILAST|17 mg/kg|4 d|CMC .5 %|ORAL GAVAGE"
fc_z <- cbind(fc_z[,1:578],fc_z[,580:646])
# check to see that all values are numbers
fc_num <- apply(fc_z, 2, is.nan)
# if all values are numbers, should return TRUE
all(!is.nan(fc_num))
# which(fc_num, arr.ind = TRUE)
# don't need this matrix anymore, so remove it.
rm(fc_num)
# only need this if environment is loaded from DM_Liver_Preprocess_18Apr16.RData
# fc_z <- data.matrix(fc_z)
# rowSds is dependent on genefilter package
library(genefilter)
# normalize log ratio values and calculate z score
z_score <- (fc_z-rowMeans(fc_z))/rowSds(fc_z)

# activation score
