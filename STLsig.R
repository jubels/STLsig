

library(caret)
library(survival)
library(fields)
library(igraph)


######################################  
##### evaluate genepairs  ############
#####################################

source("evaluateGenePairs.r")
#### surv = survivaldata, should contain in order survival time, status and treatmetn variable 
#### gepAll = gene expression data. Column names should match rownames survival data
#### nbrs = how many neighbours should be taken into account when calculating zPFS
#### sets = matrix with one genepair per row. Names genes should match rownames gene expression
#### toi = which treatment in the treatment variable should benefit be predicted for
#### repeats = how many repeats should be used for the consensus network
#### threshold = in which fraction of the top should genes be to be considered synergistic
#### pathways = list containing indexes which genepairs include a certain gene. Name of the list should match genenames in sets object

network = evaluateGenePairs(surv, gepAll , nbrs , sets , toi , repeats, threshold, pathways)

#######################################################
####### select gene networks for signature ############
######################################################

source("makeSignature.r")

##### surv = same survival data as in previous function
##### gepAll = same gene expression data as previous function
##### surv.val = other survival data (fold C in paper) to select gene sets for signature 
#### gep.val = gene epxression data for the patients in surv.val
#### nbrs = how many neighbours should be taken into account when calculating zPFS
#### toi = which treatment in the treatment variable should benefit be predicted for
#### network = output from evaluateGenePairs function
signature = makeSignature(surv, gepAll, surv.val, gep.val, nbrs, toi, network)

#######################################################
############ predict zPFS in hold out data ############
######################################################

source("predictZPFS.r")

###### surv.train = all survival training data (combined surv and surv.val)
###### gep.train = gene expression data for patients in surv.train
###### surv.test = survivaldata for the test set, to return validation performance. Can be missing, zPFS will still be returned
###### gep.test = gene expression for the patients in surv.test
###### signature = output of makeSignature function
#### nbrs = how many neighbours should be taken into account when calculating zPFS
#### toi = which treatment in the treatment variable should benefit be predicted for

predictedZPFS = predictZPFS(surv.train, gep.train, surv.test, gep.test, signature, nbrs, toi)




