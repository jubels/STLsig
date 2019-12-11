makeSignature = function(surv, surv.val, gep, gep.val, network, nbrs,toi){
  graphBestPairs = graph_from_edgelist(network[[1]][,1:2], directed = FALSE)
  geneNetworks = components(graphBestPairs)
  clusterSets = matrix(NA, length(geneNetworks[[2]]), max(geneNetworks[[2]]))
  
  for(p in 1:nrow(clusterSets)){
    clusterSets[p,1:length(which(geneNetworks[[1]] == p))] = names(geneNetworks[[1]])[which(geneNetworks[[1]] == p)]
  }
  source("calcZPFS.r")
  zpfsTrain =  calcZPFS(surv.z = surv, gep =  gep, sets = clusterSets, nbrs = nbrs, toi = toi)
  
  
  
  HRsTrain = matrix(NA, 4, nrow(zpfsTrain))
  
  for(i in 1:nrow(zpfsTrain)){
    surv$solution = rep(0, nrow(surv))
    surv$solution[which(zpfsTrain[i,] > summary(zpfsTrain[i,])[5])] = 1
    HRsTrain[1,i] = summary(coxph(Surv(surv[,1],surv[,2])~surv[,3], data = surv, subset = solution == 1))$coef[1]
    HRsTrain[2,i] = summary(coxph(Surv(surv[,1],surv[,2])~surv[,3], data = surv, subset = solution == 1))$coef[5]
    HRsTrain[3,i] = summary(coxph(Surv(surv[,1],surv[,2])~surv[,3], data = surv, subset = solution == 0))$coef[1]
    HRsTrain[4,i] = summary(coxph(Surv(surv[,1],surv[,2])~surv[,3], data = surv, subset = solution == 0))$coef[5]
  }
  

  zpfsTop = calcZPFS(surv.z = surv.val, gep =  gep.val, sets = clusterSets, nbrs = nbrs, toi = toi)
  
  orderSets = sort((HRsTrain[3,]-HRsTrain[1,]), decreasing = T, index = T)$ix
  HRsFF = matrix(NA, 4, length(orderSets))
  for(i in 2:length(orderSets)){
    surv.val$score = apply(zpfsTop[orderSets[1:i],],2,sum)
    surv.val$solution = rep(0, nrow(surv.val))
    surv.val$solution[which(surv.val$score > summary(surv.val$score)[5])] = 1
    HRsFF[1,i] = summary(coxph(Surv(surv.val[,1], surv.val[,2])~surv.val[,3], data = surv.val, subset = solution == 1))$coef[1]
    HRsFF[2,i] = summary(coxph(Surv(surv.val[,1], surv.val[,2])~surv.val[,3], data = surv.val, subset = solution == 1))$coef[5]
    HRsFF[3,i] = summary(coxph(Surv(surv.val[,1], surv.val[,2])~surv.val[,3], data = surv.val, subset = solution == 0))$coef[1]
    HRsFF[4,i] = summary(coxph(Surv(surv.val[,1], surv.val[,2])~surv.val[,3], data = surv.val, subset = solution == 0))$coef[5]
    
  }
  use = orderSets[1:which((HRsFF[3,] - HRsFF[1,]) == max((HRsFF[3,] - HRsFF[1,]), na.rm = T))]
  signature = clusterSets[use,]
  return(signature)
  
}