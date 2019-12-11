predictZPFS = function(surv.train, gep.train, surv.test =NA,  gep.test, signature, nbrs, toi, minSize){
  source("calcZPFS.r")
  zpfsVal = calcZPFS(surv.z = surv.train, gep =  gep.train, sets = signature, nbrs = nbrs, toi = toi)
  
  zpfsValInh = matrix(NA, nrow(signature), ncol(gep.test))
  for(i in 1:nrow(signature)){
    genes = signature[i,which(!is.na(signature[i,]))]
    distRef = rdist(t(cbind(gep.train[genes,], gep.test[genes,])))
    rownames(distRef) = c(colnames(gep.train), colnames(gep.test))
    colnames(distRef) = c(colnames(gep.train), colnames(gep.test))
    distRef = distRef[colnames(gep.test), rownames(surv.train)]
    zPFSUse = zpfsVal[i, which(!is.na(zpfsVal[i,]))]
    distRefUse = distRef[,which(!is.na(zpfsVal[i,]))]
    for(p in 1:ncol(gep.test)){
      weightsUse = (1/(sort(distRefUse[p,], decreasing = F)))
      zpfsValInh[i,p] = sum(zPFSUse[colnames(distRefUse)[sort(distRefUse[p,], decreasing = F, index = T)$ix]]*weightsUse,na.rm = T)/sum(weightsUse, na.rm = T)
      
    }
  }
  
  predZPFS = apply(zpfsValInh, 2, sum)
  names(predZPFS) = colnames(gep.test)
  
  results = predZPFS
  
  if(!is.na(surv.test)){
    surv.test$score = predZPFS
    HRsT = matrix(NA, 6,round((1-minSize)*nrow(surv.test)))
    count = 0
    for(i in sort(unique(surv.test$score), decreasing = F)[round(minSize*nrow(surv.test)):round((1-minSize)*nrow(surv.test))]){
      count = count +1
      surv.test$solution = rep(0, nrow(surv.test))
      surv.test$solution[which(surv.test$score > i)] = 1
      HRsT[1,count] = summary(coxph(Surv(surv.test[,1], surv.test[,2])~surv.test[,3], data = surv.test, subset = solution == 1))$coef[1]
      HRsT[2,count] = summary(coxph(Surv(surv.test[,1], surv.test[,2])~surv.test[,3], data = surv.test, subset = solution == 1))$coef[5]
      HRsT[3,count] = summary(coxph(Surv(surv.test[,1], surv.test[,2])~surv.test[,3], data = surv.test, subset = solution == 0))$coef[1]
      HRsT[4,count] = summary(coxph(Surv(surv.test[,1], surv.test[,2])~surv.test[,3], data = surv.test, subset = solution == 0))$coef[5]
      HRsT[5,count] = mean(surv.test$solution)
      HRsT[6,count] = i
    }
    
    diffHR = HRsT[3,] - HRsT[1,]
    threshold = HRsT[6,which(diffHR == max(diffHR, na.rm = T))]
    solution = rep("no benefit", nrow(surv.test))
    solution[which(predZPFS > threshold)] = "benefit"
    
    results = list(predZPFS, threshold, solution)
    names(results) = c("predictedZPFS", "threshold", "classification")
  }
  return(results)
  
  
}