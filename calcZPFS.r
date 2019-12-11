
calcZPFS = function(surv.z, gep, nbrs = 10, sets, toi = "YES"){
  

  treatment <- rownames(surv.z)[which(surv.z[,3] == toi)]
  other <- rownames(surv.z)[which(surv.z[,3] != toi)]


  ### calculating deltaPFS
  surv.diff <- matrix(NA, length(other), length(treatment))
  for (i in 1:length(treatment)){
    for (j in 1:length(other)){
      if(surv.z[treatment[i],2] == 1 & surv.z[other[j],1] > surv.z[treatment[i],1]){
        surv.diff[j,i] <- surv.z[treatment[i],1] - surv.z[other[j],1] 
      } 
      if(surv.z[other[j],2] == 1 & surv.z[other[j],1] < surv.z[treatment[i],1]){
        surv.diff[j,i] <- surv.z[treatment[i],1] - surv.z[other[j],1] 
      }
    }
  }
  
  
  rownames(surv.diff) <- other
  colnames(surv.diff) <- treatment
  rand.nbr <- array(NA, c(nbrs, length(treatment), 1000)) 
  rand.dpfs <- matrix(NA, 1000, length(treatment))
  #### calculating RPFS
  for (r in 1:1000){
    for (i in 1:length(treatment)){
      rand.nbr[,i,r] <- other[sample(1:length(other), nbrs, replace = FALSE)] 
      rand.dpfs[r,i] <- mean(surv.diff[rand.nbr[,i,r],i], na.rm = TRUE)
    }
  } 
  colnames(rand.dpfs) <- treatment
  
  ##########################################
  ### calculate RPFS from other pov ###
  ##########################################
  treatment.2 <- rownames(surv.z)[which(surv.z[,3] != toi)]
  other.2 <- rownames(surv.z)[which(surv.z[,3] == toi)]
  
  ### calculating deltaPFS
  surv.diff.2 <- matrix(NA, length(other.2), length(treatment.2))
  for (i in 1:length(treatment.2)){
    for (j in 1:length(other.2)){
      if(surv.z[treatment.2[i],2] == 1 & surv.z[other.2[j],1] > surv.z[treatment.2[i],1]){
        surv.diff.2[j,i] <- surv.z[treatment.2[i],1] - surv.z[other.2[j],1] 
      } 
      if(surv.z[other.2[j],2] == 1 & surv.z[other.2[j],1] < surv.z[treatment.2[i],1]){
        surv.diff.2[j,i] <- surv.z[treatment.2[i],1] - surv.z[other.2[j],1] 
      }
    }
  }
  
  
  rownames(surv.diff.2) <- other.2
  colnames(surv.diff.2) <- treatment.2
  rand.nbr.2 <- array(NA, c(nbrs, length(treatment.2), 1000)) 
  rand.dpfs.2 <- matrix(NA, 1000, length(treatment.2))
  #### calculating RPFS
  for (r in 1:1000){
    for (i in 1:length(treatment.2)){
      rand.nbr.2[,i,r] <- other.2[sample(1:length(other.2), nbrs, replace = FALSE)] 
      rand.dpfs.2[r,i] <- mean(surv.diff.2[rand.nbr.2[,i,r],i], na.rm = TRUE)
    }
  } 
  colnames(rand.dpfs.2) <- treatment.2
  ### calculate zPFS for all sets 
  zpfs.all = matrix(NA, nrow(sets), nrow(surv.z))
  colnames(zpfs.all) = rownames(surv.z)
  for(l in 1:nrow(sets)){
    model = rdist(t(gep[sets[l,which(!is.na(sets[l,]))],]))
    rownames(model) <- rownames(surv.z)
    colnames(model) <- rownames(surv.z) 
    
    
    
    distanceA <- model[other,treatment] 
    
    
    
    true.nbr <- matrix(NA, nbrs, length(treatment)) 
    true.dpfs <- rep(NA, length(treatment))
    
    zpfs <- rep(NA, length(treatment))
    
    colnames(true.nbr) <- treatment
    names(true.dpfs) <- treatment
    names(zpfs) <- treatment
    
    for(i in 1:length(treatment)){
      true.nbr[,i] <- rownames(distanceA)[sort(distanceA[,i], decreasing = FALSE, index = TRUE)$ix[2:(nbrs+1)]]
      true.dpfs[i] <- mean(surv.diff[true.nbr[,i],i], na.rm = TRUE)
    } 
    
    #### calculating zPFS
    for (i in 1:length(treatment)){
      zpfs[i] <- (true.dpfs[i] - mean(rand.dpfs[,i], na.rm = TRUE)) / sd(rand.dpfs[,i], na.rm = TRUE) 
    }
    
    zpfs.all[l,names(zpfs)] = zpfs 
    
    
    distanceB <- model[other.2,treatment.2] 
    
    ### calculating deltaPFS
    
    
    true.nbr <- matrix(NA, nbrs, length(treatment.2)) 
    true.dpfs <- rep(NA, length(treatment.2))
    
    zpfs <- rep(NA, length(treatment.2))
    
    colnames(true.nbr) <- treatment.2
    
    names(true.dpfs) <- treatment.2
    names(zpfs) <- treatment.2
    
    for(i in 1:length(treatment.2)){
      true.nbr[,i] <- rownames(distanceB)[sort(distanceB[,i], decreasing = FALSE, index = TRUE)$ix[2:(nbrs+1)]]
      true.dpfs[i] <- mean(surv.diff.2[true.nbr[,i],i], na.rm = TRUE)
    } 
    
    
    #### calculating zPFS
    for (i in 1:length(treatment.2)){
      zpfs[i] <- (true.dpfs[i] - mean(rand.dpfs.2[,i], na.rm = TRUE)) / sd(rand.dpfs.2[,i], na.rm = TRUE) 
    }
    
    zpfs.all[l, names(zpfs)] = -zpfs

  }
 
  return(zpfs.all)
}