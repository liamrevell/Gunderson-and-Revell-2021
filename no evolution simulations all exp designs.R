################ Phenotypic states do not evolve. Experimental Designs 1 and 3 ################
################ so all variation is measurement error.  #############################
pcor<-NULL #phylogenetic correlation coefficients. Experimental design 1
rcor<-NULL #raw correlation coefficients. Experimental design 1
rslp<-NULL #raw slopes. Experimental design 1
pslp<-NULL #phylo slopes. Experimental design 1

cpcor<-NULL #phylogenetic correlation coefficients with data corrected by Kelly and Price method. Experimental design 1
crcor<-NULL #raw correlation coefficients with data corrected by Kelly and Price method. Experimental design 1
crslp<-NULL #raw slopes with data corrected by Kelly and Price method. Experimental design 1
cpslp<-NULL #phylo slopes with data corrected by Kelly and Price method. Experimental design 1

npcor<-NULL #phylogenetic correlation coefficients w independent baseline. Experimental design 3
nrcor<-NULL #raw correlation coefficients  w independent baseline. Experimental design 3
nrslp<-NULL #raw slopes w independent baseline. Experimental design 3
npslp<-NULL #phylo slopes w independent baseline. Experimental design 3
repeat{
  tree<-pbtree(n=50,scale=1) #make a tree
  ctmax<-fastBM(tree, sig2 = 0) # evolution of basal trait along tree
  #ctmaxp<-fastBM(tree,a=1, sig2 = 0.1) # evolution of induced plastic trait across tree 
  ctmaxp<-fastBM(tree,a=1, sig2 = 0) # use this for no plasticity in the trait 
  CTMAX<-lapply(ctmax,sampleFrom,xvar=1,n=20) # For each species randomly sample 20 individual values for baseline and induced pheno
  CTMAXp<-lapply(ctmaxp,sampleFrom,xvar=1,n=20) 
  RN<-mapply(function(x,y) y-x,x=CTMAX,y=CTMAXp,SIMPLIFY=FALSE) ## Compute individual-level reaction norms (difference between CTmax’ and CTmax).
  ctmax.hat<-sapply(CTMAX,mean) # Compute mean basline and induced pheno and the mean reaction norm) for each species to use in interspecific analyses
  ctmaxp.hat<-sapply(CTMAXp, mean)
  rn.hat<-sapply(RN,mean)
  r.slp<-lm(rn.hat~ctmax.hat)$coefficients[2] # raw slope
  rslp<-c(rslp, as.data.frame(r.slp)[1,1])
  p.slp<-lm(pic(rn.hat,tree)~pic(ctmax.hat,tree)+0)$coefficients # pic slope
  pslp<-c(pslp, as.data.frame(p.slp)[1,1])
  alldat<-cbind(as.data.frame(rn.hat), as.data.frame(ctmax.hat)[,1])
  names(alldat)<-c("rn.hat", "ctmax.hat")
  obj<-phyl.vcv(as.matrix(alldat),vcv(tree),1) ## phylogenetic covariance among traits
  p.cor<-cov2cor(obj$R)["rn.hat","ctmax.hat"]## phylo correlation between ct & diff
  pcor<-c(pcor, p.cor)
  r.cor<-cor.test(ctmax.hat, rn.hat, method = "pearson")$estimate # raw correlation
  rcor<-c(rcor, as.data.frame(r.cor)[1,1])
  
  # Adjusting values by Kelly and Price Method
  n<-Ntip(tree) # number of tips on the tree
  X<-cbind(ctmax.hat[tree$tip.label],
           ctmaxp.hat[tree$tip.label]) # mean tip values for each trait for each species
  colnames(X)<-c("ctmax1","ctmax2")
  object<-phyl.vcv(X,vcv(tree),lambda=1) # phylogenetic variance/covariance matrix
  s1<-sqrt(object$R[1,1]) # square root of phylogenetic variance in ctmax1
  s2<-sqrt(object$R[2,2]) # square root of phylogenetic variance in ctmax2
  r<-cov2cor(object$R)[1,2] # phylogenetic correlation coefficeint bewteen ctmax1 and ctmax2
  T<-sqrt(n-2)*(s1/s2-s2/s1)/(2*(1-r^2)) # T statistic for  Pittman's  test of homogeneity of variances
  Pval<-2*pt(abs(T),df=n-2,lower.tail=FALSE) # P value for T statistic
  rho.hat<-if(Pval>=0.05) r else 2*r*s1*s2/(s1^2+s2^2) # rho value for phenotype correction based on whether or not variances are equal 
  Xhat<-object$alpha[,1] # phylopgenetic means for ctmax1 and ctmax2
  Dstar<-rho.hat*(X[,1]-Xhat[1])-(X[,2]-Xhat[2]) # corrected reaction norm values
  cr.slp<-lm(Dstar~ctmax.hat)$coefficients[2] # raw slope corrected rn versus baseline
  crslp<-c(crslp, as.data.frame(cr.slp)[1,1])
  cp.slp<-lm(pic(Dstar,tree)~pic(ctmax.hat,tree)+0)$coefficients # pic slope with corrected rn values
  cpslp<-c(cpslp, as.data.frame(cp.slp)[1,1])
  c.alldat<-cbind(as.data.frame(Dstar), as.data.frame(ctmax.hat)[,1]) # corrected data frame for correlation analysis
  names(c.alldat)<-c("Dstar", "ctmax.hat")
  c.obj<-phyl.vcv(as.matrix(c.alldat),vcv(tree),1) ## phylogenetic covariance among traits
  cp.cor<-cov2cor(c.obj$R)["Dstar","ctmax.hat"]## phylo correlation between ct & diff corrected
  cpcor<-c(cpcor, cp.cor)
  cr.cor<-cor.test(ctmax.hat, Dstar, method = "pearson")$estimate # raw correlation
  crcor<-c(crcor, as.data.frame(cr.cor)[1,1])
  
  # same analysis, but sample different individuals for baseline pheno
  nCTMAX<-lapply(ctmax,sampleFrom,xvar=1,n=20) # sample new individuals for independent baseline
  nctmax.hat<-sapply(nCTMAX,mean) # mean of new values for each spp
  nr.slp<-lm(rn.hat~nctmax.hat)$coefficients[2] # raw slope of rxn norm vs new baseline
  nrslp<-c(nrslp, as.data.frame(nr.slp)[1,1])
  np.slp<-lm(pic(rn.hat,tree)~pic(nctmax.hat,tree)+0)$coefficients # pic slope of rxn norm versus new baseline
  npslp<-c(npslp, as.data.frame(np.slp)[1,1])
  nalldat<-cbind(as.data.frame(rn.hat), as.data.frame(nctmax.hat)[,1])
  names(nalldat)<-c("rn.hat", "nctmax.hat")
  nobj<-phyl.vcv(as.matrix(nalldat),vcv(tree),1) ## phylogenetic covariance among traits w new baseline
  np.cor<-cov2cor(nobj$R)["rn.hat","nctmax.hat"]## phylo correlation between ct & diff w new baseline
  npcor<-c(npcor, np.cor)
  nr.cor<-cor.test(nctmax.hat, rn.hat, method = "pearson")$estimate # raw correlation w new baseline
  nrcor<-c(nrcor, as.data.frame(nr.cor)[1,1])
  if (length(nrcor) == 1000) {break}
}
sumdat<-as.data.frame(cbind(pcor, pslp,rcor, rslp, cpcor, cpslp, crcor, crslp, npcor, npslp, nrcor, nrslp))
head(sumdat)

################ Phenotypic states do not evolve.Experimental Design 2 ################
################ so all variation is measurement error. Different individ measured in each env  #############################
pcor<-NULL #phylogenetic correlation coefficients
rcor<-NULL #raw correlation coefficients
rslp<-NULL #raw slopes
pslp<-NULL #phylo slopes
cpcor<-NULL #phylogenetic correlation coefficients with data corrected by Kelly and Price method
crcor<-NULL #raw correlation coefficients with data corrected by Kelly and Price method
crslp<-NULL #raw slopes with data corrected by Kelly and Price method
cpslp<-NULL #phylo slopes with data corrected by Kelly and Price method
repeat{
  tree<-pbtree(n=50,scale=1) #make a tree
  ctmax<-fastBM(tree, sig2 = 0) # evolution of basal trait along tree
  #ctmaxp<-fastBM(tree,a=1, sig2 = 0) # evolution of induced plastic trait across tree 
  ctmaxp<-fastBM(tree, sig2 = 0) # use if there is no plasticity 
  CTMAX<-lapply(ctmax,sampleFrom,xvar=0.01,n=20) # For each species randomly sample 20 individual values for baseline and induced pheno
  CTMAXp<-lapply(ctmaxp,sampleFrom,xvar=0.01,n=20) 
  ctmax.hat<-sapply(CTMAX,mean) # Compute mean basline and induced pheno and the mean reaction norm) for each species to use in interspecific analyses
  ctmaxp.hat<-sapply(CTMAXp, mean)
  rn.hat<-ctmaxp.hat-ctmax.hat
  #rn.hat<-sapply(RN,mean)
  r.slp<-lm(rn.hat~ctmax.hat)$coefficients[2] # raw slope
  rslp<-c(rslp, as.data.frame(r.slp)[1,1])
  p.slp<-lm(pic(rn.hat,tree)~pic(ctmax.hat,tree)+0)$coefficients # pic slope
  pslp<-c(pslp, as.data.frame(p.slp)[1,1])
  alldat<-cbind(as.data.frame(rn.hat), as.data.frame(ctmax.hat)[,1])
  names(alldat)<-c("rn.hat", "ctmax.hat")
  obj<-phyl.vcv(as.matrix(alldat),vcv(tree),1) ## phylogenetic covariance among traits
  p.cor<-cov2cor(obj$R)["rn.hat","ctmax.hat"]## phylo correlation between ct & diff
  pcor<-c(pcor, p.cor)
  r.cor<-cor.test(ctmax.hat, rn.hat, method = "pearson")$estimate # raw correlation
  rcor<-c(rcor, as.data.frame(r.cor)[1,1])
  
  # Adjusting values by Kelly and Price Method
  n<-Ntip(tree) # number of tips on the tree
  X<-cbind(ctmax.hat[tree$tip.label],
           ctmaxp.hat[tree$tip.label]) # mean tip values for each trait for each species
  colnames(X)<-c("ctmax1","ctmax2")
  object<-phyl.vcv(X,vcv(tree),lambda=1) # phylogenetic variance/covariance matrix
  s1<-sqrt(object$R[1,1]) # square root of phylogenetic variance in ctmax1
  s2<-sqrt(object$R[2,2]) # square root of phylogenetic variance in ctmax2
  r<-cov2cor(object$R)[1,2] # phylogenetic correlation coefficeint bewteen ctmax1 and ctmax2
  T<-sqrt(n-2)*(s1/s2-s2/s1)/(2*(1-r^2)) # T statistic for  Pittman's  test of homogeneity of variances
  Pval<-2*pt(abs(T),df=n-2,lower.tail=FALSE) # P value for T statistic
  rho.hat<-if(Pval>=0.05) r else 2*r*s1*s2/(s1^2+s2^2) # rho value for phenotype correction based on whether or not variances are equal 
  Xhat<-object$alpha[,1] # phylopgenetic means for ctmax1 and ctmax2
  Dstar<-rho.hat*(X[,1]-Xhat[1])-(X[,2]-Xhat[2]) # corrected reaction norm values
  cr.slp<-lm(Dstar~ctmax.hat)$coefficients[2] # raw slope corrected rn versus baseline
  crslp<-c(crslp, as.data.frame(cr.slp)[1,1])
  cp.slp<-lm(pic(Dstar,tree)~pic(ctmax.hat,tree)+0)$coefficients # pic slope with corrected rn values
  cpslp<-c(cpslp, as.data.frame(cp.slp)[1,1])
  c.alldat<-cbind(as.data.frame(Dstar), as.data.frame(ctmax.hat)[,1]) # corrected data frame for correlation analysis
  names(c.alldat)<-c("Dstar", "ctmax.hat")
  c.obj<-phyl.vcv(as.matrix(c.alldat),vcv(tree),1) ## phylogenetic covariance among traits
  cp.cor<-cov2cor(c.obj$R)["Dstar","ctmax.hat"]## phylo correlation between ct & diff corrected
  cpcor<-c(cpcor, cp.cor)
  cr.cor<-cor.test(ctmax.hat, Dstar, method = "pearson")$estimate # raw correlation
  crcor<-c(crcor, as.data.frame(cr.cor)[1,1])
  if (length(crcor) == 1000) {break}
}

sumdat<-as.data.frame(cbind(pcor, pslp,rcor, rslp, cpcor, cpslp, crcor, crslp))



