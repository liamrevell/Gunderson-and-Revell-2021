########################### alternative model with evolution. Experimental Design 2 ###################
######################### phenotypes measured with and without error, and Kelly and Price correction ###########################################################
tpcor<-NULL #phylogenetic correlation coefficients without error
trcor<-NULL #raw correlation coefficients without error
trslp<-NULL #raw slopes without error
tpslp<-NULL #phylo slopes without error
pcor<-NULL #phylogenetic correlation coefficients
rcor<-NULL #raw correlation coefficients
rslp<-NULL #raw slopes
pslp<-NULL #phylo slopes
cpcor<-NULL #phylogenetic correlation coefficients with data corrected by Kelly and Price method
crcor<-NULL #raw correlation coefficients with data corrected by Kelly and Price method
crslp<-NULL #raw slopes with data corrected by Kelly and Price method
cpslp<-NULL #phylo slopes with data corrected by Kelly and Price method
npcor<-NULL # phylo cor different individ with error
nrcor<-NULL # raw cor different individ with error
nrslp<-NULL # raw slope different indiv with error
npslp<-NULL # phylo slope different indiv with error

repeat{
  tree<-pbtree(n=50,scale=1) #make a tree
  ctmax<-fastBM(tree, sig2 = 1) # evolution of basal trait along tree
  ctmaxp<-fastBM(tree,a=1, sig2 = 1) # evolution of induced plastic trait across tree 
  
  # analysis with true phenotypic values
  ctmax.hat<-sapply(ctmax,mean) # Compute mean basline and induced pheno and the mean reaction norm) for each species to use in interspecific analyses
  ctmaxp.hat<-sapply(ctmaxp, mean)
  rn.hat<-ctmaxp.hat-ctmax.hat
  tr.slp<-lm(rn.hat~ctmax.hat)$coefficients[2] # raw slope
  trslp<-c(trslp, as.data.frame(tr.slp)[1,1])
  tp.slp<-lm(pic(rn.hat,tree)~pic(ctmax.hat,tree)+0)$coefficients # pic slope
  tpslp<-c(tpslp, as.data.frame(tp.slp)[1,1])
  talldat<-cbind(as.data.frame(rn.hat), as.data.frame(ctmax.hat))
  names(talldat)<-c("rn.hat", "ctmax.hat")
  tobj<-phyl.vcv(as.matrix(talldat),vcv(tree),1) ## phylogenetic covariance among traits
  tp.cor<-cov2cor(tobj$R)["ctmax.hat","rn.hat"]## phylo correlation between ct & diff
  tpcor<-c(tpcor, tp.cor)
  tr.cor<-cor.test(ctmax.hat, rn.hat, method = "pearson")$estimate # raw correlation
  trcor<-c(trcor, as.data.frame(tr.cor)[1,1])
  
  # analysis with phenotypiic values measured with error
  CTMAX<-lapply(ctmax,sampleFrom,xvar=1,n=20) # For each species randomly sample 20 individual values for baseline and induced pheno
  CTMAXp<-lapply(ctmaxp,sampleFrom,xvar=1,n=20) 
  CTMAX.hat<-sapply(CTMAX,mean) # Compute mean basline and induced pheno and the mean reaction norm) for each species to use in interspecific analyses
  CTMAXp.hat<-sapply(CTMAXp, mean)
  RN.hat<-CTMAXp.hat-CTMAX.hat
  r.slp<-lm(RN.hat~CTMAX.hat)$coefficients[2] # raw slope
  rslp<-c(rslp, as.data.frame(r.slp)[1,1])
  p.slp<-lm(pic(RN.hat,tree)~pic(CTMAX.hat,tree)+0)$coefficients # pic slope
  pslp<-c(pslp, as.data.frame(p.slp)[1,1])
  alldat<-cbind(as.data.frame(RN.hat), as.data.frame(CTMAX.hat)[,1])
  names(alldat)<-c("rn.hat", "ctmax.hat")
  obj<-phyl.vcv(as.matrix(alldat),vcv(tree),1) ## phylogenetic covariance among traits
  p.cor<-cov2cor(obj$R)["rn.hat","ctmax.hat"]## phylo correlation between ct & diff
  pcor<-c(pcor, p.cor)
  r.cor<-cor.test(CTMAX.hat, RN.hat, method = "pearson")$estimate # raw correlation
  rcor<-c(rcor, as.data.frame(r.cor)[1,1])
  
  # Adjusting values with error by Kelly and Price Method
  n<-Ntip(tree) # number of tips on the tree
  X<-cbind(CTMAX.hat[tree$tip.label],
           CTMAXp.hat[tree$tip.label]) # mean tip values for each trait for each species
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
  cr.slp<-lm(Dstar~CTMAX.hat)$coefficients[2] # raw slope corrected rn versus baseline
  crslp<-c(crslp, as.data.frame(cr.slp)[1,1])
  cp.slp<-lm(pic(Dstar,tree)~pic(CTMAX.hat,tree)+0)$coefficients # pic slope with corrected rn values
  cpslp<-c(cpslp, as.data.frame(cp.slp)[1,1])
  c.alldat<-cbind(as.data.frame(Dstar), as.data.frame(CTMAX.hat)[,1]) # corrected data frame for correlation analysis
  names(c.alldat)<-c("Dstar", "ctmax.hat")
  c.obj<-phyl.vcv(as.matrix(c.alldat),vcv(tree),1) ## phylogenetic covariance among traits
  cp.cor<-cov2cor(c.obj$R)["Dstar","ctmax.hat"]## phylo correlation between ct & diff corrected
  cpcor<-c(cpcor, cp.cor)
  cr.cor<-cor.test(CTMAX.hat, Dstar, method = "pearson")$estimate # raw correlation
  crcor<-c(crcor, as.data.frame(cr.cor)[1,1])
  if (length(crcor) == 1000) {break}
}
sumdat<-as.data.frame(cbind(tpcor, tpslp, trcor, trslp, pcor, pslp,rcor, rslp, cpcor, cpslp, crcor, crslp))
