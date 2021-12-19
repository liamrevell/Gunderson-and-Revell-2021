########### What happens when modelled in the "classical" sense of slope & intercept, experimental design 2 ######
########### Individual phenotypes measured with error ################
### With evolution in this example
tpcor<-NULL #phylogenetic correlation coefficients
trcor<-NULL #raw correlation coefficients
trslp<-NULL #raw slopes
tpslp<-NULL #phylo slopes
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
  intercept<-fastBM(tree, sig2 = 1) # evolution of basal trait along tree
  slope<-fastBM(tree,a=1, sig2 = 1) # evolution of induced plastic trait across tree 
  e1.est<-lapply(0*slope+intercept,sampleFrom,xvar=1,n=20) # estimated pheno in env 1 w error 
  e2.est<-lapply(1*slope+intercept,sampleFrom,xvar=1,n=20) # estimated pheno in env 2 w error
  e1.est.hat<-sapply(e1.est,mean) # Compute mean basline and induced pheno and the mean reaction norm) for each species to use in interspecific analyses
  e2.est.hat<-sapply(e2.est, mean)
  rn.hat<-e2.est.hat-e1.est.hat
  tr.slp<-lm(slope~intercept)$coefficients[2] # true raw slope
  trslp<-c(trslp, as.data.frame(tr.slp)[1,1])
  tp.slp<-lm(pic(slope,tree)~pic(intercept,tree)+0)$coefficients # true pic slope
  tpslp<-c(tpslp, as.data.frame(tp.slp)[1,1])
  r.slp<-lm(rn.hat~e1.est.hat)$coefficients[2] # raw slope
  rslp<-c(rslp, as.data.frame(r.slp)[1,1])
  p.slp<-lm(pic(rn.hat,tree)~pic(e1.est.hat,tree)+0)$coefficients # pic slope
  pslp<-c(pslp, as.data.frame(p.slp)[1,1])
  alldat<-cbind(as.data.frame(rn.hat), as.data.frame(e1.est.hat))
  names(alldat)<-c("slope", "intercept")
  obj<-phyl.vcv(as.matrix(alldat),vcv(tree),1) ## phylogenetic covariance among traits
  p.cor<-cov2cor(obj$R)["intercept","slope"]## phylo correlation between ct & diff
  pcor<-c(pcor, p.cor)
  r.cor<-cor.test(e1.est.hat, rn.hat, method = "pearson")$estimate # raw correlation
  rcor<-c(rcor, as.data.frame(r.cor)[1,1])
  
  talldat<-cbind(as.data.frame(slope), as.data.frame(intercept))
  names(talldat)<-c("slope", "intercept")
  tobj<-phyl.vcv(as.matrix(talldat),vcv(tree),1) ## phylogenetic covariance among traits
  tp.cor<-cov2cor(tobj$R)["intercept","slope"]## phylo correlation between ct & diff
  tpcor<-c(tpcor, tp.cor)
  tr.cor<-cor.test(intercept, slope, method = "pearson")$estimate # raw correlation
  trcor<-c(trcor, as.data.frame(tr.cor)[1,1])
  
  # Adjusting values by Kelly and Price Method
  n<-Ntip(tree) # number of tips on the tree
  X<-cbind(e1.est.hat[tree$tip.label],
           e2.est.hat[tree$tip.label]) # mean tip values for each trait for each species
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
  cr.slp<-lm(Dstar~e1.est.hat)$coefficients[2] # raw slope corrected rn versus baseline
  crslp<-c(crslp, as.data.frame(cr.slp)[1,1])
  cp.slp<-lm(pic(Dstar,tree)~pic(e1.est.hat,tree)+0)$coefficients # pic slope with corrected rn values
  cpslp<-c(cpslp, as.data.frame(cp.slp)[1,1])
  c.alldat<-cbind(as.data.frame(Dstar), as.data.frame(e1.est.hat)[,1]) # corrected data frame for correlation analysis
  names(c.alldat)<-c("Dstar", "ctmax.hat")
  c.obj<-phyl.vcv(as.matrix(c.alldat),vcv(tree),1) ## phylogenetic covariance among traits
  cp.cor<-cov2cor(c.obj$R)["Dstar","ctmax.hat"]## phylo correlation between ct & diff corrected
  cpcor<-c(cpcor, cp.cor)
  cr.cor<-cor.test(e1.est.hat, Dstar, method = "pearson")$estimate # raw correlation
  crcor<-c(crcor, as.data.frame(cr.cor)[1,1])
  if (length(rcor) == 1000) {break}
}

sumdat<-as.data.frame(cbind(tpcor, tpslp, trcor, trslp, pcor, pslp, rcor, rslp,  cpcor, cpslp, crcor, crslp)) 
