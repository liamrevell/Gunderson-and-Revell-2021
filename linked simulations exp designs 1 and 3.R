library(phytools)
########### What happens when modelled in the "classical" sense of slope & intercept, Experimental designs 1 and 3 ######
########### Individual phenotypes measured with error ################
tpcor<-NULL #phylogenetic correlation coefficients without error
trcor<-NULL #raw correlation coefficients without error
trslp<-NULL #raw slopes without error
tpslp<-NULL #phylo slopes without error
pcor<-NULL #phylogenetic correlation coefficients with error
rcor<-NULL #raw correlation coefficients with error
rslp<-NULL #raw slopes with error
pslp<-NULL #phylo slopes with error
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
  intercept<-fastBM(tree, sig2 = 0.1) # evolution of intercept along tree
  slope<-fastBM(tree,a=1, sig2 = 0.1) # evolution of reacton norm (slope) across tree 
  
  # analysis with true phenotypes
  tr.slp<-lm(slope~intercept)$coefficients[2] # raw slope
  trslp<-c(trslp, as.data.frame(tr.slp)[1,1])
  tp.slp<-lm(pic(slope,tree)~pic(intercept,tree)+0)$coefficients # pic slope
  tpslp<-c(tpslp, as.data.frame(tp.slp)[1,1])
  talldat<-cbind(as.data.frame(slope), as.data.frame(intercept))
  names(talldat)<-c("slope", "intercept")
  tobj<-phyl.vcv(as.matrix(talldat),vcv(tree),1) ## phylogenetic covariance among traits
  tp.cor<-cov2cor(tobj$R)["intercept","slope"]## phylo correlation between ct & diff
  tpcor<-c(tpcor, tp.cor)
  tr.cor<-cor.test(intercept, slope, method = "pearson")$estimate # raw correlation
  trcor<-c(trcor, as.data.frame(tr.cor)[1,1])
  
  # analysis with phenotypes measured with error
  e1.est<-lapply(0*slope+intercept,sampleFrom,xvar=0.01,n=20) # estimated pheno in env 1 w error 
  e2.est<-lapply(1*slope+intercept,sampleFrom,xvar=0.01,n=20) # estimated pheno in env 2 w error
  ind.slp<-mapply(function(x,y) y-x,x=e1.est,y=e2.est,SIMPLIFY=FALSE) ## Compute individual-level reaction norms (difference between CTmax’ and CTmax).
  e1.est.hat<-sapply(e1.est,mean) # Compute mean basline and induced pheno and the mean reaction norm) for each species to use in interspecific analyses
  e2.est.hat<-sapply(e2.est, mean)
  rn.hat<-sapply(ind.slp,mean)
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
  
  # Adjusting values with error by Kelly and Price Method
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
  
  # same analysis, but sample different individuals for baseline pheno
  n.e1.est<-lapply(0*slope+intercept,sampleFrom,xvar=0.01,n=20) # estimated pheno in env 1 w error, sample new individuals 
  n.e1.est.hat<-sapply(n.e1.est,mean) # mean of new values for each spp
  nr.slp<-lm(rn.hat~n.e1.est.hat)$coefficients[2] # raw slope
  nrslp<-c(nrslp, as.data.frame(nr.slp)[1,1])
  np.slp<-lm(pic(rn.hat,tree)~pic(n.e1.est.hat,tree)+0)$coefficients # pic slope
  npslp<-c(npslp, as.data.frame(np.slp)[1,1])
  nalldat<-cbind(as.data.frame(rn.hat), as.data.frame(n.e1.est.hat)[,1])
  names(nalldat)<-c("slope", "intercept")
  nobj<-phyl.vcv(as.matrix(nalldat),vcv(tree),1) ## phylogenetic covariance among traits w new baseline
  np.cor<-cov2cor(nobj$R)["intercept","slope"]## phylo correlation between ct & diff w new baseline
  npcor<-c(npcor, np.cor)
  nr.cor<-cor.test(n.e1.est.hat, rn.hat, method = "pearson")$estimate # raw correlation w new baseline
  nrcor<-c(nrcor, as.data.frame(nr.cor)[1,1])
  if (length(nrcor) == 1000) {break}
}

sumdat<-as.data.frame(cbind(tpcor, tpslp, trcor, trslp, pcor, pslp, rcor, rslp,  cpcor, cpslp, crcor, crslp, npcor, npslp, nrcor, nrslp))
