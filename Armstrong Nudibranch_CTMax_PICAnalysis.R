# Script for Making Nudibranch Phylogeny for PhyloANOVA
# E. Armstrong 23 August 2016

## Nudibranch Phylogeny and PIC Analysis ##
###########################################################

# Author: ARMSTRONG Eric 
# Created: 23 August 2016
# Last Edited: 4 February 2021

# Modified from a tutorial by Liam Revell found here:
# http://blog.phytools.org/2012/03/contrasts-regression-with.html

# Another useful tutorial:
# http://www.phytools.org/Cordoba2017/ex/3/PICs.html

##### A - Set Defaults #####
############################

  # Local working directory
setwd("~/Dropbox/SICB 2021/Plasticity tradeoff simulations/Armstrong re-analysis")

  # Set the prefix for each output file name
    #outputPrefix <- "Nudibranch_CTMax_PICAnalysis"


##### B - Load AND PREPARE PHYSIOLOGY Data #####
################################################
#library(xlsx)
    
  # Step 1 - Load CTmax data
    data_CTMax <- read.csv("NudibranchCTMaxData.csv")
    
  # Step 2 - Calculate mean CTmax for each experimental group (i.e., each species at each acclimation temperature)
  # NOTE: At the same time, we remove all species for which we don't have full datasets (i.e., species with only one acclimation temp)
    library(Rmisc)
    data_CTMax_summ <- summarySE(subset(data_CTMax, !species %in% c("D. odhneri", "D. picta")), 
                                 measurevar="Ctmax", 
                                 groupvars = c("specacc","genus","species","acclimgroup"))
 
  # Step 3 - Calculate Acclimation Response Ratios (ARRs) for each species
    
    # 3a - Create a column for denoting paired acclimation groups within the same species
      data_CTMax_summ$pairs <- c(0, 1:(nrow(data_CTMax_summ)-1) %/% 2)
      
    # 3b - Calculate the ARR for each acclimation group pair
      # NOTE: ARR = (CTMax @ acclimtemp2 - CTMax @ acclimtemp1) / (acclimtemp2 - acclimtemp1)
      # We have two acclimation temperatures: 13?C and 17?C, thus the values in the FUN below
      ARRdata <- aggregate(list(ARR = data_CTMax_summ$Ctmax),
                           by = list(pairs = data_CTMax_summ$pairs, species = data_CTMax_summ$species),
                           FUN=function(y) diff(y)/(17 - 13))
    
    # 3c - Calculate error (95%-CI) associated with each ARR
    # NOTE: New CIs calculated using error propagation rules for addition/subtraction
      ARRCIdata <- aggregate(list(ARRCI = data_CTMax_summ$ci),
                           by = list(pairs = data_CTMax_summ$pairs, species = data_CTMax_summ$species),
                           FUN=function(y) sqrt(sum(y^2))/(17 - 13))
      
    # 3d - Merge all data relevant for downstream PIC into new data frame
      library(expss)
      data_ARR_CT <- data.frame("Species" = ARRdata$species,
                                "CTMax13" = vlookup(lookup_value = paste(ARRdata$species, 13, sep = ' '),
                                                    dict = data_CTMax_summ,
                                                    result_column = "Ctmax",
                                                    lookup_column = "specacc"),
                                "CTMax13ci" = vlookup(lookup_value = paste(ARRdata$species, 13, sep = ' '),
                                                    dict = data_CTMax_summ,
                                                    result_column = "ci",
                                                    lookup_column = "specacc"),
                                "CTMax17" = vlookup(lookup_value = paste(ARRdata$species, 17, sep = ' '),
                                                    dict = data_CTMax_summ,
                                                    result_column = "Ctmax",
                                                    lookup_column = "specacc"),
                                "CTMax17ci" = vlookup(lookup_value = paste(ARRdata$species, 17, sep = ' '),
                                                    dict = data_CTMax_summ,
                                                    result_column = "ci",
                                                    lookup_column = "specacc"),
                                "ARR" = ARRdata$ARR,
                                "ARRci" = ARRCIdata$ARRCI)
  
    
##### LOAD AND PREPARE PHYLOGENY DATA #####
###########################################
# Potentially useful information here:
# https://www.biostars.org/p/210401/
      
library(ape)
library(phytools)
      
  # Step 1 - Read in the phylogeny (created in Geneious based on MSA of Genbank COI sequences)
    tree_Nudi <- read.nexus("GeneiousTree_Bootstrap_Genbank_Nudi_COI.nex")
    
    # 1b - Delete the single nodes (i.e., with a single descendant) in the tree.
      # tree_Nudi <- collapse.singles(tree5)
    
    # Show tip labels (i.e., Genbank COI IDs)
    # tree_Nudi$tip.label
    
  # Step 2 - Remove the outgroup (Aplyisa spp.) for which we have no physiological data
    tree_Nudi <- drop.tip(tree_Nudi, tip=c("AF077759"))
    
    # 2a - Rename the tips to match species
    # NOTE: Doing this manually bc we have so few speceies, but could use vlookup for larger datasets
      tree_Nudi$tip.label = c("Hermissenda crassocornis", "Triopha catalinae", "Triopha maculata", "Okenia rosacea", "Doriopsilla albopunctata", "Doriopsiall fulva" )
    
    # 2b - Visualize the tree
      plot(tree_Nudi, edge.width = 2)
      

##### PERFORM THE PIC ANALYSIS #####
####################################

  # Step 1 - Resolve multichotomies randomly (if present) in phylogenetic tree
    rt1 <- multi2di(tree_Nudi, random=TRUE)
      
  # Step 2 - Specify X- and Y-values in PIC plot
  # NOTE: We are contrasting CTMax @ the "ambient" acclimation temp of 13?C against the CTMax ARR (Y).
  # We also remove the data for the species Phidiana hiltoni because there was no corresonding Genbank COI sequence data.
    
    # Explanatory variable - CTMax13 
    x <- subset(data_ARR_CT, !Species %in% c("P. hiltoni"))$CTMax13
    
    # Dependent variable - CTMax ARR
    y <- subset(data_ARR_CT, !Species %in% c("P. hiltoni"))$ARR
      
  # Step 3 - Compute the contrasts using the resolved phylogeny
    pic.x <- pic(x, rt1)
    pic.y <- pic(y, rt1)
      
  # Step 4 - Fit a bivariate regression model to the contrasts
    Phyloline <- lm(pic.y ~ pic.x - 1) # fits the regression model without an intercept
    
    # 4a - View the results to see if there is still a significant correlation between X and Y once phylogeny has been accounted for.
      summary(Phyloline)
    
    # 4b - Quick visualization of the PICs
      plot(pic.x, pic.y)
      abline(Phyloline)
      
#########################################      
#### Now, applying Kelly and Price method
      n<-Ntip(rt1) # number of tips on the tree
      X<-cbind(data_ARR_CT$CTMax13[data_ARR_CT$Species != "P. hiltoni"],
               data_ARR_CT$CTMax17[data_ARR_CT$Species != "P. hiltoni"]) # mean tip values for each trait for each species
      colnames(X)<-c("ctmax1","ctmax2")
      object<-phyl.vcv(X,vcv(rt1),lambda=1) # phylogenetic variance/covariance matrix
      s1<-sqrt(object$R[1,1]) # square root of phylogenetic variance in ctmax1
      s2<-sqrt(object$R[2,2]) # square root of phylogenetic variance in ctmax2
      r<-cov2cor(object$R)[1,2] # phylogenetic correlation coefficeint bewteen ctmax1 and ctmax2
      
      T<-sqrt(n-2)*(s1/s2-s2/s1)/(2*(1-r^2)) # T statistic for  Pittman's  test of homogeneity of variances
      
      Pval<-2*pt(abs(T),df=n-2,lower.tail=FALSE) # P value for T statistic
      
      rho.hat<-if(Pval>=0.05) r else 2*r*s1*s2/(s1^2+s2^2) # rho value for phenotype correction based on whether or not variances are equal 
      
      Xhat<-object$alpha[,1] # phylopgenetic means for ctmax1 and ctmax2
      
      Dstar<-rho.hat*(X[,1]-Xhat[1])-(X[,2]-Xhat[2]) # corrected eraction norm values
      
      # PIC with adjusted values
      cpic.x <- pic(x, rt1)
      cpic.y <- pic(Dstar, rt1)
      cpic.y <- cpic.y*-1 # flip sign to test trade-off orientation
      
      # Step 4 - Fit a bivariate regression model to the contrasts
      cPhyloline <- lm(cpic.y ~ cpic.x - 1) # fits the regression model without an intercept
      
      # 4a - View the results to see if there is still a significant correlation between X and Y once phylogeny has been accounted for.
      summary(cPhyloline)
      
      # 4b - Quick visualization of the PICs
      ## Manuscript figure, saved as 4 x 4
      plot(cpic.x, cpic.y, pch = 21, col = "black", bg = "black", 
           xlab = "Heat tolerance contrasts", ylab = "Plasticity contrast")
      abline(cPhyloline, col = "black", lwd = 1.5)
      
      points(pic.x, pic.y, pch = 21, col = "grey", bg = "grey")
      abline(Phyloline, lwd = 1.5, col = "grey")
      
      # phylogenetic least squarers of corrected plasticity values vs ctmax1
      #spp<-names(Dstar)
      #corBM<-corBrownian(1,rt1,form=~spp)
      #fit<-nlme::gls(Dstar~X[,1],        
      #               data=data.frame(Dstar=Dstar,ctmax.hat=X[,1]),
      #               correlation=corBM)
    #  
    #  summary(fit)
    #  plot(ctmax.hat,Dstar)
    #  abline(fit)
    #  
    #  beta1[i]<-coef(fit)[2]
    #  
    #  }

mean(beta1)
hist(beta1)
hist(Pval)

      
##### VISUALIZATION #####
#########################
library(ggplot2)
library(gridExtra)
library(ggpubr)

##### Figure 1 - PIC Analysis #####
  # Step 1 - Create ggplot-friendly data frame and specific the formula to use for "best-fit" line
    data_PIC <- data.frame("CTMax13contr" = pic.x, 
                           "ARRcontr" = pic.y)
    my.formula <- y ~ x

  # Step 2 - Plot PIC Results
    plotPIC = ggplot(data=data_PIC, aes(x=CTMax13contr, y=ARRcontr)) +
      geom_point(size=4) +
      geom_smooth(method = "lm", se=FALSE, 
                  formula = my.formula, color = "Grey80") +
      stat_cor(label.x = 3, label.y = 1.8,
               aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
      stat_regline_equation(label.x = 3, label.y = 2) +
      scale_x_continuous(limits=c(-15,10.5), breaks=seq(from=-15, to=10, by=5), expand=c(0,0)) +
      scale_y_continuous(limits=c(-2.2,2.5), breaks=seq(from=-2, to=2, by=1), expand=c(0.05,0)) +
      ylab(expression(paste(" ", CT[max], "  ARR Contrasts ", sep=""))) +
      xlab(expression(paste(" ", CT[max], " Contrasts", sep=""))) +
      theme_bw() +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            text = element_text(size=16, face="bold"))
    plotPIC


    pdf(file = paste0(outputPrefix, "_ContrastFigure_CTMaxARRvsCTMax13.pdf"), width = 12, height = 8.5)
    plotPIC
    dev.off()

    png(file = paste0(outputPrefix, "_ContrastFigure_CTMaxARRvsCTMax13.png"), width = 12, height = 8.5, units = 'in', res = 300)
    plotPIC
    dev.off()
    
    
   
