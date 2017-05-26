###Simulation of direct and competition genotypic effect and analysis of simulated data##
##Ani A. Elias ## June, 2016##

library(nlme) # to fit the first base model
library(LMERConvenienceFunctions) # remove outliers 
library(regress) #to fit models 
library(reshape2) #to melt matrix
library(rrBLUP) #to get K matrix
library(MASS) # for ginv
library(data.table) # to convert row names to a column


#madII is a matrix from real data with at least two replicates for each genotypes
  #- contain genotype name and field coordinates
#use simulated data on using madII 
#wrapper function to retrieve effects from simulated data
#remove outliers (2.5 times the sd) before proceeding to analysis
#this functon assumes that field coordinates are used for spatial analysis. It is required to change the 
  #method of distance calculation etc. if geocoordinates are used. 
       
#assuming unit variance for the direct effect, a variance for competiton (compG.var), and some correlation(G.cor)
#E.cor is the correlation between competition and residual errors
#dirh2 is the heritability of direct genetic effect
#fra.CompE is the fraction of total error which accounts for competiton error
#plotSize is the dimension of plots - plotSize[1] for the longer side and plotSize[2] for shorter

#test by ignoring the genetic and error correlation indicated by 0 for G.cor and E.cor
#G.cor <- {0,0.4}
#E.cor <- {0,0.8}
#dirh2 <- {0.7,0.3}
#because compG.var < (1-dirh2)/dirh2,
#when dirh2 = 0.7; compG.var < {0.4,0.2,0.1} 
  #when dirh2 = 0.3; compG.var <- {1,0.5,0.2,0.1}
#fra.CompE <- {0.9,0.6,0.3}
#G.cov is the covariance between direct and competition effect
#repl indicates number of minimum replication per genotype in the dataset

#NOTE: base model is to check the impact of avoiding the competition effect present in the data 
      # does this decrease genetic variance and heritability? because competition can decrease genetic resemblence among sibs

geno_compet_simulate_snp_SD <- function(madII, G.cor, E.cor, dirh2, compG.var, fra.CompE, plotSize, repl){

entry.cor <- matrix(NA, 20, 13)
colnames(entry.cor) <- c("Base", "Direct", "Direct_2", "Compet", "Compet_2","CS_dir","CS_comp","CS2_dir","CS2_comp",
                      "CF_dir","CF_comp","CF2_dir","CF2_comp")

entry.mse <- matrix(NA, 20, 13)
colnames(entry.mse) <- c("Base", "Direct", "Direct_2", "Compet", "Compet_2","CS_dir","CS_comp","CS2_dir","CS2_comp",
                      "CF_dir","CF_comp","CF2_dir","CF2_comp")

error.cor <- matrix(NA,20,3)
colnames(error.cor) <- c("Model_2","CS2","CF2")

error.mse <- matrix(NA,20,3)
colnames(error.mse) <- c("Model_2","CS2","CF2")

herit <- matrix(NA,20,15)
colnames(herit) <- c("dirh2", "comph2", "baseh2", "Direct_h2", "Direct2_h2", "Comp_h2","Comp2_h2","CS_dirh","CS_comph",
                      "CS2_dir","CS2_comp","CF_dirh","CF_comph","CF2_dirh","CF2_comph")

fra.Er <- matrix(NA,20,8)
colnames(fra.Er) <- c("fraSp", "Model_fraSp","E.cor","modelC2_Ecor","CS2_Evar","CS2_Ecor","CF2_Evar","CF2_Ecor")

varGs <- matrix(NA,20,22)
colnames(varGs) <- c("dirVar","compVar","G.cor","base_dirVar","modelC_dirVar","modelC_compVar", "modelC_genoCor",
      "modelC2_dirVar","modelC2_compVar","modelC2_genoCor","CS_dirVar","CS_compVar","CS_genoCor","CS2_dirVar",
      "CS2_compVar","CS2_genoCor","CF_dirVar","CF_compVar","CF_genoCor","CF2_dirVar","CF2_compVar","CF2_genoCor")
  

base.fail <- modelC.fail <- modelC2.fail <- modelCS.fail <- modelCS2.fail <- modelCF.fail <- modelCF2.fail <- vector()

for(s in c(1:20)){
#create data1
  Comp.list <- simulateComp_snp_SD(madII, G.cor, E.cor, dirh2, compG.var, fra.CompE, plotSize)
  data1 <- Comp.list$fieldObs
 
  fixed <- formula("phenoVal ~ 1")
  base <- lme(fixed = fixed, random = ~ 1|Entry, data = na.omit(data1), method = "ML")
  data1a <- (romr.fnc(base, na.omit(data1), trim = 2.5))$data #outlier removal

  #K matrix aligning with data1a
  #subsetting K matrix  - helps to remove clones that are not genotyped but phenotyped  
  table(data1a[,"Entry"] %in% rownames(snps)) 
  snp2.names <- intersect(data1a[,"Entry"], rownames(snps))
  snps2 <- subset(snps, rownames(snps) %in% snp2.names)

  #calculate K matrix
  K <- A.mat(snps2-1)
  diag(K)<-diag(K)+1e-6

  data2 <- subset(data1a,data1a[,"Entry"] %in% snp2.names)
    write.csv(data2, "data2.csv")
    data2 <- read.csv("data2.csv", h=T, sep=",") 

  entryEff <- Comp.list$entryEffects

  entryComp <- unique(data2[,c("Entry","competEff")])
  compError <- data2[,"competErr"]

#incidence matrix 
  entries <- sort(unique(data2$Entry))

  id <- factor(as.character(data2[,"Entry"]), levels = rownames(K))
  Z.base <- model.matrix(~id - 1)

  #genotypic design matrix
  idg <- factor(as.character(data2[,'Entry']), levels = rownames(K))
   Z.geno <- model.matrix(~idg -1) 
   geno <- Z.geno%*%K%*%t(Z.geno) 
  
#calculate the distance matrices
  #modified field coordinates and distance matrices
  rowVec<- data2$Range
  colVec <- data2$Column

  rowVec2 <- rowVec * plotSize[1] 
  rowDistMat2 <- sapply(rowVec2, function(rowNum) abs(rowNum - rowVec2))

  colVec2 <- colVec * plotSize[2]
  colDistMat2 <- sapply(colVec2, function(colNum) abs(colNum - colVec2))

  distMat <- sqrt(rowDistMat2^2 + colDistMat2^2)
  rownames(distMat) <- colnames(distMat) <- 1:nrow(distMat)

#incidence matrix for competition 
  # NN on longer edge only
  distMat.4 <- distMat
  distMat.4[distMat.4 > plotSize[1]] <- 0
  distMat.4[distMat.4 == plotSize[1]] <- 1
     Z.comp.4 <- distMat.4 %*% Z.base

  G.comp.4 <- Z.comp.4 %*% K %*% t(Z.comp.4)

  #slow decay ending at plotSize[1] *3
  distMat.2 <- distMat
  distMat.inv.2 <- 1/distMat.2
  distMat.inv.2[distMat.inv.2 == "Inf"] <- 0
  distMat.1 <- 1/(0.999 + 0.001^distMat.inv.2)
  a <- 1/(3*plotSize[1])
  distMat.1[distMat.1 < (1/(0.999 + 0.001^a))] <- 0
  distMat.1[is.na(distMat.1)] <- 0
  diag(distMat.1) <- 0
      Z.comp.1 <- distMat.1 %*% Z.base

   G.comp.1 <- Z.comp.1 %*% K %*% t(Z.comp.1)

  #fast decay with k=0.4
  distMat.inv <- 1/distMat
  diag(distMat.inv) <- 0 

  distMat.3 <- 0.4*(distMat.inv)/(0.4-(distMat.inv) + 1)
  distMat.3[distMat.3 == "NaN"] <- 0

      Z.comp.3 <- distMat.3 %*% Z.base
      G.comp.3 <- Z.comp.3 %*% K %*% t(Z.comp.3) 

  #competition error
  Identity <- diag(nrow(data2))
   #distMat.4 %*% Identiti matrix is distMat.4; so ignoring that line   
  G.compE <- distMat.4 %*% Identity %*% t(distMat.4) #can ignore multiplying with Identity but just ZZ'

  
  
#models  
  base <- try(regress(phenoVal ~ 1, ~geno , pos=rep(TRUE,2), tol = 1e-4, data=data2), silent = TRUE)
  if (class(base) != "try-error"){
    #blup
    base.blup <- BLUP(base)$Mean

    #entry blup
    base.EB <- base.blup[grep("geno", names(base.blup), fixed=TRUE)]
    data2a <- cbind(data2,base.EB)

    #heritability
    base.dirh2 <- base$sigma[[1]]/(base$sigma[[1]] + base$sigma[[2]])

    #correlation and mse
    entry.cor[s,1] <- cor(data2a$directG, data2a$base.EB) 
    entry.mse[s,1] <- mean(abs(data2a$directG - data2a$base.EB))
    herit[s,3] <- base.dirh2
    varGs[s,4] <- base$sigma[[1]]
  } else{
    base.fail[s] <- 1
  }  
  modelC <- try(regress(phenoVal ~ 1, ~geno + G.comp.4, pos=rep(TRUE,3), tol = 1e-4, data=data2), silent = TRUE)
  if(class(modelC) != "try-error"){
    #blup
    modelC.blup <- BLUP(modelC)$Mean

    #competition blup
    modelC.comp.a <- modelC.blup[grep("G.comp.4", names(modelC.blup), fixed=TRUE)] 
    modelC.comp <- ginv(Z.comp.4) %*% modelC.comp.a
    comp.var <- var(modelC.comp)
    rownames(modelC.comp) <- colnames(Z.comp.4)
    modelC.comp <- as.data.frame(modelC.comp)
    modelC.comp <- setDT(modelC.comp, keep.rownames = TRUE)[]
    colnames(modelC.comp) <- c("Entry","modelC.comp")
    modelC.comp$Entry <- gsub("id","",modelC.comp$Entry)
    data2.comp <- merge(entryComp,modelC.comp)

    #entry blup
    modelC.EB <- modelC.blup[grep("geno", names(modelC.blup), fixed=TRUE)]
    data2a <- cbind(data2a,modelC.EB)

    #heritability
    modelC.comph2 <- comp.var/(modelC$sigma[[1]] + comp.var + modelC$sigma[[3]])
    modelC.dirh2 <- modelC$sigma[[1]]/(modelC$sigma[[1]] + comp.var + modelC$sigma[[3]])

    #correlation and mse
    entry.cor[s,2] <- cor(data2a$directG, data2a$modelC.EB)
    entry.cor[s,4] <- cor(data2.comp$competEff, data2.comp$modelC.comp)
    entry.mse[s,2] <- mean(abs(data2a$directG - data2a$modelC.EB))
    entry.mse[s,4] <- mean(abs(data2.comp$competEff - data2.comp$modelC.comp))
    herit[s,4] <- modelC.dirh2
    herit[s,6] <- modelC.comph2
    varGs[s,5] <- modelC$sigma[[1]]
    varGs[s,6] <- comp.var
   #varGs[s,7] <- cor(modelC.EB, modelC.comp)    
  } else {
    modelC.fail[s] <- 1
  }
  modelC.2 <- try(regress(phenoVal ~ 1, ~geno + G.comp.4 + G.compE, pos=rep(TRUE,4), tol = 1e-4, data=data2), silent = TRUE)
  if(class(modelC.2) != "try-error") {
    #blup
    modelC2.blup <- BLUP(modelC.2)$Mean  

    #competition blup
    modelC2.comp.a <- modelC2.blup[grep("G.comp.4", names(modelC2.blup), fixed=TRUE)]
    modelC2.comp <- ginv(Z.comp.4) %*% modelC2.comp.a
    comp2.var <- var(modelC2.comp)
    rownames(modelC2.comp) <- colnames(Z.comp.4)
    modelC2.comp <- as.data.frame(modelC2.comp)
    modelC2.comp <- setDT(modelC2.comp, keep.rownames = TRUE)[]
    colnames(modelC2.comp) <- c("Entry","modelC2.comp")
    modelC2.comp$Entry <- gsub("id","",modelC2.comp$Entry)
    data2.comp <- merge(data2.comp,modelC2.comp)

    #entry blup
    modelC2.EB <- modelC2.blup[grep("geno", names(modelC2.blup), fixed=TRUE)]
    data2a <- cbind(data2a,modelC2.EB)
    
    #competition error
    modelC2.compE.a <- modelC2.blup[grep("G.compE", names(modelC2.blup), fixed=TRUE)] 
    modelC2.compE <- ginv(distMat.4) %*% modelC2.compE.a
    data2a <- cbind(data2a, modelC2.compE)

    #heritability
    modelC2.comph2 <- comp2.var/(modelC.2$sigma[[1]] + comp2.var + modelC.2$sigma[[3]] + modelC.2$sigma[[4]])
    modelC2.dirh2 <- modelC.2$sigma[[1]]/(modelC.2$sigma[[1]] + comp2.var + modelC.2$sigma[[3]] + modelC.2$sigma[[4]])

    #correlation and mse
    entry.cor[s,3] <- cor(data2$directG, modelC2.EB)
    entry.cor[s,5] <- cor(data2.comp$competEff, data2.comp$modelC2.comp)
    entry.mse[s,3] <- mean(abs(data2$directG - modelC2.EB))
    entry.mse[s,5] <- mean(abs(data2.comp$competEff - data2.comp$modelC2.comp)) 
    error.cor[s,1] <- cor(compError, modelC2.compE)
    error.mse[s,1] <- mean(abs(compError - modelC2.compE))
    herit[s,5] <- modelC2.dirh2
    herit[s,7] <- modelC2.comph2
    fra.Er[s,2] <- modelC.2$sigma[[3]]/modelC.2$sigma[[3]] + modelC.2$sigma[[4]]
    varGs[s,8] <- modelC.2$sigma[[1]]
    varGs[s,9] <- comp2.var
    #varGs[s,10] <- cor(modelC2.EB, modelC2.comp)
    fra.Er[s,4] <- cor((data2[,"phenoVal"] - modelC.2$predicted),modelC2.compE)
  } else{
    modelC2.fail[s] <- 1
  }

 modelCS <- try(regress(phenoVal ~ 1, ~geno + G.comp.1, pos=rep(TRUE,3), tol = 1e-4, data=data2), silent = TRUE)
  if(class(modelCS) != "try-error"){
    #blup
    modelCS.blup <- BLUP(modelCS)$Mean

    #competition blup
    modelCS.comp.a <- modelCS.blup[grep("G.comp.1", names(modelCS.blup), fixed=TRUE)] 
    modelCS.comp <- ginv(Z.comp.1) %*% modelCS.comp.a
    compCS.var <- var(modelCS.comp)
    rownames(modelCS.comp) <- colnames(Z.comp.1)
    modelCS.comp <- as.data.frame(modelCS.comp)
    modelCS.comp <- setDT(modelCS.comp, keep.rownames = TRUE)[]
    colnames(modelCS.comp) <- c("Entry","modelCS.comp")
    modelCS.comp$Entry <- gsub("id","",modelCS.comp$Entry)
    data2.comp <- merge(data2.comp,modelCS.comp)

    #entry blup
    modelCS.EB <- modelCS.blup[grep("geno", names(modelCS.blup), fixed=TRUE)]
    data2a <- cbind(data2a,modelCS.EB)

    #heritability
    modelCS.comph2 <- compCS.var/(modelCS$sigma[[1]] + compCS.var + modelCS$sigma[[3]])
    modelCS.dirh2 <- modelCS$sigma[[1]]/(modelCS$sigma[[1]] + compCS.var + modelCS$sigma[[3]])

    #correlation and mse
    entry.cor[s,6] <- cor(data2a$directG, data2a$modelCS.EB)
    entry.cor[s,7] <- cor(data2.comp$competEff, data2.comp$modelCS.comp)
    entry.mse[s,6] <- mean(abs(data2a$directG - data2a$modelCS.EB))
    entry.mse[s,7] <- mean(abs(data2.comp$competEff - data2.comp$modelCS.comp))
    herit[s,8] <- modelCS.dirh2
    herit[s,9] <- modelCS.comph2
    varGs[s,11] <- modelCS$sigma[[1]]
    varGs[s,12] <- compCS.var
    #varGs[s,13] <- cor(modelCS.EB, modelCS.comp)    
  } else {
    modelCS.fail[s] <- 1
  }
  modelCS.2 <- try(regress(phenoVal ~ 1, ~geno + G.comp.1 + G.compE, pos=rep(TRUE,4), tol = 1e-4, data=data2), silent = TRUE)
  if(class(modelCS.2) != "try-error") {
    #blup
    modelCS2.blup <- BLUP(modelCS.2)$Mean  

    #competition blup
    modelCS2.comp.a <- modelCS2.blup[grep("G.comp.1", names(modelCS2.blup), fixed=TRUE)]
    modelCS2.comp <- ginv(Z.comp.1) %*% modelCS2.comp.a
    compCS2.var <- var(modelCS2.comp)
    rownames(modelCS2.comp) <- colnames(Z.comp.1)
    modelCS2.comp <- as.data.frame(modelCS2.comp)
    modelCS2.comp <- setDT(modelCS2.comp, keep.rownames = TRUE)[]
    colnames(modelCS2.comp) <- c("Entry","modelCS2.comp")
    modelCS2.comp$Entry <- gsub("id","",modelCS2.comp$Entry)
    data2.comp <- merge(data2.comp,modelCS2.comp)

    #entry blup
    modelCS2.EB <- modelCS2.blup[grep("geno", names(modelCS2.blup), fixed=TRUE)]
    data2a <- cbind(data2a,modelCS2.EB)
    
    #competition error
    modelCS2.compE.a <- modelCS2.blup[grep("G.compE", names(modelCS2.blup), fixed=TRUE)] 
    modelCS2.compE <- ginv(distMat.4) %*% modelCS2.compE.a
    data2a <- cbind(data2a,modelCS2.compE)

    #heritability
    modelCS2.comph2 <- compCS2.var/(modelCS.2$sigma[[1]] + compCS2.var + modelCS.2$sigma[[3]] + modelCS.2$sigma[[4]])
    modelCS2.dirh2 <- modelCS.2$sigma[[1]]/(modelCS.2$sigma[[1]] + compCS2.var + modelCS.2$sigma[[3]] + modelCS.2$sigma[[4]])

    #correlation and mse
    entry.cor[s,8] <- cor(data2$directG, modelCS2.EB)
    entry.cor[s,9] <- cor(data2.comp$competEff, data2.comp$modelCS2.comp)
    entry.mse[s,8] <- mean(abs(data2$directG - modelCS2.EB))
    entry.mse[s,9] <- mean(abs(data2.comp$competEff - data2.comp$modelCS2.comp)) 
    error.cor[s,2] <- cor(compError, modelCS2.compE)
    error.mse[s,2] <- mean(abs(compError - modelCS2.compE))
    herit[s,10] <- modelCS2.dirh2
    herit[s,11] <- modelCS2.comph2
    fra.Er[s,5] <- modelCS.2$sigma[[3]]/modelCS.2$sigma[[3]] + modelCS.2$sigma[[4]]
    varGs[s,14] <- modelCS.2$sigma[[1]]
    varGs[s,15] <- compCS2.var
    #varGs[s,16] <- cor(modelCS2.EB, modelCS2.comp)
    fra.Er[s,6] <- cor((data2[,"phenoVal"] - modelCS.2$predicted),modelCS2.compE)
  } else{
    modelCS2.fail[s] <- 1
  }

  modelCF <- try(regress(phenoVal ~ 1, ~geno + G.comp.3, pos=rep(TRUE,3), tol = 1e-4, data=data2), silent = TRUE)
  if(class(modelCF) != "try-error"){
    #blup
    modelCF.blup <- BLUP(modelCF)$Mean

    #competition blup
    modelCF.comp.a <- modelCF.blup[grep("G.comp.3", names(modelCF.blup), fixed=TRUE)] 
    modelCF.comp <- ginv(Z.comp.3) %*% modelCF.comp.a
    compCF.var <- var(modelCF.comp)
    rownames(modelCF.comp) <- colnames(Z.comp.3)
    modelCF.comp <- as.data.frame(modelCF.comp)
    modelCF.comp <- setDT(modelCF.comp, keep.rownames = TRUE)[]
    colnames(modelCF.comp) <- c("Entry","modelCF.comp")
    modelCF.comp$Entry <- gsub("id","",modelCF.comp$Entry)
    data2.comp <- merge(data2.comp,modelCF.comp)

    #entry blup
    modelCF.EB <- modelCF.blup[grep("geno", names(modelCF.blup), fixed=TRUE)]
    data2a <- cbind(data2a,modelCF.EB)

    #heritability
    modelCF.comph2 <- compCF.var/(modelCF$sigma[[1]] + compCF.var + modelCF$sigma[[3]])
    modelCF.dirh2 <- modelCF$sigma[[1]]/(modelCF$sigma[[1]] + compCF.var + modelCF$sigma[[3]])

    #correlation and mse
    entry.cor[s,10] <- cor(data2a$directG, data2a$modelCF.EB)
    entry.cor[s,11] <- cor(data2.comp$competEff, data2.comp$modelCF.comp)
    entry.mse[s,10] <- mean(abs(data2a$directG - data2a$modelCF.EB))
    entry.mse[s,11] <- mean(abs(data2.comp$competEff - data2.comp$modelCF.comp))
    herit[s,12] <- modelCF.dirh2
    herit[s,13] <- modelCF.comph2
    varGs[s,17] <- modelCF$sigma[[1]]
    varGs[s,18] <- compCF.var
    #varGs[s,19] <- cor(modelCF.EB, modelCF.comp)    
  } else {
    modelCF.fail[s] <- 1
  }
  modelCF.2 <- try(regress(phenoVal ~ 1, ~geno + G.comp.3 + G.compE, pos=rep(TRUE,4), tol = 1e-4, data=data2), silent = TRUE)
  if(class(modelCF.2) != "try-error") {
    #blup
    modelCF2.blup <- BLUP(modelCF.2)$Mean  

    #competition blup
    modelCF2.comp.a <- modelCF2.blup[grep("G.comp.3", names(modelCF2.blup), fixed=TRUE)]
    modelCF2.comp <- ginv(Z.comp.3) %*% modelCF2.comp.a
    compCF2.var <- var(modelCF2.comp)
    rownames(modelCF2.comp) <- colnames(Z.comp.3)
    modelCF2.comp <- as.data.frame(modelCF2.comp)
    modelCF2.comp <- setDT(modelCF2.comp, keep.rownames = TRUE)[]
    colnames(modelCF2.comp) <- c("Entry","modelCF2.comp")
    modelCF2.comp$Entry <- gsub("id","",modelCF2.comp$Entry)
    data2.comp <- merge(data2.comp,modelCF2.comp)

    #entry blup
    modelCF2.EB <- modelCF2.blup[grep("geno", names(modelCF2.blup), fixed=TRUE)]
    data2a <- cbind(data2a,modelCF2.EB)
    
    #competition error
    modelCF2.compE.a <- modelCF2.blup[grep("G.compE", names(modelCF2.blup), fixed=TRUE)] 
    modelCF2.compE <- ginv(distMat.4) %*% modelCF2.compE.a
    data2a <- cbind(data2a,modelCF2.compE)

    #heritability
    modelCF2.comph2 <- compCF2.var/(modelCF.2$sigma[[1]] + compCF2.var + modelCF.2$sigma[[3]] + modelCF.2$sigma[[4]])
    modelCF2.dirh2 <- modelCF.2$sigma[[1]]/(modelCF.2$sigma[[1]] + compCF2.var + modelCF.2$sigma[[3]] + modelCF.2$sigma[[4]])

    #correlation and mse
    entry.cor[s,12] <- cor(data2$directG, modelCF2.EB)
    entry.cor[s,13] <- cor(data2.comp$competEff, data2.comp$modelCF2.comp)
    entry.mse[s,12] <- mean(abs(data2$directG - modelCF2.EB))
    entry.mse[s,13] <- mean(abs(data2.comp$competEff - data2.comp$modelCF2.comp)) 
    error.cor[s,3] <- cor(compError, modelCF2.compE)
    error.mse[s,3] <- mean(abs(compError - modelCF2.compE))
    herit[s,14] <- modelCF2.dirh2
    herit[s,15] <- modelCF2.comph2
    fra.Er[s,7] <- modelCF.2$sigma[[3]]/modelCF.2$sigma[[3]] + modelCF.2$sigma[[4]]
    varGs[s,20] <- modelCF.2$sigma[[1]]
    varGs[s,21] <- compCF2.var
    #varGs[s,22] <- cor(modelCF2.EB, modelCF2.comp)
    fra.Er[s,8] <- cor((data2[,"phenoVal"] - modelCF.2$predicted),modelCF2.compE)
  } else{
    modelCF2.fail[s] <- 1
  }
    

#heritability
  comph2 <- Comp.list$comph2
  
  herit[s,1] <- dirh2
  herit[s,2] <- comph2  

  #fraction of competition error variance
  fra.Er[s,1] <- fra.CompE  
  fra.Er[s,3] <- E.cor

  #variance components
  varGs[s,1] <- Comp.list$dirG.var.2
  varGs[s,2] <- Comp.list$compG.var.2
  varGs[s,3] <- Comp.list$G.cor.2
    
} # end of simulation

#for visualization of cor and mse value 
 entry.cor.melt <- melt(entry.cor)
 entry.mse.melt <- melt(entry.mse)
 error.cor.melt <- melt(error.cor)
 error.mse.melt <- melt(error.mse)
 herit.melt <- melt(herit)
 fraEr.melt <- melt(fra.Er)
 varGs.melt <- melt(varGs)

#output the files simulation for one phi
 write.csv(entry.cor.melt, paste0("CompSD","-G.cor", G.cor,"-E.cor",E.cor,"-dirh2", dirh2, "-compGvar", compG.var, "-fraE",fra.CompE,"-rep_",repl,"_entry_Effect_cor_melt.csv"))
 write.csv(entry.mse.melt, paste0("CompSD","-G.cor", G.cor,"-E.cor",E.cor,"-dirh2", dirh2, "-compGvar", compG.var, "-fraE",fra.CompE,"-rep_",repl,"_entry_Effect_rmse_melt.csv"))
 write.csv(error.cor.melt, paste0("CompSD","-G.cor", G.cor,"-E.cor",E.cor,"-dirh2", dirh2, "-compGvar", compG.var, "-fraE",fra.CompE,"-rep_",repl,"_error_Effect_cor_melt.csv"))
 write.csv(error.mse.melt, paste0("CompSD","-G.cor", G.cor,"-E.cor",E.cor,"-dirh2", dirh2, "-compGvar", compG.var, "-fraE",fra.CompE,"-rep_",repl,"_error_Effect_rmse_melt.csv"))
 write.csv(herit.melt, paste0("CompSD","-G.cor", G.cor,"-E.cor",E.cor,"-dirh2", dirh2, "-compGvar", compG.var, "-fraE",fra.CompE,"-rep_",repl,"_herit_melt.csv"))
 write.csv(fraEr.melt, paste0("CompSD","-G.cor", G.cor,"-E.cor",E.cor,"-dirh2", dirh2, "-compGvar", compG.var, "-fraE",fra.CompE,"-rep_",repl,"_fraSp_melt.csv"))
 write.csv(varGs.melt, paste0("CompSD","-G.cor", G.cor,"-E.cor",E.cor,"-dirh2", dirh2, "-compGvar", compG.var, "-fraE",fra.CompE,"-rep_",repl,"_varGs_melt.csv"))

 #mean of correlation and mse   
 entry.cor.2 <- apply(entry.cor,2,mean, na.rm = TRUE) 
 entry.mse.2 <- apply(entry.mse,2,mean, na.rm = TRUE) 
 error.cor.2 <- apply(error.cor,2,mean, na.rm = TRUE)
 error.mse.2 <- apply(error.mse,2,mean, na.rm = TRUE)
 herit.2 <- apply(herit,2,mean, na.rm = TRUE)
 fra.Er.2 <- apply(fra.Er,2,mean, na.rm = TRUE) 
 fail.sum <- apply(cbind(base.fail, modelC.fail, modelC2.fail, modelCS.fail, modelCS2.fail, modelCF.fail, modelCF2.fail),2,sum, na.rm = TRUE)
 varGs.2 <- apply(varGs,2,mean, na.rm = TRUE)

 write.csv(entry.cor.2, paste0("CompSD","-G.cor", G.cor,"-E.cor",E.cor,"-dirh2", dirh2, "-compGvar", compG.var, "-fraE",fra.CompE,"-rep_",repl,"_entry_Effect_cor_final.csv"))
 write.csv(entry.mse.2, paste0("CompSD","-G.cor", G.cor,"-E.cor",E.cor,"-dirh2", dirh2, "-compGvar", compG.var, "-fraE",fra.CompE,"-rep_",repl,"_entry_Effect_mse_final.csv"))
 write.csv(error.cor.2, paste0("CompSD","-G.cor", G.cor,"-E.cor",E.cor,"-dirh2", dirh2, "-compGvar", compG.var, "-fraE",fra.CompE,"-rep_",repl,"_error_Effect_cor_final.csv"))
 write.csv(error.mse.2, paste0("CompSD","-G.cor", G.cor,"-E.cor",E.cor,"-dirh2", dirh2, "-compGvar", compG.var, "-fraE",fra.CompE,"-rep_",repl,"_error_Effect_mse_final.csv"))
 write.csv(herit.2, paste0("CompSD","-G.cor", G.cor,"-E.cor",E.cor,"-dirh2", dirh2, "-compGvar", compG.var, "-fraE",fra.CompE,"-rep_",repl,"_herit_final.csv"))
 write.csv(fra.Er.2, paste0("CompSD","-G.cor", G.cor,"-E.cor",E.cor,"-dirh2", dirh2, "-compGvar", compG.var, "-fraE",fra.CompE,"-rep_",repl,"_fraCompE_final.csv"))
 write.csv(varGs.2, paste0("CompSD","-G.cor", G.cor,"-E.cor",E.cor,"-dirh2", dirh2, "-compGvar", compG.var, "-fraE",fra.CompE,"-rep_",repl,"_varGs_final.csv"))
} #end of function




##function to call analysis - demonstration
#Gcor <- c(0,0.4)
#Ecor <- c(0,0.8)
#dirh2 <- 0.7
#compGvar < c(0.4,0.2,0.1)  
#fraCompE <- c(0.9,0.6,0.3)
#plotSize <- c(1,2)
#repl <- 1

#for (i in c(1:length(Gcor))){
#  for (j in c(1:length(Ecor))){
#    for (k in c(1:length(compGvar))){
#      for (l in c(1:length(fraCompE))){
#        G.cor <- Gcor[i]
#        E.cor <- Ecor[j]
#        compG.var <- compGvar[k]
#        fra.CompE <- fraCompE[l]
#        geno_compet_simulate_snp(madII,G.cor,E.cor, 0.7, compG.var, fra.CompE, c(1,2), 1))
#      }
#    }
#  }
#}


