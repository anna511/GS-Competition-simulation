#Simulation of genotypic direct and competition effect##
#Ani A. Elias ## June, 2016###

makeMADIIfromFieldMap <- function(fieldMap){
  madII <- data.frame(Entry=data2$CLONE, Range=data2$Range, Column=data2$Column)
  return(madII)
}

#assuming madII has at least two replicates per genotype
#Direct and competition genotypic effects are related - have a correlation matrix 
#start the simulation function with the correlation matrix
#assuming unit variance for the direct effect, a variance for competiton (compG.var), and some correlation(G.cor)
#E.cor is the correlation between competition and residual errors
#dirh2 is the heritability of direct genetic effect
#fra.CompE is the fraction of total error which accounts for competiton error
#plotSize is the dimension of plots - plotSize[1] the width and plotSize[2] the length as they are arranged
 #with their longer edge sharing between ranges making the distance between adjascent ranges the smallest.

#test by ignoring the genetic and error correlation indicated by 0 for G.cor and E.cor
#G.cor <- {0,0.4}
#E.cor <- {0,0.8}
#dirh2 <- {0.7,0.3}
#because compG.var < (1-dirh2)/dirh2,
#when dirh2 = 0.7; compG.var < {0.4,0.2,0.1} 
  #when dirh2 = 0.3; compG.var <- {1,0.5,0.2,0.1}
#fra.CompE <- {0.9,0.6,0.3}
      
simulateComp_snp_SD <- function(madII, G.cor, E.cor, dirh2, compG.var, fra.CompE, plotSize){
  madII <- na.omit(madII)

##Genotypic covariance matrix - for direct effect
  G.cov <- G.cor * (sqrt(1 * compG.var))
  G <- matrix(c(1,G.cov,G.cov,compG.var), nrow=2,ncol=2)

  table(madII[,"Entry"] %in% rownames(snps)) 
  snp2.names <- intersect(madII[,"Entry"], rownames(snps))
  snps2 <- subset(snps, rownames(snps) %in% snp2.names)

  #calculate K matrix
  K <- A.mat(snps2-1)
  diag(K)<-diag(K)+1e-6

  #Genotypic effects
  sqrtCov <- chol(G)
  sqrtK <- chol(K)

    entries <- sort(unique(madII$Entry))
    entryEffects.G <- t(t(sqrtCov) %*% matrix(rnorm(2*length(entries)), nrow=2, ncol = length(entries)))
    direct <-  t(sqrtK) %*% entryEffects.G[,1]
    names(direct) <- entries
    compet <- t(sqrtK) %*% entryEffects.G[,2]     

##Total error variance and effects
  errVar <- ((1-dirh2) /dirh2) - compG.var
  compE <- fra.CompE * errVar
  residE <- (1-fra.CompE) * errVar

  E.cov <- E.cor * (sqrt(compE * residE))

  #Error effects
  E <- matrix(c(compE,E.cov,E.cov,residE), nrow=2,ncol=2)
  sqrtCovE <- chol(E)
  plotErrorEffects <- t(t(sqrtCovE) %*% matrix(rnorm(2*nrow(madII)), nrow=2, ncol = nrow(madII)))
  compet.Error <- plotErrorEffects[,1]
  residual.Error <- plotErrorEffects[,2]

#competition heritability
  comph2 <- compG.var/(1 + compG.var + errVar)

#to get the effects for the dataset
  idg <- factor(as.character(madII[,'Entry']), levels = entries)
   Z.base <- model.matrix(~idg -1) 
   colnames(Z.base) <- entries

  #distance matrix 
  rowVec<- madII$Range
  colVec <- madII$Column

  rowVec2 <- rowVec * plotSize[1] 
  rowDistMat2 <- sapply(rowVec2, function(rowNum) abs(rowNum - rowVec2))

  colVec2 <- colVec * plotSize[2]
  colDistMat2 <- sapply(colVec2, function(colNum) abs(colNum - colVec2))

  distMat <- sqrt(rowDistMat2^2 + colDistMat2^2)
  rownames(distMat) <- colnames(distMat) <- 1:nrow(distMat)

  #competition incidence matrix 
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

   
  #genotypic effects for the dataset
  directEffect <- Z.base %*% direct
  competinfluence <- Z.comp.1 %*% compet  
  competEffect <- Z.base %*% compet 

  #competition error effect w.r.t incidence matrix
  distMat.4 <- distMat
  distMat.4[distMat.4 > plotSize[1]] <- 0
  distMat.4[distMat.4 == plotSize[1]] <- 1
  compet.Errorinfluence <- distMat.4 %*% compet.Error  

  #final dataset
  fieldObs.1 <- cbind(directEffect, competinfluence, compet.Errorinfluence, residual.Error)
  rownames(fieldObs.1) <- NULL
  fieldObs <- data.frame(madII$Entry, madII$Range, madII$Column, competEffect, compet.Error, fieldObs.1, apply(fieldObs.1, 1, sum))
  colnames(fieldObs) <- c("Entry", "Range", "Column", "competEff", "competErr", "directG", "competG", "competE", "residE", "phenoVal")
  return(list(fieldObs=fieldObs, entryEffects=directEffect, competEffect = competEffect, compErrEff = compet.Error,
    comph2 = comph2, dirG.var.2 = var(direct), compG.var.2 = var(compet), 
    G.cor.2 = cor(direct,compet)))
}




  



