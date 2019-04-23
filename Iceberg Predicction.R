# Iceberg Prediction ####

library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(MCMCglmm)

tFullSTMatrix <- 1 - (FullSTMatrix - min(FullSTMatrix))/max(FullSTMatrix)

NonEutherians <- c("Diprotodontia",
                   "Dasyuromorphia",
                   "Paucituberculata",
                   "Didelphimorphia",
                   "Microbiotheria",
                   "Peramelemorphia", 
                   "Notoryctemorphia",
                   "Monotremata")

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

NonEutherianSp <- Panth1[Panth1$hOrder%in%NonEutherians,"Sp"]

IcebergAdjList <- list(CurrentAdj, FutureAdj1) #, FutureAdj2, FutureAdj3, FutureAdj4)
PredReps <- c("Currents", paste0("Futures", 1:4))
names(IcebergAdjList) <- PredReps[1:2]

CORES = 1

for(x in 1:length(IcebergAdjList)){
  
  FileLoc <- paste0("Iceberg Output Files", PredReps[1])
  
  FullRangeAdj <- IcebergAdjList[[x]]
  
  AllMammals <- intersect(rownames(FullSTMatrix), rownames(FullRangeAdj)) %>% setdiff(NonEutherianSp)
  
  AllMammals <- sort(AllMammals)
  
  AllMammalMatrix <- data.frame(
    Sp = as.character(rep(AllMammals,each = length(AllMammals))),
    Sp2 = as.character(rep(AllMammals,length(AllMammals))),
    Space = c(FullRangeAdj[AllMammals,AllMammals]),
    Phylo = c(tFullSTMatrix[AllMammals,AllMammals])
  ) %>% 
    mutate(Phylo = (Phylo - min(Phylo))/max(Phylo - min(Phylo)),
           Gz = as.numeric(Space==0)) %>% droplevels
  
  UpperMammals <- which(upper.tri(FullSTMatrix[AllMammals, AllMammals], diag = T))
  
  AllMammaldf <- AllMammalMatrix[-UpperMammals,]
  
  N = nrow(AllMammaldf)
  
  load("Output Files/BAMList.Rdata")
  
  SpCoefNames <- names(BAMList[[1]]$coef)[substr(names(BAMList[[1]]$coef),1,5)=="SppSp"]
  SpCoef <- BAMList[[1]]$coef[SpCoefNames]
  
  Divisions = round(seq(0, nrow(AllMammaldf), length = 21))
  
  print("Prediction Effects!")
  
  if(file.exists(paste0(FileLoc,"/AllPredictions1b.Rdata"))) load(paste0(FileLoc,"/AllPredictions1b.Rdata")) else{
    
    FakeSpp <- matrix(0 , nrow = N, ncol = length(SpCoef))# %>% as("dgCMatrix")
    
    AllMammaldf$Spp <- FakeSpp; remove(FakeSpp)
    AllMammaldf$MinCites <- mean(c(log(FinalHostMatrix$MinCites+1), log(FinalHostMatrix$MinCites.Sp2+1)))
    AllMammaldf$Domestic <- 0
    
    AllPredictions1b <- predict.bam(BAMList[[1]], 
                                    newdata = AllMammaldf, # %>% select(-Spp),
                                    type = "terms",
                                    exclude = "Spp")
    
    AllPredictions1b <- mclapply(2:21, function(i){
      predict.bam(BAMList[[1]], 
                  newdata = AllMammaldf[(Divisions[i-1]+1):Divisions[i],], 
                  type = "terms",
                  exclude = "Spp")
    }, mc.cores = CORES)
    
    save(AllPredictions1b, file = paste0(FileLoc,"/AllPredictions1b.Rdata"))
    
  }
  
  print("Predicting All Links!")
  
  if(file.exists(paste0(FileLoc, "/AllPredList.Rdata"))) load(paste0(FileLoc, "/AllPredList.Rdata")) else{
    
    AllIntercept <- attr(AllPredictions1b[[1]], "constant")
    
    AllPredictions <- lapply(AllPredictions1b, as.data.frame) %>% bind_rows
    
    AllPredictions[,"Intercept"] <- AllIntercept
    
    AllPredList <- parallel::mclapply(1:100, function(x){ # to do something non-specific
      
      AllPredictions[,"Spp"] <- sample(SpCoef, N, replace = T) + 
        sample(SpCoef, N, replace = T)
      
      BinPred <- rbinom(n = N,
                        prob = logistic(rowSums(AllPredictions)),
                        size  = 1)
      
      BinPred
      
    }, mc.cores = CORES)
    
    save(AllPredList, file = paste0(FileLoc, "/AllPredList.Rdata"))
    
  }
  
  PredDF1 <- data.frame(AllPredList)
  
  print("Simulating All Networks!")
  
  # Simulating the network #####
  
  if(file.exists(paste0(FileLoc, "/AllSims.Rdata"))) load(paste0(FileLoc,"/AllSims.Rdata")) else{
    
    AllSims <- parallel::mclapply(1:length(AllPredList), function(i){
      
      AssMat <- matrix(NA, 
                       nrow = length(AllMammals), 
                       ncol = length(AllMammals))
      
      AssMat[lower.tri(AssMat)] <- round(AllPredList[[i]])
      AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
      diag(AssMat) <- 0
      dimnames(AssMat) <- list(union(AllMammaldf$Sp,AllMammaldf$Sp2),
                               union(AllMammaldf$Sp,AllMammaldf$Sp2))
      
      as(AssMat, "dgCMatrix")
      
    }, mc.cores = CORES)
    
    if(length(which(sapply(AllSims, is.null)))>0){
      AllSims <- AllSims[-which(sapply(AllSims, is.null))]
      print("Something went wrong UGH")
      print(paste("New Length = ", length(AllSims)))
    }
    
    save(AllSims, file = paste0(FileLoc, "/AllSims.Rdata"))
    
  }
  
  # Making summed matrix ####
  
  print("Summing Matrix!")
  
  if(file.exists(paste0(FileLoc,"/AllSums.Rdata"))) load(paste0(FileLoc,"/AllSums.Rdata")) else{
    
    AllPredDF <- AllPredList %>% as.data.frame()
    
    #AllPredSums <- apply(AllPredDF,1,sum)
    AllPredSums <- logistic(rowSums(AllPredictions))
    
    AssMat <- matrix(NA, 
                     nrow = length(AllMammals), #length(union(AllMammaldf$Sp,AllMammaldf$Sp2)), 
                     ncol = length(AllMammals)) #length(union(AllMammaldf$Sp,AllMammaldf$Sp2)))
    
    AssMat[lower.tri(AssMat)] <- AllPredSums
    AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
    diag(AssMat) <- 0
    
    dimnames(AssMat) <- list(AllMammals,
                             AllMammals)
    
    AllSums <- as(AssMat, "dgCMatrix")
    
    save(AllSums, file = paste0(FileLoc,"/AllSums.Rdata"))
    
  }
  
}
