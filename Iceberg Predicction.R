# Iceberg Prediction ####

library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(MCMCglmm); library(SpRanger)

PredNetworkList <- list()

load("Output Files/BAMList.Rdata")

NonEutherians <- c("Diprotodontia",
                   "Dasyuromorphia",
                   "Paucituberculata",
                   "Didelphimorphia",
                   "Microbiotheria",
                   "Peramelemorphia", 
                   "Notoryctemorphia",
                   "Monotremata")

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

NonEutherianSp <- Panth1[Panth1$hOrder%in%NonEutherians,"Sp"]

tFullSTMatrix <- 1 - (FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp] - 
                        min(FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp]))/
  max(FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp])


# IcebergAdjList <- list(CurrentAdj) #, FutureAdj1) #, FutureAdj2, FutureAdj3, FutureAdj4)
# PredReps <- c("Currents", paste0("Futures", 1:4))
# names(IcebergAdjList) <- PredReps[1]

rownames(IcebergAdjList[[2]]) <- colnames(IcebergAdjList[[2]]) <- rownames(IcebergAdjList[[2]]) %>% str_replace("[.]", "_")

CORES = 1

for(x in 1:length(IcebergAdjList)){
  
  FileLoc <- paste0("Iceberg Output Files/", PredReps[x])
  
  FullRangeAdj <- IcebergAdjList[[x]]
  
  AllMammals <- intersect(rownames(FullSTMatrix), rownames(FullRangeAdj)) %>% setdiff(NonEutherianSp)
  
  AllMammals <- sort(AllMammals)
  
  AllMammalMatrix <- data.frame(
    Sp = as.character(rep(AllMammals,each = length(AllMammals))),
    Sp2 = as.character(rep(AllMammals,length(AllMammals))),
    Space = c(FullRangeAdj[AllMammals,AllMammals]),
    Phylo = c(tFullSTMatrix[AllMammals,AllMammals])
  ) %>% 
    mutate(#Phylo = (Phylo - min(Phylo))/max(Phylo - min(Phylo)),
      Gz = as.numeric(Space==0)) %>% droplevels
  
  UpperMammals <- which(upper.tri(FullSTMatrix[AllMammals, AllMammals], diag = T))
  
  AllMammaldf <- AllMammalMatrix[-UpperMammals,]
  
  N = nrow(AllMammaldf)
  
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
    
    save(AllPredictions1b, file = paste0(FileLoc,"/AllPredictions1b.Rdata"))
    
  }
  
  print("Predicting All Links!")
  
  
  AllIntercept <- attr(AllPredictions1b, "constant")
  
  AllPredictions <- AllPredictions1b %>% as.data.frame
  
  AllPredictions[,"Intercept"] <- AllIntercept
  
  DoNetworks = F
  
  if(DoNetworks){
    
    if(file.exists(paste0(FileLoc, "/AllPredList.Rdata"))&DoNetworks) load(paste0(FileLoc, "/AllPredList.Rdata")) else{
      
      
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
    
  }
  
  # Making summed matrix ####
  
  print("Summing Matrix!")
  
  if(file.exists(paste0(FileLoc,"/AllSums.Rdata"))) load(paste0(FileLoc,"/AllSums.Rdata")) else{
    
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
    
    save(AllSums, file = paste0(FileLoc,"/", PredReps[x], "AllSums.Rdata"))
    
  }
  
  PredNetworkList[[x]] <- AllSums
  
}

save(PredNetworkList, file = "Iceberg Output Files/PredNetworkList.Rdata")

MatrixPoints(PredNetworkList[[1]], 
             PredNetworkList[[2]], 
             Names = AllMammals[1:1000], 
             Axes = c("Current", "Future")) + 
  coord_fixed()

# Comparing currents and futures ####

AllMammals2 <- reduce(list(rownames(FullSTMatrix), 
                           rownames(IcebergAdjList[[1]]), 
                           rownames(IcebergAdjList[[2]])), 
                      intersect) %>% 
  setdiff(NonEutherianSp) %>% 
  sort()

AllMammalMatrix <- data.frame(
  Sp = as.character(rep(AllMammals2,each = length(AllMammals2))),
  Sp2 = as.character(rep(AllMammals2,length(AllMammals2))),
  SpaceCurrent = c(IcebergAdjList[[1]][AllMammals2,AllMammals2]),
  SpaceFuture = c(IcebergAdjList[[2]][AllMammals2,AllMammals2]),
  Phylo = c(tFullSTMatrix[AllMammals2,AllMammals2])
) %>% 
  mutate(#Phylo = (Phylo - min(Phylo))/max(Phylo - min(Phylo)),
    Gz = as.numeric(SpaceCurrent==0&SpaceFuture==0),
    DeltaOverlap = SpaceFuture - SpaceCurrent,
    NewOverlap = as.numeric(SpaceCurrent == 0 & SpaceFuture>0)) %>% droplevels

UpperMammals <- which(upper.tri(FullSTMatrix[AllMammals2, AllMammals2], diag = T))

AllMammaldf2 <- AllMammalMatrix[-UpperMammals,]

AllMammaldf2 <- AllMammaldf2 %>% #filter(SpaceCurrent == 0 & SpaceFuture>0) %>% 
  mutate(CurrentSharing = c(as.matrix(PredNetworkList[[1]])[AllMammals2, AllMammals2])[-UpperMammals],
         FutureSharing = c(as.matrix(PredNetworkList[[2]])[AllMammals2, AllMammals2])[-UpperMammals],
         DeltaSharing = FutureSharing - CurrentSharing) %>%
  slice(order(DeltaSharing, decreasing = T)) %>%
  inner_join(Panth1[,c("Sp", "hFamily", "hOrder")], by = c("Sp" = "Sp")) %>%
  inner_join(Panth1[,c("Sp", "hFamily", "hOrder")], by = c("Sp2" = "Sp"))

NewEncounters <- AllMammaldf2 %>% filter(SpaceCurrent == 0 & SpaceFuture>0)

NewEncounters %>% group_by(hOrder.x, hOrder.y) %>% 
  summarise(DeltaSharing = mean(DeltaSharing), N = n()) %>% 
  ggplot(aes(hOrder.x, hOrder.y, fill = DeltaSharing)) + 
  geom_tile() + 
  coord_fixed() + 
  geom_text(aes(label = N)) + 
  scale_fill_continuous_sequential(palette = "reds") 

AllMammaldf2 %>% slice(order(DeltaSharing, decreasing = T)) %>% head(25)

# Making the raster overlaps ####

names(RasterBrick) <- names(RasterBrick) %>% str_replace("[.]", "_")

NewOverlapSp <- unique(c(NewEncounters$Sp, NewEncounters$Sp2))

BinaryNewOverlaps <- as.matrix(IcebergAdjList[[1]])[AllMammals2,AllMammals2]==0&as.matrix(IcebergAdjList[[2]])[AllMammals2,AllMammals2]>0

NewIntersects <- IntersectGet(RasterBrick, 
                              Names = NewOverlapSp, 
                              Predicate = BinaryNewOverlaps)
