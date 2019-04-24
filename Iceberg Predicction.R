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

load("Iceberg Output Files/Futures1/Futures1IcebergRasterBrick.Rdata")

names(RasterBrick) <- names(RasterBrick) %>% str_replace("[.]", "_")

NewIntersectsManual <- list()

for(i in 1:nrow(NewEncounters)){
  
  print(i)
  
  NewIntersectsManual[[paste(NewEncounters[i,c("Sp","Sp2")], collapse = ".")]] <- 
    raster::intersect(RasterBrick[[NewEncounters[i,"Sp"]]], 
                      RasterBrick[[NewEncounters[i,"Sp2"]]])
  
}

saveRDS(NewIntersectsManual, file = "Iceberg Output Files/Futures1/NewIntersectsManual.rds")
save(NewIntersectsManual, file = "Iceberg Output Files/Futures1/NewIntersectsManual.Rdata")

i = length(list.files("GetIntersects"))

for(i in i:length(NewIntersectsManual)){
  print(i)
  writeRaster(NewIntersectsManual[[i]], file = paste0("GetIntersects/",names(NewIntersectsManual)[i],".tif"))
}

OverlapBrick <- raster::brick(NewIntersectsManual)
OverlapBrickSave <- OverlapBrick

sum.raw <- sum(OverlapBrick, na.rm=TRUE)
sum.raw[sum.raw==0] <- NA

for(i in 1:length(OverlapBrick)){
  print(i)
  name <- names(NewIntersectsManual[i])
  twoname <- strsplit(name,'\\.')[[1]]
  OverlapBrick[[i]][!is.na(OverlapBrick[[i]])] <- AllSums[twoname[1],twoname[2]]
  
}

mean.p <- mean(OverlapBrick, na.rm = TRUE)
mean.p[mean.p==0] <- NA

sum.p <- sum(OverlapBrick, na.rm = TRUE)
sum.p[sum.p==0] <- NA

par(mfrow = c(3,1))
plot(sum.raw, main = 'Raw')
plot(sum.p, main = 'Sum viral sharing')
plot(mean.p, main = 'Mean Viral Sharing')

writeRaster(sum.raw, file = "SumOverlap.tif")
writeRaster(sum.p, file = "SumSharing.tif")
writeRaster(mean.p, file = "MeanSharing.tif")

# Getting taxonomic patterns of predictions ####

Patterndf <- data.frame(
  Sp = AllMammals2,
  CurrentOverlaps = rowSums(IcebergAdjList[[1]][AllMammals2, AllMammals2]),
  FutureOverlaps = rowSums(IcebergAdjList[[2]][AllMammals2, AllMammals2])
) %>% mutate(DeltaOverlaps = FutureOverlaps - CurrentOverlaps) %>%
  mutate(CurrentSharing = rowSums(PredNetworkList[[1]][AllMammals2, AllMammals2]),
         FutureSharing = rowSums(PredNetworkList[[2]][AllMammals2, AllMammals2])) %>% mutate(
           DeltaSharing = FutureSharing - CurrentSharing
         ) %>% full_join(Panth1)

BarGraph(Patterndf, "hOrder", "DeltaSharing", Order = T, Just = T, Text = "N")
BarGraph(Patterndf, "hOrder", "CurrentSharing", Order = T, Just = T, Text = "N")
BarGraph(Patterndf, "hOrder", "FutureSharing", Order = T, Just = T, Text = "N")

colnames(Patterndf)[2:7] %>% lapply(function(a){
  
  BarGraph(Patterndf, "hOrder", a, Order = T, Just = T, Text = "N") + 
    theme(legend.position = "none") +
    scale_fill_discrete_sequential(palette = AlberPalettes[[1]])
  
}) %>% arrange_ggplot2(ncol = 3)

# Just new encounters

NewEncountersGraph <- graph_from_edgelist(NewEncounters[,c("Sp", "Sp2")] %>% as.matrix, directed = F)
E(NewEncountersGraph)$weights <- NewEncounters$CurrentSharing
E(NewEncountersGraph)$weights2 <- NewEncounters$FutureSharing

NewEncountersAdj <- get.adjacency(NewEncountersGraph, 
                                  attr = "weights"
) %>% 
  as.matrix #%>% as.data.frame

NewEncountersAdj2 <- get.adjacency(NewEncountersGraph, 
                                   attr = "weights2"
) %>% 
  as.matrix #%>% as.data.frame

NewEncountersPatterns <- NewEncountersAdj %>% 
  reshape2::melt(value.name = "CurrentSharing") %>% 
  rename(Sp = Var1, Sp2 = Var2) %>%
  inner_join(Panth1[,c("Sp", "hFamily", "hOrder")], by = c("Sp" = "Sp")) %>%
  inner_join(Panth1[,c("Sp", "hFamily", "hOrder")], by = c("Sp2" = "Sp")) %>% 
  right_join(NewEncountersAdj2 %>% 
               reshape2::melt(value.name = "FutureSharing") %>%
               rename(Sp = Var1, Sp2 = Var2)) %>%
  filter(!FutureSharing == 0) %>%
  mutate(#DeltaOverlap = 
    DeltaSharing = FutureSharing - CurrentSharing)

NewEncountersPatterns %>% group_by(hOrder.x, hOrder.y) %>% 
  summarise(FutureSharing = mean(FutureSharing), N = n()) %>% 
  ggplot(aes(hOrder.x, hOrder.y, fill = log(FutureSharing))) + 
  geom_tile() + 
  coord_fixed() + 
  geom_text(aes(label = N)) + 
  scale_fill_continuous_sequential(palette = "reds") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

colnames(Patterndf)[2:7] %>% lapply(function(a){
  
  BarGraph(Patterndf, "hOrder", a, Order = T, Just = T, Text = "N") + 
    theme(legend.position = "none") +
    scale_fill_discrete_sequential(palette = AlberPalettes[[1]])
  
}) %>% arrange_ggplot2(ncol = 3)

# Future network #### 

Futures1Graph <- graph_from_edgelist(AllMammaldf2[,c("Sp", "Sp2")] %>% as.matrix, directed = F)
E(Futures1Graph)$weights <- AllMammaldf2$CurrentSharing
E(Futures1Graph)$weights2 <- AllMammaldf2$FutureSharing



