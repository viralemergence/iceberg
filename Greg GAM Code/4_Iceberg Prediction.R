
# 4_Iceberg Sharing Prediction ####

# Iceberg Prediction ####

library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(MCMCglmm); library(SpRanger)

PredNetworkList <- list()

load("Iceberg Output Files/BAMList.Rdata")

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

PredReps <- c("Currents", paste0("Futures", 1:4))

rownames(IcebergAdjList[[2]]) <- colnames(IcebergAdjList[[2]]) <- rownames(IcebergAdjList[[2]]) %>% 
  str_replace("[.]", "_") %>% str_replace(" ", "_")

# Making the prediction data frame ####

AllMammals <- reduce(lapply(IcebergAdjList, rownames), 
                     intersect) %>% intersect(rownames(FullSTMatrix)) %>%
  setdiff(NonEutherianSp) %>% 
  sort()

AllMammalMatrix <- data.frame(
  Sp = as.character(rep(AllMammals,each = length(AllMammals))),
  Sp2 = as.character(rep(AllMammals,length(AllMammals))),
  Phylo = c(tFullSTMatrix[AllMammals,AllMammals])
) %>%
  inner_join(Panth1[,c("Sp", "hFamily", "hOrder")], by = c("Sp" = "Sp")) %>%
  inner_join(Panth1[,c("Sp", "hFamily", "hOrder")], by = c("Sp2" = "Sp")) %>% droplevels()

SpaceVars <- paste("Space",PredReps[1:length(IcebergAdjList)], sep = ".")
SharingVars <- paste("Sharing",PredReps[1:length(IcebergAdjList)], sep = ".")

names(SpaceVars) <- names(SharingVars) <- PredReps

AllMammalMatrix[,SpaceVars] <-
  IcebergAdjList %>% lapply(function(a) c(as.matrix(a[AllMammals,AllMammals]))) %>% bind_cols

AllMammalMatrix[,paste0("Delta",SpaceVars[2:length(SpaceVars)])]<-
  apply(AllMammalMatrix[,SpaceVars[2:length(SpaceVars)]], 2, function(a) a - AllMammalMatrix$Space.Currents)

UpperMammals <- which(upper.tri(tFullSTMatrix[AllMammals, AllMammals], diag = T))

AllMammaldf <- AllMammalMatrix[-UpperMammals,]

N = nrow(AllMammaldf)

FakeSpp <- matrix(0 , nrow = N, ncol = length(SpCoef))# %>% as("dgCMatrix")

SpCoefNames <- names(BAMList[[1]]$coef)[substr(names(BAMList[[1]]$coef),1,5)=="SppSp"]
SpCoef <- BAMList[[1]]$coef[SpCoefNames]

AllMammaldf$Spp <- FakeSpp; remove(FakeSpp)
AllMammaldf$MinCites <- mean(c(log(FinalHostMatrix$MinCites+1), log(FinalHostMatrix$MinCites.Sp2+1)), na.rm = T)
AllMammaldf$Domestic <- 0

PredNetworkList <- list()

if(file.exists(paste0("Iceberg Output Files/",Method,"/PredNetworkList.Rdata"))) load(paste0("Iceberg Output Files/",Method,"/PredNetworkList.Rdata")) else{
  
  for(x in 2:length(IcebergAdjList)){
    
    print(PredReps[x])
    
    FileLoc <- paste0("Iceberg Output Files/",Method,"/", PredReps[x])
    
    AllMammaldf$Space = AllMammaldf[,SpaceVars[x]]
    AllMammaldf$Gz = as.numeric(AllMammaldf$Space==0)
    
    print("Prediction Effects!")
    
    AllPredictions1b <- predict.bam(BAMList[[1]], 
                                    newdata = AllMammaldf, # %>% select(-Spp),
                                    type = "terms",
                                    exclude = "Spp")
    
    AllIntercept <- attr(AllPredictions1b, "constant")
    
    AllPredictions <- AllPredictions1b %>% as.data.frame
    
    AllPredictions[,"Intercept"] <- AllIntercept
    
    AllPredList <- mclapply(1:1000, function(a){
      
      AllPredictions[,"Spp"] <- sample(SpCoef, N, replace = T) + 
        sample(SpCoef, N, replace = T)
      
      # Making summed matrix ####
      
      AllPredSums <- logistic(rowSums(AllPredictions)) %>% 
        rbinom(n = N,
               prob = .,
               size  = 1)
      
      return(AllPredSums)
      
    })
    
    PredNetworkList[[PredReps[x]]] <- AllPredList
    
  }
  
  save(PredNetworkList, file = paste0("Iceberg Output Files/",Method,"/PredNetworkList.Rdata"))
  
}

AllMammaldf <- AllMammaldf %>% dplyr::select(-Spp)

PredSums = list()

for(i in PredReps){
  
  print(i)
  
  PredSums[[i]] <- PredNetworkList[[i]][[1]]
  
  for(j in 2:length(PredNetworkList[[i]])){
    
    if(j %% 50) print(j) 
    
    PredSums[[i]] <- PredSums[[i]] + PredNetworkList[[i]][[j]]
    
  }
  
}

DF <- PredSums %>% bind_cols

AllMammaldf[,SharingVars] <-
  DF/length(PredNetworkList[[1]])

AllMammaldf[,paste0("Delta",SharingVars[2:length(SharingVars)])]<-
  apply(AllMammaldf[,paste0(SharingVars[2:length(SharingVars)])], 2, function(a) a - AllMammaldf$Sharing.Currents)

# Making new encounters ####

NewEncountersList <- lapply(SpaceVars[2:length(SpaceVars)], function(a) AllMammaldf[AllMammaldf$Space.Currents==0&AllMammaldf[,a]>0,])

names(NewEncountersList) <- PredReps[2:5]

save(AllMammaldf, file = paste0("Iceberg Output Files/",Method,"/AllMammaldf.Rdata"))
save(NewEncountersList, file = paste0("Iceberg Output Files/",Method,"/NewEncounters.Rdata"))

AMDF1 <- AllMammaldf %>% rename(hOrder = hOrder.x, hFamily = hFamily.x) %>% 
  dplyr::select(-ends_with(".y")) %>% 
  dplyr::select(-ends_with(".x"))

AMDF2 <- AllMammaldf %>% 
  rename(hOrder = hOrder.y, hFamily = hFamily.y,
         Sp = Sp2, Sp2 = Sp) %>% 
  dplyr::select(-ends_with(".y")) %>% 
  dplyr::select(-ends_with(".x"))

AllMammalLong <- bind_rows(AMDF1, AMDF2)





