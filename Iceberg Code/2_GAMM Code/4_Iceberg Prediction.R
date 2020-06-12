
# Rscript "Final Iceberg Code/4_Iceberg Prediction.R" ####

# Iceberg Prediction ####

CORES <- 60

library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(SpRanger)

source("~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code/2_Iceberg Data Import.R")

PipelineReps <- LETTERS[1:4]

load(paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/",
            "BAMList.Rdata"))

IcebergAdjList <- readRDS("Iceberg Output Files/IcebergAdjList.rds")

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

if(CoryClimateReps[CR] == "gf"){
  
  PredReps <- c("Currents", paste0("Futures", 1:4))[c(1, 2, 4)]
  
}

# Making the prediction data frame ####

AllMammals <- reduce(lapply(IcebergAdjList$A, rownames), 
                     intersect) %>% intersect(rownames(FullSTMatrix)) %>%
  setdiff(NonEutherianSp) %>% 
  sort()

AllMammalMatrix <- data.frame(
  Sp = as.character(rep(AllMammals, each = length(AllMammals))),
  Sp2 = as.character(rep(AllMammals, length(AllMammals))),
  Phylo = c(tFullSTMatrix[AllMammals, AllMammals])
) %>%
  inner_join(Panth1[,c("Sp", "hFamily", "hOrder")], by = c("Sp" = "Sp")) %>%
  inner_join(Panth1[,c("Sp", "hFamily", "hOrder")], by = c("Sp2" = "Sp")) %>% droplevels()

SpaceVars <- paste0(paste("Space", PredReps, sep = "."),
                    rep(PipelineReps, each = length(PredReps)))

SharingVars <- paste0(paste("Sharing",PredReps, sep = "."), 
                      rep(PipelineReps, each = length(PredReps)))

names(SpaceVars) <- names(SharingVars) <- paste0(PredReps, rep(PipelineReps, each = length(PredReps)))

AllMammalMatrix[,SpaceVars] <-
  IcebergAdjList %>% lapply(function(a){
    
    lapply(a, function(b){
      
      c(as.matrix(b[AllMammals,AllMammals]))
      
    }) %>% bind_cols()
    
  }) %>% bind_cols

UpperMammals <- which(upper.tri(tFullSTMatrix[AllMammals, AllMammals], diag = T))

AllMammaldf <- AllMammalMatrix[-UpperMammals,]

N = nrow(AllMammaldf)

Model <- BAMList[[1]]

SpCoefNames <- names(Model$coef)[substr(names(Model$coef),1,5)=="SppSp"]
SpCoef <- Model$coef[SpCoefNames]

FakeSpp <- matrix(0 , nrow = N, ncol = length(SpCoef))# %>% as("dgCMatrix")

AllMammaldf$Spp <- FakeSpp; remove(FakeSpp)

AllMammaldf$MinCites <- mean(c(log(FinalHostMatrix$MinCites+1), 
                               log(FinalHostMatrix$MinCites.Sp2+1)), na.rm = T)

AllMammaldf$Domestic <- 0

Pipeline = "A"

print("Prediction Effects!")

Random = F

for(Pipeline in PipelineReps){
  
  print(Pipeline)
  
  Model <- BAMList[[2 - which(LETTERS == Pipeline) %% 2]]
  
  SpCoefNames <- names(Model$coef)[substr(names(Model$coef),1,5)=="SppSp"]
  SpCoef <- Model$coef[SpCoefNames]
  
  for(x in 1:length(PredReps)){
    
    print(PredReps[x])
    
    AllMammaldf$Space = AllMammaldf[,SpaceVars[paste0(PredReps[x], Pipeline)]]
    AllMammaldf$Gz = as.numeric(AllMammaldf$Space==0)
    
    AllPredictions1b <- predict.bam(Model, 
                                    newdata = AllMammaldf, # %>% select(-Spp),
                                    type = "terms",
                                    exclude = "Spp")
    
    AllIntercept <- attr(AllPredictions1b, "constant")
    
    AllPredictions <- AllPredictions1b %>% as.data.frame
    
    AllPredictions[,"Intercept"] <- AllIntercept
    
    if(Random){
      
      AllPredList <- mclapply(1:100, function(a){
        
        AllPredictions[,"Spp"] <- sample(SpCoef, N, replace = T) + 
          sample(SpCoef, N, replace = T)
        
        AllPredSums <- logistic(rowSums(AllPredictions))
        
        return(AllPredSums)
        
      }, mc.cores = CORES)
      
      PredSums <- AllPredList %>% bind_cols %>% rowSums
      
      AllMammaldf[,paste0(SharingVars[paste0(PredReps[x], Pipeline)])] <- PredSums/length(AllPredList)
      
    }else{      
      
      PredSums <- logistic(rowSums(AllPredictions))
      
      AllMammaldf[,paste0(SharingVars[paste0(PredReps[x], Pipeline)])] <- PredSums
      
    }
  }
}

AllMammaldf <- AllMammaldf %>% dplyr::select(-Spp)

for(i in 1:length(PipelineReps)){
  
  AllMammaldf[, paste0("Delta", SpaceVars[2:length(PredReps) + (i-1)*length(PredReps)])] <-
    
    apply(AllMammaldf[, SpaceVars[2:length(PredReps) + (i-1)*length(PredReps)]], 2, function(a){
      a - AllMammaldf[, paste0("Space.Currents", PipelineReps[i])]
      
    })
}

for(i in 1:length(PipelineReps)){
  
  AllMammaldf[, paste0("Delta", SharingVars[2:length(PredReps) + (i-1)*length(PredReps)])] <-
    
    apply(AllMammaldf[, SharingVars[2:length(PredReps) + (i-1)*length(PredReps)]], 2, function(a){
      a - AllMammaldf[, paste0("Sharing.Currents", PipelineReps[i])]
      
    })
}

# Making new encounters ####

NewEncountersList <- 
  
  lapply(PipelineReps, function(b){
    
    l1 <- lapply(PredReps[2:length(PredReps)], function(a){
      
      AllMammaldf[AllMammaldf[,paste0("Space.Currents",b)]==0&
                    AllMammaldf[,paste0("Space.", a, b)]>0,]
      
    })
    
    names(l1) <- PredReps[2:length(PredReps)]
    
    return(l1)
    
  })

names(NewEncountersList) <- PipelineReps

# Making old encounters ####

OldEncountersList <- 
  
  lapply(PipelineReps, function(b){
    
    l1 <- lapply(PredReps[2:length(PredReps)], function(a){
      
      AllMammaldf[AllMammaldf[,paste0("Space.Currents", b)]>0&
                    AllMammaldf[,paste0("Space.", a, b)]==0,]
      
    })
    
    names(l1) <- PredReps[2:length(PredReps)]
    
    return(l1)
    
  })

names(OldEncountersList) <- PipelineReps

saveRDS(AllMammaldf, file = "Iceberg Output Files/AllMammaldf.rds")
saveRDS(NewEncountersList, file = paste0("Iceberg Output Files/NewEncounters.rds"))
saveRDS(OldEncountersList, file = paste0("Iceberg Output Files/OldEncounters.rds"))
