
# Rscript "Iceberg Code/06_Prediction.R" ####

# Iceberg Prediction ####

CORES <- 45

library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(SpRanger)

source("~/Albersnet/Iceberg Code/04_Model Data Import.R")

# PipelineReps <- LETTERS[1:4]

load(paste0("~/Albersnet/Iceberg Files/CHELSA/Output Files/",
            "BAMList.Rdata"))

setwd("Iceberg Files/CHELSA")

IcebergAdjList <- #readRDS("Iceberg Output Files/IcebergAdjList.rds")
  list(readRDS("Output Files/CurrentsLandUseRangeAdj.rds"),
       readRDS("Output Files/CurrentsRangeAdj.rds")) %>% 
  append(
    
    "Output Files/Ranges" %>% dir_ls() %>% map(readRDS)
    
  )

names(IcebergAdjList) <- 
  c("CLUCurrents", "CCurrents") %>% 
  c("Output Files/Ranges" %>% list.files %>% 
      str_remove(".rds$") %>% 
      str_remove("RangeAdj"))

CurrentSpecies <- rownames(IcebergAdjList$CCurrents)

for(x in 1:length(IcebergAdjList)){
  
  NewAdj <- IcebergAdjList[[x]]
  
  InsertSpecies <- setdiff(CurrentSpecies, rownames(NewAdj))
  
  if(length(InsertSpecies)>0){
    
    NewAdj <- NewAdj %>% data.frame()
    NewAdj[InsertSpecies,] <- 0; NewAdj[,InsertSpecies] <- 0
    NewAdj <- NewAdj %>% as.matrix
    
    IcebergAdjList[[x]] <- NewAdj[CurrentSpecies, CurrentSpecies]
  }
}

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

tFullSTMatrix <- 1 - (FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,
                                   !rownames(FullSTMatrix)%in%NonEutherianSp] - 
                        min(FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,
                                         !rownames(FullSTMatrix)%in%NonEutherianSp]))/
  max(FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,
                   !rownames(FullSTMatrix)%in%NonEutherianSp])

PredReps <- 
  paste0("~/Albersnet/Iceberg Files/", 
         "CHELSA", "/FinalRasters") %>% 
  list.files() %>% 
  setdiff("presen")

# Making the prediction data frame ####

AllMammals <- 
  IcebergAdjList %>% map(rownames) %>% reduce(intersect) %>% 
  intersect(rownames(FullSTMatrix)) %>%
  setdiff(NonEutherianSp) %>% 
  sort()

AllMammalMatrix <- data.frame(
  Sp = as.character(rep(AllMammals, each = length(AllMammals))),
  Sp2 = as.character(rep(AllMammals, length(AllMammals))),
  Phylo = c(tFullSTMatrix[AllMammals, AllMammals])
  
) %>%
  
  inner_join(Panth1[,c("Sp", "hFamily", "hOrder")], by = c("Sp" = "Sp")) %>%
  inner_join(Panth1[,c("Sp", "hFamily", "hOrder")], by = c("Sp2" = "Sp")) %>% droplevels()

SpaceVars <- paste0("Space.", names(IcebergAdjList))

SharingVars <- paste0("Sharing.", names(IcebergAdjList))

# names(SpaceVars) <- names(SharingVars) <- 
#   paste0(PredReps, rep(PipelineReps, each = length(PredReps)))

AllMammalMatrix[,SpaceVars] <-
  IcebergAdjList %>% lapply(function(a){
    
    c(as.matrix(a[AllMammals,AllMammals]))
    
  }) %>% bind_cols()

UpperMammals <- which(upper.tri(tFullSTMatrix[AllMammals, AllMammals], diag = T))

AllMammaldf <- AllMammalMatrix[-UpperMammals,]

AllMammaldf[,c("Sp", "Sp2")] %>% saveRDS("Output Files/SpeciesPairs.rds")

N = nrow(AllMammaldf)

Model <- BAMList[[1]]

SpCoefNames <- names(Model$coef)[substr(names(Model$coef),1,5)=="SppSp"]
SpCoef <- Model$coef[SpCoefNames]

FakeSpp <- matrix(0 , nrow = N, ncol = length(SpCoef))# %>% as("dgCMatrix")

AllMammaldf$Spp <- FakeSpp; remove(FakeSpp)

AllMammaldf$MinCites <- mean(c(log(FinalHostMatrix$MinCites+1), 
                               log(FinalHostMatrix$MinCites.Sp2+1)), na.rm = T)

AllMammaldf$Domestic <- 0

print("Prediction Effects!")

Random = F

# Predict the currents ####

if(file.exists("Output Files/CurrentSharing.rds")){
  
  AllMammaldf[,paste0(rep(c("Space.", "Sharing."), 2),
                      rep(c("CLUCurrents", "CCurrents"), each = 2))] <-
    
    readRDS("Output Files/CurrentSharing.rds")
  
}else{
  
  for(Rep in c("CLUCurrents", "CCurrents")){
    
    print(Rep)
    
    CurrentsVar <- paste0("Space.", Rep)
    
    if(substr(Rep, 2, 3) == "LU"){
      
      Model <- BAMList[[1]]
      
    }else{
      
      Model <- BAMList[[2]]
      
    }
    
    SpCoefNames <- names(Model$coef)[substr(names(Model$coef),1,5)=="SppSp"]
    SpCoef <- Model$coef[SpCoefNames]
    
    AllMammaldf$Space = AllMammaldf[,paste0("Space", ".", Rep)]
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
      
      AllMammaldf[,paste0("Sharing.", Rep)] <- 
        
        PredSums/length(AllPredList)
      
    }else{      
      
      PredSums <- logistic(rowSums(AllPredictions))
      
      AllMammaldf[,paste0("Sharing.", Rep)] <- 
        PredSums
      
    }
  }
  
  AllMammaldf[,paste0(rep(c("Space.", "Sharing."), 2),
                      rep(c("CLUCurrents", "CCurrents"), each = 2))] %>% 
    
    saveRDS("Output Files/Predictions/CurrentSharing.rds")
  
}

# Predict the futures ####

ToRun <-
  "Output Files/Predictions/" %>% list.files %>% 
  str_remove(".rds$") %>% 
  str_remove("^Predictions") %>% 
  setdiff("Output Files/Ranges" %>% list.files %>% 
            str_remove(".rds$") %>% 
            str_remove("RangeAdj"), .)

print(paste0("To run: ", ToRun))

for(Rep in ToRun){
  
  print(Rep)
  
  Pipeline <- Rep %>% str_split("[.]") %>% 
    map_chr(1)
  
  if(substr(Pipeline, 1, 3) == "LU"){
    
    Model <- BAMList[[1]]
    
    CurrentsVar <- "Space.CLUCurrents"
    
    CurrentsSharing <- "Sharing.CLUCurrents"
    
  }else{
    
    Model <- BAMList[[2]]
    
    CurrentsVar <- "Space.CCurrents"
    
    CurrentsSharing <- "Sharing.CCurrents"
    
  }
  
  SpCoefNames <- names(Model$coef)[substr(names(Model$coef),1,5)=="SppSp"]
  SpCoef <- Model$coef[SpCoefNames]
  
  AllMammaldf$Space = AllMammaldf[,paste0("Space", ".", Rep)]
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
    
    AllMammaldf[,paste0("Sharing.", Rep)] <- 
      PredSums/length(AllPredList)
    
  }else{      
    
    PredSums <- logistic(rowSums(AllPredictions))
    
    AllMammaldf[,paste0("Sharing.", Rep)] <- 
      PredSums
    
  }
  
  AllMammaldf[,paste0("DeltaSpace", ".", Rep)] <-
    AllMammaldf[,paste0("Space", ".", Rep)] - 
    AllMammaldf[,CurrentsVar]
  
  AllMammaldf[,paste0("DeltaSharing", ".", Rep)] <-
    AllMammaldf[,paste0("Sharing", ".", Rep)] - 
    AllMammaldf[,CurrentsSharing]
  
  AllMammaldf[,paste0(c("Space.", "Sharing.", 
                        "DeltaSpace.", "DeltaSharing."),
                      Rep)] %>% 
    saveRDS(paste0("Output Files/Predictions/Predictions", 
                   Rep, ".rds"))
  
  AllMammaldf[,paste0(c("Space.", "Sharing.", 
                        "DeltaSpace.", "DeltaSharing."),
                      Rep)] <- NULL
  
}

AllMammaldf <- AllMammaldf %>% dplyr::select(-Spp)

setwd(here::here())

stop()

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
