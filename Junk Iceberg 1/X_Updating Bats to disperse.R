
Species <- SpeciesList %>% unlist %>% sort()

PredReps <- c("Currents", paste0("Futures", 1:4))

PipelineReps <- LETTERS[1:4]

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)

Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

Panth1 %>% filter(hOrder == "Chiroptera") %>% pull(Sp) ->
  BatSpecies

BatSpecies <- intersect(BatSpecies, Species)

i = 1

for(i in i:length(BatSpecies)){
  
  Sp <- BatSpecies[i]
  
  print(Sp)
  
  GretCDF <- readRDS(file = paste0("Iceberg Input Files/GretCDF/Currents/", Sp, ".rds")) %>%
    as.matrix %>% as.data.frame()
  
  GretCDF[,paste0("Buffer",c("Climate","ClimateLandUse"))] <- 1
  
  GretCDF %>%
    as.matrix %>% as("dgCMatrix") %>%
    saveRDS(file = paste0("Iceberg Input Files/GretCDF/Currents/", Sp, ".rds"))
  
  GretCDF <- readRDS(file = paste0("Iceberg Input Files/GretCDF/Futures/", Sp, ".rds")) %>%
    as.matrix %>% as.data.frame()
  
  GretCDF[,glue::glue("BufferClimateLandUse.{PredReps[2:5]}")] <-
    GretCDF[,glue::glue("ClimateLandUse.{PredReps[2:5]}")]
  
  GretCDF[,glue::glue("BufferClimate.{PredReps[2:5]}")] <-
    GretCDF[,glue::glue("Climate.{PredReps[2:5]}")]
  
  GretCDF %>%
    as.matrix %>% as("dgCMatrix") %>%
    saveRDS(file = paste0("Iceberg Input Files/GretCDF/Futures/", Sp, ".rds"))
  
}

# Bat richness ####

AllMammaldf <- readRDS("~/Albersnet/Iceberg Output Files/AllMammaldf.rds")

Sharing <- AllMammaldf %>% 
  dplyr::select(Sp, Sp2, 
                starts_with("Sharing.Futures"))

"Iceberg Input Files/GretCDF/Currents" %>% 
  list.files(full.names = T) ->
  CurrentFiles

paste0("Iceberg Input Files/GretCDF/Currents") %>% 
  list.files() %>% str_remove(".rds$") ->
  names(CurrentFiles)

Species <- intersect(Species, names(CurrentFiles))

FocalSp = BatSpecies[1]

GridDF <- readRDS(CurrentFiles[[FocalSp]]) %>%
  as.matrix %>%
  as.data.frame() %>%
  dplyr::select(X, Y, Climate, ClimateLandUse)

SubSharing <- Sharing %>% filter(Sp%in%FocalSp|Sp2%in%FocalSp)

GridDF[,paste0("Sharing.Currents",LETTERS[1:2])] <- cbind(
  sum(SubSharing$Sharing.CurrentsA)*GridDF[,"ClimateLandUse"],
  sum(SubSharing$Sharing.CurrentsB)*GridDF[,"Climate"]
)

for(a in a:length(BatSpecies)){
  
  FocalSp = BatSpecies[a]
  
  print(FocalSp)
  
  if(a %% 1000 == 0)  print(CurrentFiles[a])
  
  SubGretCDF <- readRDS(CurrentFiles[[FocalSp]]) %>% as.matrix %>% as.data.frame()# %>% dplyr::select(Climate, ClimateLandUse)
  
  SubSharing <- Sharing %>% filter(Sp%in%FocalSp|Sp2%in%FocalSp)
  
  SubGretCDF[,paste0("Sharing.Currents",LETTERS[1:2])] <- cbind(
    sum(SubSharing$Sharing.CurrentsA)*SubGretCDF[,"ClimateLandUse"],
    sum(SubSharing$Sharing.CurrentsB)*SubGretCDF[,"Climate"]
  )
  
  GridDF[,c("Climate", "ClimateLandUse",paste0("Sharing.Currents",LETTERS[1:2]))] <-
    GridDF[,c("Climate", "ClimateLandUse",paste0("Sharing.Currents",LETTERS[1:2]))] +
    SubGretCDF[,c("Climate", "ClimateLandUse",paste0("Sharing.Currents",LETTERS[1:2]))]
  
}

saveRDS(GridDF, file = "Iceberg Output Files/BatCurrentsGridDF.rds")

GridDF %>% ggplot(aes(X,Y, fill = ClimateLandUse)) + geom_tile()

# Future distribution
print("Futures!")

FutureCDFList <- list()

Sharing <- AllMammaldf %>% 
  dplyr::select(Sp, Sp2, 
                starts_with("Sharing.Futures"))


paste0("Iceberg Input Files/GretCDF/Futures") %>% 
  list.files(full.names = T) ->
  FutureFiles

paste0("Iceberg Input Files/GretCDF/Futures") %>% 
  list.files() %>% str_remove(".rds$") ->
  names(FutureFiles)

FocalSp = BatSpecies[1]

GridDF <- readRDS(FutureFiles[[FocalSp]]) %>% 
  as.matrix %>% 
  as.data.frame()

FutureCDFList[[FocalSp]] <- GridDF %>% dplyr::select(X,Y,contains("Futures"), -starts_with("LandUse"))

SubSharing <- Sharing %>% filter(Sp%in%FocalSp|Sp2%in%FocalSp)

for(x in 2:length(PredReps)){
  
  GridDF[,paste0("Sharing.",PredReps[x],PipelineReps)] <- cbind(
    sum(SubSharing[[paste0("Sharing.",PredReps[x],"A")]])*GridDF[,paste0("BufferClimateLandUse.",PredReps[x])],
    sum(SubSharing[[paste0("Sharing.",PredReps[x],"B")]])*GridDF[,paste0("BufferClimate.",PredReps[x])],
    sum(SubSharing[[paste0("Sharing.",PredReps[x],"C")]])*GridDF[,paste0("ClimateLandUse.",PredReps[x])],
    sum(SubSharing[[paste0("Sharing.",PredReps[x],"D")]])*GridDF[,paste0("Climate.",PredReps[x])]
    
  )
}

for(a in 2:length(BatSpecies)){
  
  FocalSp = BatSpecies[a]
  
  print(FocalSp)
  
  if(a %% 1000 == 0)  print(FutureFiles[a])
  
  SubGretCDF <- readRDS(FutureFiles[[FocalSp]]) %>% 
    as.matrix %>% 
    as.data.frame()
  
  FutureCDFList[[FocalSp]] <- 
    SubGretCDF %>% 
    dplyr::select(X, Y, contains("Futures"), -starts_with("LandUse"))
  
  SubSharing <- Sharing %>% filter(Sp%in%FocalSp|Sp2%in%FocalSp)
  
  for(x in 2:length(PredReps)){
    
    SubGretCDF[,paste0("Sharing.",PredReps[x],PipelineReps)] <- cbind(
      sum(SubSharing[[paste0("Sharing.",PredReps[x],"A")]])*SubGretCDF[,paste0("BufferClimateLandUse.",PredReps[x])],
      sum(SubSharing[[paste0("Sharing.",PredReps[x],"B")]])*SubGretCDF[,paste0("BufferClimate.",PredReps[x])],
      sum(SubSharing[[paste0("Sharing.",PredReps[x],"C")]])*SubGretCDF[,paste0("ClimateLandUse.",PredReps[x])],
      sum(SubSharing[[paste0("Sharing.",PredReps[x],"D")]])*SubGretCDF[,paste0("Climate.",PredReps[x])]
    )
    
  }
  
  GridDF[,2:ncol(GridDF)] <- GridDF[,2:ncol(GridDF)] + SubGretCDF[,2:ncol(SubGretCDF)]
  
}

saveRDS(GridDF, file = "Iceberg Output Files/BatFuturesGridDF.rds")


# Removing pinnipeds ####


lapply(IcebergAdjList, function(a) lapply(a, function(b){
  
  NonMarine <- setdiff(rownames(b), MarineSp)
  
  b[NonMarine,NonMarine]
  
})) -> IcebergAdjList


lapply(PredReps, function(a) lapply(LETTERS[1:4], function(b){
  
  IcebergAdjList[[b]][[a]] %>% saveRDS(file = paste0("Iceberg Output Files/", a, "RangeAdj", b, ".rds"))
  
}))

saveRDS(IcebergAdjList, file = "Iceberg Output Files/IcebergAdjList.rds")

