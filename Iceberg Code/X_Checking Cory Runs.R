
# X_Comparing Cory Runs ####

FocalSp <- "Zygodontomys_brevicauda"

FocalClimateReps <- ClimateReps[1:9][-c(4:5)]

FocalClimateReps %>% map(~readRDS(paste0("Iceberg Files/", .x, "/Iceberg Input Files/GretCDF/Futures/",
                                         FocalSp, ".rds")) %>% as.matrix %>% as.data.frame) ->
  
  GretCDFList

GretCDFList %>% bind_rows(.id = "ClimateRep") %>%
  ggplot(aes(X, Y)) + 
  geom_tile(aes(fill = Climate.Futures1)) +
  facet_wrap(~ClimateRep)

GretCDFList %>% bind_rows(.id = "ClimateRep") %>%
  group_by(ClimateRep) %>% summarise_at("BufferClimateLandUse.Futures4", ~sum(.x, .na.rm = T))

GretCDFList %>% bind_rows(.id = "ClimateRep") %>%
  group_by(X, Y) %>% summarise(Sum = sum(BufferClimateLandUse.Futures4, na.rm = T)) %>% 
  pull(Sum) %>% table()

FocalClimateReps %>% map(~{
  
  Files <- paste0("~/Albersnet/Iceberg Files/", .x, 
                  "/Iceberg Input Files/MaxEnt/01_Raw/Futures1") %>% 
    list.files(pattern = FocalSp, full.names = T)
  
  raster::raster(Files)}) ->
  
  RasterList

plot(RasterList[[1]])

FocalClimateReps %>% 
  map(~readRDS(paste0("Iceberg Files/", .x, "/Iceberg Output Files/Futures1RangeAdjA.rds"))) ->
  
  GretCDFList

identical(GretCDFList[[1]], GretCDFList[[2]])

Sample <- sample(1:length(c(GretCDFList[[1]])), 1000)

qplot(c(GretCDFList[[1]])[Sample])
qplot(c(GretCDFList[[2]])[Sample])

qplot(c(GretCDFList[[1]])[Sample],
      c(GretCDFList[[2]])[Sample])

# Currents ####

"Climate1" %>% map(~{
  
  Files <- paste0("~/Albersnet/Iceberg Files/", .x, 
                  "/Iceberg Input Files/MaxEnt/01_Raw/Currents") %>% 
    list.files(pattern = FocalSp, full.names = T)
  
  raster::raster(Files)
  
}) ->
  
  CurrentRasterList

CurrentRasterList[[1]] %>% plot

readRDS(paste0("Iceberg Files/", "Climate1", "/Iceberg Input Files/GretCDF/Currents/",
               FocalSp, ".rds")) %>% as.matrix %>% as.data.frame ->
  
  CurrentGretCDF

CurrentGretCDF %>% 
  ggplot(aes(X, Y)) + 
  geom_tile(aes(fill = Climate))

# Summary of new encounters ####

FocalClimateReps %>% 
  map(~readRDS(paste0("Iceberg Files/", .x, "/Iceberg Output Files/IcebergAdjList.rds"))) ->
  
  GretCDFList


GretCDFList %>% map(~{
  
  .x[[1]][lower.tri(.x[[1]])&.x[[1]]==0]
  
})[[1]]

PipelineReps <- LETTERS[1:4]

PredReps <- c("Currents", paste0("Futures", 1:4))

# NewEncountersList <- 



GretCDFList %>% map(~{
  
  l2 <- lapply(PipelineReps, function(b){
    
    l1 <- lapply(PredReps[2:5], function(a){
      
      table(((.x[[b]][[1]]==0)&lower.tri(.x[[b]][[1]]))&
              ((.x[[b]][[a]]>0)&lower.tri(.x[[b]][[1]])))
      
    })
    
    names(l1) <- PredReps[2:5]
    
    return(l1)
    
  })
  
  names(l2) <- PipelineReps
  
  return(l2)
  
}) -> FullNEList


names(FullNEList) <- FocalClimateReps

FullNEList %>% unlist -> EncounterNumbers

EncounterNumbers[str_detect(names(EncounterNumbers), "TRUE")]

#FullNEList %>% map(~map(.x))

EncounterNumbers[str_detect(names(EncounterNumbers), "TRUE")] %>% 
  matrix(nrow = 4) %>% as.data.frame() -> EncounterDF

colnames(EncounterDF) <- paste0(rep(FocalClimateReps, each = 4), PipelineReps)

EncounterDF$RCP <- PredReps[2:5]

EncounterDF %>% gather("Key", "Value", -RCP) %>% 
  mutate(Pipeline = substr(Key, 9, 9),
         Climate = substr(Key, 8, 8)) -> 
  
  EncounterDF

EncounterDF %>% ggplot(aes(Pipeline, RCP, fill = log10(Value))) + 
  facet_grid(Climate~.) +
  geom_tile() +
  scale_fill_continuous_sequential(palette = "Reds")


EncounterDF %>% 
  ggplot(aes(RCP, Value, colour = Climate)) + geom_point() +
  facet_wrap(~Pipeline, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


EncounterDF %>% 
  ggplot(aes(as.numeric(as.factor(RCP)), Value, colour = Climate)) + 
  geom_point() + geom_line() + 
  facet_wrap(~Pipeline, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Getting names ####

GretCDFList %>% map(~{
  
  l2 <- lapply(PipelineReps, function(b){
    
    .x[[b]][[1]] %>% reshape2::melt() %>% 
      slice(which(lower.tri(.x[[b]][[1]]))) -> 
      
      Currents
    
    l1 <- lapply(PredReps[2:5], function(a){
      
      .x[[b]][[a]] %>% reshape2::melt() %>% 
        slice(which(lower.tri(.x[[b]][[1]]))) %>%
        slice(which(Currents$value==0)) %>%
        filter(value>0)
      
    })
    
    names(l1) <- PredReps[2:5]
    
    return(l1)
    
  })
  
  names(l2) <- PipelineReps
  
  return(l2)
  
}) -> NEList

GretCDFList %>% map(~{
  
  l2 <- lapply(PipelineReps, function(b){
    
    .x[[b]][[1]] %>% reshape2::melt() %>% 
      slice(which(lower.tri(.x[[b]][[1]]))) -> 
      
      Currents
    
    l1 <- lapply(PredReps[2:5], function(a){
      
      .x[[b]][[a]] %>% reshape2::melt() %>% 
        slice(which(lower.tri(.x[[b]][[1]]))) %>%
        slice(which(Currents$value==0)) %>%
        filter(value>0) %>% mutate(Pair = paste0(Var1, Var2)) %>%
        pull(Pair)
      
    })
    
    names(l1) <- PredReps[2:5]
    
    return(l1)
    
  })
  
  names(l2) <- PipelineReps
  
  return(l2)
  
}) -> NEListNames

NEListNames %>% unlist %>% table()
