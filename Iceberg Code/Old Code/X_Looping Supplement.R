

# X_Iceberg Looping supplement ####

"~/Albersnet/Iceberg Code/Iceberg Greg GAMM Code" %>% 
  list.files(full.names = T) ->
  
  GAMMScripts

GAMMScripts[!str_detect(GAMMScripts, "Master")] -> GAMMScripts

CR <- 1

for(CR in CR:length(ClimateReps)){
  
  print(CR)
  
  setwd(paste0("~/Albersnet/Iceberg Files/", ClimateReps[CR]))
  
  library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr); library(SpRanger); library(cowplot)
  dir_create("Iceberg Extended Data")
  theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))
  PredReps <- c("Currents", paste0("Futures", 1:4))
  PipelineReps <- LETTERS[1:4]
  SpaceVars <- paste0(paste("Space", PredReps, sep = "."),rep(PipelineReps, each = length(PredReps)))
  SharingVars <- paste0(paste("Sharing",PredReps, sep = "."), rep(PipelineReps, each = length(PredReps)))
  names(SpaceVars) <- names(SharingVars) <- paste0(PredReps,rep(PipelineReps, each = length(PredReps)))
  Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
    dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)
  Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")
  AllMammaldf <- readRDS("Iceberg Output Files/AllMammaldf.rds")
  # Import #####
  CurrentsGridDF <- readRDS("~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/CurrentsGridDF.rds")
  FuturesGridDF <- readRDS("Iceberg Output Files/FuturesGridDF.rds")
  library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(SpRanger); library(raster)
  library(sf); library(fasterize);library(ggregplot); library(igraph);library(maptools)
  PredReps <- c("Currents", paste0("Futures", 1:4))
  PipelineReps <- LETTERS[1:4]
  AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
  AreaValues <- raster::values(AreaRaster)
  IcebergAdjList <- readRDS(paste0("Iceberg Output Files/","IcebergAdjList.rds"))
  NewEncountersList <- readRDS(paste0("Iceberg Output Files/","NewEncounters.rds"))
  AllMammaldf <- readRDS(paste0("Iceberg Output Files/","AllMammaldf.rds"))
  SpaceVars <- paste0(paste("Space", PredReps, sep = "."),rep(PipelineReps, each = length(PredReps)))
  SharingVars <- paste0(paste("Sharing",PredReps, sep = "."), rep(PipelineReps, each = length(PredReps)))
  names(SpaceVars) <- names(SharingVars) <- paste0(PredReps,rep(PipelineReps, each = length(PredReps)))
  # Iceberg General Maps ####
  GridList <- NEGridList <- list()
  CORES = 65
  t1 <- Sys.time()
  # Blanks
  blank <- matrix(0,360*2,720*2) # proper resolution
  blank <- raster(blank)
  extent(blank) <- c(-180,180,-90,90)
  projection(blank) <- CRS("+proj=longlat +datum=WGS84")
  UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
  Land = which(raster::values(UniversalBlank)==0)
  Sea = which(is.na(raster::values(UniversalBlank)))
  Species <-
    paste0("Iceberg Input Files/GretCDF/Currents") %>% list.files() %>% str_remove(".rds$") %>%
    sort
  paste0("Iceberg Input Files/GretCDF/Currents") %>%
    list.files(full.names = T) ->
    CurrentFiles
  paste0("Iceberg Input Files/GretCDF/Currents") %>%
    list.files() %>% str_remove(".rds$") ->
    names(CurrentFiles)
  paste0("Iceberg Input Files/GretCDF/Futures") %>%
    list.files(full.names = T) ->
    FutureFiles
  paste0("Iceberg Input Files/GretCDF/Futures") %>%
    list.files() %>% str_remove(".rds$") ->
    names(FutureFiles)
  FutureCDFList <- list()
  PredReps[2:5] %>% rep(4) %>%
    paste0(rep(LETTERS[1:4], each = 4)) -> NEReps
  recode1 <- c(A = "Climate + Land + Dispersal",
               B = "Climate + Dispersal",
               C = "Climate + Land",
               D = "Climate Only")
  recode2 <- c(Futures1 = 'RCP 2.6',
               Futures2 = 'RCP 4.5',
               Futures3 = 'RCP 7.0',
               Futures4 = 'RCP 8.5')
  PipeLabels <- rep(recode1, each = 4) %>% paste(rep(recode2, 4))
  names(PipeLabels) <- NEReps
  NewIntersects <- readRDS("Iceberg Output Files/NewIntersects.rds")
  NewIntersects %>% dplyr::select(X, Y, starts_with("Overlap.Futures")) %>%
    dplyr::select(-ends_with("B"), -ends_with("D")) %>%
    tidyr::gather("Key", "Encounters", -c(X,Y)) %>%
    mutate(Key = Key %>% str_remove("Overlap.")) %>%
    tidyr::separate(Key, sep = 8, into = c("RCP", "Pipeline")) %>%
    mutate(Pipeline = recode1[Pipeline],
           RCP = recode2[RCP]) %>%
    ggplot(., aes(X, Y, fill = log10(Encounters + 1))) + geom_tile() +
    scale_fill_continuous_sequential(palette = AlberPalettes[[1]],
                                     labels = c(1, 11, 101, 1001, 10001)) +
    facet_grid(RCP~Pipeline) +
    coord_fixed() +
    theme_void() + labs(fill = "Encounters") +
    theme(legend.position = "bottom", legend.text = element_text(angle = 45, hjust = 1)) +
    ggsave(filename = paste0("Iceberg Extended Data/", "A_C_Encounters", ".jpg"), width = 10, height = 8)
  
  NewIntersects %>% dplyr::select(X, Y, starts_with("DeltaOverlapSharing.Futures")) %>%
    dplyr::select(-ends_with("B"), -ends_with("D")) %>%
    tidyr::gather("Key", "NewSharings", -c(X,Y)) %>%
    mutate(Key = Key %>% str_remove("DeltaOverlapSharing.")) %>%
    tidyr::separate(Key, sep = 8, into = c("RCP", "Pipeline")) %>%
    mutate(Pipeline = recode1[Pipeline],
           RCP = recode2[RCP]) %>%
    ggplot(., aes(X, Y, fill = log10(NewSharings + 1))) + geom_tile() +
    scale_fill_continuous_sequential(palette = AlberPalettes[[1]],
                                     labels = c(1, 11, 101, 1001)) +
    facet_grid(RCP~Pipeline) +
    coord_fixed() +
    theme_void() + labs(fill = "New sharings") +
    theme(legend.position = "bottom", legend.text = element_text(angle = 45, hjust = 1)) +
    ggsave(filename = paste0("Iceberg Extended Data/", "A_C_NewSharings", ".jpg"), width = 10, height = 8)
  
}
