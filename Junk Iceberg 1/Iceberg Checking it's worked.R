
# 5_Iceberg Making Map Data ####

library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(SpRanger); library(raster)
library(sf); library(fasterize);library(ggregplot); library(igraph);library(maptools)

# Rscript "Final Iceberg Code/5_Iceberg Mapping.R"

PredReps <- c("Currents", paste0("Futures", 1:4))
PipelineReps <- LETTERS[1:4]

AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

IcebergAdjList <- readRDS(paste0("Iceberg Output Files/","IcebergAdjList.rds"))
load(paste0("Iceberg Output Files/","AllMammaldf.Rdata"))
load(paste0("Iceberg Output Files/","NewEncounters.Rdata"))

SpaceVars <- paste0(paste("Space", PredReps, sep = "."),rep(PipelineReps, each = length(PredReps)))
SharingVars <- paste0(paste("Sharing", PredReps, sep = "."), rep(PipelineReps, each = length(PredReps)))

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

CurrentCDFList <- FutureCDFList <- list()

Species <- Species %>% sort %>% intersect(names(CurrentFiles)) #%>% intersect(names(FutureFiles))

CurrentFiles <- CurrentFiles[Species]
FutureFiles <- FutureFiles[Species]

PipelineReps <- LETTERS[1:4]


# Trying new encounters ####

Pipeline <- sample(PipelineReps, 1)
x <- sample(2:5, 1)

NewEncounters <- NewEncountersList[[Pipeline]][[PredReps[x]]]

FocalSp <- "Canis_lupus"
FocalSp <- "Meles_meles"
FocalSp <- "Cervus_elaphus"
FocalSp2 <- "Capra_falconeri"

CurrentCDF <- readRDS(CurrentFiles[[FocalSp]]) %>% as.matrix %>% as.data.frame()

FutureCDF <- readRDS(FutureFiles[[FocalSp]]) %>% as.matrix %>% as.data.frame()

GridCDF <- bind_cols(CurrentCDF, FutureCDF)

GridCDF %>% 
  gather("Prediction", "Value", Continent, Climate, BufferClimate, ClimateLandUse, BufferClimateLandUse, contains("Futures")) %>%
  pull(Prediction) %>% unique ->
  PrintNames

for(i in PrintNames){
  
  GridCDF %>%
    gather("Prediction", "Value", Continent, Climate, BufferClimate, ClimateLandUse, BufferClimateLandUse, contains("Futures")) %>% 
    filter(Prediction == i) %>%
    ggplot(aes(X, Y, fill = Value)) + geom_tile() +
    scale_fill_continuous_sequential(palette = AlberPalettes[[2]]) +
    coord_fixed() + ggtitle(i) +
    ggsave(paste0("Iceberg Figures/", FocalSp, ".", i, ".jpeg"), units = "mm", height = 200, width = 300)
  
}

CurrentCDF <- readRDS(CurrentFiles[[FocalSp2]]) %>% as.matrix %>% as.data.frame()

FutureCDF <- readRDS(FutureFiles[[FocalSp2]]) %>% as.matrix %>% as.data.frame()

PracticeDF2 <- bind_cols(CurrentCDF, FutureCDF)

PracticeDF2 %>% dplyr::select(X,Y,contains(paste0("Climate.Futures"))) %>% 
  mutate(BufferClimate = PracticeDF2$BufferClimate) %>%
  mutate(Climate = PracticeDF2[,"Climate"]) %>%
  mutate_at(vars(contains(paste0("Climate.Futures"))), function(a){
    
    ifelse(a==1&PracticeDF2$Climate==0, "Gain",
           ifelse(a==0&PracticeDF2$Climate==1, "Loss",
                  ifelse(a==1&PracticeDF2$Climate==1, "Kept", 0)))
    
  }) %>% 
  mutate_all(function(a) ifelse(a == 0, NA, a)) %>% gather("Key", "Value", 3:Climate) %>% 
  na.omit ->
  LongDF2

ggplot(LongDF2, aes(X, Y, fill = Value)) + geom_tile() + facet_wrap(~Key) + coord_fixed()

LongDF %>% mutate(Sp = FocalSp) %>% bind_rows(LongDF2 %>% mutate(Sp = FocalSp2)) ->
  LongerDF

WideDF <- left_join(LongDF, LongDF2, by = c("X", "Y"), suffix = c(FocalSp, FocalSp2))

WideDF %>% filter(Climate.Futures1D.Cervus_elaphus%in%c("Gain", "Kept")&Climate.Futures1D.Capra_falconeri%in%c("Gain", "Kept"))

LongerDF %>% ggplot(aes(X, Y, fill = Value)) + geom_tile() + geom_tile() + facet_wrap(~Key) + coord_fixed()

PracticeDF %>% dplyr::select(contains(paste0("Climate.Futures"))) %>% 
  mutate(BufferClimate = PracticeDF$BufferClimate) %>%
  mutate(Climate = PracticeDF[,"Climate"]) %>%
  mutate_at(vars(contains(paste0("Climate.Futures"))), function(a){
    
    ifelse(a==1&PracticeDF$Climate==0, "Gain",
           ifelse(a==0&PracticeDF$Climate==1, "Loss",
                  ifelse(a==1&PracticeDF$Climate==1, "Kept", 0)))
    
  }) %>% 
  #mutate_all(function(a) ifelse(a == 0, NA, a)) %>%
  na.omit() %>% 
  apply(2, table) -> Totals
