
# X_Checking Thresholding ####

library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr); library(SpRanger); library(cowplot)
library(fs)

t1 <- Sys.time()

#print("Dropping Species!")

#source("Iceberg Code/Iceberg Greg ENM Code/00_Iceberg Species Dropping.R")

print("Doing Currents!")

dir_ls(paste0("Iceberg Files/","ThresholdMaps"), 
       recurse = T) %>% as.character() ->
  FullFiles

FullFiles <- FullFiles[str_detect(FullFiles, "[.]tif$")]

FullFiles %>% str_split("/") %>% map_chr(last) -> 
  names(FullFiles)

Species <- 
  names(FullFiles) %>% str_split("__") %>% map_chr(2) %>% unique %>% 
  intersect("Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents" %>% 
              list.files %>% str_remove(".rds$"))

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("Iceberg Files/Climate1/Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

# Grid areas
AreaRaster <- raster("Iceberg Files/Climate1/Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

# Land Use Data
iucndat <- read.csv('Iceberg Files/Climate1/Iceberg Input Files/IucnHabitatData.csv')
convcodes <- read.csv('Iceberg Files/Climate1/Iceberg Input Files/IUCN_LUH_conversion_table.csv')

iucndat %>%
  left_join(convcodes, by = c("code" = "IUCN_hab")) %>%
  mutate(name = name %>% str_replace(" ", "_")) ->
  Habitats

lapply(Species, function(a){
  
  Habitats %>% filter(name == a) %>% pull(LUH)
  
}) -> HabitatList

names(HabitatList) <- Species

landuse2017 <- brick('Iceberg Files/Climate1/Iceberg Input Files/landuse2017.grd')

# Continents ####
print("Continents!")

ContinentRaster <- raster("Iceberg Files/Climate1/Iceberg Input Files/continents-madagascar.tif") %>%
  resample(blank, method = "ngb")

ContinentWhich <- 
  lapply(1:max(values(ContinentRaster), na.rm = T), function(a) which(values(ContinentRaster)==a))
names(ContinentWhich) <- c("Africa", "Eurasia", "Greenland", "Madagascar", "NAm", "Oceania", "SAm")

# IUCN ranges for continent clipping ####
print("IUCN!")

load("~/LargeFiles/MammalStackFullMercator.Rdata")
load("Iceberg Files/Climate1/Iceberg Input Files/IUCNBuffers.Rdata")

IUCNSp <- names(MammalStackFull) %>% intersect(Species) %>% intersect(names(IUCNBuffers))
MammalStackFull <- MammalStackFull[IUCNSp]

# Dispersals ####

Dispersals <- read.csv("Iceberg Files/Climate1/Iceberg Input Files/Data for dispersal_Corrected.csv", header = T)

Dispersals <- Dispersals %>% filter(!is.na(Scientific_name), !is.na(disp50))
Dispersals$Scientific_name <- Dispersals$Scientific_name %>% str_replace(" ", "_")

ToBuffer <- intersect(Species, Dispersals$Scientific_name)

# Adding in exception for bats ####

Panth1 %>% filter(hOrder == "Chiroptera") %>% pull(Sp) ->
  BatSpecies

ToBuffer <- setdiff(ToBuffer, BatSpecies)

BatSpecies <- intersect(BatSpecies, Species)

# 01_Processing Currents ####

# Setting raster standards ####

XMin <- blank %>% extent %>% magrittr::extract(1)
XMax <- blank %>% extent %>% magrittr::extract(2)
YMin <- blank %>% extent %>% magrittr::extract(3)
YMax <- blank %>% extent %>% magrittr::extract(4)

NCol <- ncol(blank)
NRow <- nrow(blank)

i = 1  

Processed <- paste0("Iceberg Files/CleanedThresholds") %>% 
  list.files %>% str_remove(".rds$")

ToProcess <- setdiff(Species, Processed)

# ToProcess <- Species

TPs <- c(glue::glue("TP{c('01','05','10')}"))

if(length(ToProcess)>0){
  
  print("The rasters!")
  
  lapply(1:length(ToProcess), function(i){
    
    Sp <- ToProcess[i]
    
    print(Sp)
    
    SubFiles <- FullFiles[str_detect(FullFiles, Sp)] # %>% list.files(full.names = T)
    SubFiles <- SubFiles[!str_detect(SubFiles, "_ssp._")]
    
    if(length(SubFiles)>0){
      
      RasterLista <- lapply(SubFiles, raster)
      
      # 02_Resampling rasters ####
      
      RasterLista <- map(RasterLista, ~raster::resample(.x, blank, method = 'ngb'))
      
      GretCDF <- data.frame(
        
        X = seq(from = XMin, to = XMax, length.out = NCol) %>% rep(NRow),
        Y = seq(from = YMax, to = YMin, length.out = NRow) %>% rep(each = NCol)
        
      )
      
      RasterLista %>% map_dfc(values) %>% 
        rename_all(~.x %>% str_split("__") %>% map_chr(1)) %>% 
        bind_cols(GretCDF, .) ->
        GretCDF
      
      RasterLista %>% map_dfc(values) %>% 
        rename_all(~.x %>% str_split("__") %>% map_chr(1) %>% paste0("PreClipClim", .)) %>% 
        bind_cols(GretCDF, .) ->
        GretCDF
      
      # IUCN Buffer clipping ####
      
      if(Sp%in%IUCNSp){
        
        if(nrow(IUCNBuffers[[Sp]])>0){
          
          sf1 <- st_cast(IUCNBuffers[[Sp]], "MULTIPOLYGON")
          r1 <- fasterize::fasterize(sf1, blank)
          
          IUCNValues <- values(r1)
          
          if(length(IUCNValues)>0){
            
            GretCDF$IUCN <- 0
            GretCDF[which(!is.na(IUCNValues)),"IUCN"] <- 1
            
          }else{
            
            GretCDF$IUCN <- 1
            
          }
          
        } else{
          
          GretCDF$IUCN <- 1
          
        }
        
        
      } else{
        
        GretCDF$IUCN <- 1
        
      }
      
      GretCDF[GretCDF$IUCN == 0, TPs] <- 0
      
      # Continent clipping ####
      
      if(Sp%in%IUCNSp){
        
        r1 <- MammalStackFull[[Sp]]
        SpWhich <- which(!is.na(values(r1)))
        ContinentsInhabited <- unique(values(ContinentRaster)[SpWhich])
        
        if(Sp == "Canis_lupus"){
          
          ContinentsInhabited <- 
            ContinentsInhabited %>% setdiff(1)
          
        }
        
        if(length(ContinentsInhabited)>0){
          
          GretCDF$Continent <- 0
          
          GretCDF[unlist(ContinentWhich[ContinentsInhabited]),"Continent"] <- 1
          
        }else{
          
          GretCDF$Continent <- 1
          
        }
        
      } else{
        
        GretCDF$Continent <- 1
        
      }
      
      GretCDF[GretCDF$Continent == 0, TPs] <- 0
      
      # Land Use filters #####
      
      SpHabitat <- HabitatList[[Sp]] %>% 
        as.character %>% 
        str_split("[.]") %>% 
        unlist %>% unique %>%
        na.omit
      
      if(length(SpHabitat)>0){
        
        landuse2017[[SpHabitat]] %>% getValues -> 
          
          ValueDF
        
        if(length(SpHabitat)==1){
          
          Habitable <- as.numeric(ValueDF==1)
          Habitable[is.na(Habitable)] <- 0
          
        }else{
          
          Habitable <- as.numeric((ValueDF %>% rowSums(na.rm = T))>0)
          
        }
        
      }else{
        
        Habitable <- rep(1, nrow(GretCDF))
        
      }
      
      GretCDF$LandUse <- Habitable
      
      GretCDF %>% 
        dplyr::select(all_of(TPs)) %>% 
        mutate_all(~as.numeric((.x + GretCDF$LandUse))) ->
        
        GretCDF[, c(glue::glue("{TPs}LandUse") %>% as.character)]
      
      GretCDF %>% slice(-Sea) %>% 
        as.matrix %>% as("dgCMatrix") %>% 
        saveRDS(file = paste0("Iceberg Files/CleanedThresholds/", Sp, ".rds"))
      
    }
    
  })
  
  # stop()
  
  t2 <- Sys.time()
  
  print(t2 - t1)
  
}else print("None to process!")


# Importing and Pairswiselying ####

CORES <- 10

Files <- paste0("Iceberg Files/CleanedThresholds") %>% 
  list.files(full.names = T)

names(Files) <- Species <-
  Files %>% str_split("/") %>% map_chr(last) %>% str_remove("[.]rds$")

GretCDFList <- 
  
  mclapply(1:length(Species), function(a){
    
    Sp = Species[a]
    
    print(Files[a])
    
    readRDS(Files[[Sp]]) %>% as.matrix %>% 
      as.data.frame() %>% dplyr::select(contains("TP"))
    
  }, mc.preschedule = F, mc.cores = CORES)

object.size(GretCDFList)/(10^9)

names(GretCDFList) <- Species

TPVars <- colnames(GretCDFList[[1]])

i <- 7

for(i in i:length(TPVars)){
  
  print(TPVars[i])
  
  ValueList <- GretCDFList %>% 
    map(TPVars[i]) %>% 
    map(~as.numeric(.x > 1)) %>% 
    map(function(a) a*(AreaValues[-Sea]))
  
  IsNull <- ValueList %>% map_lgl(~length(.x) == 0) %>% #map_lgl(is.null) %>% 
    which()
  
  if(length(IsNull)>0){
    
    ValueList <- 
      ValueList[-IsNull]
    
  }
  
  ValueList %>% bind_cols() %>% as.data.frame() ->
    ValueDF
  
  RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)
  
  saveRDS(RangeAdj, file = paste0("Iceberg Files/ThresholdMatrices/", TPVars[i], ".rds"))
  
}

Grimport <- function(Dir, Str = "", Function = NULL){
  
  Files <- Dir %>% list.files(full.names = T, pattern = Str)
  
  if(is.null(Function)){
    
    if(str_detect(Str, ".rds")){
      
      Files %>% map(readRDS) -> FileList  
      
    }
    
    if(str_detect(Str, ".csv")){
      
      Files %>% map(read.csv) -> FileList  
      
    }
    
  }else{
    
    Files %>% map(Function) -> FileList
    
  }
  
  names(FileList) <- Files %>% str_split("/") %>% map_chr(last)
  
  return(FileList)
  
}

"Iceberg Files/ThresholdMatrices" %>% Grimport("[.]rds") -> RangeAdjList

RangeAdjList %>% map(rownames) %>% reduce(intersect) -> 
  IncludeSpecies

RangeAdjList[[1]][IncludeSpecies, IncludeSpecies] %>% lower.tri(diag = F) %>% 
  unlist %>% which() -> Extract

RangeAdjList %>% 
  map(~.x[IncludeSpecies, IncludeSpecies] %>% reshape2::melt() %>% slice(Extract)) %>% 
  bind_cols() %>% dplyr::select(1, 2, contains("value")) -> 
  
  DF

names(DF) <- c("Sp1", "Sp2", RangeAdjList %>% names %>% str_remove(".rds$"))

# DF %>% object.size() %>% divide_by(10^9)

TPVars <- DF %>% dplyr::select(-c(Sp1, Sp2)) %>% colnames

DF %>% summarise_at(TPVars, mean) %>% 
  bind_rows(DF %>% summarise_at(TPVars, Prev)) %>% 
  mutate_all(~round(.x, 3)) %>% 
  mutate(Variable = c("Mean", "Prev")) %>% 
  dplyr::select(Variable, all_of(TPVars))

DF %>% gather("Rep", "Value", TPVars) -> LongDF

LongDF %<>% mutate(SPPair = paste0(Sp1, "_", Sp2))

library(lme4)

LM1 <- lmer(data = LongDF, Value ~ (1|Rep) + (1|SPPair))

LMSummary <- LM1 %>% summary

LMSummary %>% str

LMSummary %>% saveRDS("")

LMSummary$varcor %>% data.frame %>% 
  mutate(Repeatability = (vcov/sum(vcov)) %>% round(3)) %>% 
  dplyr::select(Var = grp, Repeatability)
