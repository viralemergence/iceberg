
# Iceberg pre-GAM spatial processing ####

# Rscript "01_Iceberg ENM Currents.R"

library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr); library(SpRanger); library(cowplot)

CORES <- 45

t1 <- Sys.time()

print("Dropping Species!")

Method = "MaxEnt"
# Method = "RangeBags"

source("00_Iceberg Species Dropping.R")

paste0("Iceberg Input Files/","MaxEnt","/01_Raw/Currents") %>% 
  list.files(full.names = T) %>% 
  append(paste0("Iceberg Input Files/","RangeBags","/01_Raw/Currents") %>% list.files(full.names = T)) ->
  FullFiles

paste0("Iceberg Input Files/","MaxEnt","/01_Raw/Currents") %>% 
  list.files() %>% str_remove(".rds$") %>%
  append(paste0("Iceberg Input Files/","RangeBags","/01_Raw/Currents") %>% list.files() %>% str_remove(".rds$")) ->
  names(FullFiles)

Species <- SpeciesList %>% unlist %>% sort()

Files <- FullFiles[Species]

PredReps <- c("Currents", paste0("Futures", 1:4))

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

# Grid areas
AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

# Land Use Data
iucndat <- read.csv('Iceberg Input Files/IucnHabitatData.csv')
convcodes <- read.csv('Iceberg Input Files/IUCN_LUH_conversion_table.csv')

iucndat %>%
  left_join(convcodes, by = c("code" = "IUCN_hab")) %>%
  mutate(name = name %>% str_replace(" ", "_")) ->
  Habitats

lapply(Species, function(a){
  
  Habitats %>% filter(name == a) %>% pull(LUH)
  
}) -> HabitatList

names(HabitatList) <- Species

landuse2017 <- brick('Iceberg Input Files/landuse2017.grd')

# Continents ####
print("Continents!")

ContinentRaster <- raster("Iceberg Input Files/continents-greenland.tif") %>%
  resample(blank, method = "ngb")

ContinentWhich <- lapply(1:6, function(a) which(values(ContinentRaster)==a))
names(ContinentWhich) <- c("Africa", "Eurasia", "Greenland", "NAm", "Oceania", "SAm")

# IUCN ranges for continent clipping ####
print("IUCN!")

load("~/LargeFiles/MammalStackFullMercator.Rdata")
load("Iceberg Input Files/IUCNBuffers.Rdata")

IUCNSp <- names(MammalStackFull) %>% intersect(Species)
MammalStackFull <- MammalStackFull[IUCNSp]

# Dispersals ####

Dispersals <- read.csv("Iceberg Input Files/Data for dispersal.csv", header = T)

Dispersals <- Dispersals %>% filter(!is.na(Scientific_name), !is.na(disp50))
Dispersals$Scientific_name <- Dispersals$Scientific_name %>% str_replace(" ", "_")

ToBuffer <- intersect(Species, Dispersals$Scientific_name)

# 01_Processing Currents ####

Root <- paste0("Iceberg Input Files/", Method,"/01_Raw/","Currents")

# Setting raster standards ####

XMin <- blank %>% extent %>% extract(1)
XMax <- blank %>% extent %>% extract(2)
YMin <- blank %>% extent %>% extract(3)
YMax <- blank %>% extent %>% extract(4)

NCol <- ncol(blank)
NRow <- nrow(blank)

i = 1  

Processed <- paste0("Iceberg Input Files/GretCDF/","Currents") %>% 
  list.files %>% str_remove(".rds$")

ToProcess <- setdiff(Species, Processed)

print("The rasters!")

mclapply(1:length(ToProcess), function(i){
  
  Sp <- ToProcess[i]
  
  print(Sp)
  
  SubFiles <- Files[[Sp]] %>% list.files(full.names = T)
  
  if(length(SubFiles)>0){
    
    RasterLista <- raster(SubFiles[1])
    
    # 02_Resampling rasters ####
    
    RasterLista <- raster::resample(RasterLista, blank, method = 'ngb')
    
    GretCDF <- data.frame(
      X = seq(from = XMin, to = XMax, length.out = NCol) %>% rep(NRow),
      Y = seq(from = YMax, to = YMin, length.out = NRow) %>% rep(each = NCol),
      Climate = as.numeric(!is.na(values(RasterLista)))
    ) 
    
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
    
    GretCDF[GretCDF$IUCN == 0, c("Climate")] <- 0
    
    # Continent clipping ####
    
    if(Sp%in%IUCNSp){
      
      r1 <- MammalStackFull[[Sp]]
      SpWhich <- which(!is.na(values(r1)))
      ContinentsInhabited <- unique(values(ContinentRaster)[SpWhich])
      
      if(length(ContinentsInhabited)>0){
        
        GretCDF$Continent <- 0
        
        GretCDF[unlist(ContinentWhich[ContinentsInhabited]),"Continent"] <- 1
        
      }else{
        
        GretCDF$Continent <- 1
        
      }
      
    } else{
      
      GretCDF$Continent <- 1
      
    }
    
    GretCDF[GretCDF$Continent == 0, c("Climate")] <- 0
    
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
    
    GretCDF$ClimateLandUse <- as.numeric(rowSums(GretCDF[,c("Climate","LandUse")])>1)
    
    # Dispersals ####
    
    if(Sp%in%ToBuffer){
      
      Dist <- Dispersals %>% filter(Scientific_name == Sp) %>% pull(disp50)*1000
      
      for(j in c("Climate","ClimateLandUse")){
        
        if(j == "ClimateLandUse"&all(GretCDF$Climate==GretCDF$ClimateLandUse)){
          
          GretCDF[,"BufferClimateLandUse"] <- GretCDF[,"BufferClimate"]
          
        }else{
          
          r1 <- blank
          
          values(r1) <- c(1)[ifelse(GretCDF[,j]==1,1,NA)]
          r2 <- buffer2(r1, Dist)
          r3 <- resample(r2, blank, method = "ngb")
          
          GretCDF[,paste0("Buffer",j)] <- values(r3)
          GretCDF[,paste0("Buffer",j)][is.na(GretCDF[,paste0("Buffer",j)])] <- 0
          
        }
        
      }
      
    } else {
      
      GretCDF[,paste0("Buffer",c("Climate","ClimateLandUse"))] <- 1
      
    }
    
    GretCDF %>% slice(-Sea) %>% 
      as.matrix %>% as("dgCMatrix") %>% 
      saveRDS(file = paste0("Iceberg Input Files/GretCDF/Currents/", Sp, ".rds"))
    
  }
  
}, mc.preschedule = F, mc.cores = CORES)

t2 <- Sys.time()

print(t2 - t1)


# Troubleshooting ####

iucndat <- read.csv('Iceberg Input Files/IucnHabitatData.csv')
convcodes <- read.csv('Iceberg Input Files/IUCN_LUH_conversion_table.csv')
iucndat %>%
  left_join(convcodes, by = c("code" = "IUCN_hab")) %>%
  mutate(name = name %>% str_replace(" ", "_")) ->
  Habitats

list.files(paste0("~/Albersnet/Iceberg Input Files/GretCDF/Currents")) %>% str_remove(".rds") -> 
  Species

lapply(Species, function(a){
  Habitats %>% filter(name == a) %>% pull(LUH)
}) -> HabitatList
names(HabitatList) <- Species
HabitatList[sapply(HabitatList, length)>0]

Sp <- "Microtus_daghestanicus"
Sp <- "Acinonyx_jubatus"
Sp <- "Lycaon_pictus"

CDF <- readRDS(paste0("~/Albersnet/Iceberg Input Files/GretCDF/Currents/",Sp,".rds")) %>%
  as.matrix %>% as.data.frame()

CDF %>% colSums()

CDF %>% gather(key = "Key", value = "Value", -X, -Y) %>% 
  mutate(Key = factor(Key, levels = unique(Key))) %>%
  ggplot(aes(X, Y, fill = Value)) + 
  coord_fixed() + theme_cowplot() +
  geom_tile() + facet_wrap(~Key)

