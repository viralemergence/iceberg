
# 00b_Assigning Dispersals to unknown generation times


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

Dispersals$Scientific_name <- Dispersals$Scientific_name %>% str_replace(" ", "_")
Dispersals <- Dispersals %>% filter(!is.na(Scientific_name))

ToBuffer <- intersect(Species, Dispersals %>% filter(!is.na(disp50)) %>% pull(Scientific_name))

ToFill <- setdiff(Species, ToBuffer)

write.csv(ToFill, file = "Iceberg Input Files/ImputedDispersals.csv")

# Assigning generation times ####

Dispersals[Dispersals$Scientific_name%in%ToFill,"Carnivore"] # All there
Dispersals[Dispersals$Scientific_name%in%ToFill,"BodyMass.Value"] # All there
Dispersals[Dispersals$Scientific_name%in%ToFill,"GenerationLength_d"] # All NA's

# Importing nearest phylogenetic neighbour ####

if(!file.exists("Iceberg Input Files/FullSTMatrix.csv")){
  
  library(geiger);library(ape);library(picante);library(dplyr)
  
  STFull <- read.nexus("data/ele_1307_sm_sa1.tre")[[1]]
  FullSTMatrix <- as.data.frame(cophenetic(STFull)) %>% as.matrix
  
  write.csv(FullSTMatrix, file = "Iceberg Input Files/FullSTMatrix.csv", row.names = F)
  
} else{FullSTMatrix <- as.matrix(read.csv("Iceberg Input Files/FullSTMatrix.csv", header = T)) }

colnames(FullSTMatrix) %>% setdiff(ToFill, .)

rownames(FullSTMatrix) <- colnames(FullSTMatrix)

ThinSTMatrix <- FullSTMatrix[,intersect(colnames(FullSTMatrix), 
                                        Dispersals %>% filter(!is.na(disp50)) %>% pull(Scientific_name))]

for(a in 1:length(ToFill)){
  
  Sp <- ToFill[a]
  
  Vector <- ThinSTMatrix[Sp,]
  
  MinVector <- which(Vector==min(Vector))
  
  if(length(MinVector)>1){
    
    MinVector %>% sample(1) %>% names ->
      Sp2
    
  }else{
    
    MinVector %>% names ->
      Sp2
    
  }
  
  FillGenerationLength <- Dispersals[Dispersals$Scientific_name == Sp2, "GenerationLength_d"] ->
    
    Dispersals[Dispersals$Scientific_name == Sp, "GenerationLength_d"]
  
}

Dispersals[Dispersals$Scientific_name%in%ToFill,"GenerationLength_d"] # it worked!

# Colin's loop ####

for(i in 1:nrow(Dispersals)){
  
  if(is.na(Dispersals$Carnivore[i])) {} else {
    
    if(Dispersals$Carnivore[i]==1) {
      
      Dispersals$disp[i] <- (40.7*(Dispersals$BodyMass.Value[i]/1000)^0.81)/(as.numeric(Dispersals$GenerationLength_d[i])/365.25)
      
    } else {
      
      Dispersals$disp[i] <- (3.31*(Dispersals$BodyMass.Value[i]/1000)^0.65)/(as.numeric(Dispersals$GenerationLength_d[i])/365.25)
    }
    
  }
}

Dispersals$disp50 <- 50*Dispersals$disp

write.csv(Dispersals,'Iceberg Input Files/Data for dispersal_Corrected.csv')

paste0("Iceberg Input Files/GretCDF/Currents/", ToFill,".rds") %>% file.remove()
paste0("Iceberg Input Files/GretCDF/Futures/", ToFill,".rds") %>% file.remove()
