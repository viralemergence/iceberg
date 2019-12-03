# Graskby Terminal Code ####
# Rscript "FastGraskBy.R"

# Dropping species ahead of the pipeline ####
library(sf); library(tidyverse); library(raster); library(parallel)

Method = "MaxEnt"
x = "Currents"

blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

# Removing Marine Hosts ####

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

Panth1 %>% filter(hOrder%in%c("Cetacea", "Sirenia")) %>% pull(Sp) ->
  MarineSp

# Land use ####

iucndat <- read.csv('Iceberg Input Files/IucnHabitatData.csv')
convcodes <- read.csv('Iceberg Input Files/IUCN_LUH_conversion_table.csv')

iucndat %>% left_join(convcodes, by = c("code" = "IUCN_hab")) %>% group_by(name) %>%
  summarise(Habitats = paste(unique(LUH), collapse = ".")) %>%
  mutate(name = name %>% str_replace(" ", "_")) ->
  AllHabitats

AllHabitats %>% 
  filter(Habitats == "MARINE"| 
           Habitats == "MARINE.NA") %>% pull(name) -> MarineSp2

MarineSp <- c(MarineSp, MarineSp2)

# Removing Non-Eutherian Hosts ####

NonEutherians <- c("Diprotodontia",
                   "Dasyuromorphia",
                   "Paucituberculata",
                   "Didelphimorphia",
                   "Microbiotheria",
                   "Peramelemorphia", 
                   "Notoryctemorphia",
                   "Monotremata")

Panth1 %>% filter(hOrder%in%NonEutherians) %>% pull(Sp) ->
  NonEutherianSp

# Keeping Hosts with Phylogenetic Data ####

if(!file.exists("Iceberg Input Files/FullSTMatrix.csv")){
  
  library(geiger);library(ape);library(picante);library(dplyr)
  
  STFull <- read.nexus("data/ele_1307_sm_sa1.tre")[[1]]
  FullSTMatrix <- as.data.frame(cophenetic(STFull)) %>% as.matrix
  
  write.csv(FullSTMatrix, file = "Iceberg Input Files/FullSTMatrix.csv", row.names = F)
  
} else{FullSTMatrix <- as.matrix(read.csv("Iceberg Input Files/FullSTMatrix.csv", header = T)) }

PhylogeneticSp <- colnames(FullSTMatrix)

# Removing null range bag commons ####

Root <- paste0("Iceberg Input Files/RangeBags/","Commons")

Files <- list.files(Root)

Species <- Files %>% str_remove(".tif")

RangeBagLista <- mclapply(Files, function(a){
  
  SubFiles <- paste0(Root,"/",a) %>% list.files
  
  if(length(SubFiles>0)) raster(paste0(Root,"/",a,"/",SubFiles[1]))
  
})

NullRangeBagCommons <- Species[sapply(RangeBagLista, is.null)]

# Removing null range bag rares ####

Root <- paste0("Iceberg Input Files/RangeBags/","Rares")

Files <- list.files(Root)

Species <- Files %>% str_remove(".tif")

RangeBagLista <- mclapply(Files, function(a){
  
  SubFiles <- paste0(Root,"/",a) %>% list.files
  
  if(length(SubFiles>0)) raster(paste0(Root,"/",a,"/",SubFiles[1]))
  
})

NullRangeBagRares <- Species[sapply(RangeBagLista, is.null)]

# Range Bag Species ####

Root <- paste0("Iceberg Input Files/Resampled RangeBags")
Files <- list.files(Root)
RangeBagSp <- Files %>% str_remove(".tif")

# MaxEnt Species ####

Root <- paste0("Iceberg Input Files/",Method,"/Raw/Currents")
Files <- list.files(Root)
MaxEntSp <- Files %>% str_remove(".tif")

# Combining them ####

PhylogeneticSp %>% 
  setdiff(MarineSp) %>% 
  setdiff(NonEutherianSp) %>% 
  intersect(MaxEntSp) -> 
  FullMaxEntSp

# Resampling rasters ####

Root <- paste0("Iceberg Input Files/", Method,"/Raw/",x)

Files <- list.files(Root)

Species1 <- Files %>% str_remove(".tif")

RangeBagLista <- lapply(Files, function(a){
  
  SubFiles <- paste0(Root,"/",a) %>% list.files
  
  if(length(SubFiles>0)) raster(paste0(Root,"/",a,"/",SubFiles[1]))
  
})

names(RangeBagLista) <- Species1

RangeBagLista <- RangeBagLista[!sapply(RangeBagLista, is.null)]

Species <- intersect(names(RangeBagLista), Species1)
SpeciesLost <- setdiff(Species1, names(RangeBagLista))
Species <- intersect(Species, FullMaxEntSp)

# Graskby ####

Mammal_Shapes <- st_read("~/ShapeFiles")

Mammal_Shapes$Binomial = str_replace(Mammal_Shapes$binomial, " ", "_")
Mammal_Shapes <- Mammal_Shapes[order(Mammal_Shapes$Binomial),]

Mammal_Shapes <- Mammal_Shapes[Mammal_Shapes$Binomial%in%Species,]

buffer2 <- function(r, dist) {
  
  if("raster"%in%class(r)|"RasterLayer"%in%class(r)){
    projR <- projectRaster(r, crs=CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    projRb <- raster::buffer(projR, dist)
    projRb <- projectRaster(projRb, crs=CRS("+proj=longlat +datum=WGS84"))
    projRb[!is.na(projRb)] <- 1
  }
  
  if("sf"%in%class(r)){
    projR <- st_transform(r, crs=CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
    projRb <- st_buffer(projR, dist)
    projRb <- st_transform(projRb, crs=CRS("+proj=longlat +datum=WGS84"))
  }
  
  return(projRb)
  
}

GraskBy <- function(ShapeFile, Column, Distance = 1000){
  
  Sp <- ShapeFile[[Column]] %>% unique %>% sort
  
  BufferedList <- mclapply(Sp, function(a){
    
    print(a)
    
    ShapeFile[ShapeFile[[Column]] == a,] %>% buffer2(dist = Distance)
    
  }, mc.preschedule = F)
  
  names(BufferedList) <- Sp
  
  return(BufferedList)
  
}

IUCNBuffers <- GraskBy(Mammal_Shapes, "Binomial", Distance = 10^6)

save(IUCNBuffers, file = paste0("Iceberg Input Files/IUCNBuffers",Method,".Rdata"))

length(IUCNBuffers)

Species = names(IUCNBuffers)

IUCNBufferRasters <- list()

for(i in 1:length(Species)){
  
  Sp = Species[i]
  print(Sp)
  
  sf1 <- st_cast(IUCNBuffers[[Sp]], "MULTIPOLYGON")
  
  r1 <- fasterize::fasterize(sf1, blank)
  IUCNBufferRasters[[Sp]] <- r1
  
}

save(IUCNBufferRasters, file = "Iceberg Input Files/IUCNBufferRasters.Rdata")





