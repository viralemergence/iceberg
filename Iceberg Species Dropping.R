
# Dropping species ahead of the pipeline ####

library(tidyverse); library(raster); library(parallel)

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

# Those with dispersal data ####

Dispersals <- read.csv("Iceberg Input Files/Data for dispersal.csv", header = T)
Dispersals <- Dispersals %>% filter(!is.na(Scientific_name), !is.na(disp50))
Dispersals$Scientific_name <- DispersalSp <- Dispersals$Scientific_name %>% str_replace(" ", "_")

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
  intersect(DispersalSp) %>%
  intersect(MaxEntSp) -> 
  FullMaxEntSp

FullMaxEntSp; length(FullMaxEntSp)

PhylogeneticSp %>% 
  setdiff(MarineSp) %>% 
  setdiff(NonEutherianSp)%>% 
  setdiff(NullRangeBagCommons) %>% 
  setdiff(NullRangeBagRares) %>% 
  intersect(DispersalSp) %>%
  intersect(RangeBagSp) -> 
  FullRangeBagSp

FullRangeBagSp; length(FullRangeBagSp)

OnlyRangeBags <- setdiff(FullRangeBagSp, FullMaxEntSp)

x = "Rares"

Root <- paste0("Iceberg Input Files/RangeBags/",x)

Rares <- Root %>% list.files %>% str_remove(".tif")

intersect(Rares, OnlyRangeBags)




