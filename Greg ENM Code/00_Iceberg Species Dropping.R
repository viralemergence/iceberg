
# Dropping species ahead of the pipeline ####

# Rscript("Iceberg Species Dropping.R")

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

# Removing null range bag rares ####

Root <- paste0("Iceberg Input Files/RangeBags/","01_Raw/Currents")

Files <- list.files(Root)

Species <- Files %>% str_remove(".tif$")

RangeBagLista <- mclapply(Files, function(a){
  
  SubFiles <- paste0(Root,"/",a) %>% list.files
  
  if(length(SubFiles)>0) "Y" else "N"
  
}, mc.preschedule = F, mc.cores = 45)

NullRangeBagRares <- Species[unlist(RangeBagLista) == "N"]

# Range Bag Species ####

Root <- paste0("Iceberg Input Files/RangeBags/01_Raw/Currents")
Files <- list.files(Root)
RangeBagSp <- Files %>% str_remove(".tif$")

# MaxEnt Species ####

Root <- paste0("Iceberg Input Files/","MaxEnt","/01_Raw/Currents")
Files <- list.files(Root)
MaxEntSp <- Files %>% str_remove(".tif$")

# Combining them ####

PhylogeneticSp %>% 
  setdiff(MarineSp) %>% 
  setdiff(NonEutherianSp) %>% 
  intersect(MaxEntSp) -> 
  FullMaxEntSp

FullMaxEntSp; length(FullMaxEntSp)

PhylogeneticSp %>% 
  setdiff(MarineSp) %>% 
  setdiff(NonEutherianSp)%>% 
  setdiff(NullRangeBagRares) %>% 
  intersect(RangeBagSp) -> 
  FullRangeBagSp

FullRangeBagSp; length(FullRangeBagSp)

OnlyRangeBags <- setdiff(FullRangeBagSp, FullMaxEntSp)

SpeciesList <- list(MaxEnt = FullMaxEntSp, 
                    RangeBags = FullRangeBagSp)

stop()

# Dispersals ####

Dispersals <- read.csv("Iceberg Input Files/Data for dispersal.csv", header = T)

Dispersals <- Dispersals %>% filter(!is.na(Scientific_name), !is.na(disp50))
DispersalSp <- Dispersals$Scientific_name <- Dispersals$Scientific_name %>% str_replace(" ", "_")

Panth1 %>% mutate(DispersalKnown = as.numeric(Sp%in%DispersalSp)) %>% 
  filter(Sp%in%unlist(SpeciesList)) %>% ggregplot::ggMMplot("hOrder", "DispersalKnown")

lapply(NewEncountersList, function(a){
  
  lapply(a, function(b){
    
    AllSp <- c(b$Sp, b$Sp2)
    
    table(AllSp%in%DispersalSp)/sum(table(AllSp%in%DispersalSp))
    
  })
  
})


lapply(NewEncountersList, function(a){
  
  lapply(a, function(b){
    
    b %>% filter(!Sp%in%DispersalSp|!Sp2%in%DispersalSp) %>% nrow %>% magrittr::divide_by(nrow(b))
    
  })
  
})



