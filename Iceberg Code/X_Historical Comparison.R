
# Rscript "X_Historical Comparison.R"

# 1_Processing rasters ####

library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr); library(SpRanger); library(cowplot)

t1 <- Sys.time()

#print("Dropping Species!")

#source("Iceberg Code/Iceberg Greg ENM Code/00_Iceberg Species Dropping.R")

print("Doing Currents!")

paste0("Iceberg Files/","Historical/HistoricalRasters") %>% 
  list.files(full.names = T) -> 
  FullFiles

paste0("Iceberg Files/","Historical/HistoricalRasters") %>% 
  list.files() %>% str_remove(".tif$") %>% 
  str_split("__") %>% map_chr(2) %>% 
  str_remove_all("_mean_noBias_1e.06_0_1e.06_0_all_all_none_all_maxnet_none_equalWeights$") ->
  names(FullFiles)

Species <- 
  names(FullFiles) %>% 
  intersect("Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents" %>% 
              list.files %>% str_remove(".rds$"))

Files <- FullFiles[Species]

PredReps <- c("Currents", paste0("Futures", 1:4))

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

IUCNSp <- names(MammalStackFull) %>% intersect(Species)
MammalStackFull <- MammalStackFull[IUCNSp]

# Dispersals ####

Dispersals <- read.csv("Iceberg Files/Climate1/Iceberg Input Files/Data for dispersal_Corrected.csv", header = T)

Dispersals <- Dispersals %>% filter(!is.na(Scientific_name), !is.na(disp50))
Dispersals$Scientific_name <- Dispersals$Scientific_name %>% str_replace(" ", "_")

ToBuffer <- intersect(Species, Dispersals$Scientific_name)

# Adding in exception for bats ####

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

Panth1 %>% filter(hOrder%in%c("Cetacea", "Sirenia")|
                    hFamily%in%c("Phocidae", "Odobenidae", "Otariidae")) %>% pull(Sp) ->
  MarineSp


Panth1 %>% filter(hOrder == "Chiroptera") %>% pull(Sp) ->
  BatSpecies

ToBuffer <- setdiff(ToBuffer, BatSpecies)

BatSpecies <- intersect(BatSpecies, Species)

# Setting raster standards ####

XMin <- blank %>% extent %>% magrittr::extract(1)
XMax <- blank %>% extent %>% magrittr::extract(2)
YMin <- blank %>% extent %>% magrittr::extract(3)
YMax <- blank %>% extent %>% magrittr::extract(4)

NCol <- ncol(blank)
NRow <- nrow(blank)

i = 1  

Processed <- paste0("Iceberg Files/Historical/GretCDF") %>% 
  list.files %>% str_remove(".rds$")

ToProcess <- setdiff(Species, Processed)

# ToProcess <- Species

if(length(ToProcess)>0){
  
  print("The rasters!")
  
  mclapply(1:length(ToProcess), function(i){
    
    Sp <- ToProcess[i]
    
    print(Sp)
    
    SubFiles <- Files[[Sp]]# %>% list.files(full.names = T)
    
    if(length(SubFiles)>0){
      
      RasterLista <- raster(SubFiles[1])
      
      # 02_Resampling rasters ####
      
      RasterLista <- raster::resample(RasterLista, blank, method = 'ngb')
      
      GretCDF <- data.frame(
        X = seq(from = XMin, to = XMax, length.out = NCol) %>% rep(NRow),
        Y = seq(from = YMax, to = YMin, length.out = NRow) %>% rep(each = NCol),
        
        Climate = as.numeric(!is.na(values(RasterLista)))
        
      ) %>% mutate(PreClipClim = Climate)
      # 
      # # IUCN Buffer clipping ####
      # 
      # if(Sp%in%IUCNSp){
      #   
      #   if(nrow(IUCNBuffers[[Sp]])>0){
      #     
      #     sf1 <- st_cast(IUCNBuffers[[Sp]], "MULTIPOLYGON")
      #     r1 <- fasterize::fasterize(sf1, blank)
      #     
      #     IUCNValues <- values(r1)
      #     
      #     if(length(IUCNValues)>0){
      #       
      #       GretCDF$IUCN <- 0
      #       GretCDF[which(!is.na(IUCNValues)),"IUCN"] <- 1
      #       
      #     }else{
      #       
      #       GretCDF$IUCN <- 1
      #       
      #     }
      #     
      #   } else{
      #     
      #     GretCDF$IUCN <- 1
      #     
      #   }
      #   
      #   
      # } else{
      #   
      #   GretCDF$IUCN <- 1
      #   
      # }
      # 
      # GretCDF[GretCDF$IUCN == 0, c("Climate")] <- 0
      # 
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
      
      GretCDF[GretCDF$Continent == 0, c("Climate")] <- 0
      # 
      # # Land Use filters #####
      # 
      # SpHabitat <- HabitatList[[Sp]] %>% 
      #   as.character %>% 
      #   str_split("[.]") %>% 
      #   unlist %>% unique %>%
      #   na.omit
      # 
      # if(length(SpHabitat)>0){
      #   
      #   landuse2017[[SpHabitat]] %>% getValues -> 
      #     
      #     ValueDF
      #   
      #   if(length(SpHabitat)==1){
      #     
      #     Habitable <- as.numeric(ValueDF==1)
      #     Habitable[is.na(Habitable)] <- 0
      #     
      #   }else{
      #     
      #     Habitable <- as.numeric((ValueDF %>% rowSums(na.rm = T))>0)
      #     
      #   }
      #   
      # }else{
      #   
      #   Habitable <- rep(1, nrow(GretCDF))
      #   
      # }
      # 
      # GretCDF$LandUse <- Habitable
      # 
      # GretCDF$ClimateLandUse <- as.numeric(rowSums(GretCDF[,c("Climate","LandUse")])>1)
      # 
      # # Dispersals ####
      
      if(Sp%in%ToBuffer){
        
        Dist <- Dispersals %>% filter(Scientific_name == Sp) %>% pull(disp50)*1000
        
        j = "Climate"
        
        # for(j in c("Climate","ClimateLandUse")){
        
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
        
        # }
        
      } else {
        
        GretCDF[,paste0("Buffer",c("Climate", "ClimateLandUse"))] <- 1
        
      }
      
      GretCDF %>% slice(-Sea) %>% 
        as.matrix %>% as("dgCMatrix") %>% 
        saveRDS(file = paste0("Iceberg Files/Historical/GretCDF/", Sp, ".rds"))
      
    }
    
  }, mc.preschedule = F, mc.cores = 10)
  
  # stop()
  
  t2 <- Sys.time()
  
  print(t2 - t1)
  
}else print("None to process!")

# 2_Making Pairwise ####

paste0("~/Albersnet/Iceberg Files/", 
       "Historical/GretCDF") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("~/Albersnet/Iceberg Files/", 
       "Historical/GretCDF") %>% 
  list.files(full.names = T) ->
  CurrentFiles

Species ->
  names(CurrentFiles)

CurrentCDFList <- list()

Species <- Species %>% sort %>% intersect(names(CurrentFiles)) #%>% intersect(names(FutureFiles))

CurrentFiles <- CurrentFiles[Species]

CurrentCDFList <- 
  
  mclapply(1:length(Species), function(a){
    
    Sp = Species[a]
    
    print(CurrentFiles[a])
    
    readRDS(CurrentFiles[[Sp]]) %>% as.matrix %>% 
      as.data.frame() %>% dplyr::select(Climate)
    
  }, mc.preschedule = F, mc.cores = CORES)

object.size(CurrentCDFList)/(10^9)

names(CurrentCDFList) <- Species

CurrentCDFList %>% map("Climate") %>% 
  map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
  ValueDF

print("Calculating overlap!")

RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)

saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", "HistoricalRangeAdj.rds"))

# Importing the non-IUCN Clip ones ####

paste0("~/Albersnet/Iceberg Files/", 
       "Historical/GretCDF") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("~/Albersnet/Iceberg Files/", 
       "Historical/GretCDF") %>% 
  list.files(full.names = T) ->
  CurrentFiles

Species ->
  names(CurrentFiles)

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/IUCNCheck") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/IUCNCheck")  %>% 
  list.files(full.names = T) ->
  CurrentFiles2

Species ->
  names(CurrentFiles2)

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/Currents") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

paste0("~/Albersnet/Iceberg Files/",
       "Climate1/Iceberg Input Files/",
       "GretCDF/Currents")  %>% 
  list.files(full.names = T) ->
  CurrentFiles3

Species ->
  names(CurrentFiles3)

CurrentCDFList <- list()

Species <- Species %>% sort %>% 
  intersect(names(CurrentFiles))  %>% 
  intersect(names(CurrentFiles2)) %>% 
  intersect(names(CurrentFiles3))

CurrentFiles2 <- CurrentFiles2[Species]
CurrentFiles3 <- CurrentFiles3[Species]

HistoricalBufferClip <- T

CurrentCDFList <- 
  
  mclapply(1:length(Species), function(a){
    
    Sp = Species[a]
    
    print(CurrentFiles2[a])
    
    DF <- readRDS(CurrentFiles2[[Sp]]) %>% as.matrix %>% 
      as.data.frame() %>% dplyr::select(PreClipClim) %>% 
      mutate(BufferPreClipClim = PreClipClim)
    
    # Clip by buffer? ####
    
    # if(HistoricalBufferClip){
    
    DF2 <- readRDS(CurrentFiles[[Sp]]) %>% as.matrix %>% 
      as.data.frame() %>% pull(BufferClimate)
    
    DF$BufferPreClipClim[!DF2==1] <- 0
    
    DF$IUCNClippedClimate <- 
      DF$DispersalIUCNClippedClimate <- 
      readRDS(CurrentFiles3[[Sp]]) %>% as.matrix %>% 
      as.data.frame() %>% pull(Climate)
    
    DF$DispersalIUCNClippedClimate[!DF2==1] <- 0
    
    
    # }
    
    return(DF)
    
  }, mc.preschedule = F, mc.cores = CORES)

object.size(CurrentCDFList)/(10^9)

names(CurrentCDFList) <- Species

CurrentCDFList %>% map("IUCNClippedClimate") %>% 
  map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
  ValueDF

print("Calculating overlap!")

RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)

saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", 
                                "IUCNClippedCurrentsRangeAdj.rds"))

CurrentCDFList %>% map("BufferPreClipClim") %>% 
  map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
  ValueDF

print("Calculating overlap!")

RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)

saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", 
                                "PreClipRangeAdjDispersalClipped.rds"))

CurrentCDFList %>% map("PreClipClim") %>% 
  map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
  ValueDF

print("Calculating overlap!")

RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)

saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", 
                                "PreClipRangeAdj.rds"))

CurrentCDFList %>% map("DispersalIUCNClippedClimate") %>% 
  map(function(a) a*(AreaValues[-Sea])) %>% bind_cols() %>% as.data.frame() ->
  ValueDF

print("Calculating overlap!")

RangeAdj <- PairsWisely(Rasterstack = ValueDF, Area = T)

saveRDS(RangeAdj, file = paste0("Iceberg Files/Historical/", 
                                "DispersalIUCNClippedClimateRangeAdj.rds"))

"Iceberg Files/Historical" %>% list.files(pattern = "RangeAdj", full.names = T) %>% 
  map(readRDS) -> RangeAdjList

names(RangeAdjList) <- 
  "Iceberg Files/Historical" %>% list.files(pattern = "RangeAdj")

RangeAdjList %>% map(c(unlist, Prev))
RangeAdjList %>% map(c(unlist, mean))

# 3_Making data frame ####
# Final Iceberg Code/Rscript "2_Iceberg Data Import.R" ####

Unlist1 <- function(x) unlist(x, recursive = F)

# Viral data import ###

library(igraph); library(magrittr); library(dplyr); library(ggplot2); require(RCurl); library(readr);
library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(cowplot)

HistoricalRangeAdj <- 
  
  readRDS(paste0("Iceberg Files/Historical/", 
                 "HistoricalRangeAdj.rds"))
PreClipRangeAdj <- 
  
  readRDS(paste0("Iceberg Files/Historical/", 
                 "PreClipRangeAdj.rds"))

AssocsBase <- read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/associations.csv") %>% data.frame()
HostTraits <- read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/hosts.csv") %>% data.frame()
VirusTraits <- read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/viruses.csv") %>% data.frame()

names(AssocsBase)[1:2] <- c("Virus", "Host")
AssocsBase <- mutate(AssocsBase, Virus = as.factor(Virus), Host = as.factor(Host))

AssocsBase2 <- AssocsBase
AssocsBase2 <- droplevels(AssocsBase[!AssocsBase$Host == "Homo_sapiens"&
                                       !AssocsBase$Virus == "Rabies_virus",])

# Making bipartite projections ####

AssocsTraits <- AssocsBase2[,1:2]

m <- table(AssocsTraits)
M <- as.matrix(m)

bipgraph <- graph.incidence(M, weighted = NULL)

Hostgraph <- bipartite.projection(bipgraph)$proj2

HostAdj <- as.matrix(get.adjacency(Hostgraph, attr = "weight"))
diag(HostAdj) <- table(AssocsBase2$Host)
Remove <- which(rowSums(HostAdj)==diag(HostAdj))
HostAdj <- HostAdj[-Remove,-Remove]

# Deriving metrics from the networks ####

Hosts <- data.frame(Sp = names(V(Hostgraph)),
                    Degree = degree(Hostgraph),
                    Eigenvector = eigen_centrality(Hostgraph)$vector,
                    Kcore = coreness(Hostgraph),
                    Between = betweenness(Hostgraph))

Hosts <- merge(Hosts, HostTraits, by.x = "Sp", by.y = "hHostNameFinal", all.x = T)
Hosts <- Hosts %>% dplyr::rename(hDom = hWildDomFAO)

Domestics <- Hosts[Hosts$hDom == "domestic", "Sp"]
Wildlife <- Hosts[Hosts$hDom == "wild", "Sp"]

AssocsTraits <- merge(AssocsTraits, HostTraits, by.x = "Host", by.y = "hHostNameFinal", all.x = T)

AssocsTraits$Domestic <- ifelse(AssocsTraits$Host%in%Domestics,1,0)
AssocsTraits$Wildlife <- ifelse(AssocsTraits$Host%in%Wildlife,1,0)

ZoonoticViruses <- AssocsBase %>% filter(Host == "Homo_sapiens") %>% dplyr::select(Virus)

Hosts <- Hosts %>% 
  mutate(
    Domestic = ifelse(Sp %in% Domestics, 1, 0),
    Wildlife = ifelse(Sp %in% Wildlife, 1, 0),
    hZoonosisCount = c(table(AssocsTraits[AssocsTraits$Virus%in%ZoonoticViruses$Virus,"Host"])),
    Records = c(table(AssocsTraits$Host))
  )

#devtools::install_github("gfalbery/ggregplot")
library(ggregplot); library(ggplot2); library(RColorBrewer)

ParasitePalettes<-c("PuRd","PuBu","BuGn","Purples","Oranges")
ParasiteColours<-c("#DD1c77","#2B8CBE","#2CA25F",brewer.pal(5,"Purples")[4],brewer.pal(5,"Oranges")[4])

AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
AlberColours[length(AlberColours)+1:2] <- RColorBrewer::brewer.pal(11, AlberPalettes[[4]])[c(2,10)]

theme_set(theme_cowplot())

# Adding in all mammal supertree ####

if(!file.exists("Iceberg Files/Climate1/Iceberg Input Files/FullSTMatrix.csv")){
  
  library(geiger);library(ape);library(picante);library(dplyr)
  
  STFull <- read.nexus("data/ele_1307_sm_sa1.tre")[[1]]
  FullSTMatrix <- as.data.frame(cophenetic(STFull)) %>% as.matrix
  
  write.csv(FullSTMatrix, file = "Iceberg Files/Climate1/Iceberg Input Files/FullSTMatrix.csv", row.names = F)
  
} else{FullSTMatrix <- as.matrix(read.csv("Iceberg Files/Climate1/Iceberg Input Files/FullSTMatrix.csv", header = T)) }

# Making Viral Associations and Polygons ####

VirusAssocs <- apply(M, 1, function(a) names(a[a>0]))

# Creating final dataset

# Replacing absent names in the full ST matrix ####

NonEutherians <- c("Diprotodontia",
                   "Dasyuromorphia",
                   "Paucituberculata",
                   "Didelphimorphia",
                   "Microbiotheria",
                   "Peramelemorphia", 
                   "Notoryctemorphia",
                   "Monotremata")

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

NonEutherianSp <- Panth1[Panth1$hOrder%in%NonEutherians,"Sp"]

FinalHostNames <- reduce(list(
  rownames(HistoricalRangeAdj), 
  rownames(PreClipRangeAdj),
  # rownames(HostAdj), 
  colnames(FullSTMatrix)), 
  intersect)

FHN <- FinalHostNames %>% setdiff(NonEutherianSp); length(FHN)

AllMammals <- intersect(colnames(FullSTMatrix),colnames(HistoricalRangeAdj))
AllMammals <- AllMammals[order(AllMammals)]
AbsentHosts <- FHN[which(!FHN%in%AllMammals)]

NameReplace <- c(
  "Micaelamys_namaquensis",
  "Akodon_paranaensis",
  "Bos_frontalis",
  "Bos_grunniens",
  "Bubalus_arnee", # Absent
  "Capra_hircus",
  "Hexaprotodon_liberiensis",
  "Equus_burchellii",
  "Oryzomys_alfaroi" ,
  "Oryzomys_laticeps",
  "Oryzomys_megacephalus",
  "Callithrix_argentata",
  "Miniopterus_schreibersii",
  "Myotis_ricketti",
  "Oryzomys_albigularis",
  "Ovis_aries",
  "Piliocolobus_badius",
  "Piliocolobus_rufomitratus" ,
  "Lycalopex_gymnocercus" ,
  "Rhinolophus_hildebrandtii",
  "Oryzomys_angouya",
  "Mops_condylurus",
  "Chaerephon_plicatus",
  "Chaerephon_pumilus",
  "Taurotragus_oryx")

names(NameReplace) <- AbsentHosts

rownames(FullSTMatrix) <- colnames(FullSTMatrix) <- sapply(colnames(FullSTMatrix), function(a) ifelse(a%in%AbsentHosts, NameReplace[a], a))

NonEutherians <- c("Diprotodontia",
                   "Dasyuromorphia",
                   "Paucituberculata",
                   "Didelphimorphia",
                   "Microbiotheria",
                   "Peramelemorphia", 
                   "Notoryctemorphia",
                   "Monotremata")

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

NonEutherianSp <- Panth1[Panth1$hOrder%in%NonEutherians,"Sp"]

tFullSTMatrix <- 1 - 
  (FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp] - 
     min(FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp]))/
  max(FullSTMatrix[!rownames(FullSTMatrix)%in%NonEutherianSp,!rownames(FullSTMatrix)%in%NonEutherianSp])

tSTMatrix <- tFullSTMatrix

# Going Ahead ####

rownames(Hosts) = Hosts$Sp

FinalHostNames <- reduce(list(
  rownames(HistoricalRangeAdj), 
  rownames(PreClipRangeAdj),
  # rownames(HostAdj), 
  colnames(FullSTMatrix)), 
  intersect)

FinalHostNames %>% setdiff(NonEutherianSp)

FHN <- FinalHostNames; length(FHN)

# UpperHosts <- # Removing diagonals, as they're uninformative
#   which(upper.tri(HostAdj[FHN,FHN], diag = T))

HostMatrixdf <- data.frame(#Virus = c(HostAdj[FHN, FHN]),
  Historical = c(HistoricalRangeAdj[FHN, FHN]),
  PreClipClim = c(PreClipRangeAdj[FHN, FHN]),
  #Phylo = c(tFullSTMatrix[FHN, FHN]),
  Sp = as.character(rep(FHN, each = length(FHN))),
  Sp2 = as.character(rep(FHN, length(FHN)))
)

HostMatrixdf$Sp <- as.character(HostMatrixdf$Sp)
HostMatrixdf$Sp2 <- as.character(HostMatrixdf$Sp2)

HostMatrixVar <- c("hOrder", "hFamily", "hDom", "hAllZACites", "hDiseaseZACites"
                   #"LongMean", "LatMean")
)

HostMatrixdf[,HostMatrixVar] <- Hosts[HostMatrixdf$Sp, HostMatrixVar]
HostMatrixdf[,paste0(HostMatrixVar,".Sp2")] <- Hosts[HostMatrixdf$Sp2, HostMatrixVar]
HostMatrixdf[HostMatrixdf$Sp == "Lynx_lynx",] <- HostMatrixdf[HostMatrixdf$Sp == "Lynx_lynx",] %>% mutate(hAllZACites = 1167, hDiseaseZACites = 115)

HostMatrixdf <- HostMatrixdf %>% mutate(
  
  hOrder = Hosts[HostMatrixdf$Sp,"hOrder"],
  hFamily = Hosts[HostMatrixdf$Sp,"hFamily"],
  hDom = Hosts[HostMatrixdf$Sp,"hDom"]
  
)

HostMatrixdf$Space0 <- ifelse(HostMatrixdf$Space == 0, "No Overlap", "Overlap")
HostMatrixdf$Cites <- log(HostMatrixdf$hAllZACites + 1)
HostMatrixdf$TotalCites <- log(HostMatrixdf$hAllZACites + HostMatrixdf$hAllZACites.Sp2 + 1)
HostMatrixdf$MinCites <- apply(HostMatrixdf[,c("hAllZACites", "hAllZACites.Sp2")],1, function(a) min(a, na.rm = T))

HostMatrixdf$DCites <- log(HostMatrixdf$hDiseaseZACites + 1)
HostMatrixdf$MinDCites <- apply(HostMatrixdf[,c("hDiseaseZACites", "hDiseaseZACites.Sp2")],1, function(a) min(a, na.rm = T))
HostMatrixdf$TotalDCites <- log(HostMatrixdf$hDiseaseZACites + HostMatrixdf$hAllZACites.Sp2 + 1)

HostMatrixdf$DomDom <- paste(HostMatrixdf$hDom, HostMatrixdf$hDom.Sp2)
HostMatrixdf$DomDom <- ifelse(HostMatrixdf$DomDom == "domestic wild", "wild domestic", HostMatrixdf$DomDom) %>%
  factor(levels = c("wild wild", "domestic domestic", "wild domestic"))

UpperHosts <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

FinalHostMatrix <- HostMatrixdf[-UpperHosts,]

FinalHostMatrix$Phylo <- FinalHostMatrix$Phylo
FinalHostMatrix$MinDCites <- log(FinalHostMatrix$MinDCites + 1)
FinalHostMatrix$VirusBinary <- ifelse(FinalHostMatrix$Virus>0, 1, 0)

Remove1 <- FinalHostMatrix %>% group_by(Sp) %>% dplyr::summarise(Mean = mean(VirusBinary)) %>% slice(order(Mean)) %>% filter(Mean==0) %>% dplyr::select(Sp)
Remove2 <- FinalHostMatrix %>% group_by(Sp2) %>% dplyr::summarise(Mean = mean(VirusBinary)) %>% slice(order(Mean)) %>% filter(Mean==0) %>% dplyr::select(Sp2)

Remove3 <- which(table(c((FinalHostMatrix %>% filter(Phylo < 0.25) %>% dplyr::select(Sp, Sp2))$Sp %>% as.character(),
                         (FinalHostMatrix %>% filter(Phylo < 0.25) %>% dplyr::select(Sp, Sp2))$Sp2 %>% as.character()))>20) %>% 
  names

RemoveSp <- intersect(Remove1$Sp, Remove2$Sp2)

FinalHostMatrix <- FinalHostMatrix %>% filter(!Sp%in%RemoveSp&!Sp2%in%RemoveSp)

FinalHostMatrix$Sp <- factor(FinalHostMatrix$Sp, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
FinalHostMatrix$Sp2 <- factor(FinalHostMatrix$Sp2, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))

FinalHostMatrix <- FinalHostMatrix %>% slice(order(Sp,Sp2))

# Let's try this a second time ####

FHN <- levels(FinalHostMatrix$Sp)

HostMatrixdf <- data.frame(Virus = c(HostAdj[FHN, FHN]),
                           Space = c(FullRangeAdj[FHN, FHN]),
                           Phylo = c(tFullSTMatrix[FHN, FHN]),
                           Sp = as.character(rep(FHN, each = length(FHN))),
                           Sp2 = as.character(rep(FHN, length(FHN)))
)

HostMatrixdf$Sp <- as.character(HostMatrixdf$Sp)
HostMatrixdf$Sp2 <- as.character(HostMatrixdf$Sp2)

HostMatrixVar <- c("hOrder", "hFamily", "hDom", "hAllZACites", "hDiseaseZACites")

HostMatrixdf[,HostMatrixVar] <- Hosts[HostMatrixdf$Sp, HostMatrixVar]
HostMatrixdf[,paste0(HostMatrixVar,".Sp2")] <- Hosts[HostMatrixdf$Sp2, HostMatrixVar]

HostMatrixdf$Cites <- log(HostMatrixdf$hAllZACites + 1)
HostMatrixdf$TotalCites <- log(HostMatrixdf$hAllZACites + HostMatrixdf$hAllZACites.Sp2 + 1)
HostMatrixdf$MinCites <- apply(HostMatrixdf[,c("hAllZACites", "hAllZACites.Sp2")],1, function(a) min(a, na.rm = T))

HostMatrixdf$DCites <- log(HostMatrixdf$hDiseaseZACites + 1)
HostMatrixdf$MinDCites <- apply(HostMatrixdf[,c("hDiseaseZACites", "hDiseaseZACites.Sp2")],1, function(a) min(a, na.rm = T))
HostMatrixdf$TotalDCites <- log(HostMatrixdf$hDiseaseZACites + HostMatrixdf$hAllZACites.Sp2 + 1)

HostMatrixdf$DomDom <- paste(HostMatrixdf$hDom, HostMatrixdf$hDom.Sp2)
HostMatrixdf$DomDom <- ifelse(HostMatrixdf$DomDom == "domestic wild", "wild domestic", HostMatrixdf$DomDom) %>%
  factor(levels = c("wild wild", "domestic domestic", "wild domestic"))

UpperHosts <- which(upper.tri(HostAdj[FHN,FHN], diag = T))

FinalHostMatrix <- HostMatrixdf[-UpperHosts,]

FinalHostMatrix$MinDCites <- log(FinalHostMatrix$MinDCites + 1)
FinalHostMatrix$VirusBinary <- ifelse(FinalHostMatrix$Virus>0, 1, 0)

FinalHostMatrix$Gz <- as.numeric(FinalHostMatrix$Space==0)

FinalHostMatrix$Sp <- factor(FinalHostMatrix$Sp, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
FinalHostMatrix$Sp2 <- factor(FinalHostMatrix$Sp2, levels = sort(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))

FinalHostMatrix <- FinalHostMatrix %>% slice(order(Sp,Sp2))

# Simulating on the full network ####

AllMammals <- intersect(rownames(FullSTMatrix), rownames(FullRangeAdj)) %>% setdiff(NonEutherianSp)

AllMammals <- sort(AllMammals)

AllMammalMatrix <- data.frame(
  Sp = as.character(rep(AllMammals,each = length(AllMammals))),
  Sp2 = as.character(rep(AllMammals,length(AllMammals))),
  Space = c(FullRangeAdj[AllMammals,AllMammals]),
  Phylo = c(tFullSTMatrix[AllMammals,AllMammals])
) %>% 
  mutate(Gz = as.numeric(Space==0)) %>% droplevels

UpperMammals <- which(upper.tri(FullSTMatrix[AllMammals, AllMammals], diag = T))

AllMammaldf <- AllMammalMatrix[-UpperMammals,]

N = nrow(AllMammaldf); N

# Adding on space 2 #### 

CurrentsRangeAdjB <- readRDS("Iceberg Output Files/CurrentsRangeAdjB.rds")

CurrentsRangeAdjB %>% reshape2::melt() -> LongBSpace

LongBSpace %>% filter(paste(Var1,Var2, sep = ".")%in%paste(FinalHostMatrix$Sp,FinalHostMatrix$Sp2, sep = ".")) %>%
  slice(order(Var1,Var2))->
  SubLongBSpace

FinalHostMatrix$SpaceA <- FinalHostMatrix$Space
FinalHostMatrix$SpaceB <- SubLongBSpace$value



