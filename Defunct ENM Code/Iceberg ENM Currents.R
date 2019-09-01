
# Iceberg pre-GAM spatial processing ####

library(tidyverse); library(raster); library(parallel); library(sf)

Method = "MaxEnt"
#Method = "RangeBags"

PredReps <- c("Currents", paste0("Futures", 1:4))

blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

# 01_Processing Currents ####

# 02_Resampling rasters ####

Root <- paste0("Iceberg Input Files/", Method,"/01_Raw/","Currents")

Files <- list.files(Root)

Species1 <- Root %>% list.files %>% str_remove(".tif")

RasterLista <- lapply(Files, function(a){
  
  SubFiles <- paste0(Root,"/",a) %>% list.files
  
  if(length(SubFiles>0)) raster(paste0(Root,"/",a,"/",SubFiles[1]))
  
})

names(RasterLista) <- Species1

RasterLista <- RasterLista[!sapply(RasterLista, is.null)]

Species <- intersect(names(RasterLista), Species1)
SpeciesLost <- setdiff(Species1, names(RasterLista))

Species <- intersect(Species, SpeciesList[[Method]])

RasterLista <- RasterLista[Species]

mclapply(1:length(Species), function(a){
  
  if(a %% 500==0) print(a)
  
  testraster <- RasterLista[[Species[a]]]
  testraster <- raster::resample(testraster, blank, method = 'ngb')
  
  writeRaster(testraster,
              file = paste0("Iceberg Input Files/",
                            Method,"/02_Resampled/","Currents","/",Species[a],".tif"), overwrite = T)
  
  return(testraster)
  
}, mc.preschedule = F)

"Iceberg Input Files/" %>% paste0(Method,"/02_Resampled/","Currents") %>% list.files(full.names = T) %>%
  mclapply(raster) -> RasterListb

names(RasterListb) <- Species <-
  "Iceberg Input Files/" %>% paste0(Method,"/02_Resampled/","Currents") %>%
  list.files() %>%
  str_remove(".tif")

# 03_Clipping by the IUCN ranges ####

if(file.exists(paste0("Iceberg Input Files/IUCNBuffers",Method,".Rdata"))){ load(paste0("Iceberg Input Files/IUCNBuffers",Method,".Rdata"))} else{
  
  Mammal_Shapes <- st_read("~/ShapeFiles")
  Mammal_Shapes$Binomial = str_replace(Mammal_Shapes$binomial, " ", "_")
  Mammal_Shapes <- Mammal_Shapes[order(Mammal_Shapes$Binomial),]
  Mammal_Shapes <- Mammal_Shapes[Mammal_Shapes$Binomial%in%Species,]
  
  IUCNBuffers <- GraskBy(Mammal_Shapes, "Binomial", Distance = 10^6)
  
  save(IUCNBuffers, file = paste0("Iceberg Input Files/IUCNBuffers",Method,".Rdata"))
  
}

IUCNBufferRasters <- list()

for(i in 1:length(Species)){
  
  Sp = Species[i]
  print(Sp)
  
  if(Sp%in%names(IUCNBuffers)){
    
    sf1 <- st_cast(IUCNBuffers[[Sp]], "MULTIPOLYGON")
    
    r1 <- fasterize::fasterize(sf1, blank)
    IUCNBufferRasters[[Sp]] <- r1
    
  } else{
    
    r1 <- blank
    
    values(r1) <- 1
    
    IUCNBufferRasters[[Sp]] <- r1
    
  }
  
}

mclapply(Species, function(a){
  
  print(a)
  
  r1 <- RasterListb[[a]]
  r2 <- IUCNBufferRasters[[a]]
  
  r1[is.na(r2)] <- NA
  
  writeRaster(r1,
              file = paste0("Iceberg Input Files/",
                            Method,"/03_MaskedCurrents/",a,".tif"), overwrite = T)
  
}, mc.preschedule = F)

"Iceberg Input Files/" %>% paste0(Method,"/03_MaskedCurrents") %>% list.files(full.names = T) %>%
  mclapply(raster) -> RasterListc

names(RasterListc) <- Species <-
  "Iceberg Input Files/" %>% paste0(Method,"/03_MaskedCurrents") %>% list.files() %>%
  str_remove(".tif")

# 04_Land Use filters #####

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

PreLUArea <- list() -> PostLUArea

RasterListd <- list()

LUPre_Post <- mclapply(1:length(Species), function(i){
  
  Sp <- Species[i]
  
  print(Sp)
  
  SpHabitat <- HabitatList[[Sp]] %>% as.character %>% str_split("[.]") %>% unlist %>% unique %>%
    na.omit # is this right??
  
  if(length(SpHabitat)>0){
    
    landuse2017[[SpHabitat]] %>% getValues -> ValueDF
    
    if(length(SpHabitat)==1){
      
      Uninhabitable <- which(ValueDF==0)
      
    }else{
      
      Uninhabitable <- which((ValueDF %>% rowSums(na.rm = T))==0)
      
    }
    
    r1 <- RasterListc[[Sp]]
    
    PreLUArea <- r1 %>% values %>% na.omit %>% length
    
    values(r1)[Uninhabitable] <- NA
    
    PostLUArea <- r1 %>% values %>% na.omit %>% length
    
  }else{
    
    r1 <- RasterListc[[Sp]]
    
    PreLUArea <- PostLUArea <- r1 %>% values %>% na.omit %>% length
    
  }
  
  writeRaster(r1, file = paste0("Iceberg Input Files/",
                                Method,"/04_LandUseClipped/Currents/",Sp,".tif"), overwrite = T)
  
  return(c(PreLUArea, PostLUArea))
  
})

names(HabitatList)[which(HabitatList %>% sapply(function(a) a %>% na.omit %>% length)==0)] ->
  NoLUSpecies

LUPre_Post %>% bind_cols() %>% t %>% as.data.frame() %>%
  mutate(Species = Species[1:n()]) %>%
  rename(PreLUArea = V1, PostLUArea = V2) %>%
  mutate(PercentLost = 1- (PostLUArea/PreLUArea)) ->
  LandUseLossDF
LandUseLossDF %>%
  mutate(NoLUInfo = as.numeric(Species %in% NoLUSpecies)) ->
  LandUseLossDF
qplot(LandUseLossDF$PercentLost)
mean(LandUseLossDF$PercentLost, na.rm = T)
LandUseLossDF %>% filter(NoLUInfo == 0) %>% pull(PercentLost) %>% mean(na.rm = T)
save(LandUseLossDF, file = paste0("Iceberg Input Files/", Method, "Currents", "LandUseLoss.Rdata"))

"Iceberg Input Files/" %>% paste0(Method,"/04_LandUseClipped/","Currents") %>% list.files(full.names = T) %>%
  mclapply(raster) -> RasterListd

names(RasterListd) <- Species <- "Iceberg Input Files/" %>%
  paste0(Method,"/04_LandUseClipped/","Currents") %>%
  list.files() %>% str_remove(".tif")

# 04b_Continent clipping and area-correcting ####

ContinentRaster <- raster("Iceberg Input Files/continents-final.tif") %>%
  resample(blank, method = "ngb")

ContinentWhich <- lapply(1:5, function(a) which(values(ContinentRaster)==a))
names(ContinentWhich) <- c("Africa", "Eurasia", "NAm", "SAm", "Oceania")

"Iceberg Input Files/" %>% paste0(Method,"/04_LandUseClipped/","Currents") %>% list.files(full.names = T) %>%
  mclapply(raster) -> RasterListd

names(RasterListd) <- Species <- "Iceberg Input Files/" %>%
  paste0(Method,"/04_LandUseClipped/","Currents") %>%
  list.files() %>% str_remove(".tif")

paste0("Iceberg Input Files/",Method,"/Final/Currents") %>% list.files %>% str_remove(".tif") %>%
  setdiff(Species, .) -> ToClip

load("~/LargeFiles/MammalStackFullMercator.Rdata")

NoIUCN <- intersect(ToClip, names(MammalStackFull))[which(sapply(intersect(ToClip, names(MammalStackFull)), function(a) MammalStackFull[[a]] %>% values %>% na.omit %>% length)==0)]
NoIUCN2 <- setdiff(ToClip, names(MammalStackFull))
NoIUCN <- union(NoIUCN, NoIUCN2)

ToClip <- setdiff(ToClip, NoIUCN)

mclapply(1:length(ToClip), function(a){
  
  Sp <- ToClip[a]
  
  r1 <- MammalStackFull[[Sp]]
  SpWhich <- which(!is.na(values(r1)))
  ContinentsInhabited <- unique(values(ContinentRaster)[SpWhich])
  
  r2 <- RasterListd[[Sp]]*AreaRaster
  values(r2)[-unlist(ContinentWhich[ContinentsInhabited])] <- NA
  
  writeRaster(r2, file = paste0("Iceberg Input Files/",
                                Method,"/Final/","Currents","/",Sp,".tif"),
              overwrite = T)
  
}, mc.preschedule = F, mc.cores = 45)

mclapply(1:length(NoIUCN), function(a){
  
  Sp <- NoIUCN[a]
  
  r2 <- RasterListd[[Sp]]*AreaRaster
  
  writeRaster(r2, file = paste0("Iceberg Input Files/",
                                Method,"/Final/","Currents","/",Sp,".tif"),
              overwrite = T)
  
}, mc.preschedule = F)

# ~~~ These are the final currents #####

# 05_Dispersals ####

Dispersals <- read.csv("Iceberg Input Files/Data for dispersal.csv", header = T)

Dispersals <- Dispersals %>% filter(!is.na(Scientific_name), !is.na(disp50))
Dispersals$Scientific_name <- Dispersals$Scientific_name %>% str_replace(" ", "_")

paste0("Iceberg Input Files/",
       Method,"/04_LandUseClipped/", "Currents") %>% list.files %>% str_remove(".tif") ->
  Species

ToBuffer <- intersect(Species, Dispersals$Scientific_name)

NotToBuffer <- setdiff(Species, ToBuffer)

mclapply(1:length(ToBuffer), function(i){
  
  Sp <- ToBuffer[i]
  
  Dist <- Dispersals %>% filter(Scientific_name == Sp) %>% pull(disp50)*1000
  
  r1 <- RasterListd[[Sp]]
  
  r2 <- buffer2(r1, Dist)
  
  writeRaster(r2, file = paste0("Iceberg Input Files/",
                                Method,"/05_DispersalBuffers/",Sp,".tif"), overwrite = T)
  
  r3 <- resample(r2, blank, method = "ngb")
  
  writeRaster(r3, file = paste0("Iceberg Input Files/",
                                Method,"/06_DispersalBuffersResampled/",Sp,".tif"),
              overwrite = T)
  
}, mc.preschedule = F)


mclapply(1:length(NotToBuffer), function(i){
  
  Sp <- NotToBuffer[i]
  
  r1 <- RasterListd[[Sp]]
  
  r3 <- resample(r1, blank, method = "ngb")
  
  values(r3) <- 1
  
  writeRaster(r3, file = paste0("Iceberg Input Files/",
                                Method,"/06_DispersalBuffersResampled/",Sp,".tif"),
              overwrite = T)
  
}, mc.preschedule = F)

paste0("Iceberg Input Files/",
       Method,"/05_DispersalBuffers/") %>% list.files(full.names = T) %>%
  mclapply(raster) -> RasterListe

# 06_Resampling Dispersals ####

mclapply(1:length(Species), function(a){
  
  Sp <- Species[a]
  
  r1 <- RasterListe[[Sp]]
  r2 <- resample(r1, blank, method = "ngb")
  
  writeRaster(r2, file = paste0("Iceberg Input Files/",
                                Method,"/06_DispersalBuffers_Resampled/",Sp,".tif"),
              overwrite = T)
  
}, mc.preschedule = F)

paste0("Iceberg Input Files/",
       Method,"/06_DispersalBuffers_Resampled/") %>% list.files(full.names = T) %>%
  mclapply(raster) -> RasterListf

paste0("Iceberg Input Files/",
       Method,"/06_DispersalBuffers_Resampled/") %>% list.files %>% str_remove(".tif") ->
  names(RasterListf)
