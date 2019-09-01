
# ~~~~~~~ Iceberg ENM Futures Processing ####

library(tidyverse); library(raster); library(parallel); library(sf)

Method = "MaxEnt"
# Method = "RangeBags"

PredReps <- c("Currents", paste0("Futures", 1:4))

blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

paste0("Iceberg Input Files/", 
       Method,"/06_DispersalBuffersResampled/") %>% list.files %>% str_remove(".tif") -> 
  Species

paste0("Iceberg Input Files/", 
       Method,"/06_DispersalBuffersResampled") %>% list.files(full.names = T) %>%
  mclapply(raster) -> DispersalRasters

names(DispersalRasters) <- Species

for(x in 2:length(PredReps)){
  
  # 01_Initial import ####
  
  x = PredReps[x]
  
  Root <- paste0("Iceberg Input Files/", Method,"/01_Raw/",x)
  
  Files <- list.files(Root)
  
  Species1 <- Root %>% list.files %>% str_remove(".tif")
  
  RasterLista <- lapply(Files, function(a){
    
    SubFiles <- paste0(Root,"/",a) %>% list.files
    
    if(length(SubFiles>0)) raster(paste0(Root,"/",a,"/",SubFiles[1]))
    
  })
  
  names(RasterLista) <- Species1# %>% intersect(Species)
  
  RasterLista <- RasterLista[!sapply(RasterLista, is.null)]
  
  RasterLista <- RasterLista[intersect(names(RasterLista), Species)]
  
  Species <- intersect(names(RasterLista), Species1)
  SpeciesLost <- setdiff(Species1, names(RasterLista))
  
  RasterLista <- RasterLista[Species]
  
  # 02_Resampling ####
  
  mclapply(1:length(Species), function(a){
    
    if(a %% 500==0) print(a)
    
    testraster <- RasterLista[[Species[a]]]
    testraster <- raster::resample(testraster, blank, method = 'ngb')#*AreaRaster
    
    writeRaster(testraster, 
                file = paste0("Iceberg Input Files/",
                              Method,"/02_Resampled/",x,"/",Species[a],".tif"), overwrite = T)
    
    return(testraster)
    
  }, mc.preschedule = F, mc.cores = 45)
  
  remove(RasterLista)
  
  "Iceberg Input Files/" %>% paste0(Method,"/02_Resampled/",x) %>% list.files(full.names = T) %>%
    mclapply(raster) -> RasterListb
  
  names(RasterListb) <- Species <- 
    "Iceberg Input Files/" %>% paste0(Method,"/02_Resampled/",x) %>% 
    list.files() %>%
    str_remove(".tif")
  
  # 03_Dispersal Clipping ####
  
  mclapply(1:length(Species), function(i){
    
    Sp = Species[i]
    
    r1 <- RasterListb[[Sp]]
    r2 <- DispersalRasters[[Sp]]
    
    r1[is.na(r2)] <- NA
    
    writeRaster(r1, 
                file = paste0("Iceberg Input Files/",
                              Method,"/07_DispersalClippedFutures/",x,"/",Species[i],".tif"), overwrite = T)
    
  }, mc.cores = 45)
  
  remove(RasterListb)
  
  paste0("Iceberg Input Files/",
         Method,"/07_DispersalClippedFutures/",x) %>% list.files(full.names = T) %>% mclapply(raster) ->
    RasterListc
  
  paste0("Iceberg Input Files/",
         Method,"/07_DispersalClippedFutures/",x) %>% list.files() %>% str_remove(".tif") ->
    names(RasterListc)
  
  # 04_Land use clipping ####
  
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
  
  LandUse <- brick(paste0("Iceberg Input Files/", x, ".grd"))
  
  PreLUArea <- list() -> PostLUArea
  
  LUPre_Post <- mclapply(1:length(Species), function(i){
    
    Sp <- Species[i]
    
    print(Sp)
    
    SpHabitat <- HabitatList[[Sp]] %>% as.character %>% str_split("[.]") %>% unlist %>% unique %>% 
      na.omit # is this right??
    
    if(length(SpHabitat)>0){
      
      LandUse[[SpHabitat]] %>% getValues -> ValueDF
      
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
                                  Method,"/04_LandUseClipped/",x,"/",Sp,".tif"), overwrite = T)
    
    return(c(PreLUArea, PostLUArea))
    
  })
  
  remove(RasterListc)
  
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
  
  saveRDS(LandUseLossDF, file = paste0("Iceberg Input Files/",Method, x, "LandUseLoss.rds"))
  
}

# 05_Continent clipping and area-correcting ####

load("~/LargeFiles/MammalStackFullMercator.Rdata")

ContinentRaster <- raster("Iceberg Input Files/continents-final.tif") %>%
  resample(blank, method = "ngb")

ContinentWhich <- lapply(1:5, function(a) which(values(ContinentRaster)==a))
names(ContinentWhich) <- c("Africa", "Eurasia", "NAm", "SAm", "Oceania")

paste0("Iceberg Input Files/",Method,"/Final/Futures1") %>% 
  list.files %>% str_remove(".tif") -> ToClip

NoIUCN <- intersect(ToClip, names(MammalStackFull))[which(sapply(intersect(ToClip, names(MammalStackFull)), function(a) MammalStackFull[[a]] %>% values %>% na.omit %>% length)==0)]
NoIUCN2 <- setdiff(ToClip, names(MammalStackFull))
NoIUCN <- union(NoIUCN, NoIUCN2)

ToClip <- setdiff(ToClip, NoIUCN)

for(x in 2:length(PredReps)){
  
  x = PredReps[x]
  
  print(x)
  
  Root <- paste0("Iceberg Input Files/", Method,"/04_LandUseClipped/",x)
  
  Files <- list.files(Root, full.names = T)
  Species <- Root %>% list.files() %>% str_remove(".tif")
  
  Files %>% mclapply(raster) -> RasterListd
  names(RasterListd) <- Species
  
  mclapply(1:length(ToClip), function(a){
    
    Sp <- ToClip[a]
    
    r1 <- MammalStackFull[[Sp]]
    SpWhich <- which(!is.na(values(r1)))
    ContinentsInhabited <- unique(values(ContinentRaster)[SpWhich])
    
    r2 <- RasterListd[[Sp]]*AreaRaster
    values(r2)[-unlist(ContinentWhich[ContinentsInhabited])] <- NA
    
    writeRaster(r2, file = paste0("Iceberg Input Files/",
                                  Method,"/Final/",x,"/",Sp,".tif"),
                overwrite = T)
    
  }, mc.preschedule = F)
  
  mclapply(1:length(NoIUCN), function(a){
    
    Sp <- NoIUCN[a]
    
    r2 <- RasterListd[[Sp]]*AreaRaster
    
    writeRaster(r2, file = paste0("Iceberg Input Files/",
                                  Method,"/Final/",x,"/",Sp,".tif"),
                overwrite = T)
    
  }, mc.preschedule = F)
  
}

# ~~~ These are the final futures #####
