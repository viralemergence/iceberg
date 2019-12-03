
# Iceberg Quick Dispersals ####
# Rscript "Iceberg Quick Continent Clipping.R"

library(tidyverse); library(raster); library(parallel)

load("~/LargeFiles/MammalStackFullMercator.Rdata")

Method = "MaxEnt"
x = "Currents"

blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

# Clipping by continents ####

ContinentRaster <- raster("Iceberg Input Files/continents-final.tif") %>%
  resample(blank, method = "ngb")

ContinentWhich <- lapply(1:5, function(a) which(values(ContinentRaster)==a))
names(ContinentWhich) <- c("Africa", "Eurasia", "NAm", "SAm", "Oceania")

"Iceberg Input Files/" %>% paste0(Method,"/04_LandUseClipped/",x) %>% list.files(full.names = T) %>%
  mclapply(raster) -> RasterListd

names(RasterListd) <- Species <- "Iceberg Input Files/" %>% 
  paste0(Method,"/04_LandUseClipped/",x) %>% 
  list.files() %>% str_remove(".tif")

"Iceberg Input Files/MaxEnt/Final/Currents" %>% list.files %>% str_remove(".tif") %>%
  setdiff(Species, .) -> ToClip

NoIUCN <- intersect(ToClip, names(MammalStackFull))[which(sapply(intersect(ToClip, names(MammalStackFull)), function(a) MammalStackFull[[a]] %>% values %>% na.omit %>% length)==0)]
NoIUCN2 <- setdiff(ToClip, names(MammalStackFull))
NoIUCN <- union(NoIUCN, NoIUCN2)

ToClip <- setdiff(ToClip, NoIUCN)

mclapply(1:length(ToClip), function(a){
  
  Sp <- ToClip[a]
  
  r1 <- MammalStackFull[[Sp]]
  SpWhich <- which(!is.na(values(r1)))
  ContinentsInhabited <- unique(values(ContinentRaster)[SpWhich])
  
  r2 <- RasterListd[[Sp]]
  values(r2)[-unlist(ContinentWhich[ContinentsInhabited])] <- NA 
  
  writeRaster(r2, file = paste0("Iceberg Input Files/", 
                                Method,"/Final/",x,"/",Sp,".tif"), 
              overwrite = T)
  
}, mc.preschedule = F, mc.cores = 45)


mclapply(1:length(NoIUCN), function(a){
  
  Sp <- NoIUCN[a]
  
  r2 <- RasterListd[[Sp]]
  
  writeRaster(r2, file = paste0("Iceberg Input Files/", 
                                Method,"/Final/",x,"/",Sp,".tif"), 
              overwrite = T)
  
}, mc.preschedule = F, mc.cores = 45)

