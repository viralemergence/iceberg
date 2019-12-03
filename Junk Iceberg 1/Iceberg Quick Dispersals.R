
# Iceberg Quick Dispersals ####
# Rscript "Iceberg Quick Dispersals.R"

library(tidyverse); library(raster); library(parallel)

Method = "MaxEnt"
x = "Currents"

"Iceberg Input Files/" %>% paste0(Method,"/04_LandUseClipped/",x) %>% list.files(full.names = T) %>%
  mclapply(raster) -> RasterListd

names(RasterListd) <- Species <- "Iceberg Input Files/" %>% 
  paste0(Method,"/04_LandUseClipped/",x) %>% 
  list.files() %>% str_remove(".tif")

# Dispersals ####

Dispersals <- read.csv("Iceberg Input Files/Data for dispersal.csv", header = T)

Dispersals <- Dispersals %>% filter(!is.na(Scientific_name), !is.na(disp50))
Dispersals$Scientific_name <- Dispersals$Scientific_name %>% str_replace(" ", "_")

Species <- intersect(Species, Dispersals$Scientific_name)

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

paste0("Iceberg Input Files/", 
       Method,"/05_DispersalBuffers/") %>% 
  list.files %>% str_remove(".tif") %>%
  setdiff(Species, .) -> ToBuffer

mclapply(1:length(ToBuffer), function(i){
  
  Sp <- ToBuffer[i]
  
  Dist <- Dispersals %>% filter(Scientific_name == Sp) %>% pull(disp50)*1000
  
  r1 <- RasterListd[[Sp]]
  r2 <- buffer2(r1, Dist)
  
  writeRaster(r2, file = paste0("Iceberg Input Files/", 
                                Method,"/05_DispersalBuffers/",Sp,".tif"), overwrite = T)
  
}, mc.preschedule = F, mc.cores = 45)
