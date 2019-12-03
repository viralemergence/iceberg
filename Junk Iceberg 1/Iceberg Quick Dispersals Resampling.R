
# Iceberg Quick Dispersals resampling ####
# Rscript "Iceberg Quick Dispersals Resampling.R"

library(tidyverse); library(raster); library(parallel)

blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

Method = "MaxEnt"
x = "Currents"

"Iceberg Input Files/" %>% paste0(Method,"/05_DispersalBuffers") %>% list.files(full.names = T) %>%
  mclapply(raster) -> RasterListe

names(RasterListe) <- Species <- "Iceberg Input Files/" %>% 
  paste0(Method,"/05_DispersalBuffers") %>% 
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
       Method,"/06_DispersalBuffersResampled/") %>% 
  list.files %>% str_remove(".tif") %>%
  setdiff(Species, .) -> ToBuffer

mclapply(1:length(ToBuffer), function(i){
  
  Sp <- ToBuffer[i]
  
  r1 <- RasterListe[[Sp]]
  r2 <- resample(r1, blank, method = "ngb")
  
  writeRaster(r2, file = paste0("Iceberg Input Files/", 
                                Method,"/06_DispersalBuffersResampled/",Sp,".tif"), 
              overwrite = T)
  
}, mc.preschedule = F, mc.cores = 45)
