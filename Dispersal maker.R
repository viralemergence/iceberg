# Dispersal ####

library(ncdf4)
library(raster)
library(parallel)

### fix blank

blank <- matrix(0,720,360)
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

# Dispersal script

Dispersals <- read.csv("data/Data for dispersal.csv", header = T)
disp <- Dispersals %>% filter(!is.na(Scientific_name))

RasterSp <- list.files("IceMaps") %>% str_split(" ") %>% map(1) %>% sapply(function(a) a[1]) %>% sort

all(RasterSp%in%disp$Scientific_name)
all(disp$Scientific_name %in% RasterSp)

t1 <- Sys.time()

i = 1
j = length(RasterSp)

mclapply(RasterSp[i:j], function(a){
  
  testraster <- raster(paste0('IceMaps/',a,' .tif')) # CAREFUL OF SPACE BEFORE PERIOD
  
  testraster <- raster::resample(testraster, blank, method = 'ngb')
  dk <- disp[disp$Scientific_name == a, 'disp50']
  buff <- buffer(testraster, dk)
  writeRaster(buff, filename = paste0("PostDispersal/",a,'.tif'))
  
}, mc.cores = 20)

t2 <- Sys.time()

t2 - t1
