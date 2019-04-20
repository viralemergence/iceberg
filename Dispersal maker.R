library(ncdf4)
library(raster)
library(parallel)

### fix blank

blank <- matrix(0,720,360)
blank <- blanker(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

# Dispersal script

disp <- na.omit(read.csv('C:/Users/cjcar/Dropbox/ViralChatter/Data for dispersal.csv'))

i = 1
j = nrow(disp)

t1 <- Sys.time()

mclapply(disp$Scientific_name[i:j], function(a){
  
  testraster <- raster(paste0('D:/ICEBERG/RawENMs/PPM/BinaryLU/',a,'.tif'))

  testraster <- raster::resample(testraster, blank, method = 'ngb')
  dk <- disp[disp$Scientific_name == a, 'disp50']
  buff <- buffer(testraster, dk)
  writeRaster(buff, filename = paste0("PostDispersal/",a,'.tif'))
  
}, mc.cores = 20)

t2 <- Sys.time()

t2 - t1