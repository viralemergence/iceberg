library(ncdf4)
library(raster)


### fix blank

blank <- matrix(0,720,360)
blank <- blanker(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")
  
# Dispersal script


disp <- na.omit(read.csv('C:/Users/cjcar/Dropbox/ViralChatter/Data for dispersal.csv'))
testraster <- raster('D:/ICEBERG/RawENMs/PPM/BinaryLU/Abrocoma_bennettii .tif')
testname <- "Abrocoma bennettii" # this needs to be part of the parallel

testraster <- raster::resample(testraster,blank,method='ngb')
dk <- disp[disp$Scientific_name==testname,'disp50']
buff <- buffer(testraster, dk)
writeRaster(buff, 'filename') # needs to be fixed

