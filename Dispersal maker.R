library(ncdf4)
library(raster)

# Blank raster

setwd("D:/ICEBERG")
ncin <- nc_open("states.nc")

get <- function(varname) {
  soilm_array <- ncvar_get(ncin,varname,start=c(1,1,1166))
  slice <- raster(t(soilm_array))
  raster::extent(slice) <- c(-180,180,-90,90)
  return(slice)
}
lutypes <- names(ncin$var)[1]
blank <- lapply(lutypes,get)[[1]]


# Dispersal script

disp <- na.omit(read.csv('C:/Users/cjcar/Dropbox/ViralChatter/Data for dispersal.csv'))
testraster <- raster('D:/ICEBERG/RawENMs/PPM/BinaryLU/Abrocoma_bennettii .tif')
testname <- "Abrocoma bennettii"

testraster <- raster::resample(testraster,blank,method='ngb')
dk <- disp[disp$Scientific_name==testname,'disp50']
buff <- buffer(testraster, dk)

