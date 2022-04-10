
library(dismo)
library(raster)

temp <- raster::brick("C:/Users/cjcar/Downloads/adaptor.mars.internal-1605386567.6335182-13770-31-377b4dc6-4585-4c7e-88fb-57b7de0a8859.nc")
temp1 <- temp[[1:180]]
temp2 <- temp[[181:360]]

prep <- raster::brick("C:/Users/cjcar/Downloads/adaptor.mars.internal-1605386711.0461376-16515-23-dfb714e2-04ed-42b3-b4c5-805def808674.nc")
prep1 <- prep[[1:180]]
prep2 <- prep[[181:360]]

########

# Past 

bio1a <- raster::calc(temp1, fun = mean)
bio4a <- 100 * raster::calc(temp1, fun = sd)

bio1b <- raster::calc(temp2, fun = mean)
bio4b <- 100 * raster::calc(temp2, fun = sd)

setwd("C:/Users/cjcar/Desktop")

writeRaster(bio1a, 'bio1-19811995.tif')
writeRaster(bio4a, 'bio4-19811995.tif')
writeRaster(bio1b, 'bio1-20052019.tif')
writeRaster(bio4b, 'bio4-20052019.tif')