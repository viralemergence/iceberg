library(raster)
library(ncdf4)

tas <- brick("C:/Users/cjcar/Downloads/adaptor.mars.internal-1605386567.6335182-13770-31-377b4dc6-4585-4c7e-88fb-57b7de0a8859.nc")
pre <- brick("C:/Users/cjcar/Downloads/adaptor.mars.internal-1605386711.0461376-16515-23-dfb714e2-04ed-42b3-b4c5-805def808674.nc")

blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180, 180, -90, 90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

pre <- resample(pre, blank)
tas <- resample(tas, blank)

tas <- tas[[1:360]]
writeRaster(pre, "C:/Users/cjcar/Dropbox/ERA5 Iceberg/pre.nc")
writeRaster(tas, "C:/Users/cjcar/Dropbox/ERA5 Iceberg/tas.nc")

