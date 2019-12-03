
# Getting Colin Mouse 

buffer2 <- function(r, dist) {
  projR <- projectRaster(r, crs=CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
  projRb <- raster::buffer(projR, dist)
  projRb <- projectRaster(projRb, crs=CRS("+proj=longlat +datum=WGS84"))
  projRb[!is.na(projRb)] <- 1
  return(projRb)
}

mammal_shapes <- st_read("~/ShapeFiles")

mammal_shapes <- st_transform(mammal_shapes, 
                              CRS("+proj=longlat +datum=WGS84"))

mammal_shapes$binomial = str_replace(mammal_shapes$binomial, " ", "_")
mammal_shapes <- mammal_shapes[mammal_shapes$binomial=="Mus_musculus",]

random <- raster(list.files('~/IcebergRasters/Abrocoma_bennettii', full.names = TRUE))
random <- extend(random, c(-180,180,-90,90))
random[is.na(random)] <- 1
mammal_raster_full <- fasterize(mammal_shapes, random)

print("Fasterising!")

FullMammalRanges <- fasterize(mammal_shapes, mammal_raster_full, by = "binomial")

BufferMouse <- buffer2(FullMammalRanges, 1000000)

writeRaster(BufferMouse, "Mus_musculus.tif")
