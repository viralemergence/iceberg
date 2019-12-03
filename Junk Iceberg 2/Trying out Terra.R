
# Trying out Terra ####

library(raster); library(terra)


i = "Abrocoma_bennettii"

r1 <- raster(paste0("Iceberg Input Files/Currents/",i,".tif"))                         
length(na.omit(values(r1)))

t1 = Sys.time()

r2 <- terra::buffer(r1, 100000)

t2 = Sys.time()

r3 <- raster::buffer(r1, 100000)

t3 = Sys.time()

t2-t1 # Terra Time
t3-t2 # Raster Time


# Trying a bigger one ####

i = Fail[1000]

i

r1 <- raster(paste0("Iceberg Input Files/Currents/",i,".tif"))                         
length(na.omit(values(r1)))

t1 = Sys.time()

r2 <- terra::buffer(r1, 100000)

t2 = Sys.time()

r3 <- raster::buffer(r1, 100000)

t3 = Sys.time()

t2-t1 # Terra Time
t3-t2 # Raster Time

