# Iceberg Data Parsing ####

a = "Bison_bonasus"

testraster1 <- raster(paste0('Iceberg Input Files/Currents/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD
testraster1 <- raster::resample(testraster1, blank, method = 'ngb')

testraster2 <- raster(paste0('Iceberg Input Files/Futures1/',a %>% str_replace("_", " "),'.tif')) # CAREFUL OF SPACE BEFORE PERIOD
testraster3 <- raster(paste0('Iceberg Input Files/Futures2/',a %>% str_replace("_", " "),'.tif')) # CAREFUL OF SPACE BEFORE PERIOD

par(mfrow = c(1,3))

plot(testraster1)
plot(testraster2)
plot(testraster3)

sapply(list(testraster1, testraster2, testraster3), function(a) length(na.omit(values(a))))
