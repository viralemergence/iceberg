
# Making raster figures for Colin ####

library(raster)

blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

# Figure 1 ####

Figure1DF = data.frame(
  
  A = GridList$Currents$Richness,
  
  B = GridList$Futures1$Richness - GridList$Currents$Richness,
  D = NEGridList$Futures1$OverlapSum,
  F = NEGridList$Futures1$SharingSum,
  
  C = GridList$Futures4$Richness - GridList$Currents$Richness,
  E = NEGridList$Futures4$OverlapSum,
  G = NEGridList$Futures4$SharingSum
  
)

Figure1Rasters <- lapply(Figure1DF, function(a){
  
  values(blank)[-Sea] <- a
  
  return(blank)
  
})

for(i in 1:length(Figure1Rasters)){
  
  writeRaster(Figure1Rasters[[i]], file = paste0("Iceberg Figures/Figure1",LETTERS[i],".tif"))
  
}

# Figure 3 ####

Figure3DF = data.frame(
  
  A = EbolaGridList$Currents$Richness,
  
  B = EbolaGridList$Futures1$Richness,
  
  C = EbolaNEGridList$Futures1$OverlapSum,
  
  D = EbolaNEGridList$Futures1$SharingSum,
  
  E = BPNEGridList$Futures1$OverlapSum
  
)

Figure3Rasters <- lapply(Figure3DF, function(a){
  
  values(blank)[-Sea] <- a
  
  return(blank)
  
})

for(i in 1:length(Figure3Rasters)){
  
  writeRaster(Figure3Rasters[[i]], file = paste0("Iceberg Figures/Figure3",LETTERS[i],".tif"), overwrite = T)
  
}

