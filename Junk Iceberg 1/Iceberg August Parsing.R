
# August parsing iceberg ####

library(ncdf4)
library(raster)
library(parallel)
library(tidyverse)

Root <- "Pre-Dispersal/Currents"

Files <- list.files(Root)

Species <- Files %>% str_remove(".tif")

Primates <- intersect(Species, Panth1 %>% filter(hOrder == "Primates") %>% pull(Sp))

RasterLista <- lapply(Files, function(a) raster(paste0(Root,"/",a)))

blank <- matrix(0, 360*2, 720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180, 180, -90, 90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

RasterListb <- mclapply(1:length(RasterLista), function(a){
  
  if(a %% 500==0) print(a)
  
  testraster <- RasterLista[[a]]
  testraster <- raster::resample(testraster, blank, method = 'ngb')#*AreaRaster
  
  return(testraster)
  
}, mc.cores = 45)

names(RasterListb) <- Species

print("Getting richness")

OverlapSums <- rep(0, ncol(RasterListb[[1]])*nrow(RasterListb[[1]]))

for(i in 1:length(Primates)){
  
  if(i %% 500 == 0) print(i)
  SubSums <- raster::getValues(RasterListb[[Primates[i]]])
  SubSums[is.na(SubSums)] <- 0
  OverlapSums <- OverlapSums + SubSums
  
}

i = 1

for(i in i:length(RasterListb)){
  
  if(i %% 500 == 0) print(i)
  SubSums <- raster::getValues(RasterListb[[i]])
  SubSums[is.na(SubSums)] <- 0
  OverlapSums <- OverlapSums + SubSums
  
}

GridDF <- data.frame(
  Richness = OverlapSums,
  X = rep(1:ncol(RasterListb[[1]]), nrow(RasterListb[[1]])),
  Y = rep(nrow(RasterListb[[1]]):1, each = ncol(RasterListb[[1]]))
)

UniversalBlank <- raster("Old Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

GridDF <- GridDF[-Sea,]

ggplot(GridDF, aes(X, Y)) + geom_tile(aes(fill = Richness, colour = Richness)) + theme_void() + 
  # lims(x = c(800,950), y = c(300,500)) +
  coord_fixed()

GridDF <- GridDF %>% mutate(DalliIUCN = values(MammalStackFull[["Ovis_dalli"]]),
                            DalliENM = values(RasterListb[["Ovis_dalli"]])) %>%
  mutate(DalliSum = as.numeric(!is.na(DalliIUCN))+as.numeric(!is.na(DalliENM)))

GridDF %>% slice(-Sea) %>%
  ggplot(aes(X, Y)) + 
  geom_tile(aes(fill = DalliSum)) +
  geom_vline(xintercept = 170) +
  coord_fixed()

GridDF %>% mutate(Dalli = values(MammalStackFull[["Ovis_dalli"]])) %>%
  mutate(X = 180*(X - mean(na.omit(X)))/(max(X)/2)) %>%
  #na.omit() %>% pull(X) %>% min
  ggplot(aes(X, Y)) + geom_tile(aes(fill = Richness, colour = Richness)) + # theme_void() + 
  # lims(x = c(800,950), y = c(300,500)) +
  # coord_fixed() + 
  facet_wrap(~is.na(Dalli)) +
  geom_vline(aes(xintercept = -115))



# Pivoting to rangebags because of Cory ####

for(x in c("Commons", "Rares")){
  
  Root <- paste0("Iceberg Input Files/RangeBags/",x)
  
  Files <- list.files(Root)
  
  Species1 <- Files %>% str_remove(".tif")
  
  RangeBagLista <- lapply(Files, function(a){
    
    SubFiles <- paste0(Root,"/",a) %>% list.files
    
    if(length(SubFiles>0)) raster(paste0(Root,"/",a,"/",SubFiles[1]))
    
  })
  
  names(RangeBagLista) <- Species1
  
  RangeBagLista <- RangeBagLista[!sapply(RangeBagLista, is.null)]
  
  Species <- intersect(names(RangeBagLista), Species1)
  SpeciesLost <- setdiff(Species1, names(RangeBagLista))
  
  RangeBagListb <- mclapply(1:length(Species), function(a){
    
    if(a %% 500==0) print(a)
    
    testraster <- RangeBagLista[[a]]
    testraster <- raster::resample(testraster, blank, method = 'ngb')#*AreaRaster
    
    return(testraster)
    
  }, mc.cores = 45)
  
  names(RangeBagListb) <- Species
  
  i = 1
  
  for(i in i:length(Species)) writeRaster(RangeBagListb[[Species[i]]], 
                                          file = paste0("Iceberg Input Files/Resampled RangeBags/",Species[i],".tif"))
  
}

CommonsFiles <- "Iceberg Input Files/RangeBags/Commons" %>% list.files %>% paste0(".tif")#%>% str_remove(".tif")
RaresFiles <- "Iceberg Input Files/RangeBags/Rares" %>% list.files %>% paste0(".tif") #%>% str_remove(".tif")

Fail <- setdiff(CommonsFiles, list.files("Iceberg Input Files/Resampled RangeBags")) %>%
  intersect(names(RangeBagListb) %>% paste0(".tif"))

mclapply(Fail %>% str_remove(".tif"), function(a){ 
  
  writeRaster(RangeBagListb[[a]], 
              file = paste0("Iceberg Input Files/Resampled RangeBags/",a,".tif"))
  
})

# Doing some commons checking ####

Root <- "Iceberg Input Files/RangeBags/Commons"

Files <- list.files(Root)

RangeBagSpecies <- Files %>% str_remove(".tif")

Root <- "Pre-Dispersal/Currents"

Files <- list.files(Root)

MaxEntSpecies <- Files %>% str_remove(".tif")
