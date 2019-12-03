
# Making the Amazon Hole Map for Colin ####

load("~/LargeFiles/MammalStackFullMercator.Rdata")

ToBagSpecies <- read_csv("Iceberg Input Files/ToBagSpecies.csv")$hull %>%
  str_replace(" ", "_")

NotBagSpecies <- read_csv("Iceberg Input Files/NotBagSpecies.csv")$point %>%
  str_replace(" ", "_")

# Iceberg General Maps ####

ToBagList <- list()

RasterListb <- lapply(1:length(ToBagSpecies), function(a){
  
  if(a %% 500==0) print(a)
  
  MammalStackFull[[ToBagSpecies[a]]]
  
})

ToSkip <- which(sapply(RasterListb, is.null))

names(RasterListb) <- ToBagSpecies

RasterListb <- RasterListb[-ToSkip]
ToBagSpecies <- ToBagSpecies[-ToSkip]

OverlapSums <- rep(0, ncol(RasterListb[[1]])*nrow(RasterListb[[1]]))

for(i in 1:length(RasterListb)){
  
  if(i %% 500 == 0) print(i)
  SubSums <- raster::getValues(RasterListb[[i]]) > 0 %>% as.numeric
  SubSums[is.na(SubSums)] <- 0
  OverlapSums <- OverlapSums + SubSums
  
}

ToBagDF <- data.frame(
  Richness = OverlapSums,
  X = rep(1:ncol(RasterListb[[1]]), nrow(RasterListb[[1]])),
  Y = rep(nrow(RasterListb[[1]]):1, each = ncol(RasterListb[[1]]))
)

remove(RasterListb)

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

ToBagDF <- ToBagDF[-Sea,]

ggplot(ToBagDF, aes(X, Y, fill = Richness)) + geom_tile() + 
  coord_fixed()

# Iceberg General Maps ####

RasterListb <- lapply(1:length(NotBagSpecies), function(a){
  
  if(a %% 500==0) print(a)
  
  MammalStackFull[[NotBagSpecies[a]]]
  
})

ToSkip <- which(sapply(RasterListb, is.null))

names(RasterListb) <- NotBagSpecies

RasterListb <- RasterListb[-ToSkip]
NotBagSpecies <- NotBagSpecies[-ToSkip]

OverlapSums <- rep(0, ncol(RasterListb[[1]])*nrow(RasterListb[[1]]))

for(i in 1:length(RasterListb)){
  
  if(i %% 500 == 0) print(i)
  SubSums <- raster::getValues(RasterListb[[i]]) > 0 %>% as.numeric
  SubSums[is.na(SubSums)] <- 0
  OverlapSums <- OverlapSums + SubSums
  
}

NotBagDF <- data.frame(
  Richness = OverlapSums,
  X = rep(1:ncol(RasterListb[[1]]), nrow(RasterListb[[1]])),
  Y = rep(nrow(RasterListb[[1]]):1, each = ncol(RasterListb[[1]]))
)

remove(RasterListb)

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

NotBagDF <- NotBagDF[-Sea,]

ggplot(NotBagDF, aes(X, Y, fill = Richness)) + geom_tile() + 
  coord_fixed()



for(i in 1:length(names(MammalStackFull))){
  
  if(i %% 500 == 0) print(i)
  SubSums <- raster::getValues(MammalStackFull[[i]]) > 0 %>% as.numeric
  SubSums[is.na(SubSums)] <- 0
  OverlapSums <- OverlapSums + SubSums
  
}

AllDF <- data.frame(
  X = rep(1:ncol(MammalStackFull[[1]]), nrow(MammalStackFull[[1]])),
  Y = rep(nrow(MammalStackFull[[1]]):1, each = ncol(MammalStackFull[[1]])),
  Richness = OverlapSums
)

AllDF <- AllDF[-Sea,]

ggplot(AllDF, aes(X, Y, fill = Richness)) + geom_tile() + 
  coord_fixed()

AllDF <- AllDF %>% mutate(ToBag = ToBagDF$Richness,
                          NotBag = NotBagDF$Richness) %>%
  mutate(ENM = Richness - (ToBag + NotBag))


AllDFLong <- gather(AllDF, key = "key", value = "value", Richness, ToBag, NotBag, ENM) %>%
  group_by(key) %>% mutate(value = value/max(value))

ggplot(AllDFLong, aes(X, Y, fill = value)) + facet_wrap(~key) + geom_tile() + coord_fixed() +
  theme(strip.background = element_rect(fill = "white"))

