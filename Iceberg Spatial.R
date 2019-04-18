
# Parsing Cory's maps ####

library(velox);
library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger)

Files <- list.files("~/IcebergRasters")

table(table(Files))

UpperFiles <- Files
substr(UpperFiles, 1,1) <- substr(UpperFiles, 1,1) %>% toupper()

write.csv(Files %>% str_replace("_", " "), file = "SpNames_Unedited.csv")
write.csv(UpperFiles %>% str_replace("_", " "), file = "SpNames_Upper.csv")

ToRemove <- Files[which(sapply(Files, function(a) length(list.files(paste("~/IcebergRasters", a, sep = '/'))))==0)]

file.remove(paste("~/IcebergRasters", ToRemove, sep = '/'))

Files <- setdiff(Files, ToRemove)

#Files <- unique(Files)

VeloxList <- lapply(Files, function(a){
  print(which(Files==a))
  print(a)
  r1 <- velox(list.files(paste("~/IcebergRasters", a, sep = '/'), full.names = TRUE)[1])
})

RasterLista <- lapply(1:length(VeloxList), function(a){
  
  print(Files[a])
  
  VeloxList[[a]]$as.RasterLayer(band = 1) #%>% rasterToPolygons(dissolve = T)
  
})

RasterListb <- lapply(1:length(RasterLista), function(a){
  
  print(Files[a])
  
  RasterLista[[a]] %>% rasterToPolygons(dissolve = T)
  
})

Files <- Files[RasterListb[-which(sapply(RasterListb, length)==0)]]
RasterListb <- RasterListb[-which(sapply(RasterListb, length)==0)]

RasterListc <- st_as_sf(bind(RasterListb))

RasterListc <- st_transform(RasterListc, 
                              "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs") # Mollweide projection 

# plot(RasterListc[,"Binomial"])

UpperFiles <- Files
substr(UpperFiles, 1,1) <- substr(UpperFiles, 1,1) %>% toupper()

RasterListc$Binomial <- UpperFiles

IcebergRaster <- raster(RasterListc, res = 25000)

RasterBrick <- fasterize(RasterListc, IcebergRaster, by = "Binomial")

IcebergAdj <- PairsWisely(RasterBrick)

save(IcebergAdj, file = "IcebergRangeAdj.Rdata")
save(RasterBrick, file = "IcebergRasterBrick.Rdata")

load("data/FullRangeOverlap.Rdata")

IUCN_ENM_Sp <- intersect(rownames(FullRangeAdj), UpperFiles)

SuboverlapIUCN <- FullRangeAdj[IUCN_ENM_Sp,IUCN_ENM_Sp]
SuboverlapENM <- IcebergAdj[IUCN_ENM_Sp,IUCN_ENM_Sp]

IUCN_ENM_DF <- data.frame(
  Sp = rep(IUCN_ENM_Sp, each = length(IUCN_ENM_Sp)),
  Sp2 = rep(IUCN_ENM_Sp, length(IUCN_ENM_Sp)),
  IUCN = c(SuboverlapIUCN),
  ENM = c(SuboverlapENM)
) %>% slice(which(lower.tri(SuboverlapIUCN)))

ggplot(IUCN_ENM_DF, aes(IUCN, ENM)) + geom_point()

ggplot(IUCN_ENM_DF, aes(IUCN, ENM)) + 
  geom_point(alpha = 0.01) + 
  coord_fixed() + 
  geom_smooth() + AlberTheme + 
  ggsave("Iceberg_Correlations.jpeg", 
         units = "mm", width = 100, height = 100, dpi = 300)

# Whales??? 

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hFamily = MSW05_Family, hOrder = MSW05_Order)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

LeftOut <- setdiff(names(RasterBrick), rownames(FullRangeAdj))

Panth1 %>% filter(Sp%in%LeftOut) %>% select(Sp, hFamily, hOrder) %>% write.csv("uwu.csv")

# Replacing dispersal data ####

Dispersals <- read.csv("data/Data for dispersal.csv", header = T)
Dispersals$Scientific_name <- Dispersals$Scientific_name %>% str_replace(" ", "_")

NeedDispersals <- read.csv("data/species without dispersal.csv", header = T)

setdiff(NeedDispersals$names, Dispersals$Scientific_name)

library(geiger);library(ape);library(picante);library(dplyr)

STFull <- read.nexus("data/ele_1307_sm_sa1.tre")[[1]]
FullSTMatrix <- as.data.frame(cophenetic(STFull)) %>% as.matrix

NeedDispersals$names <- NeedDispersals$names %>% str_replace(" ", "_")

NewSTNames <- list()
NewSTNames[1:nrow(FullSTMatrix)] <- NA

i = 1

for(i in i:length(NewSTNames)){
  print(i)
  temp <- gnr_resolve(names = rownames(FullSTMatrix)[i])[1,] %>% as.data.frame()
  temp <- temp$matched_name %>% str_split(" ")
  if(substr(temp[[1]][2],1,1) == "(") temp[[1]] <- temp[[1]][-2]
  NewSTNames[[i]] <- temp[[1]][1:2] %>% paste(collapse = "_")
}

NewSTNames2 <- NewSTNames %>% unlist()

setdiff(NeedDispersals$names, NewSTNames2)
setdiff(NewSTNames2, rownames(FullSTMatrix))

# Trying with dispersal names

NewDispNames <- list()
NewDispNames[1:length(NeedDispersals$names)] <- NA

i = 1

for(i in i:length(NewDispNames)){
  print(i)
  temp <- gnr_resolve(names = NeedDispersals$names[i])[1,] %>% as.data.frame()
  temp <- temp$matched_name %>% str_split(" ")
  if(substr(temp[[1]][2],1,1) == "(") temp[[1]] <- temp[[1]][-2]
  NewDispNames[[i]] <- temp[[1]][1:2] %>% paste(collapse = "_")
}

NewDispNames2 <- NewDispNames %>% unlist()

cbind(NeedDispersals$names, NewDispNames2)

setdiff(NeedDispersals$names, NewDispNames2)
setdiff(NeedDispersals$names, NewDispNames2) %in% NewSTNames2

# Trying with all Dispersal species

abm <- read.csv("data/all filtered binary maps.csv", header = T)
abm$names <- str_replace(abm$names, " ", "_")
NoPhylo <- setdiff(abm$names, rownames(FullSTMatrix))

NoPhylo %>% length
NeedDispersals$names %>% length

FinalDispNames <- list()
FinalDispNames[1:length(NoPhylo)] <- NA

i = 1

for(i in i:length(FinalDispNames)){
  print(i)
  temp <- gnr_resolve(names = NoPhylo[i])[1,] %>% as.data.frame()
  temp <- temp$matched_name %>% str_split(" ")
  if(substr(temp[[1]][2],1,1) == "(") temp[[1]] <- temp[[1]][-2]
  FinalDispNames[[i]] <- temp[[1]][1:2] %>% paste(collapse = "_")
}

FinalDispNames2 <- FinalDispNames %>% unlist()

cbind(NoPhylo, FinalDispNames2)

setdiff(NoPhylo, FinalDispNames2)
setdiff(NoPhylo, FinalDispNames2) %in% NewSTNames2
setdiff(FinalDispNames2, rownames(FullSTMatrix))
intersect(FinalDispNames2, rownames(FullSTMatrix))

FinalSpeciesList <- reduce(list(rownames(FullSTMatrix),
                                abm$names,
                                Dispersals$Scientific_name),
                           intersect)

length(FinalSpeciesList)

write.csv(FinalSpeciesList, file = "OwO.csv")
