# Rscript "Iceberg Code/Iceberg Continent Clipping.R" ####

library(tidyverse); library(raster); library(parallel)

ContinentFiles <- setdiff(list.files("Iceberg Input Files/Continents"), "Buffered")
Continents <- ContinentFiles %>% lapply(function(a) raster(paste0("Iceberg Input Files/Continents/",a)))

ContinentsBuffer <- list()
ContinentsBuffer <- Continents

i = 1

mclapply(1:length(Continents), function(i){
  
  print(ContinentFiles[[i]])
  
  SubBuffer <- buffer(Continents[[i]], 50000)
  
  ContinentsBuffer[[i]] <- SubBuffer
  
  writeRaster(SubBuffer, file = paste0("Iceberg Input Files/Continents/Buffered/",ContinentFiles[[i]]))
  
}, mc.cores = 5)

load("~/LargeFiles/MammalStackFullMercator.Rdata")

ContinentValues <- lapply(ContinentsBuffer, values) %>% as.data.frame() 
names(ContinentValues) <- c("Africa", "Australia", "Eurasia", "N_America", "S_America")

ContinentWhich <- apply(ContinentValues, 2, function(a) which(!is.na(a)))

ContinentsInhabited <- names(MammalStackFull) %>% lapply(function(a){
  
  print(a)
  
  Valuedf <- ContinentValues %>%
    mutate(New = values(MammalStackFull[[a]]))
  
  Valuedf <- Valuedf %>% filter(!is.na(New)) %>% dplyr::select(-New)
  
  names(Valuedf)[which(colSums(Valuedf, na.rm = T)>0)] %>% return
  
})

names(ContinentsInhabited) <- names(MammalStackFull)

NoContinents <- which(ContinentsInhabited %>% sapply(function(a) length(a)==0))
NoContinentValues <- MammalStackFull[NoContinents] %>% sapply(function(a) length(na.omit(values(a))))
KeepContinents <- names(NoContinentValues)

# Bit of troubleshooting ####

RasterListc <- RasterListb
names(RasterListc) <- Files %>% str_remove(".tif$") %>% str_replace(" ", "_")
RasterListc <- RasterListc[setdiff(names(RasterListc), DropSpecies)]

i = "Gorilla_gorilla"
WhichFill <- unlist(ContinentWhich[ContinentsInhabited[[i]]])
if(length(table(values(RasterListc[[i]])[-WhichFill]))>0) print(length(table(values(RasterListc[[i]])[-WhichFill]))/(length(na.omit(values(RasterListc[[i]])[WhichFill]))+length(na.omit(values(RasterListc[[i]])[-WhichFill]))))
values(RasterListc[[i]])[WhichFill][!is.na(values(RasterListc[[i]])[WhichFill])] <- 10
values(RasterListc[[i]])[-WhichFill][!is.na(values(RasterListc[[i]])[-WhichFill])] <- 20
plot(RasterListc[[i]])


i = "Canis_lupus"

NotInhabited <- ContinentWhich[setdiff(names(ContinentValues), ContinentsInhabited[[i]])] %>% unlist

values(RasterListc[[i]])[NotInhabited][!is.na(values(RasterListc[[i]])[NotInhabited])] <- 20
values(RasterListc[[i]])[-NotInhabited][!is.na(values(RasterListc[[i]])[-NotInhabited])] <- 10
plot(RasterListc[[i]])





