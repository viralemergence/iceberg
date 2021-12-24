
# Iceberg pre-GAM spatial processing ####

# Rscript "02_Iceberg ENM Futures.R"

t1 <- Sys.time()

library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr)

paste0("~/Albersnet/Iceberg Files/", 
       "Climate1/Iceberg Input Files/GretCDF/Currents") %>% 
  list.files() %>% 
  str_remove(".rds$") %>% sort ->
  Species

PredReps <- c("Currents", paste0("Futures", 1:4))

# Blanks
blank <- matrix(0,360*2, 720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

# Grid areas
AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

# Land Use Data
iucndat <- read.csv('Iceberg Input Files/IucnHabitatData.csv')
convcodes <- read.csv('Iceberg Input Files/IUCN_LUH_conversion_table.csv')

iucndat %>%
  left_join(convcodes, by = c("code" = "IUCN_hab")) %>%
  mutate(name = name %>% str_replace(" ", "_")) ->
  Habitats

lapply(Species, function(a){
  
  Habitats %>% filter(name == a) %>% pull(LUH)
  
}) -> HabitatList

names(HabitatList) <- Species

LandUseList <- lapply(PredReps[2:length(PredReps)], function(a){
  
  brick(paste0("Iceberg Input Files/LandUse", a, ".grd"))
  
})

names(LandUseList) <- PredReps[2:length(PredReps)]

# Continents ####
print("Continents!")

ContinentRaster <- raster("Iceberg Input Files/continents-madagascar.tif") %>%
  resample(blank, method = "ngb")

ContinentWhich <- lapply(1:max(values(ContinentRaster), na.rm = T), function(a) which(values(ContinentRaster)==a))
names(ContinentWhich) <- c("Africa", "Eurasia", "Greenland", "Madagascar", "NAm", "Oceania", "SAm")

# 01_Processing Currents ####

paste0("~/Albersnet/Iceberg Files/", "Climate1/Iceberg Input Files/MaxEnt/01_Raw/Currents") %>% 
  list.files(full.names = T) %>% 
  append(paste0("~/Albersnet/Iceberg Files/", "Climate1/Iceberg Input Files/RangeBags/01_Raw/Currents") %>% list.files(full.names = T)) ->
  FullFiles

# NFiles <- FullFiles %>% map(~.x %>% list.files %>% length)

paste0("Iceberg Input Files/","MaxEnt","/01_Raw/Currents") %>% 
  list.files() %>% #str_remove(".rds$") %>% str_split("__") %>% map_chr(2) %>%
  append(paste0("Iceberg Input Files/","RangeBags","/01_Raw/Currents") %>% 
           list.files()) ->
  names(FullFiles)

FullFiles <- FullFiles[!duplicated(names(FullFiles))]

Files <- FullFiles[Species]

Files %>% str_detect("RangeBag") %>% as.numeric %>% add(1) %>%
  
  c("MaxEnt", "RangeBag")[.] ->
  Methods

names(Methods) <- Species

# Setting raster standards ####

XMin <- blank %>% extent %>% magrittr::extract(1)
XMax <- blank %>% extent %>% magrittr::extract(2)
YMin <- blank %>% extent %>% magrittr::extract(3)
YMax <- blank %>% extent %>% magrittr::extract(4)

NCol <- ncol(blank)
NRow <- nrow(blank)

Processed <- paste0("Iceberg Input Files/GretCDF/","Futures") %>% 
  list.files %>% str_remove(".rds$")

(ToProcess <- setdiff(Species, Processed)) %>% length

# ToProcess <- Species

if(CoryClimateReps[CR] == "gf"){
  
  PredReps <- PredReps[c(1, 2, 4)]
  
}

i = 1

print(paste0("To process:", length(ToProcess)))

mclapply(i:length(ToProcess), function(i){
  
  Sp <- ToProcess[i]
  
  print(Sp)
  
  # 02_Resampling rasters ####
  
  RasterLista <- lapply(PredReps[2:length(PredReps)], function(a){
    
    # print(a)
    
    Files[str_detect(Files, paste0("__", Sp, "__"))] %>% 
      str_replace("Climate1", ClimateReps[CR])  %>% 
      str_replace("Currents", a) ->
      
      SubFiles
    
    r1 <- raster(SubFiles)
    
    raster::resample(r1, blank, method = 'ngb')
    
  })
  
  names(RasterLista) <- PredReps[2:length(PredReps)]
  
  GretCDF <- data.frame(
    X = seq(from = XMin, to = XMax, length.out = NCol) %>% rep(NRow),
    Y = seq(from = YMax, to = YMin, length.out = NRow) %>% rep(each = NCol)
  ) 
  
  GretCDF[,paste0("Climate.",PredReps[2:length(PredReps)])] <- 
    RasterLista %>% 
    lapply(function(a) as.numeric(!is.na(values(a)))) %>% 
    bind_cols
  
  # Land Use filters #####
  
  SpHabitat <- HabitatList[[Sp]] %>% 
    as.character %>% 
    str_split("[.]") %>% 
    unlist %>% unique %>%
    na.omit
  
  lapply(PredReps[2:length(PredReps)], function(a){
    
    if(length(SpHabitat)>0){
      
      LandUseList[[a]][[SpHabitat]] %>% getValues -> 
        
        ValueDF
      
      if(length(SpHabitat)==1){
        
        Habitable <- as.numeric(ValueDF==1)
        Habitable[is.na(Habitable)] <- 0
        
      }else{
        
        Habitable <- as.numeric((ValueDF %>% rowSums(na.rm = T))>0)
        
      }
      
    }else{
      
      Habitable <- rep(1, nrow(GretCDF))
      
    }
    
  }) %>% bind_cols -> 
    GretCDF[,paste0("LandUse.",PredReps[2:length(PredReps)])]
  
  # Importing currents grid ####
  
  CurrentsGretCDF <- 
    readRDS(paste0("~/Albersnet/Iceberg Files/", 
                   "Climate1/Iceberg Input Files/GretCDF/Currents/", 
                   Sp, ".rds"))
  
  GretCDF <- GretCDF %>% slice(-Sea)
  
  GretCDF[,c("Continent",paste0("Buffer",c("Climate","ClimateLandUse")))] <- 
    
    CurrentsGretCDF[,c("Continent",paste0("Buffer",c("Climate","ClimateLandUse")))]
  
  GretCDF[GretCDF$Continent == 0, paste0("Climate.",PredReps[2:length(PredReps)])] <- 0
  
  GretCDF[,paste0("ClimateLandUse.",PredReps[2:length(PredReps)])] <- 
    
    lapply(PredReps[2:length(PredReps)], function(a){
      
      as.numeric(rowSums(GretCDF[,paste0(c("Climate.","LandUse."),a)])>1)
      
    }) %>% bind_cols
  
  
  PredReps[2:length(PredReps)] %>% lapply(function(a){
    
    List1 <- lapply(c("BufferClimate", "BufferClimateLandUse"), function(b){
      
      as.numeric(rowSums(GretCDF[,c(b,paste0(substr(b, 7, nchar(b)),".",a))])>1)
      
    })
    
    names(List1) <- paste0(c("BufferClimate", "BufferClimateLandUse"),".",a)
    
    return(List1)
    
  }) %>% bind_cols() -> FillDF
  
  GretCDF %>% bind_cols(FillDF) ->
    GretCDF
  
  GretCDF %>% as.matrix %>% as("dgCMatrix") %>% 
    
    saveRDS(file = paste0("Iceberg Input Files/GretCDF/Futures/", Sp, ".rds"))
  
}, mc.preschedule = F, mc.cores = CORES)

t2 = Sys.time()

t2 - t1
