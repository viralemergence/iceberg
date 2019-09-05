
# Iceberg pre-GAM spatial processing ####

# Rscript "02_Iceberg ENM Futures.R"

CORES <- 25

library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr)

Method = "MaxEnt"
# Method = "RangeBags"

# for(Method in c("RangeBags", "MaxEnt")){

paste0("Iceberg Input Files/", Method, "/GretCDF/Currents") %>% list.files() %>% 
  str_remove(".rds$") ->
  Species

PredReps <- c("Currents", paste0("Futures", 1:4))

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
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

LandUseList <- lapply(PredReps[2:5], function(a){
  
  brick(paste0("Iceberg Input Files/LandUse", a, ".grd"))
  
})

names(LandUseList) <- PredReps[2:5]

# Continents ####
print("Continents!")

ContinentRaster <- raster("Iceberg Input Files/continents-final.tif") %>%
  resample(blank, method = "ngb")

ContinentWhich <- lapply(1:5, function(a) which(values(ContinentRaster)==a))
names(ContinentWhich) <- c("Africa", "Eurasia", "NAm", "SAm", "Oceania")

# 01_Processing Currents ####

Root <- paste0("Iceberg Input Files/", Method,"/01_Raw")

# Setting raster standards ####

XMin <- blank %>% extent %>% magrittr::extract(1)
XMax <- blank %>% extent %>% magrittr::extract(2)
YMin <- blank %>% extent %>% magrittr::extract(3)
YMax <- blank %>% extent %>% magrittr::extract(4)

NCol <- ncol(blank)
NRow <- nrow(blank)

i = 1  

Processed <- paste0("Iceberg Input Files/", Method,"/GretCDF/","Futures") %>% 
  list.files %>% str_remove(".rds$")

(ToProcess <- setdiff(Species, Processed)) %>% length

mclapply(1:length(ToProcess), function(i){
  
  Sp <- ToProcess[i]
  
  print(Sp)
  
  # 02_Resampling rasters ####
  
  RasterLista <- lapply(PredReps[2:5], function(a){
    
    SubFiles <- paste0(Root,"/",a,"/",Sp) %>% list.files
    r1 <- raster(paste0(Root,"/",a,"/",Sp,"/",SubFiles[1]))
    raster::resample(r1, blank, method = 'ngb')
    
  })
  
  names(RasterLista) <- PredReps[2:5]
  
  GretCDF <- data.frame(
    X = seq(from = XMin, to = XMax, length.out = NCol) %>% rep(NRow),
    Y = seq(from = YMax, to = YMin, length.out = NRow) %>% rep(each = NCol)
  ) 
  
  GretCDF[,paste0("Climate.",PredReps[2:5])] <- 
    RasterLista %>% 
    lapply(function(a) as.numeric(!is.na(values(a)))) %>% 
    bind_cols
  
  # Land Use filters #####
  
  SpHabitat <- HabitatList[[Sp]] %>% 
    as.character %>% 
    str_split("[.]") %>% 
    unlist %>% unique %>%
    na.omit
  
  lapply(PredReps[2:5], function(a){
    
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
    GretCDF[,paste0("LandUse.",PredReps[2:5])]
  
  # Importing currents grid ####
  
  CurrentsGretCDF <- readRDS(paste0("Iceberg Input Files/", Method, "/GretCDF/Currents/", Sp, ".rds"))
  
  GretCDF <- GretCDF %>% slice(-Sea)
  
  GretCDF[,c("Continent",paste0("Buffer",c("Climate","ClimateLandUse")))] <- 
    
    CurrentsGretCDF[,c("Continent",paste0("Buffer",c("Climate","ClimateLandUse")))]
  
  GretCDF[GretCDF$Continent == 0, paste0("Climate.",PredReps[2:5])] <- 0
  
  GretCDF[,paste0("ClimateLandUse.",PredReps[2:5])] <- 
    lapply(PredReps[2:5], function(a){
      
      as.numeric(rowSums(GretCDF[,paste0(c("Climate.","LandUse."),a)])>1)
      
    }) %>% bind_cols
  
  
  PredReps[2:5] %>% lapply(function(a){
    
    List1 <- lapply(c("BufferClimate", "BufferClimateLandUse"), function(b){
      
      as.numeric(rowSums(GretCDF[,c(b,paste0(substr(b, 7, nchar(b)),".",a))])>1)
      
    })
    
    names(List1) <- paste0(c("BufferClimate", "BufferClimateLandUse"),".",a)
    
    return(List1)
    
  }) %>% bind_cols() -> FillDF
  
  GretCDF %>% bind_cols(FillDF) ->
    GretCDF
  
  GretCDF %>% as.matrix %>% as("dgCMatrix") %>% saveRDS(file = paste0("Iceberg Input Files/", Method,"/GretCDF/Futures/",Sp,".rds"))
  
}, mc.preschedule = F, mc.cores = CORES)

t2 = Sys.time()

t2 - t1

# }

Sp <- "Bos_javanicus"
Sp <- "Antilocapra_americana"

CDF <- readRDS(paste0("~/Albersnet/Iceberg Input Files/MaxEnt/GretCDF/Futures/",Sp,".rds")) %>%
  as.matrix %>% as.data.frame()

CDF %>% colSums()

CDF %>% dplyr::select(X,Y,ends_with("Futures1"),
                      "BufferClimate", "BufferClimateLandUse") %>%
  gather(key = "Key", value = "Value", -X, -Y) %>% 
  filter(Value>0) %>% dplyr::select(X, Y) %>% lapply(range) ->
  Lims 

CDF %>% 
  dplyr::select(X,Y,
                "BufferClimate", "BufferClimateLandUse",
                ends_with("Futures1")) %>%
  gather(key = "Key", value = "Value", -X, -Y) %>% 
  mutate(Key = factor(Key, levels = unique(Key))) %>%
  ggplot(aes(X, Y, fill = Value)) + 
  coord_fixed() + theme_cowplot() +
  lims(x = Lims$X, y = Lims$Y) +
  geom_tile() + facet_wrap(~Key)


