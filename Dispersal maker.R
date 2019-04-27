# Dispersal ####

library(ncdf4)
library(raster)
library(parallel)
library(tidyverse)

### fix blank

blank <- matrix(0,360*2,720*2)
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

# Dispersal script

Dispersals <- read.csv("Iceberg Input Files/Data for dispersal.csv", header = T)
disp <- Dispersals %>% filter(!is.na(Scientific_name), !is.na(disp50))
disp$Scientific_name <- disp$Scientific_name %>% str_replace(" ", "_")

# file.rename(paste0("Iceberg Input Files/Currents/",list.files("Iceberg Input Files/Currents")), paste0("Iceberg Input Files/Currents/",list.files("Iceberg Input Files/Currents") %>% str_replace(" ","")))

#PreDispList <- mclapply(list.files("Iceberg Input Files/Currents"), function(a) raster(paste0("Iceberg Input Files/Currents/",a)), mc.cores = 30)

#NoValues <- sapply(PreDispList, function(a){
#  
#  a %>% values %>% na.omit %>% length
#  
#})>0

RasterSp <- list.files("Iceberg Input Files/Currents") %>% #[-NoValues] %>% 
  str_split(".tif$") %>% map(1) %>% 
  sapply(function(a) a[1]) %>% sort %>% 
  intersect(disp$Scientific_name)

#NoValueSp <- list.files("Iceberg Input Files/Currents")[-NoValues] %>% 
#  str_split(".tif$") %>% map(1) %>% 
#  sapply(function(a) a[1]) %>% sort

NoDispSp <- list.files("Iceberg Input Files/Currents") %>% 
  str_split(".tif$") %>% map(1) %>% 
  sapply(function(a) a[1]) %>% sort %>% 
  setdiff(disp$Scientific_name)

all(RasterSp%in%disp$Scientific_name)
all(disp$Scientific_name %in% RasterSp)

t1 <- Sys.time()

i = 1
j = length(RasterSp)

if(length(list.files("PostDispersal"))==0){
  
  mclapply(RasterSp[i:j], function(a){
    
    testraster <- raster(paste0('Iceberg Input Files/Currents/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD
    
    testraster <- raster::resample(testraster, blank, method = 'ngb')
    if(1%in%unique(values(testraster))){
      dk <- disp[disp$Scientific_name == a, 'disp50']*1000
      buff <- buffer(testraster, dk)
      writeRaster(buff, filename = paste0("PostDispersal/",a,'.tif'))
    }
  }, mc.cores = 60)
  
}

# If not all of them worked, this will rerun them quickly! 

Success <- list.files("PostDispersal") %>% str_split(".tif$") %>% purrr::map(1) %>% unlist
Fail <- setdiff(RasterSp, Success)

mclapply(Fail, function(a){
  
  testraster <- raster(paste0('Iceberg Input Files/Currents/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD
  
  testraster <- raster::resample(testraster, blank, method = 'ngb')
  if(1%in%unique(values(testraster))){
    dk <- disp[disp$Scientific_name == a, 'disp50']*1000
    buff <- buffer(testraster, dk)
    writeRaster(buff, filename = paste0("PostDispersal/",a,'.tif'))
  }
}, mc.cores = 60)

stop()

a = sample(Success, 1)
r1 <- raster(paste0('Iceberg Input Files/Currents/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD
r2 <- raster(paste0('PostDispersal/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD

plot(r2)
plot(r1, add = T, fill = "red")

# Checking it's all worked ####
# Ironically this doesn't work lol ####

PreDispList <- lapply(list.files("Iceberg Input Files/Currents") %>% sort, function(a) raster(paste0("Iceberg Input Files/Currents/",a)))
PostDispList <- lapply(list.files("PostDispersal") %>% sort, function(a) raster(paste0("PostDispersal/",a)))

PreSizes <- sapply(PreDispList[which(list.files("Iceberg Input Files/Currents")%in%list.files("PostDispersal"))][1:500], function(a) length(na.omit(values(a))))
PostSizes <- sapply(PostDispList[1:500], function(a) length(na.omit(values(a))))

PreSizes == PostSizes
table(PreSizes == PostSizes)

# Parsing failures #####

FailList <- lapply(Fail, function(a){
  
  raster(paste0('Iceberg Input Files/Currents/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD
  
})

ShouldntFail <- sapply(FailList, function(a){
  
  a %>% values %>% na.omit %>% length
  
})


sapply(Fail[ShouldntFail>0], function(a){
  
  testraster <- raster(paste0('Iceberg Input Files/Currents/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD
  
  testraster <- raster::resample(testraster, blank, method = 'ngb')
  
  if(testraster %>% values %>% na.omit %>% length == 0) print("AAAAAAAAAAAA") else{
    
    dk <- disp[disp$Scientific_name == a, 'disp50']
    buff <- buffer(testraster, dk)
    writeRaster(buff, filename = paste0("PostDispersal/",a,'.tif'))
  }
})

paste0("Iceberg Input Files/Currents/",
       setdiff(list.files("Iceberg Input Files/Currents"), list.files("PostDispersal")) %>%
         intersect(paste0(RasterSp, ".tif"))
) %>% file.remove()

# Compare pre- and post-dispersal for a given species ####

Success <- list.files("PostDispersal") %>% str_split(".tif$") %>% purrr::map(1) %>% unlist

a = last(Success)
a

testraster1 <- raster(paste0('Iceberg Input Files/Currents/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD
testraster1 <- raster::resample(testraster1, blank, method = 'ngb')

testraster2 <- raster(paste0('PostDispersal/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD
testraster2 <- raster::resample(testraster2, blank, method = 'ngb')

length(na.omit(values(testraster1)))
length(na.omit(values(testraster2)))

par(mfrow = c(2,1))

plot(testraster1)
plot(testraster2)





