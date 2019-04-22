# Dispersal ####

library(ncdf4)
library(raster)
library(parallel)

### fix blank

blank <- matrix(0,360,720)
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

# Dispersal script

Dispersals <- read.csv("data/Data for dispersal.csv", header = T)
disp <- Dispersals %>% filter(!is.na(Scientific_name), !is.na(disp50))
disp$Scientific_name <- disp$Scientific_name %>% str_replace(" ", "_")

file.rename(paste0("IceMaps/",list.files("IceMaps")), paste0("IceMaps/",list.files("IceMaps") %>% str_replace(" ","")))

RasterSp <- list.files("IceMaps") %>% 
  str_split(".tif$") %>% map(1) %>% 
  sapply(function(a) a[1]) %>% sort %>% 
  intersect(disp$Scientific_name)

NoDispSp <- list.files("IceMaps") %>% 
  str_split(".tif$") %>% map(1) %>% 
  sapply(function(a) a[1]) %>% sort %>% 
  setdiff(disp$Scientific_name)

all(RasterSp%in%disp$Scientific_name)
all(disp$Scientific_name %in% RasterSp)

t1 <- Sys.time()

i = 1
j = length(RasterSp)

mclapply(RasterSp[i:j], function(a){
  
  testraster <- raster(paste0('IceMaps/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD
  
  testraster <- raster::resample(testraster, blank, method = 'ngb')
  dk <- disp[disp$Scientific_name == a, 'disp50']
  buff <- buffer(testraster, dk)
  writeRaster(buff, filename = paste0("PostDispersal/",a,'.tif'))
  
}, mc.cores = 70)

t2 <- Sys.time()

t2 - t1

# If not all of them worked, this will rerun them quickly! 

Success <- list.files("PostDispersal") %>% str_split(".tif$") %>% map(1) %>% unlist
Fail <- setdiff(RasterSp, Success)

mclapply(Fail, function(a){
  
  testraster <- raster(paste0('IceMaps/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD
  
  testraster <- raster::resample(testraster, blank, method = 'ngb')
  dk <- disp[disp$Scientific_name == a, 'disp50']
  buff <- buffer(testraster, dk)
  writeRaster(buff, filename = paste0("PostDispersal/",a,'.tif'))
  
}, mc.cores = 15)

a = sample(Success, 1)
r1 <- raster(paste0('IceMaps/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD
r2 <- raster(paste0('PostDispersal/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD

plot(r2)
plot(r1, add = T, fill = "red")

# Checking it's all worked ####
# Ironically this doesn't work lol ####

PreDispList <- lapply(list.files("IceMaps")[1:50], function(a) raster(paste0("IceMaps/",a)))
PostDispList <- lapply(list.files("PostDispersal"), function(a) raster(paste0("PostDispersal/",a)))

PreSizes <- sapply(PreDispList[1:50], function(a) length(na.omit(values(a))))
PostSizes <- sapply(PostDispList[1:50], function(a) length(na.omit(values(a))))

PreSizes>=PostSizes

# Parsing failures #####

FailList <- lapply(Fail, function(a){
  
  raster(paste0('IceMaps/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD
  
})

ShouldntFail <- sapply(FailList, function(a){
  
  a %>% values %>% na.omit %>% length
  
})


sapply(Fail[ShouldntFail>0], function(a){
  
  testraster <- raster(paste0('IceMaps/',a,'.tif')) # CAREFUL OF SPACE BEFORE PERIOD
  
  testraster <- raster::resample(testraster, blank, method = 'ngb')
  
  if(testraster %>% values %>% na.omit %>% length == 0) print("AAAAAAAAAAAA") else{
    
    dk <- disp[disp$Scientific_name == a, 'disp50']
    buff <- buffer(testraster, dk)
    writeRaster(buff, filename = paste0("PostDispersal/",a,'.tif'))
  }
})

paste0("IceMaps/",
       setdiff(list.files("IceMaps"), list.files("PostDispersal")) %>%
         intersect(paste0(RasterSp, ".tif"))
) %>% file.remove()



