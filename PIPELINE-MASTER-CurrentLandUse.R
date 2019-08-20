
library(ncdf4)
library(raster)
setwd("D:/ICEBERG")

ncin <- nc_open("states.nc")

get <- function(varname) {
  soilm_array <- ncvar_get(ncin,varname,start=c(1,1,1166))
  slice <- raster(t(soilm_array))
  raster::extent(slice) <- c(-180,180,-90,90)
  return(slice)
}

lutypes <- names(ncin$var)[1:12]

landuse2017 <- stack(lapply(lutypes,get))
names(landuse2017) <- lutypes
landuse2017 <- (landuse2017>0)
plot(landuse2017)

###################################

#IUCN grabber: run back half of iucn habitat
#####################

# 1. read in the list of species we're including
setwd('C:/Users/cjcar/Dropbox/ViralChatter')
spname <- gsub('_',' ',read.csv('FinalSpeciesList.csv')[,2])

# 2. read in their habitat types document

iucndat <- read.csv('Land cover data/IucnHabitatData.csv')
LUHdat<-as.list(rep(NA,length(spname))) 
convcodes <- read.csv('Land cover data/IUCN_LUH_conversion_table.csv')

for (i in 1:length(LUHdat)){
  
  #get iucn habitat codes for the species  
  iucncodes<-iucndat[iucndat$name %in% spname[i],"code"]
  
  if(length(iucncodes)==0) {
    LUHdat[[i]]<-NULL
    names(LUHdat)[[i]]<-as.character(spname[i])
    next
  }
  #get matching LUH land use types  
  LUH<-unique(convcodes$LUH[match(as.numeric(as.character(iucncodes)),convcodes$IUCN_hab)])
  
  # when there is no match make LUH = NA
  # A single NA value means the species lives only in IUCN categories that have no 
  # mapping to LUH data (e.g., caves or ocean habitats)
  suppressWarnings(if (is.na(LUH) & length(LUH)==1) {
    LUHdat[[i]]<-NA
    names(LUHdat)[[i]]<-as.character(spname[i])
    next
  })
  # when there is a match    
  LUHvector<-unlist(strsplit(as.character(LUH),".", fixed=TRUE))
  
  LUHdat[[i]]<-LUHvector
  names(LUHdat)[[i]]<-as.character(spname[i])
}

########################################################

# 3. create a for loop

propred <- c()
drops <- c()

#length(spname)
for (i in 1:length(spname)) {
  
  spnamei <- spname[i]
  wdir <- paste('D:/ICEBERG/RawENMs/PPM/BinaryMaps/',
                gsub(' ','_',spnamei),sep='')
  
  if(length(list.files(wdir))==0) { next } else {
  cur <- raster(paste(wdir,list.files(wdir),sep='/'))
# 4. inside the for loop do the following
# 4a. test for NA, skip those
  
LUHi <- LUHdat[spnamei]

if (all(is.na(LUHi))) { 
  curlu <- cur
  writeRaster(curlu, paste(paste('D:/ICEBERG/RawENMs/PPM/BinaryLU/',
                                 gsub(' ','_',spnamei),sep=''),'.tif',''))
  print(i)
} else { 
  LUHi[[1]] <- LUHi[[1]][!is.na(LUHi[[1]])]
  if (all(LUHi[[1]]=='MARINE')) {
    drops <- c(drops,spnamei)
  } else {
    
# 4b. pull habitat uses if not
# 5. add together those binaries from landuse2017
    
    allhabs <- max(landuse2017[[LUHi[[1]]]])
    
    # note: max(landuse2017[[1]],landuse2017[[2]])
    # 6. filter
    
    allhabs=raster::resample(allhabs,cur,method='ngb')
    suitable.cont <- (sum(cur,allhabs,na.rm=TRUE)==2)
    suitable.cont[!(suitable.cont==1)] = NA
    
    propred <- c(propred,ncell(suitable.cont[suitable.cont==1])/ncell(cur[cur==1]))
    
    curlu <- suitable.cont
  }

# 7. write new raster "curlu" to hard drive

writeRaster(curlu, paste(paste('D:/ICEBERG/RawENMs/PPM/BinaryLU/',
              gsub(' ','_',spnamei),sep=''),'.tif',''))
print(i)
  }
  }
}

hist(propred, main="Current range remaining after land use filter",
     xlab='Proportion of range')
1-mean(propred)
