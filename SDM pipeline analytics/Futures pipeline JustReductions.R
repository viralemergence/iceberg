setwd('C:/Users/cjcar/Dropbox/ViralChatter')
spname <- read.csv('FinalSpeciesList.csv')[,2] #2445


#### 
filenames <- list.files(path='D:/ICEBERG/RawENMs/PPM/BinaryMaps/', 
                    recursive= TRUE, 
                    include.dirs=FALSE,
                    full.names=FALSE)
splitter <- function(x) {  strsplit(x,'/')[[1]][1]  }
filenames <- unlist(lapply(filenames, splitter))

drop.manual <- c(as.character(spname[!(spname %in% filenames)])) # These two don't have niche models 
drop.manual <- gsub('_',' ',drop.manual)

source('lulcdrops.R')
drop.manual <- c(drop.manual,drop.marine) # this adds the two all-marine species to drops

noL.auto <- noL.auto # this came from the source

noL.manual <- c("Ardops_nichollsi","Calomys_boliviae","Cryptotis_mayensis",
                "Dasypus_sabanicola","Falsistrellus_mackenziei","Mesembriomys_macrurus",
                "Myotis_peninsularis","Nelsonia_neotomodon","Neotoma_phenax",
                "Nyctalus_azoreum","Peromyscus_simulus","Plecotus_teneriffae",
                "Pteropus_mariannus","Reithrodontomys_zacatecae","Sylvisorex_howelli") # These ones resampled weird or had no suitable area after LULC
noL.manual <- gsub("_"," ",noL.manual)

noL.all <- unique(c(noL.auto, noL.manual)) # Final list for no LULC conversion, just use currents

disp.names <- gsub('\\.tif','',gsub('_',' ',list.files('D:/ICEBERG/RawENMs/PPM/PostDispersal')))

noD.auto <- c(spname[!(spname %in% disp.names)]) # Species with no dispersal

###########################################################

# SO here's how this species pipeline needs to work
# We read in every species (2445) from spname
# If they're in drop.manual (4) we drop them
# If they're in noL.all (38) we don't land use filter the futures
# If they're in noD.auto (43) we don't clip them by dispersal
# Final outputs are written out into a folder for their pathway

##########################################################

ncin <- nc_open("D:/ssp5rcp85.nc")

get <- function(varname) {
  soilm_array <- ncvar_get(ncin,varname,start=c(1,1,56)) # 56 is what sets it to 2070
  slice <- raster(t(soilm_array[,,1]))
  raster::extent(slice) <- c(-180,180,-90,90)
  return(slice)
}

lutypes <- names(ncin$var)[1:12]

landuse2070 <- stack(lapply(lutypes,get))
names(landuse2070) <- lutypes
landuse2070 <- (landuse2070>0)
plot(landuse2070)

################# NEED TO REPLACE landuse2070 WITH A NEW THING

future.dir <- 'D:/ICEBERG/RawENMs/PPM/OtherEnvPred/he8570/BinaryMaps'
outdir <- 'D:/ICEBERG/RawENMs/PPM/OtherEnvPred/he8570 PROCESSED'

raws.long <- list.files(future.dir, 
                   recursive= TRUE, 
                   include.dirs=FALSE,
                   full.names=TRUE)
raw.names <- list.files(future.dir, 
                        recursive= TRUE, 
                        include.dirs=FALSE,
                        full.names=FALSE)
raw.names <- gsub('_',' ',unlist(lapply(raw.names, splitter)))

lu.red <- c()
disp.red <- c()

for (i in 1:length(raw.names)) {
  
  spnamei <- gsub('_',' ',raw.names[i])
  
  # FIRST CHECK FOR MANUAL DROP
  
  if(spnamei %in% drop.manual) {
    next
  } else { if (spnamei %in% spname) {
    
    input <- raster(raws.long[i])
    input <- raster::resample(input,landuse2070[[1]],method='ngb')
    
    # SECOND CHECK LAND USE
    
    if(spnamei %in% noL.all) {} else {
      lu.init <- ncell(input[!is.na(input)])
      LUHdati <- LUHdat[spnamei]
      allhabs <- max(landuse2070[[LUHdati[[1]]]])
      input <- (sum(input,allhabs,na.rm=TRUE)==2)
      input[!(input==1)] = NA
      lu.red <- c(lu.red, ncell(input[!is.na(input)])/lu.init)
    }
    
    # THIRD CHECK DISPERSAL
    
    if(spnamei %in% noD.auto) {} else {
      lu.init <- ncell(input[!is.na(input)])
      disp <- raster(paste(paste('D:/ICEBERG/RawENMs/PPM/PostDispersal',
                           gsub(' ','_',spnamei),sep='/'),'tif',sep='.')) 
      disp <- raster::resample(disp,landuse2070[[1]],method='ngb')
      input <- (sum(input,disp,na.rm=TRUE)==2)
      input[!(input==1)] = NA
      disp.red <- c(disp.red, ncell(input[!is.na(input)])/lu.init)
    }
  
  }}
  print(i)
}

mean(lu.red[!is.nan(lu.red)])
mean(disp.red[!is.nan(disp.red)])

reductions <- data.frame(lu.red, disp.red)
write.csv(reductions, 'reductions 8.5.csv')

1-mean(na.omit(reductions[,1]))
1-mean(na.omit(reductions[,2]))

sum(na.omit(reductions[,1])==0)
sum(na.omit(reductions[,2])==0)
