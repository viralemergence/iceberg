
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

clim.ex <- c()

for (i in 1:length(raw.names)) {
  
  spnamei <- gsub('_',' ',raw.names[i])
  
  # FIRST CHECK FOR MANUAL DROP
  
  if(spnamei %in% drop.manual) {
    next
  } else { if (spnamei %in% spname) {
    
    input <- raster(raws.long[i])
    if (ncell(input[!is.na(input)])==0) {
      clim.ex <- c(clim.ex, spnamei)
    }
    print(i)
  }
  }
}

length(clim.ex)
