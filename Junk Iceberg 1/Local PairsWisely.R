# Local PairsWisely

PairsWisely <- function(Rasterstack, Species = "All", Area = F){
  
  library(raster)
  library(tidyverse)
  library(Matrix)
  
  t1 <- Sys.time()
  print("Getting the grid values")
  
  if(class(Rasterstack)=="RasterBrick"){
    
    Valuedf <- data.frame(raster::getValues(Rasterstack)) %>% as.matrix
    
  }
  
  if(class(Rasterstack)=="list"){
    
    Valuedf <- lapply(1:length(Rasterstack), function(a){
      
      getValues(Rasterstack[[a]])
      
    }) %>% bind_cols %>% as.data.frame()
    
  }
  
  if(class(Rasterstack)=="data.frame"){
    
    Valuedf <- Rasterstack
    
  }
  
  if(Area){AreaFun <- function(a){sum(a, na.rm = T)} } else {AreaFun <- function(a){length(which(a>0))}}
  
  colnames(Valuedf) <- names(Rasterstack)
  
  remove(Rasterstack)
  
  Valuedf[is.na(Valuedf)] <- 0
  Valuedf <- Valuedf %>% as.matrix() %>% as("dgCMatrix")
  
  print(paste0("Data frame size = ", dim(Valuedf)))
  
  if (Species != "All"){
    Valuedf <- Valuedf[, Species]
  }
  
  if (any(Matrix::colSums(Valuedf) == 0)) {
    print("Removing some species with no ranging data :(")
    Valuedf <- Valuedf[, -which(Matrix::colSums(Valuedf) == 0)]
  }
  
  print(paste0("Data frame size = ", dim(Valuedf)))
  
  RangeOverlap <- matrix(NA, nrow = ncol(Valuedf), ncol = ncol(Valuedf))
  dimnames(RangeOverlap) <- list(colnames(Valuedf), colnames(Valuedf))
  print("Calculating Overlap")
  for (x in 1:(ncol(Valuedf) - 1)) {
    print(colnames(Valuedf)[x])
    TrainGrids <- Valuedf[, x]
    SubRangedf <- Valuedf[which(TrainGrids > 0), x:ncol(Valuedf)]
    if (!is.null(dim(SubRangedf))) {
      RangeOverlap[x, x:ncol(Valuedf)] <- apply(SubRangedf, 2, AreaFun)
    }
    else RangeOverlap[x, x:ncol(Valuedf)] <- sapply(SubRangedf, AreaFun)
  }
  
  x = ncol(Valuedf)
  print(colnames(Valuedf)[x])
  TrainGrids <- Valuedf[, x]
  SubRangedf <- Valuedf[which(TrainGrids>0), x:ncol(Valuedf)]
  RangeOverlap[x, x:ncol(Valuedf)] <- AreaFun(SubRangedf)
  FullRangeA = matrix(rep(diag(RangeOverlap), nrow(RangeOverlap)),
                      nrow(RangeOverlap))
  FullRangeB = matrix(rep(diag(RangeOverlap), each = nrow(RangeOverlap)),
                      nrow(RangeOverlap))
  RangeAdj <- RangeOverlap/(FullRangeA + FullRangeB - RangeOverlap)
  TimeTaken = Sys.time() - t1
  print(TimeTaken)
  diag(RangeAdj) <- NA
  RangeAdj[lower.tri(RangeAdj)] <- t(RangeAdj)[!is.na(t(RangeAdj))]
  diag(RangeAdj) <- 1
  
  RangeAdj
}
