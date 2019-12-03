
# Comparing different projections ####

library(velox);
library(sf); library(fasterize); library(Matrix);library(ggplot2);
library(ggregplot); library(raster); library(tidyverse); library(igraph); 
library(maptools); library(SpRanger); library(parallel)

x = 1

Mercator = readRDS("Iceberg Output Files/IcebergAdjList_MercatorProj.rds")
Mercator_NoArea = readRDS("Iceberg Output Files/IcebergAdjList_MercatorProj_NoArea.rds")
MercatorUnc = readRDS("Iceberg Output Files/IcebergAdjList_MercatorProjUnc.rds")
MercatorUnc_NoArea = readRDS("Iceberg Output Files/IcebergAdjList_MercatorProjUnc_NoArea.rds")
Mollweide = readRDS("Iceberg Output Files/IcebergAdjList_MollweideProj.rds")
Mollweide_NoArea = readRDS("Iceberg Output Files/IcebergAdjList_MollweideProj_NoArea.rds")

ProjCompdf <- data.frame(
  
  Mercator = c(Mercator[[1]][lower.tri(Mercator[[1]])]),
  Mercator_NoArea = c(Mercator_NoArea[[1]][lower.tri(Mercator[[1]])]),
  MercatorUnc = c(MercatorUnc[[1]][lower.tri(Mercator[[1]])]),
  MercatorUnc_NoArea = c(MercatorUnc_NoArea[[1]][lower.tri(Mercator[[1]])]),
  Mollweide = c(Mollweide[[1]][lower.tri(Mollweide[[1]])]),
  Mollweide_NoArea = c(Mollweide_NoArea[[1]][lower.tri(Mollweide_NoArea[[1]])])
  
)

Namesdf = reshape2::melt(Mercator[[1]])[lower.tri(Mercator[[1]]),c("Var2","Var1")]
Names = paste(Namesdf$Var1, Namesdf$Var2, sep = ".")
ProjCompdf$Spp = Names
ProjCompdf$Discrepancy = with(ProjCompdf, Mollweide - Mollweide_NoArea)

qplot(ProjCompdf$Discrepancy)

ProjCompdf %>% filter(Discrepancy == max(Discrepancy))

a = "Dipodomys_ingens"
b = "Callorhinus_ursinus"
x = 1
y = 1

ProjReps <- c("MollweideProj","MercatorProj","MercatorProjUnc")

DropSpecies <- c("Lobodon carcinophaga","Leptonychotes weddellii","Enhydra lutris","Histriophoca fasciata") %>% 
  str_replace(" ","_")

PredReps <- c("Currents", paste0("Futures", 1:4))

Files <- list.files(paste0("Iceberg Output Files/",ProjReps[y],"/",PredReps[x])) %>% setdiff(DropSpecies)

Files <- Files[which((Files %>% str_remove(".tif"))%in%c(a,b))]

VeloxList <- lapply(Files, function(a){
  if(which(Files==a) %% 500==0) print(a)
  r1 <- velox(paste(paste0("Iceberg Output Files/",ProjReps[y],"/",PredReps[x]), a, sep = '/'))
})

RasterLista <- lapply(1:length(VeloxList), function(a){
  
  if(a %% 500==0) print(Files[a])
  
  VeloxList[[a]]$as.RasterLayer(band = 1) #%>% rasterToPolygons(dissolve = T)
  
})


AreaFun <- function(a) {
  sum(a, na.rm = T)
}

Valuedf = lapply(RasterLista, function(a) raster::getValues(a)) %>% bind_cols

colnames(Valuedf) <- c(a,b)
Valuedf[is.na(Valuedf)] <- 0
Valuedf <- Valuedf %>% as.matrix() %>% as("dgCMatrix")
print(paste0("Data frame size = ", dim(Valuedf)))
print(paste0("Data frame size = ", dim(Valuedf)))
RangeOverlap <- matrix(NA, nrow = ncol(Valuedf), ncol = ncol(Valuedf))
dimnames(RangeOverlap) <- list(colnames(Valuedf), colnames(Valuedf))
print("Calculating Overlap")
for (x in 1:(ncol(Valuedf) - 1)) {
  print(colnames(Valuedf)[x])
  TrainGrids <- Valuedf[, x]
  SubRangedf <- Valuedf[which(TrainGrids > 0), x:ncol(Valuedf)]
  if (!is.null(dim(SubRangedf))) {
    RangeOverlap[x, x:ncol(Valuedf)] <- apply(SubRangedf, 
                                              2, AreaFun)
  }
  else RangeOverlap[x, x:ncol(Valuedf)] <- sapply(SubRangedf, 
                                                  AreaFun)
}
x = ncol(Valuedf)
print(colnames(Valuedf)[x])
TrainGrids <- Valuedf[, x]
SubRangedf <- Valuedf[which(TrainGrids > 0), x:ncol(Valuedf)]
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

# Trying without Area

AreaFun <- function(a) {
  length(which(a > 0))
}

Valuedf = lapply(RasterLista, function(a) raster::getValues(a)) %>% bind_cols

colnames(Valuedf) <- c(a,b)
Valuedf[is.na(Valuedf)] <- 0
Valuedf <- Valuedf %>% as.matrix() %>% as("dgCMatrix")
print(paste0("Data frame size = ", dim(Valuedf)))
print(paste0("Data frame size = ", dim(Valuedf)))
RangeOverlap <- matrix(NA, nrow = ncol(Valuedf), ncol = ncol(Valuedf))
dimnames(RangeOverlap) <- list(colnames(Valuedf), colnames(Valuedf))
print("Calculating Overlap")
for (x in 1:(ncol(Valuedf) - 1)) {
  print(colnames(Valuedf)[x])
  TrainGrids <- Valuedf[, x]
  SubRangedf <- Valuedf[which(TrainGrids > 0), x:ncol(Valuedf)]
  if (!is.null(dim(SubRangedf))) {
    RangeOverlap[x, x:ncol(Valuedf)] <- apply(SubRangedf, 
                                              2, AreaFun)
  }
  else RangeOverlap[x, x:ncol(Valuedf)] <- sapply(SubRangedf, 
                                                  AreaFun)
}
x = ncol(Valuedf)
print(colnames(Valuedf)[x])
TrainGrids <- Valuedf[, x]
SubRangedf <- Valuedf[which(TrainGrids > 0), x:ncol(Valuedf)]
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


ProjCompdf <- data.frame(
  
  Mercator = lapply(Mercator, function(a) lapply(a, function(b) c(b[lower.tri(b)]))),
  Mercator_NoArea = lapply(Mercator_NoArea, function(a) lapply(a, function(b) c(b[lower.tri(b)]))),
  MercatorUnc = lapply(MercatorUnc, function(a) lapply(a, function(b) c(b[lower.tri(b)]))),
  MercatorUnc_NoArea = lapply(MercatorUnc_NoArea, function(a) lapply(a, function(b) c(b[lower.tri(b)]))),
  Mollweide = lapply(Mollweide, function(a) lapply(a, function(b) c(b[lower.tri(b)]))),
  Mollweide_NoArea = lapply(Mollweide_NoArea, function(a) lapply(a, function(b) c(b[lower.tri(b)])))
  
)

ProjCompdf = ProjCompdf[which(lower.tri((readRDS("Iceberg Output Files/IcebergAdjList_MercatorProj.rds")[[x]]))),]

ProjCompdf %>% select(-Spp) %>% slice(1:10000) %>% GGally::ggpairs(lower = list(continuous = "smooth"))
