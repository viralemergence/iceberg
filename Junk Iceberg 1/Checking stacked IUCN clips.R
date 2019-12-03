
# Trying a stack 


paste0("Iceberg Input Files/", Method, "/02_Resampled/Currents") %>% list.files(full.names = T) %>%
  lapply(raster) -> 
  IUCNClippedRasters

paste0("Iceberg Input Files/", Method, "/02_Resampled/Currents") %>% list.files() %>% str_remove(".tif") ->
  names(IUCNClippedRasters)

OverlapSums <- OverlapSharingSums <- rep(0, ncol(IUCNClippedRasters[[1]])*nrow(IUCNClippedRasters[[1]]))

print("Getting richness")

for(i in i:length(IUCNClippedRasters)){
  
  if(i %% 500 == 0) print(i)
  SubSums <- raster::getValues(IUCNClippedRasters[[i]])
  SubSums[is.na(SubSums)] <- 0
  OverlapSums <- OverlapSums + SubSums
  
}


GridDF <- data.frame(
  Richness = OverlapSums,
  X = rep(1:ncol(IUCNClippedRasters[[1]]), nrow(IUCNClippedRasters[[1]])),
  Y = rep(nrow(IUCNClippedRasters[[1]]):1, each = ncol(IUCNClippedRasters[[1]]))
)

ggplot(GridDF, aes(X, Y, fill = Richness)) + geom_tile()

saveRDS(GridDF, file = "02_ResampledCurrents.rds")


# Trying a stack 


paste0("Iceberg Input Files/", Method, "/03_MaskedCurrents") %>% list.files(full.names = T) %>%
  lapply(raster) -> 
  IUCNClippedRasters

paste0("Iceberg Input Files/", Method, "/03_MaskedCurrents") %>% list.files() %>% str_remove(".tif") ->
  names(IUCNClippedRasters)

OverlapSums <- OverlapSharingSums <- rep(0, ncol(IUCNClippedRasters[[1]])*nrow(IUCNClippedRasters[[1]]))

print("Getting richness")

for(i in i:length(IUCNClippedRasters)){
  
  if(i %% 500 == 0) print(i)
  SubSums <- raster::getValues(IUCNClippedRasters[[i]])
  SubSums[is.na(SubSums)] <- 0
  OverlapSums <- OverlapSums + SubSums
  
}


GridDF <- data.frame(
  Richness = OverlapSums,
  X = rep(1:ncol(IUCNClippedRasters[[1]]), nrow(IUCNClippedRasters[[1]])),
  Y = rep(nrow(IUCNClippedRasters[[1]]):1, each = ncol(IUCNClippedRasters[[1]]))
)

ggplot(GridDF, aes(X, Y, fill = Richness)) + geom_tile()

