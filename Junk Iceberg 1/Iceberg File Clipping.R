
# Removing raster files I don't need #### 

library(tidyverse)

Method = "MaxEnt"

FolderNames <- c("03_MaskedCurrents","04_LandUseClipped","05_DispersalBuffers","06_DispersalBuffersResampled","07_DispersalClippedFutures")

for(Folder in FolderNames){
  
  Pass1 <- "Iceberg Input Files/" %>% paste0(Method, "/", Folder) %>% list.files(full.names = T)
  
  if(str_ends(Pass1, ".tif") %>% all){
    
    file.remove(Pass1)
    
  }else{
    
    lapply(Pass1, function(a){
      
      a %>% list.files(full.names = T) %>% file.remove()
      
    })
    
  }

}
