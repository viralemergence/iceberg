
# Cleaning the final ip Data ####

library(tidyverse); library(fs)

setwd(here::here())

Files <- "IP" %>% list.files(full.names = T)# %>% substr(1,2)

FilesDF <- data.frame(
  
  FileName = Files,
  Name = Files %>% str_split("/") %>% map_chr(last) %>% str_remove(".zip")
  
) %>% 
  mutate(GCM = substr(Name, 1, 2),
         RCP = substr(Name, 3, 4),
         Year = substr(Name, 5, 6))

FilesDF

dir_create("Iceberg Files/CHELSA/IP Rasters")

library(zip)

i <- 1

for(i in 1:length(Files)){
  
  print(i)
  
  unzip(Files[i], 
        exdir = paste0("IP Rasters/", 
                       FilesDF[i, "Name"]))
  
}

# Directories <- "New Iceberg Files" %>% dir_ls(recurse = T)

FullFiles <- "IP Rasters" %>% dir_ls(recurse = 4, regex = ".tif$")

dir_create("Final IP Rasters")

From <- FullFiles[str_detect(FullFiles, "TP05|X0.165")]

To <-
  FullFiles[str_detect(FullFiles, "TP05|X0.165")] %>% 
    str_split("/") %>% 
    map_chr(last) %>% 
    str_split("__") %>% 
    map_chr(2) %>% 
    paste0("Final IP Rasters/", ., ".tif")

file_copy(From, To)

library(magrittr)

# FilesDF %<>% mutate_at("Method", ~str_replace_all(.x, "PP", "PPM"))

# FilesDF %<>% 
#   mutate(FilePath = paste0("New Iceberg Files/", Name, 
#                            "/Users/ctg/Documents/SDMs/ViralChatter/Revision1/VC_4_22_20_outputs/", Method, "/OtherEnvPred/", GCM, RCP, "/BinaryMaps"))

FilesDF %<>% 
  mutate(FilePath = paste0("Iceberg Files/CHELSA/Intermediate Rasters/", Name, 
                           "/PPM/OtherEnvPred/", Name, "/BinaryMaps"))

FilesDF$FilePath %<>%
  str_replace_all("present/PPM/OtherEnvPred/present/BinaryMaps", "present/PPM/BinaryMaps")

FilesDF$FilePath %>% 
  dir.exists()

FilesDF$FilePath %>% 
  str_replace_all("PPM", "RangeBag") %>%
  dir.exists()

# FilesDF %<>% 
#   mutate_at("FilePath", ~str_replace_all(.x, "/presen/", paste0("/", Name, "/")))

dir_create("Iceberg Files/CHELSA/Rasters")

FilesDF %<>% 
  mutate(NewPath = paste0("Iceberg Files/CHELSA/Rasters/", GCM, RCP, Year))

# 1:(nrow(FilesDF) - 2) %>% map(function(a){

a <- 1

for(a in 1:(nrow(FilesDF))){
  
  print(a)
  
  file.rename(FilesDF[a, "FilePath"], 
              FilesDF[a, "NewPath"] %>% paste0(., "_PPM"))
  
  file.rename(FilesDF[a, "FilePath"] %>% str_replace_all("PPM", "RangeBag"), 
              FilesDF[a, "NewPath"] %>% str_replace_all("PPM", "RangeBag") %>% paste0(., "_RB"))
  
  
  
}

# Did something manually there ####

CoryClimateReps <- paste0("Climate", 1:5)

names(CoryClimateReps) <- FilesDF$GCM %>% unique %>% sort %>% setdiff("pr")

dir_create("Iceberg Files/CHELSA/FinalRasters")

CopyFolders <- "Iceberg Files/CHELSA/Rasters" %>% 
  list.files(full.names = T) %>% 
  str_split("_") %>% map_chr(1) %>% unique %>% sort

i <- 1

for(i in (1:(length(CopyFolders) - 1)) %>% setdiff(36)){
  
  print(i)
  
  CopyFolders[i] %>% str_replace_all("Rasters", "FinalRasters") %>% dir_create
  
  CopyFolders[i] %>% paste0("_PPM") %>% list.files(full.names = T) %>% 
    map(function(a){
      
      print(a)
      
      b <- a %>% list.files(full.names = T) %>% extract2(2)
      
      file.rename(b, a %>% str_replace_all("Rasters", "FinalRasters") %>% 
                    str_remove("_PPM") %>% paste0(".tif"))
      
    })
  
  CopyFolders[i] %>% paste0("_RB") %>% list.files(full.names = T) %>% 
    map(function(a){
      
      b <- a %>% list.files(full.names = T) %>% extract2(2)
      
      file.rename(b, a %>% str_replace_all("Rasters", "FinalRasters") %>% 
                    str_remove("_RB") %>% paste0(".tif"))
      
    })
}

i <- 36

CopyFolders[i] %>% str_replace_all("Rasters", "FinalRasters") %>% dir_create

ToDelete <- CopyFolders[i] %>% paste0("_PPM") %>% list.files(full.names = T) %>% 
  map_dbl(function(a){
    
    a %>% list.files %>% length
    
  })

CopyFolders[i] %>% paste0("_PPM") %>% list.files(full.names = T) %>% 
  extract(which(ToDelete < 3)) %>% 
  file_delete()

CopyFolders[i] %>% paste0("_PPM") %>% list.files(full.names = T) %>% 
  map(function(a){
    
    print(a)
    
    b <- a %>% list.files(full.names = T) %>% extract2(2)
    
    file.rename(b, a %>% str_replace_all("Rasters", "FinalRasters") %>% 
                  str_remove("_PPM") %>% paste0(".tif"))
    
  })

ToDelete <- CopyFolders[i] %>% paste0("_RB") %>% list.files(full.names = T) %>% 
  map_dbl(function(a){
    
    a %>% list.files %>% length
    
  })

CopyFolders[i] %>% paste0("_RB") %>% list.files(full.names = T) %>% 
  extract(which(ToDelete < 3)) %>% 
  file_delete()

CopyFolders[i] %>% paste0("_RB") %>% list.files(full.names = T) %>% 
  map(function(a){
    
    b <- a %>% list.files(full.names = T) %>% extract2(2)
    
    file.rename(b, a %>% str_replace_all("Rasters", "FinalRasters") %>% 
                  str_remove("_RB") %>% paste0(".tif"))
    
  })

# Not sure what all this is ####

# "Iceberg Files/MaxEnt" %>% 
#   dir_ls %>% 
#   map(function(a){
#     
#     file.rename(a, a %>% str_remove_all("PPM$"))
#     
#   })
# 
# "Iceberg Files/RangeBag" %>% 
#   dir_ls %>% 
#   map(function(a){
#     
#     file.rename(a, a %>% str_remove_all("RB$"))
#     
#   })
# 
# "Iceberg Files/MaxEnt/" %>% paste0(CoryClimateReps)
# 
# i <- 1
# 
# for(i in 2:length(CoryClimateReps)){
#   
#   print(i)
#   
#   OldFiles <- "Iceberg Files/MaxEnt" %>% 
#     dir_ls(regexp = paste0("/", names(CoryClimateReps[i])))
#   
#   NewFiles <- OldFiles %>% str_replace_all("/MaxEnt/", paste0("/MaxEnt/", CoryClimateReps[i], "/"))
#   
#   dir_create(paste0("Iceberg Files/MaxEnt/", CoryClimateReps[i]))
#   
#   file.rename(OldFiles %>% as.character, NewFiles)
#   
# }
# 
# for(i in 1:length(CoryClimateReps)){
#   
#   print(i)
#   
#   OldFiles <- "Iceberg Files/RangeBag" %>% 
#     dir_ls(regexp = paste0("/", names(CoryClimateReps[i])))
#   
#   NewFiles <- OldFiles %>% str_replace_all("/RangeBag/", paste0("/RangeBag/", CoryClimateReps[i], "/"))
#   
#   dir_create(paste0("Iceberg Files/RangeBag/", CoryClimateReps[i]))
#   
#   file.rename(OldFiles %>% as.character, NewFiles)
#   
# }
# 
# Replace <- paste0("/Futures", 1:4)
# 
# names(Replace) <- c("2670", "4570", "7070", "8570")
# 
# "Iceberg Files" %>% dir_ls(recurse = 3) %>% 
#   str_remove_all(paste0(paste0("/", names(CoryClimateReps)), collapse = "|")) %>% 
#   str_replace_all(Replace) -> NewFiles
# 
# "Iceberg Files" %>% dir_ls(recurse = 3) %>% 
#   file.rename(NewFiles)
# 
# "Iceberg Files" %>% dir_ls(recurse = 1)