
# Rscript "Iceberg Code/07b_Summarising Maps.R"

library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(SpRanger); library(raster)
library(sf); library(fasterize);library(ggregplot); library(igraph);library(maptools); library(fs)
library(magrittr)

setwd(here::here())
setwd("Iceberg Files/CHELSA")

dir_create("Output Files/SummaryDFs")

if(!file.exists("Output Files/SummaryDFs/GridDF.rds")){
  
  GridDFList <- 
    "Output Files/GridDFs" %>%
    dir_ls %>% 
    map(readRDS) %>% 
    map(~.x %>% rename(Overlap = 3, Sharing = 4))
  
  names(GridDFList) <- 
    "Output Files/GridDFs" %>%
    list.files %>% str_remove(".rds$")
  
  UltraGridDF <- 
    GridDFList %>% 
    bind_rows(.id = "Key") %>%
    separate("Key", sep = "[.]",
             into = c("Pipeline", 
                      "Rep", 
                      "Bats")) %>% 
    mutate(GCM = substr(Rep, 1, 2),
           RCP = substr(Rep, 3, 4),
           Year = substr(Rep, 5, 6))
  
  UltraGridDF %<>% 
    group_by(Pipeline, RCP, Year, Bats, X, Y) %>% 
    summarise_at(c("Overlap", "Sharing"), ~mean(.x, na.rm = T))
  
  UltraGridDF %>% saveRDS("Output Files/SummaryDFs/GridDF.rds")
  
}

# Intersects ####

Locs <- 
  "Output Files/GridDFs" %>%
  dir_ls %>% 
  extract2(1) %>% 
  readRDS %>% 
  dplyr::select(X, Y)

if(!file.exists("Output Files/SummaryDFs/IntersectDF.rds")){
  
  IntersectDFList <- 
    "Output Files/IntersectDFs" %>%
    dir_ls %>% 
    map(readRDS) %>% 
    map(~.x %>% rename(Overlap = 1, Sharing = 2, 
                       BatOverlap = 3, BatSharing = 4)) %>% 
    map(~bind_cols(.x, Locs))
  
  names(IntersectDFList) <- 
    "Output Files/IntersectDFs" %>%
    list.files %>% str_remove(".rds$")
  
  UltraIntersectDF <- 
    IntersectDFList %>% 
    bind_rows(.id = "Key") %>%
    separate("Key", sep = "[.]",
             into = c("Pipeline", 
                      "Rep", 
                      "Lose")) %>% 
    mutate(GCM = substr(Rep, 1, 2),
           RCP = substr(Rep, 3, 4),
           Year = substr(Rep, 5, 6))
  
  UltraIntersectDF %<>% 
    group_by(Pipeline, RCP, Year, X, Y) %>% 
    summarise_at(c("Overlap", "Sharing", "BatOverlap", "BatSharing"), 
                 ~mean(.x, na.rm = T))
  
  UltraIntersectDF %>% saveRDS("Output Files/SummaryDFs/IntersectDF.rds")
  
}