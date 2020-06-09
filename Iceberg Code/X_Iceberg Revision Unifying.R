
# X_Iceberg Organisation ####

library(fs); library(tidyverse); library(magrittr)

setwd("~/Albersnet/")

ClimateReps <- glue::glue("Climate{1:9}")

"Iceberg Revision" %>% list.files(full.names = T) %>% extract2(1) %>%
  
  paste0("/PPM") %>% list.files(full.names = T) %>% extract2(2) %>% 
  list.files %>%
  substr(1, 2) %>% 
  unique -> 
  
  CoryClimateReps

names(CoryClimateReps) <- ClimateReps

paste0("~/Albersnet/Iceberg Files/", 
       "Climate1/Iceberg Input Files/GretCDF/Currents") %>% 
  list.files() %>% 
  str_remove(".rds$") ->
  Species

PredReps <- c("Currents", paste0("Futures", 1:4))

# Combining Mammal DF's ####

glue("Iceberg Files/{ClimateReps}/Iceberg Output Files/AllMammaldf.rds") %>% 
  map(readRDS) -> 
  AllMammalDFList

# SpJoin <- function(a, b) full_join(a, b, by = c("Sp", "Sp2"))

AllMammalDFList %>% 
  #reduce(SpJoin) -> 
  bind_rows(.id = "ClimateRep") ->
  AllAllMammalDF

AllAllMammalDF %>% 
  group_by(ClimateRep) %>% 
  summarise_at(vars(matches("^Sharing")), ~mean(.x, na.rm = T)) %>% 
  saveRDS("Iceberg Files/Summary/Iceberg Output Files/SummaryNetworkChanges2.rds")

saveRDS(AllAllMammalDF, file = "Iceberg Files/Summary/Iceberg Output Files/AllAllMammalDF.rds")

# Combining predictions ####

# setwd("~/Albersnet/")

SummariseFiles <- c(
  
  "NewIntersects.rds", 
  "FuturesGridDF.rds",
  "NoBatFuturesGridDF.rds",
  
  "BPNewIntersects.rds", 
  "EbolaNewIntersects.rds"
  
)[5]

# i <- "NewIntersects.rds"

for(i in SummariseFiles){
  
  print(i)
  
  glue::glue("Iceberg Files/{ClimateReps[1:9]}/Iceberg Output Files") %>% 
    map(function(a){
      
      Files <- list.files(a, full.names = T)
      
      Files <- Files[str_detect(Files, paste0("/", i))]
      
    }) %>%
    map(~readRDS(.x)) -> 
    List
  
  names(List) <- ClimateReps
  
  List %>% 
    bind_rows(.id = "ClimateRep") %>%
    mutate_at("ClimateRep", ~paste0("Climate", .x)) -> 
    
    LongList
  
  ClimateReps[1:9] %>% 
    map(function(a){
      
      RenameVars <- 
        List[[a]] %>% names %>% setdiff(c("X", "Y"))
      
      List[[a]] %>% 
        rename_at(RenameVars, ~paste0(a, ".", .x))
      
    }) %>% bind_cols -> 
    
    WideList
  
  LongList %>%
    saveRDS(paste0("Iceberg Files/Summary/Long", i))
  
  WideList %>%
    saveRDS(paste0("Iceberg Files/Summary/Wide", i))
  
  LongList %>%
    group_by(X, Y) %>% 
    summarise_all(~mean(.x, na.rm = T)) -> 
    
    MeanDF
  
  MeanDF %>%
    saveRDS(paste0("Iceberg Files/Summary/", i))
  
  LongList %>%
    group_by(X, Y) %>% 
    summarise_all(~min(.x, na.rm = T)) -> 
    
    MinDF
  
  MinDF %>%
    saveRDS(paste0("Iceberg Files/Summary/Lower", i))
  
  LongList %>%
    group_by(X, Y) %>% 
    summarise_all(~min(.x, na.rm = T)) -> 
    
    MaxDF
  
  MaxDF %>%
    saveRDS(paste0("Iceberg Files/Summary/Upper", i))
  
}


# Summarising the ebola list argh ####

SummariseFiles <- c(
  
  "EbolaGridList.rds"
  
)

# i <- "NewIntersects.rds"

for(i in SummariseFiles){
  
  print(i)
  
  glue::glue("Iceberg Files/{ClimateReps[1:9]}/Iceberg Output Files") %>% 
    map(function(a){
      
      Files <- list.files(a, full.names = T)
      
      Files <- Files[str_detect(Files, paste0("/", i))]
      
    }) %>%
    map(~readRDS(.x)) -> 
    List
  
  names(List) <- ClimateReps
  
  List %>% 
    bind_rows(.id = "ClimateRep") %>%
    mutate_at("ClimateRep", ~paste0("Climate", .x)) -> 
    
    LongList
  
  ClimateReps[1:9] %>% 
    map(function(a){
      
      RenameVars <- 
        List[[a]] %>% names %>% setdiff(c("X", "Y"))
      
      List[[a]] %>% 
        rename_at(RenameVars, ~paste0(a, ".", .x))
      
    }) %>% bind_cols -> 
    
    WideList
  
  LongList %>%
    saveRDS(paste0("Iceberg Files/Summary/Long", i))
  
  WideList %>%
    saveRDS(paste0("Iceberg Files/Summary/Wide", i))
  
  LongList %>%
    group_by(X, Y) %>% 
    summarise_all(~mean(.x, na.rm = T)) -> 
    
    MeanDF
  
  MeanDF %>%
    saveRDS(paste0("Iceberg Files/Summary/", i))
  
  remove(MeanDF)
  
  LongList %>%
    group_by(X, Y) %>% 
    summarise_all(~min(.x, na.rm = T)) -> 
    
    MinDF
  
  MinDF %>%
    saveRDS(paste0("Iceberg Files/Summary/Lower", i))
  
  LongList %>%
    group_by(X, Y) %>% 
    summarise_all(~min(.x, na.rm = T)) -> 
    
    MaxDF
  
  MaxDF %>%
    saveRDS(paste0("Iceberg Files/Summary/Upper", i))
  
}


# LongList %>% nrow -> Length
# split(1:Length, 1:length(ClimateReps))

# Summarising range changes ####

ClimateReps %>% map(~glue("Iceberg Files/{.x}/Iceberg Output Files/RangeChangeList.rds") %>% 
                      readRDS) ->
  RangeChangeListList

names(RangeChangeListList) <- ClimateReps

RangeChangeListList %>% unlist(recursive = F) %>% 
  bind_rows(.id = "ClimateRep") -> 
  
  FullRangeChangeDF

FullRangeChangeDF %<>% mutate(ClimateRep = ClimateRep %>% str_remove(paste0(".", Sp)))

FullRangeChangeDF %>% saveRDS(file = "Iceberg Files/Summary/Iceberg Output Files/FullRangeChangeDF.rds")
