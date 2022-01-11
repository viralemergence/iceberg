
# 08_Main Figures ####

library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(SpRanger); library(raster)
library(sf); library(fasterize);library(ggregplot); library(igraph);library(maptools); library(fs)
library(magrittr); library(cowplot); library(patchwork); library(colorspace)

setwd(here::here())

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

NewEncounterList <- readRDS("~/Albersnet/Iceberg Files/CHELSA/Output Files/NewEncounterList.rds")

# Figure ####

# EncounterDF <- 
#   NewEncounterList %>% 
#   reduce(~full_join(.x, .y, by = c("Var1", "Var2")))

EncounterDF <-
  NewEncounterList %>% 
  map(nrow) %>% 
  as.data.frame %>% gather("Rep", "N") %>% 
  separate(Rep, into = c("Pipeline", "Rep"), sep = "[.]") %>% 
  mutate(GCM = substr(Rep, 1, 2),
         RCP = substr(Rep, 3, 4),
         Year = substr(Rep, 5, 6))

Pipelines <- 
  NewEncounterList %>% 
  names %>% str_split("[.]") %>% 
  map_chr(1) 

SubReps <- 
  NewEncounterList %>% 
  names %>% str_split("[.]") %>% 
  map_chr(2) 

substr(SubReps, 5, 6) <- ".."

SubReps <- paste0(Pipelines, ".", SubReps) %>% unique

CumEncounterDF <- 
  SubReps %>% 
  map(function(a){
    
    print(names(NewEncounterList[str_detect(names(NewEncounterList), a)]))
    
    SubList <- NewEncounterList[str_detect(names(NewEncounterList), a)]
    
    SubList %<>% map(~.x %>% dplyr::select(1:2))
    
    names(SubList) <- c("New11", "New41", "New71")
    
    data.frame(
      
      New91 = 0,
      New11 = SubList$New11 %>% nrow,
      New41 = SubList$New11 %>% bind_rows(SubList$New41) %>% unique %>% 
        nrow,
      New71 = SubList$New11 %>% bind_rows(SubList$New41) %>% bind_rows(SubList$New71) %>% unique %>% 
        nrow
      
    )
    
  }) %>% bind_rows(.id = "Rep")

CumEncounterDF$Rep <- SubReps

CumEncounterDF %<>% 
  # dplyr::select(-New91) %>% 
  gather("Year", "CumN", -Rep) %>% 
  mutate_at("Year", ~str_remove(.x, "New")) %>% 
  mutate_at("Rep", ~str_remove(.x, "..$")) %>% 
  separate(Rep, into = c("Pipeline", "Rep"), sep = "[.]") %>% 
  mutate(GCM = substr(Rep, 1, 2),
         RCP = substr(Rep, 3, 4))

CumEncounterDF %<>% 
  mutate_at("Year", ~.x %>% as.numeric %>% add(2000)) %>%
  mutate_at("Year", ~ifelse(.x > 2090, .x - 100, .x))

(CumEncounterDF %>% 
    filter(Year > 2000) %>% 
    mutate_at("Pipeline", ~str_replace_all(.x, c("LU" = " + LU",
                                                 "D" = " + D"))) %>% 
    ggplot(aes(Year, CumN, colour = RCP)) + 
    geom_line(aes(group = paste0(GCM, RCP))) +
    geom_point() +
    facet_wrap(~Pipeline, scales = "free")) /
  
  (EncounterDF %>% 
     ggplot(aes(Year, N, colour = RCP)) + 
     geom_line(aes(group = paste0(GCM, RCP))) +
     geom_point() +
     facet_wrap(~Pipeline, scales = "free"))

# Range changes ####

RangeChangeList <- 
  "Iceberg Files/CHELSA/Output Files/Areas" %>% dir_ls() %>% 
  map(readRDS)

RangeChangeDF <- 
  RangeChangeList %>% bind_rows(.id = "Rep") %>% 
  mutate_at("Rep", ~str_split(.x, "/") %>% map_chr(last) %>% 
              str_remove(".rds$") %>% 
              str_remove("Areas"))

Species <- RangeChangeDF$Sp %>% unique %>% sort

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

Panth1 %>% filter(hOrder == "Chiroptera") %>% pull(Sp) %>%
  intersect(Species) ->
  BatSpecies

CurrentClimateLandUseAreas <- 
  readRDS("~/Albersnet/Iceberg Files/CHELSA/Output Files/CurrentClimateLandUseAreas.rds") %>% 
  rename(CurrentArea = Area)

CurrentClimateAreas <- readRDS("~/Albersnet/Iceberg Files/CHELSA/Output Files/CurrentClimateAreas.rds") %>% 
  rename(CurrentArea = Area)

CurrentRangeDF <- 
  CurrentClimateLandUseAreas %>%
  mutate(Rep = "CLUCurrents") %>%
  bind_rows(CurrentClimateAreas %>%
              mutate(Rep = "CCurrents"))

RangeChangeDF %<>% 
  separate(Rep, into = c("Pipeline", "Rep"), sep = "[.]") %>% 
  mutate(GCM = substr(Rep, 1, 2),
         RCP = substr(Rep, 3, 4),
         Year = substr(Rep, 5, 6))

RangeChangeSummary <- 
  RangeChangeDF %>% 
  filter(!Sp %in% BatSpecies) %>%
  left_join(CurrentClimateAreas, by = "Sp") %>% 
  mutate(AbsoluteChange = Area - CurrentArea,
         PercentChange = (Area - CurrentArea)/CurrentArea) %>% 
  filter(!PercentChange == Inf) %>% 
  group_by(GCM, RCP, Year, Pipeline) %>% 
  summarise_at(c("AbsoluteChange", "PercentChange"), ~mean(.x, na.rm = T))

RangeChangeSummary %>% 
  ggplot(aes(Year, PercentChange, colour = RCP)) + 
  geom_hline(yintercept = 0, lty = 2, alpha = 0.3) +
  geom_line(aes(group = paste0(GCM, RCP))) +
  geom_point() +
  facet_wrap(~Pipeline)

# New Intersects plotting ####

paste0("~/Albersnet/Iceberg Files/CHELSA/Iceberg Input Files/GretCDF/Currents") %>% 
  list.files(full.names = T) ->
  CurrentFiles

paste0("~/Albersnet/Iceberg Files/CHELSA/Iceberg Input Files/GretCDF/Currents") %>% 
  list.files() %>% str_remove(".rds$") ->
  names(CurrentFiles)

FocalRep <- "CLUD.uk2671"
FocalRep <- "CLUD.uk8571"
FocalRep <- "CLU.gf2611"

FocalIntersects <- readRDS(paste0("~/Albersnet/Iceberg Files/CHELSA/Output Files/IntersectDFs/",
                                  FocalRep, ".NewIntersects.rds")) 

FocalIntersects[,"Fill"] <- FocalIntersects[,paste0("Overlap.", FocalRep)]

readRDS(CurrentFiles[[1]]) %>% 
  as.matrix %>% 
  as.data.frame() %>% 
  dplyr::select(X, Y) %>% 
  bind_cols(FocalIntersects) %>% 
  ggplot(aes(X, Y, fill = Fill)) + 
  geom_tile() + 
  coord_sf() + 
  scale_fill_continuous_sequential(palette = "Terrain") +
  labs(Fill = FocalRep)

"Iceberg Files/CHELSA/Output Files/IntersectDFs" %>% 
  dir_ls(regex = FocalRep %>% substr(1, nchar(.) - 2)) %>% 
  map(readRDS) %>% 
  # map(~.x %>% rename())
  # bind_rows(.id = "Rep") %>% 
  bind_cols() %>% 
  bind_cols(readRDS(CurrentFiles[[1]]) %>% 
              as.matrix %>% 
              as.data.frame() %>% 
              dplyr::select(X, Y)) %>% 
  gather("Rep", "Value", -c(X, Y)) %>% 
  filter(str_detect(Rep, "^NoBatOverlap[.]")) %>%
  ggplot(aes(X, Y, fill = Value)) + 
  geom_tile() +
  facet_grid(Rep~.) + 
  coord_fixed() +
  scale_fill_continuous_sequential(palette = "Terrain")



