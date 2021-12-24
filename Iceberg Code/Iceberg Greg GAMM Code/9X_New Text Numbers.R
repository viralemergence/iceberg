
# 9X_New Text Numbers

library(tidyverse); library(fs); library(magrittr); library(glue)
library(ggregplot)

setwd("~/Albersnet/")

"Iceberg Files/Summary" %>% list.files(full.names = T)

# Setting up ####

Unlist1 <- function(a) unlist(a, recursive = F)

MeanSD <- function(a) list(Mean = mean(a, na.rm = T), 
                           SD = sd(a, na.rm = T))

glue::glue("Climate{1:9}") ->
  
  ClimateReps

# Importing New Encounters ####

glue("~/Albersnet/Iceberg Files/{ClimateReps}/Iceberg Output Files/NewEncounters.rds") %>% 
  map(readRDS) -> 
  NewEncountersListList

names(NewEncountersListList) <- ClimateReps

NewEncountersListList %<>% 
  map(Unlist1) %>% 
  Unlist1

NewEncountersListList %>% names 

# Importing Ebola ####

glue("Iceberg Files/{ClimateReps}/Iceberg Output Files/EbolaEncounters.rds") %>% 
  map(readRDS) -> 
  EbolaEncountersListList

names(EbolaEncountersListList) <- ClimateReps

EbolaEncountersListList %<>% 
  #map(Unlist1) %>% 
  Unlist1

EbolaEncountersListList %>% names 

# Importing BP ####

glue("Iceberg Files/{ClimateReps}/Iceberg Output Files/BPEncounters.rds") %>% 
  map(readRDS) -> 
  BPEncountersListList

names(BPEncountersListList) <- ClimateReps

BPEncountersListList %<>% 
  #map(Unlist1) %>% 
  Unlist1

BPEncountersListList %>% names 

# Loading All All Mammal DF

#AllAllMammalDF <- 
#  readRDS("~/Albersnet/Iceberg Files/Summary/Iceberg Output Files/AllAllMammalDF.rds")

SummaryNetworkChanges <- 
  readRDS("~/Albersnet/Iceberg Files/Summary/Iceberg Output Files/SummaryNetworkChanges2.rds")

# sink("Results.txt")

# Colin's Email ####

# ~~~~ Abstract ####

# Ballpark average number of first encounters with and without dispersal 
# (analogous to 3,000 to 13,000)

# New encounters ####
# With bats ####

NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures")] %>% 
  map_dbl(nrow) %>% 
  mean

NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures")] %>% 
  map_dbl(nrow) %>% 
  mean

# Without bats ####

NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures")] %>% 
  map_dbl(~.x %>% filter(!hOrder.x == "Chiroptera", 
                         !hOrder.y == "Chiroptera") %>% nrow) %>% 
  mean

NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures")] %>% 
  map_dbl(~.x %>% filter(!hOrder.x == "Chiroptera", 
                         !hOrder.y == "Chiroptera") %>% nrow) %>% 
  mean

# New sharings ####
# With bats ####

NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures")] %>% 
  map(~.x %>% summarise_at(vars(matches("DeltaSharing.Futures[1-4]A")), 
                           function(a) sum(a, na.rm = T))) %>% 
  bind_rows() %>% unlist %>% 
  mean(na.rm = T)

NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures")] %>% 
  map(~.x %>% summarise_at(vars(matches("DeltaSharing.Futures[1-4]C")), 
                           function(a) sum(a, na.rm = T))) %>% 
  bind_rows() %>% unlist %>% 
  mean(na.rm = T)

# Without bats ####

NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures")] %>% 
  map(~.x %>% filter(!hOrder.x == "Chiroptera", !hOrder.y == "Chiroptera") %>% 
        summarise_at(vars(matches("DeltaSharing.Futures[1-4]A")), 
                     function(a) sum(a, na.rm = T))) %>% 
  bind_rows() %>% unlist %>% 
  mean(na.rm = T)

NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures")] %>% 
  map(~.x %>% filter(!hOrder.x == "Chiroptera", !hOrder.y == "Chiroptera") %>% 
        summarise_at(vars(matches("DeltaSharing.Futures[1-4]C")), 
                     function(a) sum(a, na.rm = T))) %>% 
  bind_rows() %>% unlist %>% 
  mean(na.rm = T)

# Page 4

# % of all species with at least one first encounter in RCP 2.6 vs. 8.5, 
# climate + land, gcm mean +/- gcm s.d.

# 2.6 

NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures1")] %>% 
  map(~.x %>% dplyr::select(Sp, Sp2) %>% unlist) ->
  
  NESpeciesList1

NESpeciesList1 %>% 
  unlist %>% table() -> 
  
  NESpeciesCounts1

# Across all scenarios

table(NESpeciesCounts1>0) # qplot(NESpeciesCounts)

# Mean + SD

NESpeciesList1 %>% map_dbl(nunique) %>% mean

NESpeciesList1 %>% map_dbl(nunique) %>% sd

# 8.5

NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures4")] %>% 
  map(~.x %>% dplyr::select(Sp, Sp2) %>% unlist) ->
  
  NESpeciesList

NESpeciesList %>% 
  unlist %>% table() -> 
  
  NESpeciesCounts

# Across all scenarios

table(NESpeciesCounts>0) # qplot(NESpeciesCounts)

# Mean + SD

NESpeciesList %>% map_dbl(nunique) %>% mean

NESpeciesList %>% map_dbl(nunique) %>% sd

# of first encounters total, climate+land , rcp AND gcm mean, +/- gcm s.d. 
# With bats ####

NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures")] %>% 
  map_dbl(nrow) %>% 
  MeanSD

NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures1")] %>% 
  map_dbl(nrow) %>% 
  MeanSD

NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures4")] %>% 
  map_dbl(nrow) %>% 
  MeanSD


# Without bats ####

NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures")] %>% 
  map_dbl(~.x %>% filter(!hOrder.x == "Chiroptera", 
                         !hOrder.y == "Chiroptera") %>% nrow) %>% 
  MeanSD

# of first encounters RCP 2.6 and 8.5, climate+land , gcm mean +/- gcm .s.d. (I haven't decided which to use so let's get both)

NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures1")] %>% 
  map_dbl(~.x %>% filter(!hOrder.x == "Chiroptera", 
                         !hOrder.y == "Chiroptera") %>% nrow) %>% 
  MeanSD

NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures4")] %>% 
  map_dbl(~.x %>% filter(!hOrder.x == "Chiroptera", 
                         !hOrder.y == "Chiroptera") %>% nrow) %>% 
  MeanSD

# of new sharing events (sum delta sharing) RCP 2.6 and 8.5, climate+land, gcm mean +/- gcm .s.d. (I haven't decided which to use so let's get both)

NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures1")] %>% 
  map(~.x %>% summarise_at(vars(starts_with("DeltaSharing.Futures1C")), 
                           function(a) sum(a, na.rm = T))) %>% 
  bind_rows() %>% unlist %>% 
  MeanSD

NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures4")] %>% 
  map(~.x %>% summarise_at(vars(starts_with("DeltaSharing.Futures4C")), 
                           function(a) sum(a, na.rm = T))) %>% 
  bind_rows() %>% unlist %>% 
  MeanSD

# climate + land + dispersal RCP 2.6 for page 5:

NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures1")] %>% 
  map(~.x %>% summarise_at(vars(starts_with("DeltaSharing.Futures1A")), 
                           function(a) sum(a, na.rm = T))) %>% 
  bind_rows() %>% unlist


# ~~~~ Page 5 ####

# % reduction of first encounters AND viral sharing in (climate+land+dispersal) compared to (climate+land), RCP 2.6, gcm mean +/- gcm .s.d.

(1 - NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures1")] %>% 
   map(~.x %>% nrow) %>% unlist/(NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures1")] %>% 
                                   map(~.x %>% nrow) %>% unlist)) %>% 
  MeanSD

(1 - NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures1")] %>% 
    map(~.x %>% summarise_at(vars(starts_with("DeltaSharing.Futures1A")), 
                             function(a) sum(a, na.rm = T))) %>% 
    bind_rows() %>% unlist/(NewEncountersListList[str_detect(names(NewEncountersListList), "C.Futures1")] %>% 
                              map(~.x %>% summarise_at(vars(starts_with("DeltaSharing.Futures1C")), 
                                                       function(a) sum(a, na.rm = T))) %>% 
                              bind_rows() %>% unlist)) %>% 
  MeanSD


# # of first encounters, climate + land + dispersal, RCP 2.6, gcm mean +/- gcm .s.d.

NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures1")] %>% 
  map(~.x %>% nrow) %>% unlist %>% 
  MeanSD

# % of first encounters with a bat in them (sum of bat-bat and bat-nonbat), RCP 2.6 and 8.5, gcm mean +/- gcm .s.d.

(NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures")] %>% 
    map(~.x %>% filter(hOrder.x == "Chiroptera"|hOrder.y == "Chiroptera") %>% nrow) %>% 
    unlist/(
      NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures")] %>% 
        map(~.x %>% nrow) %>% unlist 
      
    )) %>% MeanSD

(NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures1")] %>% 
    map(~.x %>% filter(hOrder.x == "Chiroptera"|hOrder.y == "Chiroptera") %>% nrow) %>% 
    unlist/(
      NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures1")] %>% 
        map(~.x %>% nrow) %>% unlist 
      
    )) %>% MeanSD

(NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures4")] %>% 
    map(~.x %>% filter(hOrder.x == "Chiroptera"|hOrder.y == "Chiroptera") %>% nrow) %>% 
    unlist/(
      NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures4")] %>% 
        map(~.x %>% nrow) %>% unlist 
      
    )) %>% MeanSD

# ~~~~ Page 6 #####

#Ebola-nonEbola first encounters in 2.6 and 8.5, CL and CLD, gcm mean +/- gcm .s.d.

EbolaEncountersListList[str_detect(names(EbolaEncountersListList), "Futures1A")] %>% 
  map(~.x %>% nrow) %>% unlist %>% 
  MeanSD

EbolaEncountersListList[str_detect(names(EbolaEncountersListList), "Futures4A")] %>% 
  map(~.x %>% nrow) %>% unlist %>% 
  MeanSD

EbolaEncountersListList[str_detect(names(EbolaEncountersListList), "Futures1C")] %>% 
  map(~.x %>% nrow) %>% unlist %>% 
  MeanSD

EbolaEncountersListList[str_detect(names(EbolaEncountersListList), "Futures4C")] %>% 
  map(~.x %>% nrow) %>% unlist %>% 
  MeanSD

#Ebola-nonEbola viral sharing events (sum delta) in 2.6 and 8.5, CL and CLD, gcm mean +/- gcm .s.d.

EbolaEncountersListList[str_detect(names(EbolaEncountersListList), "Futures1A")] %>% 
  map(~.x %>% summarise_at("DeltaSharing.Futures1A", function(a) sum(a, na.rm = T))) %>% unlist %>% 
  MeanSD

EbolaEncountersListList[str_detect(names(EbolaEncountersListList), "Futures4A")] %>% 
  map(~.x %>% summarise_at("DeltaSharing.Futures4A", function(a) sum(a, na.rm = T))) %>% unlist %>% 
  MeanSD

EbolaEncountersListList[str_detect(names(EbolaEncountersListList), "Futures1C")] %>% 
  map(~.x %>% summarise_at("DeltaSharing.Futures1C", function(a) sum(a, na.rm = T))) %>% unlist %>% 
  MeanSD

EbolaEncountersListList[str_detect(names(EbolaEncountersListList), "Futures4C")] %>% 
  map(~.x %>% summarise_at("DeltaSharing.Futures4C", function(a) sum(a, na.rm = T))) %>% unlist %>% 
  MeanSD

#Bat-primate first encounters in 2.6 and 8.5, CL and CLD, gcm mean +/- gcm .s.d.

BPEncountersListList[str_detect(names(BPEncountersListList), "Futures1A")] %>% 
  map(~.x %>% nrow) %>% unlist %>% 
  MeanSD

BPEncountersListList[str_detect(names(BPEncountersListList), "Futures4A")] %>% 
  map(~.x %>% nrow) %>% unlist %>% 
  MeanSD

BPEncountersListList[str_detect(names(BPEncountersListList), "Futures1C")] %>% 
  map(~.x %>% nrow) %>% unlist %>% 
  MeanSD

BPEncountersListList[str_detect(names(BPEncountersListList), "Futures4C")] %>% 
  map(~.x %>% nrow) %>% unlist %>% 
  MeanSD

#Bat-primate viral sharing events (sum delta) in 2.6 and 8.5, CL and CLD, gcm mean +/- gcm .s.d.

BPEncountersListList[str_detect(names(BPEncountersListList), "Futures1A")] %>% 
  map(~.x %>% summarise_at("DeltaSharing.Futures1A", function(a) sum(a, na.rm = T))) %>% unlist %>% 
  MeanSD

BPEncountersListList[str_detect(names(BPEncountersListList), "Futures4A")] %>% 
  map(~.x %>% summarise_at("DeltaSharing.Futures4A", function(a) sum(a, na.rm = T))) %>% unlist %>% 
  MeanSD

BPEncountersListList[str_detect(names(BPEncountersListList), "Futures1C")] %>% 
  map(~.x %>% summarise_at("DeltaSharing.Futures1C", function(a) sum(a, na.rm = T))) %>% unlist %>% 
  MeanSD

BPEncountersListList[str_detect(names(BPEncountersListList), "Futures4C")] %>% 
  map(~.x %>% summarise_at("DeltaSharing.Futures4C", function(a) sum(a, na.rm = T))) %>% unlist %>% 
  MeanSD

# ~~~~ Page 7  ####

# % of species that experienced a net range gain, RCP 2.6 and 8.5, CLD, gcm mean +/- gcm .s.d.

FullRangeChangeDF %>% filter(str_detect(Rep, "BufferClimateLandUse.Futures1$")) %>% 
  mutate_at("OverallChange", ~ifelse(.x == Inf, NA, .x)) %>% 
  group_by(ClimateRep) %>% 
  summarise_at("OverallChange", ~.x %>% subtract(1) %>% Prev) %>% 
  pull(OverallChange) %>% MeanSD

FullRangeChangeDF %>% filter(str_detect(Rep, "BufferClimateLandUse.Futures4$")) %>% 
  mutate_at("OverallChange", ~ifelse(.x == Inf, NA, .x)) %>% 
  group_by(ClimateRep) %>% 
  summarise_at("OverallChange", ~.x %>% subtract(1) %>% Prev) %>% 
  pull(OverallChange) %>% MeanSD

# Average % range gain or change, RCP 2.6 and 8.5, CLD,  gcm mean +/- gcm .s.d.

# AllMammalDFList %>% map(~.x %>% dplyr::select(Sp, Sp2) %>% unlist %>% unique) -> AllMammalsList

# AllMammalsList %>% reduce(union) -> AllMammals

# saveRDS(AllMammals, file = "~/Albersnet/Iceberg Files/Summary/AllMammals.rds)

FullRangeChangeDF %>% 
  filter(str_detect(Rep, "BufferClimateLandUse.Futures1$")) %>% 
  filter(Sp %in% AllMammals) %>% 
  group_by(ClimateRep) %>% 
  summarise_at("OverallChange", ~median(.x, na.rm = T)) %>% 
  pull(OverallChange) %>% subtract(1) %>% multiply_by(100) %>% 
  MeanSD

FullRangeChangeDF %>% 
  filter(str_detect(Rep, "BufferClimateLandUse.Futures1$")) %>% 
  filter(Sp %in% AllMammals) %>% 
  group_by(ClimateRep) %>% 
  summarise_at("OverallChange", ~Prev(.x>1)) %>% 
  pull(OverallChange) %>% multiply_by(100) %>% 
  MeanSD

FullRangeChangeDF %>% filter(str_detect(Rep, "BufferClimateLandUse.Futures4$")) %>% 
  filter(Sp %in% AllMammals) %>% 
  #filter(OverallChange<100) %>% 
  #mutate_at("OverallChange", ~ifelse(.x == Inf, NA, .x)) %>% 
  group_by(ClimateRep) %>% 
  summarise_at("OverallChange", ~median(.x, na.rm = T)) %>% 
  pull(OverallChange) %>% subtract(1) %>% multiply_by(100) %>% 
  MeanSD

FullRangeChangeDF %>% filter(str_detect(Rep, "BufferClimateLandUse.Futures4$")) %>% 
  filter(Sp %in% AllMammals) %>% 
  group_by(ClimateRep) %>% 
  summarise_at("OverallChange", ~Prev(.x>1)) %>% 
  pull(OverallChange) %>% multiply_by(100) %>% 
  MeanSD

## of species that lost their entire range in CL and CLD, in 8.5, gcm mean +/- gcm .s.d.

FullRangeChangeDF %>% 
  filter(Sp %in% AllMammals) %>% 
  filter(str_detect(Rep, "^ClimateLandUse.Futures4$")) %>% 
  mutate_at("OverallChange", ~as.numeric(.x==0)) %>% 
  group_by(ClimateRep) %>% 
  summarise_at("OverallChange", ~sum(.x, na.rm = T)) %>% 
  pull(OverallChange) %>% 
  MeanSD

FullRangeChangeDF %>% 
  filter(Sp %in% AllMammals) %>% 
  filter(str_detect(Rep, "BufferClimateLandUse.Futures4$")) %>% 
  mutate_at("OverallChange", ~as.numeric(.x==0)) %>% 
  group_by(ClimateRep) %>% 
  summarise_at("OverallChange", ~sum(.x, na.rm = T)) %>% 
  pull(OverallChange) %>% 
  MeanSD

FullRangeChangeDF %>% 
  filter(Sp %in% AllMammals) %>% 
  filter(str_detect(Rep, "^ClimateLandUse.Futures4$")) %>% 
  mutate_at("OverallChange", ~as.numeric(.x==0)) %>% 
  group_by(ClimateRep) %>% 
  summarise_at("OverallChange", ~sum(.x, na.rm = T)) %>% 
  pull(OverallChange) %>% subtract(

FullRangeChangeDF %>% 
  filter(Sp %in% AllMammals) %>% 
  filter(str_detect(Rep, "BufferClimateLandUse.Futures4$")) %>% 
  mutate_at("OverallChange", ~as.numeric(.x==0)) %>% 
  group_by(ClimateRep) %>% 
  summarise_at("OverallChange", ~sum(.x, na.rm = T)) %>% 
  pull(OverallChange), .) %>% MeanSD

# % fewer first encounters in 8.5 compared to 2.6 (CLD),  gcm mean +/- gcm .s.d.
# NB EXCLUDES Climate 5 because no 8.5

(1 - NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures1")][-5] %>% 
    map(~.x %>% nrow) %>% unlist/(NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures4")] %>% 
                                    map(~.x %>% nrow) %>% unlist)) %>% multiply_by(100) %>% 
  MeanSD

# % lower connectivity in the network in 8.5 compared to 2.6 (CLD),  gcm mean +/- gcm .s.d.  

SummaryNetworkChanges %>%
  tidyr::gather("Key", "Value", -ClimateRep) %>% 
  filter(str_detect(Key, "1A|4A")) %>% na.omit %>% 
  group_by(Key) %>% 
  summarise(Mean = mean(Value), 
            SD = sd(Value))

SummaryNetworkChanges %>%
  tidyr::gather("Key", "Value", -ClimateRep) %>% 
  filter(!ClimateRep == 5) %>% 
  filter(str_detect(Key, "1A|4A")) %>% na.omit %>% 
  data.frame -> LongSummaryChanges

LongSummaryChanges %>% 
  group_by(ClimateRep) %>%
  summarise(Diff = diff(Value)) %>% 
  pull(Diff) %>% 
  divide_by(LongSummaryChanges[1:8,"Value"]) %>% MeanSD

sink()
