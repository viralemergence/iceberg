
# Betacov ####

library(tidyverse); library(fs); library(magrittr); library(glue)
library(ggregplot)

setwd("~/Albersnet/")

# Setting up ####

Unlist1 <- function(a) unlist(a, recursive = F)

glue::glue("Climate{1:9}") ->
  
  ClimateReps

Betacov <- read_csv("https://raw.githubusercontent.com/viralemergence/virionette/master/03_interaction_data/virionette.csv")

Betacov %<>% filter(virus_genus == "Betacoronavirus") %>% 
  mutate_at("host_species", ~.x %>% str_replace_all(" ", "_"))

BetacovHosts <- Betacov$host_species

# Importing New Encounters ####

print("Getting new encounters!")

glue("~/Albersnet/Iceberg Files/{ClimateReps}/Iceberg Output Files/NewEncounters.rds") %>% 
  map(~readRDS(.x) %>% Unlist1 %>% 
        map(function(a) a %>% filter(Sp%in%BetacovHosts|Sp2%in%BetacovHosts))) -> 
  NewEncountersListList

names(NewEncountersListList) <- ClimateReps

NewEncountersListList %<>% Unlist1

print("Saving!")

saveRDS(NewEncountersListList, file = "BetacovEncounters.rds")

