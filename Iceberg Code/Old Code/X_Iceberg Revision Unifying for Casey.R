
# X_Summarising files for Casey ####

library(tidyverse)

AssocsBase <- read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/associations.csv") %>% data.frame()
HP3EbolaHosts <- AssocsBase %>% filter(vVirusNameCorrected == "Zaire_ebolavirus")
LauraEbolaHosts <- read.csv("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/ZEBOV hosts.csv")

EbolaHosts <- union(HP3EbolaHosts$hHostNameFinal, LauraEbolaHosts$bat_species) %>% 
  intersect(list.files("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents") %>% 
              str_remove(".rds$"))

#EbolaHosts %>% setdiff(c("Miniopterus_schreibersii","Pipistrellus_pipistrellus", "Hipposideros_pomona",
#                         "Cynopterus_sphinx","Acerodon_jubatus", "Rousettus_leschenaultii")) ->
#  
 # EbolaSpecies

EbolaSpecies <- EbolaHosts

AllAllMammalDF <- readRDS("~/Albersnet/Iceberg Files/Summary/Iceberg Output Files/AllAllMammalDF.rds")

AllAllMammalDF %>% 
  filter(Space.CurrentsA == 0, Space.Futures1A>0) %>% 
  bind_rows(
    
    AllAllMammalDF %>% 
      filter(Space.CurrentsA == 0, Space.Futures1A>0) %>% 
      rename(hOrder.x = hOrder.y, hOrder.y = hOrder.x)
    
  ) %>% 
  
  filter(Sp %in% EbolaSpecies | Sp2 %in% EbolaSpecies, 
         !(Sp %in% EbolaSpecies & Sp2 %in% EbolaSpecies)) %>% 
  
  group_by(ClimateRep, hOrder.x, hOrder.y) %>% 
  
  summarise(DeltaSharing.Futures1A.Positive = mean(DeltaSharing.Futures1A[DeltaSharing.Futures1A>0]), 
            DeltaSharing.Futures1A.Full = mean(DeltaSharing.Futures1A), 
            Sharing.CurrentsA = mean(Sharing.CurrentsA),
            N = n()) -> SummaryLinks

SummaryLinks %>% 
  ungroup %>% 
  select(-ClimateRep) %>% 
  group_by(hOrder.x, hOrder.y) %>% 
  summarise_all(~mean(.x, na.rm = T)) %>% 
  bind_cols(SummaryLinks %>% 
              ungroup %>% 
              select(-ClimateRep) %>% 
              group_by(hOrder.x, hOrder.y) %>% 
              summarise_all(~sd(.x, na.rm = T)) %>% 
              ungroup %>% 
              select(-c(hOrder.x, hOrder.y)) %>% 
              rename_all(~paste0(.x, ".sd"))) ->
  
  MeanLinks

MeanLinks

# Overall (not just new encounters) ####

AllAllMammalDF %>% 
  bind_rows(
    
    AllAllMammalDF %>% 
      rename(hOrder.x = hOrder.y, hOrder.y = hOrder.x)
    
  ) %>% 
  
  filter(Sp %in% EbolaSpecies | Sp2 %in% EbolaSpecies, 
         !(Sp %in% EbolaSpecies & Sp2 %in% EbolaSpecies)) %>% 
  
  group_by(ClimateRep, hOrder.x, hOrder.y) %>% 
  
  summarise(DeltaSharing.Futures1A.Positive = mean(DeltaSharing.Futures1A[DeltaSharing.Futures1A>0]), 
            DeltaSharing.Futures1A.Full = mean(DeltaSharing.Futures1A), 
            Sharing.CurrentsA = mean(Sharing.CurrentsA),
            N = n()) -> 
  SummaryLinks2

# Just making a new allmammaldf for Casey ####

AllAllMammalDF %>% 
  group_by(Sp, Sp2, hOrder.x, hOrder.y, hFamily.x, hFamily.y) %>% 
  summarise_if(is.numeric, mean) -> 
  AllMammalNumericMeans

AllMammalNumericMeans %>% dim

AllMammalNumericMeans %>% 
  bind_rows(AllMammalNumericMeans %>% 
              rename(Sp = Sp2, Sp2 = Sp, 
                     hOrder.x = hOrder.y, hOrder.y = hOrder.x,
                     hFamily.x = hFamily.y, hFamily.y = hFamily.x)) -> 
  
  AllMammalMatrixMeans

AllMammalNumericMeans %>% write.csv("OrderEdges.csv")
AllMammalMatrixMeans %>% write.csv("OrderEdges2.csv")

AllMammalMatrixMeans %>% saveRDS("OrderEdges2.rds")

EbolaSpecies
