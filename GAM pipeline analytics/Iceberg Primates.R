# Iceberg just Primate encounters ####

PrimateEncounters <- lapply(NewEncountersList, function(a){ a %>%
    mutate(Primate = as.numeric(hOrder.x%in%c("Primates")|hOrder.y%in%c("Primates"))) %>% 
    filter(Primate == 1)
  
})

lapply(PrimateEncounters, nrow)

