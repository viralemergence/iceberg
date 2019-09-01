
# Iceberg SI Table for Scenarios

# New Encounters




# New Encounter Numbers
sapply(NewEncountersList, nrow)

# New Encounter Sharing Probabilies
PredReps[2:5] %>% sapply(function(a){
  
  NewEncountersList[[a]][,SharingVars[a]] %>% mean
  
})

# Expected New Viral Contacts Emerging
PredReps[2:5] %>% sapply(function(a){
  
  NewEncountersList[[a]] %>% nrow
  
})


# Mean pairwise viral sharing probability

sapply(PredReps[2:5], function(a){
  
  AllMammaldf[,SharingVars[a]] %>% mean
  
})

# Delta pairwise viral sharing probability (vs 2019)

PredReps[2:5] %>% sapply(function(a){
  
  AllMammaldf[,paste0("Delta",SharingVars[a])] %>% mean
  
})

# Number of first contacts

sapply(NewEncountersList, nrow)

# Mean increase in sharing probability of first contacts

sapply(PredReps[2:5], function(a){
  
  NewEncountersList[[a]][,paste0("Delta",SharingVars[a])] %>% mean
  
})

# Expected ``sharing events'' from first contacts

sapply(PredReps[2:5], function(a){
  
  NewEncountersList[[a]][,paste0("Delta",SharingVars[a])] %>% sum
  
})


# Mean raw sharing probability of first contacts

sapply(PredReps[2:5], function(a){
  
  NewEncountersList[[a]][,SharingVars[a]] %>% mean
  
})

# Mean new viral sharing events of first bat-primate contacts

sapply(PredReps[2:5], function(a){
  
  BPEncounters[[a]][,paste0("Delta",SharingVars[a])] %>% sum
  
})
