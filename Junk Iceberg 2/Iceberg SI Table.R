
# Iceberg SI Table for Scenarios

# New Encounters

SITable <- data.frame(
  
  RCP = c("2.6", "3.5", "4.6", "8.5"),
  
  # New Encounter Numbers
  NE_Number = sapply(NewEncountersList, nrow),
  
  # New Encounter Sharing Probabilies
  NE_MeanSharing = PredReps[2:5] %>% sapply(function(a){
    
    NewEncountersList[[a]][,SharingVars[a]] %>% mean
    
  }),
  
  # New Encounter Sharing Probabilies
  NE_SumSharing = PredReps[2:5] %>% sapply(function(a){
    
    NewEncountersList[[a]][,SharingVars[a]] %>% sum
    
  }),
  
  # New Encounter Sharing Probabilies
  NE_DeltaSharing = PredReps[2:5] %>% sapply(function(a){
    
    NewEncountersList[[a]][,paste0("Delta",SharingVars[a])] %>% sum
    
  }),
  
  # Mean pairwise viral sharing probability
  
  MeanSharing = sapply(PredReps[2:5], function(a){
    
    AllMammaldf[,SharingVars[a]] %>% mean
    
  }),
  
  # Delta pairwise viral sharing probability (vs 2019)
  
  DeltaSharing = PredReps[2:5] %>% sapply(function(a){
    
    AllMammaldf[,paste0("Delta",SharingVars[a])] %>% mean
    
  })
  
)

SISubTable <- data.frame(
  
  RCP = c("2.6", "3.5", "4.6", "8.5"),
  
  # Mean new viral sharing events of first bat-primate contacts
  
  BPNE_Number = sapply(PredReps[2:5], function(a){
    
    BPEncounters[[a]] %>% nrow
    
  }),
  
  # Mean new viral sharing events of first bat-primate contacts
  
  BPNE_SumSharing = sapply(PredReps[2:5], function(a){
    
    BPEncounters[[a]][,paste0("Delta",SharingVars[a])] %>% sum
    
  })

)
