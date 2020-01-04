
head(ChangeDF)

ChangeDF %>% filter(Rep == 'Climate.Futures1') %>% 
  pull(OverallChange) %>% na.omit() %>% mean()


ChangeDF %>% filter(Rep == 'BufferClimateLandUse.Futures1') %>% 
  pull(OverallChange) %>% na.omit() %>% mean()

# Species with new encounters ####
NewEncountersList %>% unlist(recursive = F) %>% 
  map(`[`, c("Sp", "Sp2")) %>%
  map_dbl(~.x %>% unlist %>% nunique) %>% 
  range %>% divide_by(length(AllMammals)) %>%
  round(2)

# Deviance contributions in second GAMs ####

# Modularity #### 

# Bat new encounters percentage ####

NewEncountersList[["A"]][["Futures1"]] %>%
  filter(hOrder.x == "Chiroptera"|hOrder.y == "Chiroptera") %>%
  nrow %>% 
  divide_by(nrow(NewEncountersList[["A"]][["Futures1"]]))


# Ebola new encounters ####

EbolaEncounters$Futures1C %>% nrow

EbolaEncounters$Futures1A %>% nrow

EbolaEncounters$Futures1A$DeltaSharing.Futures1A %>% 
  sum %>% round

# Bat primate encounters ####

BPEncounters %>% map_dbl(nrow) %>% range

PredReps[2:5] %>% rep(4) %>% 
  paste0(rep(LETTERS[1:4], each = 4)) %>%
  map_dbl(~BPEncounters[[.x]][,paste0("DeltaSharing.", .x)] %>% sum) %>%
  range %>% round

# Dispersal ####

# Range loss ####

# Import gretcdfs and compare how much range is lost in certain scenarios



# Range loss to zero ####

IcebergAdjList %>% 
  unlist(recursive = F) ->
  RangeSubList

RangeSubList[-c(1,6,11,16)] %>% 
  map_dbl(~sum(rowSums(.x)==0)) ->
  RangeLost

RangeLost[["A.Futures4"]] - 
  RangeLost[["C.Futures4"]] 

NewEncountersList$A$Futures4 %>% nrow %>%
  divide_by(nrow(NewEncountersList$A$Futures1))

sum(NewEncountersList$A$Futures4$DeltaSharing.Futures4A) -
  sum(NewEncountersList$A$Futures1$DeltaSharing.Futures1A)

mean(AllMammaldf$Sharing.Futures4A)/
  mean(AllMammaldf$Sharing.Futures1A)

# Methods ####

paste0("Iceberg Input Files/","MaxEnt","/01_Raw/Currents") %>% 
  list.files(full.names = T) %>% 
  append(paste0("Iceberg Input Files/","RangeBags","/01_Raw/Currents") %>% list.files(full.names = T)) ->
  FullFiles

paste0("Iceberg Input Files/","MaxEnt","/01_Raw/Currents") %>% 
  list.files() %>% str_remove(".rds$") %>%
  append(paste0("Iceberg Input Files/","RangeBags","/01_Raw/Currents") %>% list.files() %>% str_remove(".rds$")) ->
  names(FullFiles)

Species <- SpeciesList %>% unlist %>% sort()

Files <- FullFiles[Species]

names(FullFiles) %>% intersect(MarineSp)

names(FullFiles) %>% setdiff(MarineSp) -> NonMarine

NonMarine %>% setdiff(NonEutherianSp) -> Eutherian

Eutherian %>% setdiff(colnames(IcebergAdjList[["D"]][[1]]))

length(NullRangeBagRares) + dim()

HabitatList %>% map_dbl(~all(is.na(.x))) %>% sum

sum(is.na(Dispersals$disp50))

colnames(IcebergAdjList[[1]][[1]]) -> WorkingSp

names(FullFiles) %>% setdiff(NonEutherianSp) %>% length

length(FHN)

# GAMM validation ####


# Mapping methods ####

length(AllMammals)



