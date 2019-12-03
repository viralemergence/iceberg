
# Rscript "Comparing Different Projection Outputs.R"

# Running different iceberg adj projections ####

# Mercator area corrected ####

rm(list = ls())

Mercator = readRDS("Iceberg Output Files/IcebergAdjList_MercatorProj.rds")
IcebergAdjList = Mercator

CodeRoot <- "R Code/0_Data Import"

StartTime <- Sys.time()
print("Import")
source(paste0(CodeRoot,"/","0a_EHA Data Import.R"))
print("Phylo")
source(paste0(CodeRoot,"/","0b_Phylogenetic Data Import.R" ))
print("Space")
source(paste0(CodeRoot,"/","0c2_Exhaustive Spatial Data Import.R"))
source(paste0(CodeRoot,"/","0c_Kludging Spatial Data Import.R"))
#print("Ecology")
#source(paste0(CodeRoot,"/","0d_Ecological Data Import.R"))
print("Final Dataset")
source(paste0(CodeRoot,"/","0e_Creating Final Host Dataset.R"))
#print("Subsets")

CodeRoot <- "Iceberg Code"

print("GAMs")
source(paste0(CodeRoot,"/","Iceberg GAMs.R" ))

print("Predicting")
source(paste0(CodeRoot,"/","Iceberg Prediction.R"))

saveRDS(PredNetworkList, file = "Iceberg Output Files/MercatorPredNetworkList.rds")
saveRDS(AllMammaldf, file = "Iceberg Output Files/MercatorAllMammaldf.rds")
saveRDS(NewEncountersList, file = "Iceberg Output Files/MercatorNewEncounters.rds")
save(BAMList, file = "Iceberg Output Files/MercatorBAMList.Rdata")

# Mollweide Projection ####

rm(list = ls())

Mollweide_NoArea = readRDS("Iceberg Output Files/IcebergAdjList_MollweideProj_NoArea.rds")

IcebergAdjList = Mollweide_NoArea

CodeRoot <- "R Code/0_Data Import"

StartTime <- Sys.time()
print("Import")
source(paste0(CodeRoot,"/","0a_EHA Data Import.R"))
print("Phylo")
source(paste0(CodeRoot,"/","0b_Phylogenetic Data Import.R" ))
print("Space")
source(paste0(CodeRoot,"/","0c2_Exhaustive Spatial Data Import.R"))
source(paste0(CodeRoot,"/","0c_Kludging Spatial Data Import.R"))
#print("Ecology")
#source(paste0(CodeRoot,"/","0d_Ecological Data Import.R"))
print("Final Dataset")
source(paste0(CodeRoot,"/","0e_Creating Final Host Dataset.R"))
#print("Subsets")

CodeRoot <- "Iceberg Code"

print("GAMs")
source(paste0(CodeRoot,"/","Iceberg GAMs.R" ))

print("Predicting")
source(paste0(CodeRoot,"/","Iceberg Prediction.R"))

saveRDS(PredNetworkList, file = "Iceberg Output Files/MollweidePredNetworkList.rds")
saveRDS(AllMammaldf, file = "Iceberg Output Files/MollweideAllMammaldf.rds")
saveRDS(NewEncountersList, file = "Iceberg Output Files/MollweideNewEncounters.rds")
save(BAMList, file = "Iceberg Output Files/MollweideBAMList.Rdata")

# Uncorrected Mercator ####

rm(list = ls())

MercatorUnc = readRDS("Iceberg Output Files/IcebergAdjList_MercatorProjUnc.rds")

IcebergAdjList = MercatorUnc

CodeRoot <- "R Code/0_Data Import"

StartTime <- Sys.time()
print("Import")
source(paste0(CodeRoot,"/","0a_EHA Data Import.R"))
print("Phylo")
source(paste0(CodeRoot,"/","0b_Phylogenetic Data Import.R" ))
print("Space")
source(paste0(CodeRoot,"/","0c2_Exhaustive Spatial Data Import.R"))
source(paste0(CodeRoot,"/","0c_Kludging Spatial Data Import.R"))
#print("Ecology")
#source(paste0(CodeRoot,"/","0d_Ecological Data Import.R"))
print("Final Dataset")
source(paste0(CodeRoot,"/","0e_Creating Final Host Dataset.R"))
#print("Subsets")

CodeRoot <- "Iceberg Code"

print("GAMs")
source(paste0(CodeRoot,"/","Iceberg GAMs.R" ))

print("Predicting")
source(paste0(CodeRoot,"/","Iceberg Prediction.R"))

saveRDS(PredNetworkList, file = "Iceberg Output Files/MercatorUncPredNetworkList.rds")
saveRDS(AllMammaldf, file = "Iceberg Output Files/MercatorUncAllMammaldf.rds")
saveRDS(NewEncountersList, file = "Iceberg Output Files/MercatorUncNewEncounters.rds")
save(BAMList, file = "Iceberg Output Files/MercatorUncBAMList.Rdata")
