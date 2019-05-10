# Master Source Code for Iceberg Viral Chatter Analyses ####

library(dplyr)

rm(list = ls())

# Running data setup scripts ####

CodeRoot <- "~/Albersnet"

StartTime <- Sys.time()

print("Import")
source(paste0(CodeRoot,"/","Iceberg Spatial.R"))

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

CodeRoot <- "~/Albersnet"

print("Phylo")
source(paste0(CodeRoot,"/","Iceberg GAMs.R" ))

print("Predicting")
source(paste0(CodeRoot,"/","Iceberg Prediction.R"))

print("New Encounters")
source(paste0(CodeRoot,"/","Iceberg New Encounters.R"))

#print("Ecology")

EndTime <- Sys.time()

EndTime - StartTime
