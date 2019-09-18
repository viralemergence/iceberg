
# Rscript "Final Iceberg Code/0_Master Iceberg.R"

# rm(list = ls())

# Iceberg Spatial.R
# Iceberg Data Import.R
# Iceberg Gams.R
# 4_Iceberg Sharing Prediction ####

#print("Currents!")
#source("01_Iceberg ENM Currents.R")

#print("Futures!")
#source("02_Iceberg ENM Futures.R")

print("Spatial!")
source("Final Iceberg Code/1_Iceberg Spatial.R")

print("Data Import!")
source("Final Iceberg Code/2_Iceberg Data Import.R")

print("GAMs!")
source("Final Iceberg Code/3_Iceberg GAMs.R")

print("Prediction!")
source("Final Iceberg Code/4_Iceberg Prediction.R")

print("Mapping!")
source("Final Iceberg Code/5_Iceberg Mapping.R")
