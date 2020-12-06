![banner](https://github.com/cjcarlson/iceberg/blob/master/Banner.png)

# Iceberg

## Carlson, Albery, et al.: Climate change will drive novel cross-species viral transmission

This repo relates to this preprint: https://www.biorxiv.org/content/10.1101/2020.01.24.918755v1

We ask the question: how might climate change drive mammal species to shift their ranges, giving viruses opportunities to jump into new species?

The study can be separated into two conceptual and practical stages: generating and processing species distribution models (SDMs) for the present and the future, and then using them to predict how they will alter viral ecology as a result of climate change.

1. The study uses species distribution models (SDMs) to predict the movement of mammal species' ranges between 2020 and 2070, predicts how these movements are likely to result in novel encounters between mammal species that currently have no overlapping ranges.

2. After projecting species' movements, we used a published model of interspecific viral sharing (Albery et al., 2020: https://www.nature.com/articles/s41467-020-16153-4) to estimate how these movements are likely to alter viral sharing patterns in the future.

This repo uses the following datasets:

- The IUCN species ranges (https://www.iucnredlist.org/resources/spatial-data-download).  
- The PanTheria mammalian phenotypic dataset (http://esapubs.org/archive/ecol/E090/184/).  
- The Olival et al. (2017) host-virus association dataset (https://www.nature.com/articles/nature22975).  
- The Fritz et al (2009) mammalian phylogenetic supertree (https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1461-0248.2009.01307.x).
- It takes in rasters generated by an SDM protocol run by [Cory Merow](https://github.com/cmerow), across >3000 species, first for the present day and then for four future representative concentration pathways (RCPs). This was repeated using nine different climate models, totalling 1 current prediction and 4x9 = 36 future predictions per species. 37x3000 = 111,000 rasters before our processing begins.

### NB this code was run on a 72-core supercluster. It involves extremely large lists of files and memory-hungry processes. Running this on a desktop or laptop may not be possible.

This repo is organised as follows:

## For each climate model, the code in the folder `Iceberg Code/1_ENM Code` will:

### For the currents:
-- Load the rasters for each species and resample them to a global Mercator projection.  
-- Clip them to within 1000km of their IUCN ranges and ensuring that they remain within the same continent.  
-- Clip each species by their known habitat use.  
-- Buffer each species' current range by its known dispersal ability, to allow us to control for their ability to move towards their future ranges.  
-- write out data frames (`GretCDFs`) of X and Y coordinates alongside presence/absence for each of the stages of the processing pipeline.  

### For the futures:
-- Load and resample the rasters, as with the currents.  
-- Clip the futures by their habitat use, as with the currents.  
-- Use the dispersal buffer generated by the currents to prevent species from moving further than is biologically plausible.  
-- Write out data frames (`GretCDFs`), as with the currents.  

## For each climate model, the code in the folder `Iceberg Code/2_GAMM Code` will:

-- `0_Master Iceberg.R`: Loop through the rest of the code.  
-- `1_Iceberg Spatial.R`: Import the `GretCDFs` and use them to calculate current and future overlap between species' ranges. This makes use of the `PairsWisely` function from the `SpRanger` package, per Albery et al. (2020).  
-- `2_Iceberg Data Import.R`: Import the accessory phenotypic data used for the GAMMs, including PanTheria and the host-virus association dataset. Will make a pairwise viral sharing dataset, where every row is a species pair, and each column is a different pairwise trait (e.g. phylogenetic similarity and geographic overlap).  
-- `3_Iceberg GAMs.R`: Run Generalised Additive Mixed Models (GAMMs) exaining how these pairwise traits determine the probability of sharing at least one virus in the Olival et al. (2017) dataset.  
-- `4_Iceberg Prediction.R`: Use the GAMM to predict viral sharing probabilities for species in the present day and in all of the future scenarios. The difference between each future and current scenario comprises the overall change in sharing, and therefore the number of new sharing events.  
-- `5_Iceberg Mapping.R`: Map the novel encounters and sharing events, showing how viruses will be shared newly in 50 years' time.  
-- `6_Iceberg Submaps.R`: Repeat the mapping process, but with subsets of data: bat-primate encounters and known Ebola hosts.  
-- `6b_Order New Encounters.R`  
-- `6c_Range Changes.R`: Import the GretCDFs and examine how different climate scenarios alter change in species' range sizes.  
-- `7_New Encounter GAMs.R`: Run a model designed to test the determinants of spatial patterns of novel encounters, at the raster pixel level.  

### The remaining code will generate display outputs and results:

-- `8_Individual Reps Iceberg Display Outputs.R`: Generate main figures, for each of the 9 climate models.  
-- `8b_Individual Reps Extended Data Figures.R`: Generate extended data figures, for each of the 9 climate models.  
-- `8c_Summary Iceberg Display Outputs.R`: Generate main figures, for mean +/- standard deviation of the results across all 9 models.   
-- `9_Text Numbers.R`: Generate results to put in the main text.  

### The file structure is as follows:

- Iceberg Files
-- Climate1  
-- Climate2  
-- Climate3  
-- Climate4  
-- Climate5  
-- Climate6  
-- Climate7  
-- Climate8  
-- Climate9  
-- Summary  

### Each of the Climate Rep folders includes:
-- Data: Phenotypic data etc.; does not vary across Climate reps.  
-- Input files: The manipulated and resampled rasters and GretCDFs.  
-- Output Files: The GAMs, summarised data frames, etc.  
-- Figures: The maps, figures etc, for the focal climate rep.  

### The summary folder includes summarised versions of all of these files.


For questions about the pipeline, email Greg at gfalbery@gmail.com and I'll get back to you!


