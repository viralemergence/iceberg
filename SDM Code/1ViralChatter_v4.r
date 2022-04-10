# !note that i never automated the copying of future layers tthat needed to be processed
#====================================================================
#== 1. Install packages
#====================================================================

rm(list=ls())
library(devtools)
if(Sys.info()['user']=='ctg') 	library(colorout)

library(BIENWorkflow)
library(cmsdm)
library(ggplot2)
library(gridExtra)
library(textTinyR)
library(parallel)

mc.cores=1

#====================================================================
#== 2. Create directories
#====================================================================

runName='VC_12_23_21'

#isRares=length(grep('Rare',runName))>0

redoSort=FALSE; whichSpecies='all'

if(Sys.info()['user']=='ctg') {
 #baseDir=paste0('/Users/ctg/Documents/SDMs/ViralChatter/Revision2/', runName,'_inputs')
 baseDir=paste0('/Volumes/cm2/ViralChatter/', runName,'_inputs')
}

# these use the demo data included in the package
mySpeciesDir='/Users/ctg/Dropbox/Projects/Viral_Chatter/Data/Masked_Presence_Data_By_Expert_Maps_R2'
#'/Users/ctg/Dropbox/Projects/Viral_Chatter/Data/PresenceData'
#myEnvDir= '/Volumes/cm2/ViralChatter/Revision1/CurrentEnv'
myEnvDir='/Volumes/cm2/viralChatterEnv_R2/Present'
otherDirs=list(
  rdistDir=paste0(baseDir,'/rDists'),
  domains=paste0(baseDir,'/domains'))

#sortDirNames=c('PPM','Points','RangeBag')#c('Samples1','Samples3','Samples10')
#myAlgorithmNames=c('PPM','Points','RangeBag')
myAlgorithmNames=c('PPM','Points','RangeBag')
sortDirNames=c('Points','RangeBag','SDM','SDM10_20')
offsetDataDirNames=c('Expert_Maps')
myCustomBackgroundDir=NULL
myCustomBackgroundTable=NULL
# these should be stored as a single stack, in the same order as the env
#myOtherEnvDir=  '/Volumes/cm2/ViralChatter/Revision1/OtherEnv'
myOtherEnvDir=  '/Volumes/cm2/ViralChatter/fix_ip2611'
#'/Volumes/cm2/viralChatterEnv_R2/Future'
myCustomDomainRegionRaster='/Users/ctg/Dropbox/Projects/Viral_Chatter/Data/continents-final.tif'
#myCustomDomainRegionRaster='/Users/ctg/Dropbox/Projects/Viral_Chatter/Data/EcoregionRasterDinnerstein2017_QD_global.tif'
mySamplingModelDir=NULL
mySamplingDataDir=NULL
myOffsetDirs=NULL
samplingDataDirNames=NULL
customBackgroundTablePath=NULL
# run once then move other env over and rerun
allBaseDirs=sdmDirectories(baseDir,
                           myAlgorithmNames=myAlgorithmNames,
                           warn=FALSE,
                           mySpeciesDir=mySpeciesDir,
                           myEnvDir=myEnvDir,
                           myOffsetDirs=myOffsetDirs,
                           mySamplingDataDir=mySamplingDataDir,
                           mySamplingModelDir=mySamplingModelDir,
                           myCustomBackgroundDir=myCustomBackgroundDir,
                           myCustomBackgroundTable= myCustomBackgroundTable,	
                           myCustomDomainRegionRaster= myCustomDomainRegionRaster,												 
                           myOtherEnvDir=myOtherEnvDir,
                           offsetDataDirNames=offsetDataDirNames,
                           samplingDataDirNames=samplingDataDirNames,
                           overwrite=FALSE,
                           sortDirNames=sortDirNames,
                           otherDirs=otherDirs)

#====================================================================
#== 3. Prep environmental layers
#====================================================================

if(Sys.info()['user']=='ctg') { # for cory only
 #source('/Users/ctg/Dropbox/Projects/Viral_Chatter/Code/2ViralChatter_v2.r')\
 	myprj = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
	env=stack(paste0(allBaseDirs$envDir,'/AllEnv.tif'))
	layerNames=try(read.csv(list.files(allBaseDirs$envDir,'.csv',full.names=T), stringsAsFactors=F,header=F))
	world.shp=readOGR( '/Users/ctg/Dropbox/Projects/Climate_Horizons/SpeciesHorizons/HorizonsMaster/rangeHorizonsAnalysis/Model_Analysis/TM_WORLD_BORDERS_SIMPL-0.3/TM_WORLD_BORDERS_SIMPL-0.3.shp')
	#world.shp=spTransform(world.shp,projection(env))
	pointsProj= CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
	
	## aggregate env layer for subsampling very large species (rather than running the spThin algorithm)
	print('prepped coarse grid for thinning')
	coarse.grid.file=paste0(allBaseDirs$miscDir, '/coarse_grid_for_thinning.tif')
	if(!file.exists(coarse.grid.file)){
		env.template=env[[1]]
		coarse.grid=aggregate(env.template,fact=4)
		writeRaster(coarse.grid,file=coarse.grid.file,overwrite=T)
	}
} 
save(allBaseDirs,file=paste0(baseDir,'/allBaseDirs.rdata'))

#====================================================================
#== 4. Sort species
#====================================================================

#--------------------------------------------------------------------
### 4.a sort species by algorithm
#--------------------------------------------------------------------
sortDone=checkSortDone(allBaseDirs, deleteSorted=FALSE)
# here i manually copy the files from the last run to save time
sortSpeciesBulk(allBaseDirs,pointsProj,myprj=myprj,doThin=FALSE, thinCutoff=NULL,verbose=TRUE,overwrite=FALSE,doClean=TRUE,doParallel=T, mc.cores=mc.cores)

#cm: 10/23/21 ditched this to see if thinning will address the oversampling issues
if(redoSort){ sortSpeciesFast(allBaseDirs,
										pointsProj=CRS('+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'),
										myprj=myprj,
										verbose=TRUE,
										overwrite=TRUE,
										mc.cores=10,
										sampleCutoffs=c(1,3,10),
										whichSpecies=whichSpecies)}

if(length(sortDone$notSorted)>1 & !all(is.na(sortDone$notSorted))) sortSpeciesBulk(allBaseDirs,pointsProj=pointsProj,myprj=myprj,doThin=TRUE,thinCutoff=25,maxCells=20000,verbose=TRUE,overwrite=FALSE,doClean=TRUE,mc.cores=mc.cores)
#--------------------------------------------------------------------
### 4.b specify which species to run
#--------------------------------------------------------------------
speciesList=c(list.files(file.path(allBaseDirs$speciesDir,'SDM'), full.names=T),list.files(file.path(allBaseDirs$speciesDir,'SDM10_20'), full.names=T))
done=tools::file_path_sans_ext(basename(gsub('__Maps__v1','', list.files(allBaseDirs$algorithmDirs$PPM$figureDir, full.names=TRUE))))
keep=!tools::file_path_sans_ext(basename(speciesList)) %in% done 
#done=checkSDMDone(allBaseDirs,speciesList,algorithm='PPM') 
#str(done)
speciesList=speciesList[keep]
# speciesList=done$problems


#====================================================================
#== 5. Model settings
#====================================================================
cbs=commonBIENSettings(myprj,world.shp,runName=runName)
toDo=cbs$toDo
toSave=cbs$toSave
toOverwrite=cbs$toOverwrite

toDo$background$customBG=FALSE
toDo$background$makeBiasedBackground=FALSE

toDo$misc$algorithm='PPM' # must match name specified in myAlgorithmNames
toDo$sampling$samplingSettings=c('noBias') # note that the target group may give wierd predictions for the demo data, since its just a few species
toDo$fitting$expertSettings=list(prob=1e-6,rate=0,skew=1e-6,shift=0)
toDo$fitting$predictorSettings=NULL
toDo$fitting$priorSettings=NULL
toDo$fitting$algorithmSettings=c('maxnet')

# these are changes to the defaults
toDo$transfer$whichFutures='all' 
toDo$transfer$otherEnvProjections=T
toDo$plotting$plotOtherEnvPred=T
toDo$plotting$whichFuturesToPlot='ip'

toSave$domain=FALSE
toDo$domain$domainTrimMethod='ecoregionRaster'
toDo$domain$coarseCrop=FALSE
toDo$domain$fixWeirdDisjunctEcoregions=FALSE
toDo$trimDomainMask=FALSE
toDo$background$maskEvalBGByBufferedPresences=FALSE
toDo$projection$maskToOccEcoNeighbors=FALSE
toDo$evaluation$evaluationModel='fullDomain'
toDo$metadata$bienMetadata=FALSE
toDo$fitting$maxTime=2200
toDo$presences$pvalSet=1e-9 # guessing at a very conservative approach to keep data.

# write out run setting used for metadata
runSettings=list(toDo,toSave,toOverwrite)
# 11/1/20 I overwrote the orignial run stuff so new BIEN defaults may have replaced the real settings. probably just about thresholding
save(runSettings,file=paste0(dirname(allBaseDirs$envDir),'/runSettings.rdata'))

#--------------------------------------------------------------------
### 6. Run SDM Models
#--------------------------------------------------------------------
# 200-300 s per species usually
j=13; (speciesCSV=speciesList[j]);
mclapply(length(speciesList):1, function(j){ #length(speciesList)
  print(j)
  out=sdmBIEN(speciesCSV=speciesList[j],
              allBaseDirs,
              toDo=toDo,
              toSave=toSave,
              toOverwrite=toOverwrite) 
  if(class(out)=='try-error') print(out)													 
  gc()
},mc.cores=mc.cores)

#------------------------------------------------------------
### 7. Run range bagging for all species 
#------------------------------------------------------------
myCustomDomainRegionRaster=='/Users/ctg/Dropbox/Projects/Viral_Chatter/Data/EcoregionRasterDinnerstein2017_QD_global.tif'
allBaseDirs2=sdmDirectories(baseDir,
                           myAlgorithmNames=myAlgorithmNames,
                           warn=FALSE,
                           mySpeciesDir=mySpeciesDir,
                           myEnvDir=myEnvDir,
                           myOffsetDirs=myOffsetDirs,
                           mySamplingDataDir=mySamplingDataDir,
                           mySamplingModelDir=mySamplingModelDir,
                           myCustomBackgroundDir=myCustomBackgroundDir,
                           myCustomBackgroundTable= myCustomBackgroundTable,	
                           myCustomDomainRegionRaster= myCustomDomainRegionRaster,												 
                           myOtherEnvDir=myOtherEnvDir,
                           offsetDataDirNames=offsetDataDirNames,
                           samplingDataDirNames=samplingDataDirNames,
                           overwrite=FALSE,
                           sortDirNames=sortDirNames,
                           otherDirs=otherDirs)

isRares=T
# if(isRares) 
speciesListHull=speciesList=c(list.files(file.path(allBaseDirs$speciesDir,'RangeBag'),full.names=T))
# if(!isRares) 
# speciesListHull=speciesList=c(list.files(file.path(allBaseDirs$speciesDir,'Samples10'),full.names=T))
future=T

doneHull=checkSDMDone(allBaseDirs2,speciesList=speciesListHull, algorithm='RangeBag',checkFigs=F) # check whether any already done
str(doneHull) 
speciesListHull=doneHull$notRun

# to run species that failed as ppm
# # speciesList=list.files(file.path(allBaseDirs$speciesDir,'Samples10'), full.names=T)
# # done=tools::file_path_sans_ext(basename(gsub('__Maps__v1','', list.files(allBaseDirs$algorithmDirs$PPM$figureDir, full.names=TRUE))))
# # keep=!tools::file_path_sans_ext(basename(speciesList)) %in% done 
# # #done=checkSDMDone(allBaseDirs,speciesList,algorithm='PPM') 
# # #str(done)
speciesList=c(list.files(file.path(allBaseDirs$speciesDir,'SDM'), full.names=T),list.files(file.path(allBaseDirs$speciesDir,'SDM10_20'), full.names=T))
done=tools::file_path_sans_ext(basename(gsub('__Maps__v1','', list.files(allBaseDirs$algorithmDirs$PPM$figureDir, full.names=TRUE))))
keep=!tools::file_path_sans_ext(basename(speciesList)) %in% done 
speciesListHull=speciesList[keep]


#Settings
csrb=commonRangeBaggingSettings(myprj,world.shp,runName=runName)
toDo=csrb$toDo
toSave=csrb$toSave
toOverwrite=csrb$toOverwrite

if(future){
	toDo$transfer$otherEnvProjections=TRUE
	toSave$shapeFileTransfer=TRUE
	toDo$transfer$whichFutures='all'
	toDo$plotting$plotOtherEnvPred=TRUE
	toDo$plotting$whichFuturesToPlot=c('ip')
	toDo$plotting$nPlotCol=4
}
if(!future) toDo$plotting$nPlotCol=2

toDo$domain$domainTrimMethod='ecoregionRaster'
toDo$fitting$maxTime=40
toDo$misc$bienMetadata=FALSE
toDo$domain$coarseCrop=FALSE
toDo$domain$fixWeirdDisjunctEcoregions=FALSE
toDo$trimDomainMask=FALSE
toDo$background$maskEvalBGByBufferedPresences=FALSE
toDo$projection$maskToOccEcoNeighbors=FALSE
toDo$evaluation$evaluationModel='fullDomain'

if(!isRares){ 
	toDo$presences$removeSpatialOutliers=TRUE
	toDo$presences$removeEnvOutliers=TRUE
} 

# write out run setting used for metadata
runSettingsRB=list(toDo,toSave,toOverwrite)
save(runSettingsRB,file=paste0(dirname(allBaseDirs2$envDir),'/runSettingsRangeBag.rdata'))


# Run Models
j=2; (speciesCSV=speciesListHull[j]);
mclapply(1:length(speciesListHull), function(j){ #length(speciesListHull)
  print(j)
  out=sdmRangeBag(speciesCSV=speciesListHull[j],
									allBaseDirs2,
									toDo=toDo,
									toSave=toSave,
									toOverwrite=toOverwrite) 
  if(class(out)=='try-error') print(out)													 
  gc()
},mc.cores=mc.cores)


