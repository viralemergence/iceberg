
sdmRangeBag=function(speciesCSV,
										 allBaseDirs,
										 toDo=toDo,
										 toSave=toSave,
										 toOverwrite=toOverwrite){

  out=try({

	verbose=toDo$misc$verbose
	algorithm='RangeBag'
	
	start.time=proc.time()
	species=unlist(strsplit(basename(speciesCSV),'.csv'))
	if(verbose>0) print(paste0(species,': ','Process started at ', round(start.time[3])))

		#-- generate random seed to be sure to get the same results for each run
	let=strsplit(strsplit(species,'_')[[1]][2],'')[[1]][1:4]
	seed=as.integer(paste(which(letters %in% let),collapse=''))
	set.seed(seed)
	
	#-------------------------------------------------------------------
	## 1. Important Objects
	#-------------------------------------------------------------------
		#-- these are the 3 main objects referenced throughout the workflow
		#-- status monitors workflow progress and catches errors

	status=.statusDF(species)

	dirs=speciesDirectories(species=species,
													allBaseDirs=allBaseDirs,
													algorithm=toDo$misc$algorithm,
													deletePastModels=toDo$misc$deletePastModels,
													writeSDPred=toSave$SDPred,
													whichFutures=toDo$transfer$whichFutures,
													saveModelObj=toSave$modelObj)
	if(class(dirs)=='try-error') {print(paste0(species,': ','dirs are messed up')) ; return()}
	# set up temp raster path to delete later
	rasterOptions(tmpdir=dirs$temp.path)
	if(verbose>2) print(paste0(species,': ','directories created '))

	#-- specify which models to fit 
		# this is for old compatibility, might be smarter to integrate this with stats storage	
	modelSettings=list(samplingSettings=toDo$sampling$samplingSettings,
										 expertSettings=toDo$fitting$expertSettings,
										 predictorSettings=toDo$fitting$predictorSettings,	
										 algorithmSettings= toDo$fitting$algorithmSettings,
										 formulaMaker=toDo$fitting$formulaMaker)

	stats=.statsStorage(species,modelSettings,type=algorithm,toDo)
	#!! only works for a single model at a time 
	stats$rbV=toDo$fitting$rbV
	stats$rbD=toDo$fitting$rbD
	stats$rbP=toDo$fitting$rbP
	stats$uniqueID=paste0(species,'__',stats$modelNames,'__', toDo$misc$runName)
	#-------------------------------------------------------------------
	## 2. Environment
	#-------------------------------------------------------------------

	### 2.a read/clean data
	#-------------------------------------------------------------------
	envFile=list.files(dirs$envDir,'.tif',full.names=TRUE)
	env=suppressWarnings(stack(envFile))

	#### 2.a.i different layers for pres and sampling
	#-- none yet

	### 2.b preselect layers
	#-------------------------------------------------------------------
	layerNames=try(read.csv(list.files(dirs$envDir,'.csv', full.names=TRUE),header=FALSE, stringsAsFactors=FALSE))
	if(!class(layerNames)=='try-error'){
		best.var=layerNames[,1] #names(env)
		names(env)=best.var
	} else {
		print('add a file with your layer names called layerNames.csv in the Env folder or bad shit could happen')
		best.var=names(env)
	}

	#-------------------------------------------------------------------
	## 3. Occurrence
	#-------------------------------------------------------------------

	### 3.a read/clean data and optionally thin
	#-------------------------------------------------------------------
	if(verbose>2) print(paste0(species,': ','cleaning occurrences '))
	prog.points=cleanOcc(speciesCSV=speciesCSV,
											 env=env,
											 doClean=FALSE,
											 writeResult=FALSE,
											 dirs=dirs)

	status$prog.points=.statusTracker(status$prog.points,prog.points, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.points==TRUE) return(status)
	pres=prog.points$pres
	# do the extract once at beginning
	pres@data=cbind(pres@data,raster::extract(env,pres,ID=FALSE))
	if(nrow(unique(pres@data))<3){ 
		message(paste0(species,':	 there are less than 3 unique points in environmental space, so skipping'))
		return()
	}

	if(any(toDo$presences$removeSpatialOutliers,
	toDo$presences$removeEnvOutliers)){
		prog.outliers=findOutlyingPoints(pres,
																		 spOutliers= toDo$presences$removeSpatialOutliers,
																		 envOutliers= toDo$presences$removeEnvOutliers,
																		 bestVar=best.var,
																		 doPlot=TRUE,
																		 plotFile=paste0(dirs$sp.diag.path, '/OccurrenceOutliers.png'),
																		 pval=toDo$presences$pvalSet,
																		 env=env,
																		 species=species)
	
		pres=prog.outliers$pres
		#--  record number tossed
		stats$n.spatial.outliers=length(prog.outliers$spatialToss)
		stats$n.env.outliers=length(prog.outliers$envToss)	
	}
	status$has.points=ifelse(nrow(pres)>0,TRUE, FALSE)
	status$n.points=nrow(pres)
	stats$n.pres=length(pres)

	# # 	status$prog.outliers=.statusTracker(status$prog.outliers,prog.outliers, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	# # 	if(!status$prog.outliers==TRUE) return(status)

	if(verbose>2) print(paste0(species,': ',length(pres),' presence points'))

	#-------------------------------------------------------------------
	## 4. Sampling
	#-------------------------------------------------------------------

	if(verbose>0) print(paste0(species,': ','trimming domain....'))

	### 4.b determine modeling domain
	#-------------------------------------------------------------------
	#add json string to store which ecoregions were used
	env=makeDomain(dirs=dirs,
								 biased.bg=NULL,
								 env=env,
								 pres=pres,
								 method=toDo$domain$domainTrimMethod,
								 writeResult=toSave$domain,
								 overwrite=toOverwrite$domain,
								 overwriteBiasMask=toOverwrite$biasMask,
								 trimBufferkm=toDo$domain$trimBufferkm,
								 maxTime=toDo$domain$domainMaxTime,
								 domainClumpDist=toDo$domain$domainClumpDist,
								 coarseCrop=toDo$domain$coarseCrop,
								 gdalCompress=toDo$misc$gdalCompress)
	if(class(env)=='try-error') stop(paste0(species,': could not trim domain, so skipping. Maybe this species did not fall within an ecoregion.')	)		 						 
		# add plot of correlations and what was tossed!
	
		#-- mask to be used later for trimming the domain for different apps	
	if(toDo$trimDomainMask){	
		ecoMask=makeDomain(dirs=dirs,
											 biased.bg=NULL,
											 env=env,
											 pres=pres,
											 method='ecoregionAndNeighborsRaster',
											 writeResult=toSave$domain,
											 overwrite=toOverwrite$domain,
											 overwriteBiasMask=toOverwrite$biasMask,
											 trimBufferkm=toDo$domain$trimBufferkm,
											 maxTime=toDo$domain$maxTrimDomainTime,
											 domainClumpDist=toDo$domain$domainClumpDist,
											 coarseCrop=TRUE,
											 fixWeirdDisjunctEcoregions= toDo$domain$fixWeirdDisjunctEcoregions,
											 gdalCompress=toDo$misc$gdalCompress,
											 saveMasks=TRUE) 
		env$ecoMask=raster::extend(ecoMask$trim.domain.2,env)	# save mask to make bg						 					 
	} else {
		env$ecoMask=1 #hack to avoid changing any later code. leads to useless calculations below but not useless cory time
	}
	
	if(toDo$presences$removeCorrelatedPredictors){
		rcp=removeCorrelatedPredictors(env,
																	 best.var=best.var,
																	 predictorsToKeepNoMatterWhat=
																	 toDo$presences$predictorsToKeepNoMatterWhat,
																	 species)
		env=rcp[['env']]; best.var=rcp[['best.var']]
	}

	status$prog.domain=.statusTracker(status$prog.domain,env, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.domain==TRUE) return(status)
	if(verbose>0) print(paste0(species,': ','done trimming domain'))

	### 4.d sample background
	#-------------------------------------------------------------------
		#-- based on the trimmed domain
	# the biased background is going to the wrong folder and woverwriting the unbaised. turned off saving for now.
	prog.bg=makeBackground(env=env,
												 dirs=dirs,
												 stats=stats,
												 biased.bg=NULL,
												 maxBGFit=toDo$background$maxBGFit,
												 writeResult=toSave$background,
												 overwrite=toOverwrite$background)
	if(class(prog.bg$bg)=="character"){
		if(prog.bg$bg=='noData'){
			print(paste0(species, ': bailing because <10 background found '))
			return(status)
		}
	}
	bg=prog.bg$bg 

	status$prog.bg=.statusTracker(status$prog.bg,prog.bg, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.bg==TRUE) return(status)

	#-------------------------------------------------------------------
	## 6. SDM
	#-------------------------------------------------------------------

	### 6.a full model
	#-------------------------------------------------------------------
	if(verbose>0) print(paste0(species,': ','fitting full model ... '))
	#env=env[[-c(grep('targetBG.mask',names(env)), grep('trim.domain',names(env)))]]	
	# TODO record number of models failed	
	prog.fit.full=sdmRangeBag2(dirs,
														 pres.fit=pres,
														 k=NULL,
														 name='full',
														 stats=stats,
														 maxTime=toDo$fitting$maxTime,
														 best.var=best.var,
														 saveModelObj=toSave$modelObj,
														 runModel=TRUE,
														 evalFork=toDo$fitting$evalFork)
	if(prog.fit.full$mods=='noData') return(status)	
	
	status$prog.fit.full=.statusTracker(status$prog.fit.full,	prog.fit.full,toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.fit.full==TRUE) return(status)
	if(is.na(prog.fit.full$mods)) {full.mods=NULL} else {	full.mods=prog.fit.full$mods }
	if(toDo$fitting$runFullModel) stats=prog.fit.full$stats
	if(verbose>0) print(paste0(species,': ','done fitting full model '))

	#-------------------------------------------------------------------
	## 7. Projection
	#-------------------------------------------------------------------

	if(verbose>0) print(paste0(species,': ','projecting ...'))

	### 7.a. Project full model
	#-------------------------------------------------------------------
	prog.proj.full=projectSDM(dirs=dirs,
														mod.out=full.mods,
														stats=stats,
														projEnv=env,
														name=paste0(stats$species[1],'__full'),
														fileNameSuffix='_fullDomain',
														best.var=best.var)
	status$prog.proj.full=.statusTracker(status$prog.proj.full,stats, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.proj.full==TRUE) return(status)
	stats=prog.proj.full$stats
	if(toDo$misc$rmObjs) rm(prog.proj.full)

	status$time.project=round((proc.time() - start.time)[3]- status$time.input.prep - status$time.fit)
	if(verbose>0) print(paste0(species,': ','done projecting '))
	if(verbose>0) print(paste0(species,': ',"time: ",status$time.project))
	
	### 7.b. Mask continuous preds to relevant ecoregions (occ+neighbors)
	#-------------------------------------------------------------------
	# basically we won't use the continuous predictions based on buffered points for anything; they'll just hang around in case we need them in the future. 
		   			
  fullDomainFiles=list.files(dirs$sp.pred.path,pattern='full', full.names=TRUE)	
  # 10/12/21: i think this can just be plugged in from the ppm function to replace gthe commented stuff below
	if(toDo$projection$maskToOccEcoNeighbors){
		myMasks=raster::stack(list.files(dirs$sp.mask.path,full.names=T))[[1]]  
		occEcoNeighFiles=sub('_fullDomain.tif','_occEcoNeigh.tif', fullDomainFiles,fixed=T)    	
			#-- occupied and neighbors
		thisMask=raster::reclassify(myMasks,c(.99,2.01,1))
		a=lapply(seq_along(fullDomainFiles),function(x){
			fold.r=suppressWarnings(raster(fullDomainFiles[x]))
			fold.r2=raster::crop(fold.r,thisMask)
			fold.r3=raster::mask(fold.r2,thisMask)
			.writeRasterCM(fold.r3,occEcoNeighFiles[x],toDo$misc$gdalCompress)
		})
	}
	#Maybe this should happen below after i've taken the ensemble means. i think i only really need the folds for the fullDomain stored and for occEcoNeigh to evaluate
		#-- occupied only
	if(toDo$trimDomainMask){
		myMasks=raster::stack(list.files(dirs$sp.mask.path,full.names=T))[[1]]  
		occEcoFiles=sub('_fullDomain.tif','_occEco.tif',fullDomainFiles,fixed=T)
		thisMask=raster::reclassify(myMasks,c(1.01,Inf,NA))
		a=lapply(seq_along(fullDomainFiles),function(x){
			fold.r=suppressWarnings(raster(fullDomainFiles[x]))
			fold.r2=raster::crop(fold.r,thisMask)
			fold.r3=raster::mask(fold.r2,thisMask)
			.writeRasterCM(fold.r3,occEcoFiles[x],toDo$misc$gdalCompress)
		})
	}
# 	occEcoNeighFiles=sub('_fullDomain.tif','_occEcoNeigh.tif', fullDomainFiles,fixed=TRUE)
#  	
# 		#-- by convention, the first layer is the ecoregion one
# 	myMasks=suppressWarnings(raster::stack(list.files( dirs$sp.mask.path,full.names=TRUE)))[[1]]      	
# 		#-- occupied and neighbors
# 	thisMask=raster::reclassify(myMasks,c(.99,2.01,1))
# 	lapply(seq_along(fullDomainFiles),function(x){
# 		fold.r=raster(fullDomainFiles[x])
# 		fold.r2=raster::crop(fold.r,thisMask)
# 		fold.r3=raster::mask(fold.r2,thisMask)
# 		.writeRasterCM(fold.r3,occEcoNeighFiles[x],toDo$misc$gdalCompress)
# 	})
# 		#-- occupied only
# 	occEcoFiles=sub('_fullDomain.tif','_occEco.tif',fullDomainFiles,fixed=TRUE)
# 	thisMask=raster::reclassify(myMasks,c(1.01,Inf,NA))
# 	lapply(seq_along(fullDomainFiles),function(x){
# 		fold.r=raster(fullDomainFiles[x])
# 		fold.r2=raster::crop(fold.r,thisMask)
# 		fold.r3=raster::mask(fold.r2,thisMask)
# 		.writeRasterCM(fold.r3,occEcoFiles[x],toDo$misc$gdalCompress)
# 	})
	#-------------------------------------------------------------------
	## 8. Evaluate
	#-------------------------------------------------------------------

	if(verbose>0) print(paste0(species,': ','evaluating ...'))

	### 8.a Evaluate continuous predictions
	#-------------------------------------------------------------------
	prog.eval.cont=evaluateContFullCV(pres=pres,
																		bg=bg,
																		dirs=dirs,
																		stats=stats,
																		full.mods=full.mods,
																		maxBGEval=toDo$background$maxBGEval,
																		whichModel='full',
																		modelNameGrep= toDo$evaluation$evaluationModel)

	status$prog.eval.cont=.statusTracker(status$prog.eval.cont, prog.eval.cont,  toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.eval.cont==TRUE) return(status)
	stats=prog.eval.cont$stats

	### 8.b Evaluate binary predictions
	#-------------------------------------------------------------------
	prog.eval.bin=evaluateBinFullCV(pres=pres,
																	bg=bg,
																	dirs=dirs,
																	stats=stats,
																	thresholds=toDo$threshold$thresholds,
																	whichModel='full',
																	modelNameGrep= toDo$evaluation$evaluationModel)
	stats=prog.eval.bin
	status$prog.eval.bin=.statusTracker(status$prog.eval.bin,stats,toDo)
	if(!status$prog.eval.bin==TRUE) return(status)

	### 8.c Find best models
	#-------------------------------------------------------------------
	# CM 3/31/21 - seems to work but not used for anything
	# 	prog.best=bestModelIndicator(stats=stats)
	# 	stats=prog.best
	# 	status$prog.find.best=.statusTracker(status$prog.find.best,stats, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	# 	if(!status$prog.find.best==TRUE) return(status)

	status$time.eval=round((proc.time() - start.time)[3]- status$time.input.prep - status$time.fit- status$time.project)
	if(verbose>0) print(paste0(species,': ','done evaluating'))
	if(verbose>2) print(paste0(species,': ',"time: ",status$time.eval))

	#-------------------------------------------------------------------
	## 9.  Transfer model to new conditions
	#-------------------------------------------------------------------
	### 9.a. Project CV models
	#-------------------------------------------------------------------
		#-- this only project to the full domain
	prog.transfer=NULL
	if(toDo$transfer$otherEnvProjections){
		if(verbose>0) print(paste0(species,': ','transferring ... '))
		# environment for transfer
		env.new=NULL
		prog.transfer=try({
		 	new.envs=list.files(dirs$otherEnvDir,full.names=TRUE,pattern='tif')
		 	if(any(toDo$transfer$whichFutures=='all')){
				keep=seq_along(basename(new.envs))
			} else {
		 	 	keep=c(mapply(function(x) grep(x,basename(new.envs)),toDo$transfer$whichFutures))
		 	}
		 	new.envs=new.envs[keep]
		 	
			for(ii in 1:length(new.envs)){
				print(paste0(species,': ','  projecting ',tools::file_path_sans_ext(basename(new.envs[ii]))))
				env.new=raster::stack(new.envs[ii])
				layerNames=try(read.csv(list.files(dirs$envDir,'.csv', full.names=TRUE),header=FALSE, stringsAsFactors=FALSE))
				if(!class(layerNames)=='try-error'){ 
					names(env.new)=layerNames[,1]
				} else {
					names(env.new)=best.var 
				}
				# subset just the best variables used in modeling for speed
				env.new=env.new[[best.var]] 
				
				aa=.cropcrop(env.new,env[['trim.domain']])
				env.new=raster::stack(aa[[1]],aa[[2]])
				#cv.stats=vector('list',length(unique(pres$folds)))
				#for(k in sort(unique(pres$folds))){
				transfer=tools::file_path_sans_ext(basename(new.envs[ii]))
				tmp.proj=projectSDM(dirs=dirs,
														mod.out=full.mods,
														stats=stats,
														projEnv=env.new,
														priors=NA,
														name=paste0(stats$species[1],'__full'),
														best.var=best.var,
														transfer=transfer,
														fileNameSuffix='_fullDomain',
														algorithm=toDo$misc$algorithm,
														verbose=1)
				 #}
			}
		})
		if(toDo$misc$rmObjs) rm(env.new)
	}
	
	status$prog.transfer=.statusTracker(status$prog.transfer,stats, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.transfer==TRUE) return(status)

	status$time.transfer=round((proc.time() - start.time)[3]- status$time.input.prep - status$time.fit- status$time.project- status$time.eval)
	if(verbose>0) print(paste0(species,': ','done transferring'))
	if(verbose>2) print(paste0(species,': ',"time: ",status$time.transfer))

	#-------------------------------------------------------------------
	## 10. Threshold predictions
	#-------------------------------------------------------------------
	### 10.a. under fitting conditions
	#-------------------------------------------------------------------
		#-- this is only for full domains; we mask later
		#-- note that the shapefile is for the full domain, which is dumb, but we only use it for quick plots so doesn't really matter.
	if(toDo$threshold$thresholding){
		if(verbose>0) print(paste0(species,': ','thresholding ... '))
		prog.threshold=threshold(dirs=dirs,
														 stats=stats,
														 pres=pres,
														 env=env,
														 saveShapeFile=toSave$shapeFile,
														 makeThresholdMaps= toDo$threshold$thresholds,
														 whichModel='full',
														 modelNameGrep='fullDomain',
														 algorithm=toDo$misc$algorithm,
														 makeSparse=toDo$threshold$makeSparseMatrix,
														 gdalCompress=toDo$misc$gdalCompress)

		status$prog.threshold= .statusTracker(status$prog.threshold,stats, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
		if(!status$prog.threshold==TRUE) return(status)
	}
	
	### 10.b. under new conditions 
	#-------------------------------------------------------------------
	if(toDo$transfer$otherEnvProjections & toDo$threshold$thresholding) {
		 if(verbose>0) print(paste0(species,': ','thresholding transferred models  ... '))
		 prog.threshold.new=try({
			 new.envs=list.files(dirs$otherEnvDir,full.names=TRUE,pattern='tif')
			 if(any(toDo$transfer$whichFutures=='all')){
				 keep=seq_along(basename(new.envs))
			 } else {
				 keep=c(mapply(function(x) grep(x,basename(new.envs)),toDo$transfer$whichFutures))
			 }
			 for(ii in keep){
			 	 envName=tools::file_path_sans_ext(basename(new.envs[ii]))
				 tmp.thresh=threshold(dirs=dirs, 
															stats=stats,
															pres=pres,
															env=env,
															saveShapeFile=toSave$shapeFileTransfer,
															whichModel='full',
															transfer=envName,
															makeThresholdMaps= toDo$threshold$makeThresholdMaps,
															algorithm=toDo$misc$algorithm,
															modelNameGrep='fullDomain',
															makeSparse=toDo$threshold$makeSparseMatrix,
															gdalCompress=toDo$misc$gdalCompress)
			 }
		 })
		status$prog.threshold.new=.statusTracker(status$prog.threshold.new, stats,toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
		if(!status$prog.threshold.new==TRUE) return(status)
	}

	status$time.threshold=round((proc.time() - start.time)[3]- status$time.input.prep - status$time.fit- status$time.project- status$time.eval- status$time.transfer)
	if(verbose>0) print(paste0(species,': ','done thresholding'))
	if(verbose>2) print(paste0(species,': ',"time: ",status$time.transfer))
	
	#-------------------------------------------------------------------
	## 11. Masking predictions
	#-------------------------------------------------------------------
		#-- the continuous preductions on the fitting data were already masked above to allow evaluation to occupied + neighbors
	
	# ### 11.a. Mask binary maps to Occupied + neighboring ecoregions 
	#-------------------------------------------------------------------
	if(toDo$trimDomainMask){
		myMasks=raster::stack(list.files(dirs$sp.mask.path,full.names=TRUE))[[1]]      	
		thisMask=raster::reclassify(myMasks,c(.99,2.01,1)) 
		maskBinaryMaps(thisMask,toDo,dirs,newMapName='occEcoNeigh',species)
		if(verbose>0) print(paste0(species,': done masking to occEcoNeigh'))

		### 11.b. Mask binary maps to Occupied ecoregions 
		#-------------------------------------------------------------------
		myMasks=raster::stack(list.files(dirs$sp.mask.path,full.names=TRUE))[[1]]      	
		thisMask=raster::reclassify(myMasks,c(1.01,Inf,NA)) # use occupied only
		maskBinaryMaps(thisMask,toDo,dirs,newMapName='occEco',species)
		if(verbose>0) print(paste0(species,': done masking to occEco'))
	}
	#-------------------------------------------------------------------
	## 12.  Plotting
	#-------------------------------------------------------------------
	prog.plot=compareSdmStatsPlot(stats=stats,
	                              dirs=dirs,
	                              range.poly=NULL,
	                              pres=pres,
	                              env=env,
	                              speciesInfo=FALSE,
	                              legendLoc='panel',
	                              logscaleZ=TRUE,
	                              modelNameGrep= toDo$evaluation$evaluationModel,
	                              toDo=toDo)     

	status$prog.plot=.statusTracker(status$prog.plot,prog.plot, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)

	if(!status$prog.plot==TRUE) return(status)
	status$prog.plot=.statusTracker(status$prog.plot,stats$plotPath)
	if(verbose>0) print(paste0(species,': ','done plotting'))

	#-------------------------------------------------------------------
	## 13.  Clean up / Record model metadata
	#-------------------------------------------------------------------
	# later set this up for rmm objects	
	if(toDo$misc$bienMetadata)	md=makeBIENMetadata(stats=stats,dirs=dirs,toDo=toDo)

 	if(verbose>2) print(paste0(species,': ',"writing statistics and cleaning"))
 	status$time.all=stats$runTime=round((proc.time() - start.time)[3])
 	status$finished=TRUE
 	e=exitSDMWorkflow(dirs=dirs,stats=stats,status=status,toSave=toSave, toDo=toDo)

 	# clean up tmp files
	unlink(rasterOptions(overwrite=TRUE)$tmpdir,recursive=TRUE)
		
	if(verbose>0) print(paste0(species,': ',"finished model in ", status$time.all))
	if(verbose>2) print(paste0(species,': ','-------------------------'))

	return(status)
  })
  return(out)
}

