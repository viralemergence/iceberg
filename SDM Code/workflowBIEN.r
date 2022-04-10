

sdmBIEN=function(speciesCSV,
                 allBaseDirs,
                 algorithm,
								 toDo=toDo,
								 toSave=toSave,
								 toOverwrite=toOverwrite){

  out=try({

	verbose=toDo$misc$verbose

	start.time=proc.time()
	species=unlist(strsplit(basename(speciesCSV),'.csv'))
	if(verbose>0) print(paste0(species,': Process started at ', round(start.time[3])))

	# generate random seed to be sure to get the same results for each run
	let=strsplit(strsplit(species,'_')[[1]][2],'')[[1]][1:4]
	if(all(is.na(let))) let=strsplit(species,'')[[1]][1:4]
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
	if(class(dirs)=='try-error') {
		outMessage=(paste0(species,': dirs are messed up')) 
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}
	
		#-- set up temp raster path to delete later
	rasterOptions(tmpdir=dirs$temp.path)

	#-- specify which models to fit 
		# this is for old compatibility, might be smarter to integrate this with stats storage	
	modelSettings=list(samplingSettings=toDo$sampling$samplingSettings,
										 expertSettings=toDo$fitting$expertSettings,
										 predictorSettings=toDo$fitting$predictorSettings,	
										 algorithmSettings= toDo$fitting$algorithmSettings,
										 formulaMaker=toDo$fitting$formulaMaker)
	#modelSpecs(modelSettings)
	stats=.statsStorage(species,modelSettings,type=toDo$misc$algorithm)
	stats$uniqueID=paste0(species,'__',stats$modelNames,'__', toDo$misc$runName)
	# add presence and background sample sizes
	
	#-------------------------------------------------------------------
	## 2. Environment
	#-------------------------------------------------------------------

	### 2.a read/clean data
	#-------------------------------------------------------------------
	envFile=list.files(dirs$envDir,'.tif',full.names=T)
	env=suppressWarnings(stack(envFile))

	#### 2.a.i different layers for pres and sampling
	#-- none yet

	### 2.b preselect layers
	#-------------------------------------------------------------------
	layerNames=try(read.csv(list.files(dirs$envDir,'.csv', full.names=T),header=FALSE, stringsAsFactors=FALSE))
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
	prog.points=cleanOcc(speciesCSV=speciesCSV,
											 env=env,
											 doClean=FALSE,
											 writeResult=FALSE,
											 dirs=dirs)

	status$prog.points=.statusTracker(status$prog.points,prog.points, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.points==TRUE) {
		outMessage=paste0(species,":  Skipping; cleaning occurrences failed with:",status$prog.points)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}
	pres=prog.points$pres
	# do the extract once at beginning; this is duplicated below; doesn't hurt anything but should be cleand up
	pres@data=cbind(pres@data,raster::extract(env,pres,ID=FALSE))
	keep=complete.cases(pres@data)
	pres=pres[keep,]
	if(length(pres) < toDo$presences$minSampleSize){
		outMessage=paste0(species,":  Skipping; after extracting env at presences reduced sample size to ",length(pres) ," which is less than threshold you set of ", toDo$presences$minSampleSize)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}
	# !! try rgeos::gCentroid
	if(any(toDo$presences$removeSpatialOutliers, toDo$presences$removeEnvOutliers)){
		prog.outliers=findOutlyingPoints(pres,
																		 spOutliers= toDo$presences$removeSpatialOutliers,
																		 envOutliers= toDo$presences$removeEnvOutliers,
																		 bestVar=best.var,
																		 doPlot=T,
																		 plotFile=paste0(dirs$sp.diag.path, '/', species,'__OccurrenceOutliers.png'),
																		 pval=toDo$presences$pvalSet,
																		 env=env,
																		 species=species)
		stats$n.spatial.outliers=length(prog.outliers$spatialToss)
		stats$n.env.outliers=length(prog.outliers$envToss) # record number tossed
		pres=prog.outliers$pres
		keep=complete.cases(pres@data)
		pres=pres[keep,]
			# if we tossed any, reassign folds to be sure we didn't toss an entire fold
		if(any(keep==FALSE)) pres=spatialStratify(pres,stats) 
	}
	
	if(length(pres) < toDo$presences$minSampleSize) {
		outMessage=paste0(species,":  Skipping; after extracting env at presences reduced sample size to less than ", toDo$presences$minSampleSize)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}

	stats$n.pres=length(pres)
	status$has.points=ifelse(nrow(pres)>0,TRUE, FALSE)
	status$n.points=nrow(pres)
	status$prog.outliers=.statusTracker(status$prog.outliers,prog.outliers, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.outliers==TRUE) {
		outMessage=paste0(species,":  Skipping; cleaning outliers failed with: ",status$prog.outliers)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}

	### 3.b filter/prep occurrences
	#-------------------------------------------------------------------
	#== Stratify 5 folds for CV
	pres=spatialStratify(pres,stats)

	status$prog.cv=.statusTracker(status$prog.cv,pres,toSave=toSave, toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.cv==TRUE) return(status)
	if(length(pres)==1) status$do.cv=FALSE
	if(length(pres)>1) status$do.cv=TRUE

	if(verbose>2) print(paste0(species,': ',length(pres),' presence points'))

	#-------------------------------------------------------------------
	## 4. Sampling
	#-------------------------------------------------------------------

	### 4.a read/clean data
	#-------------------------------------------------------------------
	#== define which species contribute to sampling
	biased.bg=biasedBackground(dirs=dirs,
														 env=env,
														 pointsProj=toDo$presences$pointsProj,
														 overwrite=toOverwrite$biasedBackground,
														 writeResult=toSave$biasedBackground,
														 customBG=toDo$background$customBG,
														 makeBiasedBackground= toDo$background$makeBiasedBackground)
	if(is.null(biased.bg) & toDo$background$makeBiasedBackground){
		warning(paste0(species,": I couldn't make the target group background so I'm turning those settings off and continuing with the models ignoring sampling bias"))
		toDo$background$makeBiasedBackground=F
		toDo$background$customBG=T
		toss=which(stats$samplingSet=='targetBG')
		stats=stats[-toss,]
	}
	status$prog.bias=.statusTracker(status$prog.bias,biased.bg, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
  if(!status$prog.bias==TRUE) {
  	outMessage=paste0(species,":  Skipping; biasedBackground failed with:", status$prog.bias)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
  }

	### 4.b determine modeling domain
	#-------------------------------------------------------------------
	env=makeDomain(dirs=dirs,
								 biased.bg=biased.bg,
								 env=env,
								 pres=pres,
								 method=toDo$domain$domainTrimMethod,
								 writeResult=toSave$domain,
								 overwrite=toOverwrite$domain,
								 overwriteBiasMask=toOverwrite$biasMask,
								 trimBufferkm=toDo$domain$trimBufferkm,
								 maxTime=toDo$domain$maxTrimDomainTime,
								 domainClumpDist=toDo$domain$domainClumpDist,
								 coarseCrop=toDo$domain$coarseCrop,
								 fixWeirdDisjunctEcoregions= toDo$domain$fixWeirdDisjunctEcoregions,
								 gdalCompress=toDo$misc$gdalCompress,
								 saveMasks=FALSE)
	if(class(env)=='try-error') {
		outMessage=paste0(species,': Skipping; could not trim domain. Maybe this species did not fall within an ecoregion, but just guessing.')	
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}	 						
			# this is the mask to be used later for trimming the domain for different apps	
	if(toDo$trimDomainMask){
		ecoMask=makeDomain(dirs=dirs,
											 biased.bg=biased.bg,
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
											 saveMasks=T) 
		env$ecoMask=raster::extend(ecoMask$trim.domain.2,env)	# save mask to make bg						 					 
	} else {
		env$ecoMask=1 #hack to avoid changing any later code. leads to useless calculations below but not useless cory time
	}
	# add plot of correlations?
	if(toDo$presences$removeCorrelatedPredictors){
		rcp=removeCorrelatedPredictors(env,
																	 best.var=best.var,
																	 predictorsToKeepNoMatterWhat=
																	 toDo$presences$predictorsToKeepNoMatterWhat,
																	 species)
		env=rcp[['env']]; best.var=rcp[['best.var']]
	}
	
	if(toDo$domain$standardizePredictors){
		print(paste0(species,': standardizing predictors'))
		env.masked=raster::mask(env,env$ecoMask)
		env.means=raster::cellStats(env.masked[[best.var]],'mean')
		env.sd=raster::cellStats(env.masked[[best.var]],'sd')
		out=raster::stack(lapply(seq_along(best.var),function(i) calc(env[[best.var[i]]],function(x) (x-env.means[best.var[i]])/env.sd[best.var[i]])))
		ln=names(env)
		env=raster::stack(out,env[[which(!names(env)%in%best.var)]])
		names(env)=ln
		#write these out in case we want to project somewhere else
		stats$envMeansOverDomain=rjson::toJSON(env.means)
		stats$envSDOverDomain=rjson::toJSON(env.sd)
		rm(out)
	}

  stats$predictors.kept=rjson::toJSON(best.var)
  stats$predictors.tossed=rjson::toJSON(rcp[['tossed']])

	status$prog.domain=.statusTracker(status$prog.domain,env, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.domain==TRUE) {
		outMessage=paste0(species,":  Skipping; domain failed  with: ", status$prog.domain)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)	
	}
	if(verbose>0) print(paste0(species,': done trimming domain'))

	### 4.d sample background
	#-------------------------------------------------------------------
	#== based on the trimmed domain
	prog.bg=makeBackground(env=env.masked, # borrow this from above
												 dirs=dirs,
												 stats=stats,
												 biased.bg=biased.bg,
												 maxBGFit=toDo$background$maxBGFit,
												 writeResult=toSave$background,
												 overwrite=toOverwrite$background)
	if(class(prog.bg$bg)=="character"){
		if(prog.bg$bg=='noData'){
			outMessage=paste0(species, ': bailing because <10 background found ')
			message(outMessage)
			exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
			return(outMessage)
		}
	}
	status$prog.bg=.statusTracker(status$prog.bg,prog.bg, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.bg==TRUE) {
		outMessage=paste0(species,":  Skipping; background failed with: ", status$prog.bg)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}
	
	#-- make separate background for evaluation
	if(toDo$background$maskEvalBGByBufferedPresences){
			#-- rasterize points, buffer them by 25km, sample 100x more than pres (or 10k, whichever is less)
		pres.r=pres %>% SpatialPoints %>% raster::rasterize(env,field=rep(1,length(pres)))
		# plot(pres.r); points(pres)
		pres.rb=raster::buffer(pres.r,width=2e4)
		if(	sum(values(pres.r),na.rm=T) == sum(values(pres.rb),na.rm=T) ){
			stop(paste0(species,': buffering presences didnt work'))
		}
		presMask=suppressWarnings(raster(pres.rb))
		values(presMask)=1
		values(presMask)[values(pres.rb)==1]=NA
		env.masked2=raster::mask(env.masked,presMask)
		maxBGEval=min(1e4,toDo$background$maxBGEvalMultiplier*length(pres))
		prog.pres.masked.bg=makeBackground(env=env.masked2, 
																			 dirs=dirs,
																			 stats=stats,
																			 biased.bg=biased.bg,
																			 maxBGFit=maxBGEval,
																			 writeResult=toSave$background,
																			 overwrite=toOverwrite$background)
	} else {bg.eval=prog.bg$bg}	
	#-------------------------------------------------------------------
	### 5. Expert/Prior
	#-------------------------------------------------------------------

	### 5.a read/clean data
	#-------------------------------------------------------------------
	#prog.prior=uniformPrior(env,stats)
	fx_priors=env[[1]]
	fx_not_nas=!is.na(fx_priors)
	values(fx_priors)[values(fx_not_nas)]=1/sum(values(fx_not_nas), na.rm=T)
	# center priors to avoid maxnet issues.
	fx_priors=log(fx_priors)
	fx_prior.means=try(raster::cellStats(fx_priors,mean)) # error from large rasters?
	if(class(fx_prior.means)=='try-error'){
		warning('stop - the uniform prior does not fit in memory to use with cellStats')
	}
	fx_prior.means=raster::cellStats(fx_priors,mean)
	raster::values(fx_priors)=do.call('cbind',lapply(1:length(fx_prior.means), function(x) raster::values(fx_priors[[x]])+abs(fx_prior.means[x])))
		names(fx_priors)='prior_unif'
	stats$prior.file=NULL
	prog.prior = list(priors=fx_priors,stats=stats)

	priors=prog.prior$priors
	stats=prog.prior$stats
	
	status$prog.prior=.statusTracker(status$prog.prior,priors, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.prior==TRUE) {
		outMessage=paste0(species,":  Skipping; prior failed with: ", status$prog.prior)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}
	status$has.expert=FALSE
	status$prog.expert=T 
	status$prog.prior=T
	status$prog.eval.prior=T
	status$n.models=nrow(stats)
	range.poly=NULL
	
	if(toDo$misc$rmObjs) rm(prog.prior)
	status$time.input.prep=round((proc.time() - start.time)[3])
	if(verbose>0) print(paste0(species,': done input prep in ',status$time.input.prep,'s'))

	# do extract with all the layers created. these objects should all be completely homologous except for folds in pres
	env2=raster::stack(env,priors)
	bg=SpatialPointsDataFrame(prog.bg$bg, data.frame(raster::extract(env2,prog.bg$bg),prog.bg$bg@data[,c('sp','equalWeights')]))
	if(toDo$background$maskEvalBGByBufferedPresences){
		bg.eval=SpatialPointsDataFrame(prog.pres.masked.bg$bg, data.frame(raster::extract(env2,prog.pres.masked.bg$bg),prog.pres.masked.bg$bg@data[,c('sp','equalWeights')]))
	} else {bg.eval=bg}	
	if(!is.null(biased.bg)){
		biased.bg=SpatialPointsDataFrame(prog.bg$biased.bg, data.frame(raster::extract(env2,prog.bg$biased.bg),prog.bg$biased.bg@data[,c('sp','equalWeights')]))
	}

	pres=SpatialPointsDataFrame(coordinates(pres), data.frame(raster::extract(env2,pres), pres@data[,c('sp','equalWeights','folds')]))
		# in case presences fall outside the environmental rasters (giving NA env values)
	keep=complete.cases(pres@data)
	pres=pres[keep,]
		# if we tossed any, reassign folds to be sure we didn't toss an entire fold
	if(any(keep==FALSE)) pres=spatialStratify(pres,stats)
	if(length(pres) < toDo$presences$minSampleSize){
		outMessage=paste0(species,":  Skipping; after extracting env at presences reduced sample size to less than ", toDo$presences$minSampleSize)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}
	if(toDo$misc$rmObjs) rm(prog.bg)
	if(toDo$misc$rmObjs) rm(env2)
	
	#-------------------------------------------------------------------
	## 6. SDM
	#-------------------------------------------------------------------

	### 6.a full model
	#-------------------------------------------------------------------
	full.mods=NULL
	
	### 6.b. cross validation
	#-------------------------------------------------------------------

	if(verbose>0) print(paste0(species,': fitting cv models'))
							
	cv.mods=try({
		cv.mods.k=vector('list',length(unique(pres$folds)))
		for(k in sort(unique(pres$folds))){
			if(verbose>3) print(paste0(species,':     fold ', k,':'))
			# need at least 3 fitted models to continue, so bail early if that's not going to happen so we're not waiting for slow models we won't use anyhow
			if(k>1){
				cond=(sum(sapply(cv.mods.k[1:(k-1)],function(x) is.null(x))))
				if(cond>2){
					outMessage=paste0(species,': 3 models already failed (maybe due to exceeding time limit) so I wont even bother to keep trying the others')
					message(outMessage)
					return(outMessage)
				}
			}			
			model=sdmPPMOffset.cv(dirs=dirs,
														pres.fit=pres,
														k=k,
														bg=bg,
														biased.bg=biased.bg,
														stats=stats,
														maskName=toDo$fitting$maskName,
														formulaMaker=toDo$fitting$formulaMaker,
														maxTime=toDo$fitting$maxTime,
														best.var=best.var,
														saveModelObj=toSave$modelObj,
														evalFork=toDo$fitting$evalFork,
														standardize=!toDo$domain$standardizePredictors,
														quadraticUpperLimit0= toDo$fitting$quadraticUpperLimit0,
														lambdaSeq=toDo$fitting$lambdaSeq, 
														minPresences= ceiling(.6*toDo$presences$minSampleSize))		
			# retrofitting for consistency; would be better if this were the output above
			if(is.null(model)) {model=NA; status$some.models.missing=TRUE
			} else if(model[[1]]$n.features.lambda.1se==0) { 
				model=NA; status$some.models.missing=TRUE
			}  
			cv.mods.k[[k]]=model
		}
		cv.mods.k
	})
	# add a catch in case predictions are uniform
	failed=sapply(cv.mods,function(x){ is.na(x) })
	if(all(failed)){
		outMessage=paste0(species,': Skipping; glmnet ran but no models could be fit')
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}
	if(sum(!failed)<3){
		outMessage=paste0(species,': Skipping; more than 2 folds failed')
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}

	status$prog.fit.k=.statusTracker(status$prog.fit.k,cv.mods, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.fit.k==TRUE) {
		outMessage=paste0(species,":  Skipping; fitting failed ", status$prog.fit.k)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}
	#cv.mods=prog.fit.k
	#if(toDo$misc$rmObjs) rm(prog.fit.k) 

	#-- store metrics in stats (combine metrics for each fold)
		# might be nice to use json some day...
	n.per.fold=pres$folds %>% table
	stats$eff.n.pres.fit=n.per.fold %>% mean %>% round
	stats$n.pres.by.fold=rjson::toJSON(n.per.fold %>% unname)
	stats$used.lambda.min.by.fold=rjson::toJSON(sapply(cv.mods,function(x) ifelse(is.na(x),NA,x[[1]]$used.lambda.min)))
	# only working for 1 model.extract
	stats$model.fit.ok=rjson::toJSON(which(sapply(cv.mods,function(x){ any(class(x[[1]])=='cv.glmnet')})))
	# TODO: probably doesnt work for multiple models
	stats$n.features.by.fold= rjson::toJSON(.getFromNestedLists(cv.mods,'n.features.lambda.1se') %>% unname)
	stats$formula.fit.attempt.by.fold= rjson::toJSON(.getFromNestedLists(cv.mods,'fitAttempt')%>% unname)
	#-- plot regularization curves (combine metrics for each fold)
	if(toDo$plotting$regCurve){
		plot.f=paste0(dirs$sp.diag.path,'/',species,'_RegPath.pdf')
		pdf(plot.f,w=4*nrow(stats),h=10)
			par(mfrow=c(5,nrow(stats)),mar=c(3,6,3,1))
			lapply(seq_along(cv.mods),function(x){ lapply(seq_along(cv.mods[[x]]),function(y) { 	
					if(any(class(cv.mods[[x]][[y]])=='cv.glmnet')) { 
						plot(cv.mods[[x]][[y]],cex.lab=1,cex.axis=1)
						mtext(paste0('fold ',x,' model ',y),3,line=-2) 
				}})
			})
		dev.off() # ;system(paste0('open ',plot.f))
	}

	status$time.fit=round((proc.time() - start.time)[3] - status$time.input.prep)
	if(verbose>0) print(paste0(species,': done fitting cv models in ',status$time.fit,'s'))

	#-------------------------------------------------------------------
	## 7. Projection
	#-------------------------------------------------------------------

	if(verbose>0) print(paste0(species,': projecting ...'))

	### 7.b.  Project k-fold CV
	#-------------------------------------------------------------------
	prog.proj.k=try({
		# create a mask that prevent predictions beyond N sd of the training data
					# subset just the best variables used in modeling for speed
		# don't normalize the folds - leave the predictions comparable and take their weighted avg for ensemble. then normalize, and record the normalization constant (wait to normalize so i only have 1 constant to rescale the threshold by). thresholds will be estimated, and to threshold future maps, multiply the threshold by the normalization constant to get future binary maps
		env.new=env[[best.var]] 
		extrapLimits=NULL
		if(!is.null(toDo$projection$extrapolationRule)){
			if(toDo$projection$extrapolationRule=='trainingPresSD'){
				pres.dat=raster::extract(env.new,pres)
				e.mask=extrapolationMask(env.new,
															 	 dat=pres.dat,
															 	 rule='sd',
															 	 val=toDo$projection$extrapolationValue)															 	 
				extrapLimits=e.mask$extrapLims											 	 
			}
			if(toDo$domain$standardizePredictors){
				all.masks=raster::stackApply(e.mask$masks, rep(1,length(best.var)), sum)==length(best.var)
				env.new=stack(raster::mask(env.new,all.masks,maskvalue=0), env[[which(!names(env) %in% best.var)]])
			} else {
				env.new=stack(e.mask$maskedEnv,env[[which(!names(env) %in% best.var)]])
			}
		} 
		if(toDo$projection$projectPresentInOnlyUsedEcoregions){

			if(!any(na.omit(values(env$trim.domain))==2)){ print("Note - you chose toDo$projection$projectPresentInOnlyUsedEcoregions=TRUE to avoid projecting into unoccupied (adjacent) ecoregions in with the fitting data (present) but you did not choose a domain that included adjacent ecoregions. So you're good, but just wanted to let you know that there aren't any adjacent ecoregions in your domain.")
			} else {
				myMask=env$trim.domain
				myMask[myMask==2]=NA
				env.new=raster::mask(env.new,myMask)
				if(toDo$misc$rmObjs) rm(myMask)
			}
		}

		cv.stats=vector('list',length(unique(pres$folds)))
			# TODO: might not work for multiple models

		for(k in sort(unique(pres$folds))){
			if(any(unlist(lapply(cv.mods[[k]],class))=='cv.glmnet')){
				tmp.proj=projectSDM(dirs=dirs,
														mod.out=cv.mods[[k]],
														stats=stats,
														projEnv=env.new,
														priors=priors,
														name=paste0(stats$species[1],'__fold',k),
														transfer=NULL,
														best.var=best.var,
														fileNameSuffix='_fullDomain',
														gdalCompress=toDo$misc$gdalCompress)													
				tmp.proj$stats$fold=k
				cv.stats[[k]]=tmp.proj$stats
			}	
		}
		list(cv.stats=cv.stats)
	})
	
	status$prog.proj.k=.statusTracker(status$prog.proj.k,cv.mods, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.proj.k==TRUE) {
		outMessage=paste0(species,":  Skipping; projection failed with: ", status$prog.cv.stats)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}
	stats.cv=prog.proj.k$cv.stats
	stats$normalizationConstant=rjson::toJSON(sapply(stats.cv,function(x) ifelse(is.null(x),NA,x$normalizationConstant)))
	

	status$time.project=round((proc.time() - start.time)[3]- status$time.input.prep - status$time.fit)
	if(verbose>0) print(paste0(species,': ','done projecting '))
	if(verbose>0) print(paste0(species,': ',"time: ",status$time.project))

	### 7.d.  Mask to relevant ecoregions
	#-------------------------------------------------------------------
	# basically we won't use the continuous predictions based on buffered points for anything; they'll just hang around in case we need them in the future. 
		   			
	fullDomainFiles=list.files(dirs$sp.pred.path,pattern='fold',full.names=T)	
				#-- by convention, the first layer is the ecoregion one
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

	#-------------------------------------------------------------------
	## 8. Evaluate
	#-------------------------------------------------------------------
	
	if(verbose>0) print(paste0(species,': evaluating ...'))

	### 8.a Evaluate continuous predictions
	#-------------------------------------------------------------------
	# Note- this just evaluates the models masked to occupied ecoregions
	prog.eval.cont=evaluateContFullCV(pres=pres,
																		bg1=bg.eval,
																		biased.bg=biased.bg,
																		dirs=dirs,
																		stats=stats,
																		full.mods=full.mods,
																		maxBGEval=toDo$background$maxBGEval,
																		whichModel='mean',
																		modelNameGrep= toDo$evaluation$evaluationModel,
																		NATo0=FALSE)																		
	status$prog.eval.cont=.statusTracker(status$prog.eval.cont, prog.eval.cont,  toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.eval.cont==TRUE){
		outMessage=paste0(species,":  Skipping; continuous evaluations failed with: ", status$prog.eval.cont)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}
	stats=prog.eval.cont$stats

	### 8.b Evaluate binary predictions
	#-------------------------------------------------------------------
	# consider selecting absences for eval smarter and then using SPS
		# buffer around presences to select background. distances equals x km/max.dist.between pts
	prog.eval.bin=evaluateBinFullCV(pres=pres,
																	bg=bg.eval,
																	dirs=dirs,
																	biased.bg=biased.bg,
																	stats=stats,
																	thresholds=toDo$threshold$thresholds,
																	modelNameGrep= toDo$evaluation$evaluationModel)															
	stats=prog.eval.bin
	status$prog.eval.bin=.statusTracker(status$prog.eval.bin,stats,toDo)
	if(!status$prog.eval.bin==TRUE) {
		outMessage=paste0(species,":  Skipping; binary evaluations failed with: ", status$prog.eval.bin)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}
	
	### 8.c Find best models
	#-------------------------------------------------------------------
	stats=bestModelIndicator(stats=stats)
	status$prog.find.best=.statusTracker(status$prog.find.best,stats, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.find.best==TRUE) {
		outMessage=paste0(species,":  Skipping; finding best models failed with: ", status$prog.avg)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)	
	}

	status$time.eval=round((proc.time() - start.time)[3]- status$time.input.prep - status$time.fit- status$time.project)
	if(verbose>0) print(paste0(species,': done evaluating'))
	if(verbose>2) print(paste0(species,': time: ',status$time.eval))
	
	### 8.d  Summarize k-fold CV weighting by performances
	#-------------------------------------------------------------------

	wtrans=function(y){
		y[y<0.5]=NA	# need to toss models where auc < .5 so low AUC doesnt lead to huge weights
		w=(y-.5)^2; w/sum(w,na.rm=T)}
	w1=suppressWarnings(stats$cv.byfold.test.auc %>% as.character %>% fromJSON %>% unlist %>% as.numeric %>% wtrans %>% replace_na(0)) #weights
	stats$AUCWeighted.test.AUC=suppressWarnings(stats$cv.byfold.test.auc %>% as.character %>% fromJSON %>% unlist %>% as.numeric %>% weighted.mean(w1) %>% round(3))
	w2=suppressWarnings(stats$cv.byfold.test.pAUC.sens10.8%>% as.character %>% fromJSON %>% unlist %>% as.numeric %>% wtrans %>% replace_na(0))
	stats$pAUCWeighted.test.pAUC= suppressWarnings(stats$cv.byfold.test.pAUC.sens10.8%>% as.character %>% fromJSON %>% unlist %>% as.numeric %>% weighted.mean(w2) %>% round(3))
	if(sum(w2>0)<3){
		outMessage=paste0(species,":  Skipping; only 1 or 2 folds has a useful model, which kinda means you don't have a useful model  ", status$prog.avg)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}
	# need to add the domain to the file name; think we just have occEcoNeigh
	if(!toDo$evaluation$evaluationModel=='fullDomain'){
		prog.summarize.cv=summarizeKFoldFast(dirs=dirs,
																				 stats=stats,
																				 writeSDPred=FALSE, # not implemented for weights
																				 algorithm=toDo$misc$algorithm,
																				 gdalCompress=toDo$misc$gdalCompress,
																				 ensembleWeights= toDo$projection$ensembleWeights,
																				 modelNameGrep= toDo$evaluation$evaluationModel,
																				 transfer=NULL)
	}
	#-- could also summarize only this one and mask to other domains but less code
	prog.summarize.cv2=summarizeKFoldFast(dirs=dirs,
																			 stats=stats,
																			 writeSDPred=FALSE, # not implemented for weights
																			 algorithm=toDo$misc$algorithm,
																			 gdalCompress=toDo$misc$gdalCompress,
																			 ensembleWeights= toDo$projection$ensembleWeights,
																			 modelNameGrep='fullDomain',
																			 transfer=NULL)
																			 
	status$prog.avg=.statusTracker(status$prog.avg,stats, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.avg==TRUE) {
		outMessage=paste0(species,":  Skipping; summarizing models fast failed with: ", status$prog.avg)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)	
	}

	status$time.project=round((proc.time() - start.time)[3]- status$time.input.prep - status$time.fit)
	if(verbose>0) print(paste0(species,': done projecting '))
	if(verbose>0) print(paste0(species,': time: ',status$time.project))

	#-------------------------------------------------------------------
	## 9.  Transfer model to new conditions
	#-------------------------------------------------------------------
	### 9.a. Project CV models to full domain
	#-------------------------------------------------------------------
	prog.transfer=NULL
	if(toDo$transfer$otherEnvProjections){
		if(verbose>0) print(paste0(species,': transferring ... '))
		# environment for transfer
		env.new=NULL
		prog.transfer=try({
		 	new.envs=list.files(dirs$otherEnvDir,full.names=T,pattern='tif')
		 	if(all(toDo$transfer$whichFutures=='all')){
				keep=seq_along(basename(new.envs))
			} else {
		 	 	keep=c(mapply(function(x) grep(x,basename(new.envs)),toDo$transfer$whichFutures))
		 	}
			new.envs=new.envs[keep]

			for(ii in 1:length(new.envs)){
				print(paste0(species,': ','  projecting ',tools::file_path_sans_ext(basename(new.envs[ii]))))
				env.new=raster::stack(new.envs[ii])
				prepped=.prepTransferEnv(env.new=env.new,env=env, env.means=env.means, env.sd=env.sd, priors=priors,toDo=toDo,best.var=best.var,pres=pres, dirs=dirs)
				env.new=prepped$env.new
				priors=prepped$priors
			
				for(k in sort(unique(pres$folds))){
					transfer=tools::file_path_sans_ext(basename(new.envs[ii]))
					if(any(unlist(lapply(cv.mods[[k]],class))=='cv.glmnet')){ # to avoid projecting when a model errors
						tmp.proj=projectSDM(dirs=dirs,
																mod.out=cv.mods[[k]],
																stats=stats,
																projEnv=env.new,
																priors=priors,
																name=paste0('fold',k),
																best.var=best.var,
																transfer=transfer,
																algorithm=toDo$misc$algorithm,
																fileNameSuffix='_fullDomain',
																verbose=1)
					}
				} # folds loop
			} # env loop
		}) # try
	}
	
	status$prog.transfer=.statusTracker(status$prog.transfer,stats, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(toDo$misc$rmObjs) rm(env.new)
	if(!status$prog.transfer==TRUE){
		outMessage=paste0(species,":  Skipping; transfer failed with: ", status$prog.transfer)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}

	### 9.b.  Summarize k-fold CV for transfer
	#-------------------------------------------------------------------
	prog.transfer.avg=NULL
	if(toDo$transfer$otherEnvProjections){
		prog.transfer.avg=try({
			for(ii in 1:length(new.envs)){
				transfer=tools::file_path_sans_ext(basename(new.envs[ii]))
				summarizeKFoldFast(dirs=dirs,stats=stats,writeSDPred=FALSE, 
													 transfer=transfer, 
													 algorithm=toDo$misc$algorithm,
													 gdalCompress=toDo$misc$gdalCompress,
													 ensembleWeights=toDo$projection$ensembleWeights,
													 modelNameGrep='fullDomain')
			}
		}) # try
	} # do transfer

	status$prog.transfer.avg=.statusTracker(status$prog.avg,stats, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	if(!status$prog.transfer.avg==TRUE) {
		outMessage=paste0(species,":  Skipping; summarizing transfer models failed with: ", status$prog.transfer.avg)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}

	status$time.transfer=round((proc.time() - start.time)[3]- status$time.input.prep - status$time.fit- status$time.project- status$time.eval)
	if(verbose>0) print(paste0(species,': done transferring'))
	if(verbose>2) print(paste0(species,': time: ',status$time.transfer))

	#-------------------------------------------------------------------
	## 10. Threshold predictions
	#-------------------------------------------------------------------
	### 10.a. under fitting conditions
	#-------------------------------------------------------------------
	ensemble.tr=lapply(toDo$threshold$thresholds, function(x){
			out=evaluateBin(pres=pres,bg=bg.eval,biased.bg=biased.bg,
											modelPath=intersect(list.files(dirs$sp.pred.path, 
											pattern='pAUCWeighted',full.names=T),
												list.files(dirs$sp.pred.path, 
											pattern=toDo$evaluation$evaluationModel,full.names=T)),
											tr='TP',presQuantile=as.numeric(gsub('TP','.',x)), 
											NATo0=TRUE)[-(1:3)]
			names(out)=paste0('pAUCWeighted.',x,'.',names(out))							
			out
	})
	for(x in seq_along(ensemble.tr)){
		a=ensemble.tr[[x]]
		a[-grep('threshold',names(a))]=round(a[-grep('threshold',names(a))],3)
		stats[,colnames(ensemble.tr[[x]])]=a
	}	

	## do same for continuous
	model.r=unrtrans(suppressWarnings(raster(list.files(dirs$sp.pred.path, pattern=paste0(toDo$evaluation$evaluationModel,'__pAUCWeighted'), full.names=T))))
	p.ext=raster::extract(model.r,pres)
	a.ext=raster::extract(model.r,bg.eval)

  out=evaluateCont(p=p.ext,a=a.ext, biased.bg=biased.bg,
									 modelRaster=model.r,partial.auc=list(c(1, .80)), partial.auc.focus=list("sens"))
	names(out)=paste0('pAUCWeighted.',names(out))	
	stats[,names(out)]=round(out,3)

		#-- these are applied to average predictions rather than each fold.
		#-- only do for the 1 weighted ensemble
	# this should applies to the full domain ONLY, because we apply masking below
	if(toDo$threshold$thresholding){
		if(verbose>0) print(paste0(species,': thresholding ... '))
		prog.threshold=threshold(dirs=dirs,
														 stats=stats,
														 pres=pres,
														 env=env,
														 saveShapeFile=toSave$shapeFile,
														 makeThresholdMaps=toDo$threshold$thresholds,
														 makeSparse=toDo$threshold$makeSparseMatrix,
														 modelNameGrep='fullDomain',
														 gdalCompress=toDo$misc$gdalCompress)

		status$prog.threshold= .statusTracker(status$prog.threshold,stats, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
		if(!status$prog.threshold==TRUE) {
			outMessage=paste0(species,":  Skipping; thresholding failed with: ", status$prog.threshold)
			message(outMessage)
			exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
			return(outMessage)
		}
	}
	### 10.b. threshold  transferred models

	#-------------------------------------------------------------------
	if(!all(is.null(stats.cv))){
		if(toDo$transfer$otherEnvProjections & toDo$threshold$thresholding){
			if(verbose>0) print(paste0(species,': thresholding transferred models  ... '))
			prog.threshold.new=try({
				new.envs=list.files(dirs$otherEnvDir,full.names=T,pattern='tif')
				if(all(toDo$transfer$whichFutures=='all')){
					keep=seq_along(basename(new.envs))
				} else {
					keep=c(mapply(function(x) grep(x,basename(new.envs)),toDo$transfer$whichFutures))
				}
				new.envs=new.envs[keep]
				for(ii in 1:length(new.envs)){
				  envName=tools::file_path_sans_ext(basename(new.envs[ii]))
				  if(verbose>0) print(paste0(species,': thresholding transferred models  - ',envName))
					tmp=threshold(dirs=dirs, 
												stats=stats,
												pres=pres,
												env=env, 
												saveShapeFile=toSave$shapeFileTransfer,
												whichModel='mean',
												transfer=envName,
												modelNameGrep='fullDomain',
												makeThresholdMaps=toDo$threshold$makeThresholdMaps,
												makeSparse=toDo$threshold$makeSparseMatrix,
												gdalCompress=toDo$misc$gdalCompress)
				}
			})
		}

		status$prog.threshold.new=.statusTracker(status$prog.threshold.new, stats,toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
		if(!status$prog.threshold.new==TRUE) {
			outMessage=paste0(species,":  Skipping; thresholding transfered models failed with: ", status$prog.thershold.new)
			message(outMessage)
			exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
			return(outMessage)
		}
	}

	status$time.threshold=round((proc.time() - start.time)[3]- status$time.input.prep - status$time.fit- status$time.project- status$time.eval- status$time.transfer)
	if(verbose>0) print(paste0(species,': done thresholding'))
	if(verbose>2) print(paste0(species,': time: ',status$time.transfer))


	### 10.c. Trinary Maps under fitting conditions
	#-------------------------------------------------------------------
	# pass toSave to function to save trinary (it should always be made to get upper and lower bounds on range size and pAUC)
	if(toDo$threshold$trinaryMap){
	 prog.trinary=trinaryMapSDM(dirs,
															stats,
															env,
															pres,
															bg1=bg.eval,
															openFig=FALSE,
															shapesToPlot=toDo$plotting$shapesForPlotting,
															doMapPlot=TRUE,
															doROCPlot=FALSE,
															modelNameGrep1=toDo$evaluation$evaluationModel,
															modelNameGrep2='pAUCWeighted',
															verbose=5)

	 if(!class(prog.trinary)=='try-error') {stats=prog.trinary 
	 } else {
		 outMessage=paste0(species,":  NOT Skipping the species, because its not super important, but trinary maps failed with: ", status$trinary)
		 message(outMessage)
	 }
	 status$prog.trinary=.statusTracker(status$prog.trinary,outObj=	prog.trinary, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
	 #if(!status$prog.trinary==TRUE) return(status)
	}

	#-------------------------------------------------------------------
	## 11. Masking predictions
	#-------------------------------------------------------------------
		#-- the continuous preductions on the fitting data were already masked above to allow evaluation on occupied plus neighboring. I don't think we really need
	
	#-------------------------------------------------------------------

	### 11.a. Mask binary maps to occupied (only) and Occ+neighbor ecoregions 
	#-------------------------------------------------------------------
	if(toDo$trimDomainMask){
		myMasks=raster::stack(list.files(dirs$sp.mask.path,full.names=T))[[1]]      	
		thisMask=raster::reclassify(myMasks,c(.99,2.01,1)) 
		maskBinaryMaps(thisMask,toDo,dirs,newMapName='occEcoNeigh',species)
		if(verbose>0) print(paste0(species,': done masking to occEcoNeigh'))

	### 11.b. Mask binary maps to Occupied ecoregions 
	#-------------------------------------------------------------------
		myMasks=raster::stack(list.files(dirs$sp.mask.path,full.names=T))[[1]]      	
		thisMask=raster::reclassify(myMasks,c(1.01,Inf,NA)) # use occupied only
		maskBinaryMaps(thisMask,toDo,dirs,newMapName='occEco',species)
		if(verbose>0) print(paste0(species,': done masking to occEco'))
	}


	#-------------------------------------------------------------------
	## 12.  Plotting
	#-------------------------------------------------------------------

	### 12.a  weighted ensemble maps 
	#-------------------------------------------------------------------
	prog.plot=compareSdmStatsPlot(stats=stats,
	                              dirs=dirs,
	                              range.poly=NULL,
	                              pres=pres,
	                              env=env,
	                              speciesInfo=FALSE,
	                              legendLoc='panel',
	                              logscaleZ=TRUE,
	                              modelNameGrep=toDo$evaluation$evaluationModel,
	                              toDo=toDo)
	
	status$prog.plot=.statusTracker(status$prog.plot,prog.plot, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)

	if(!status$prog.plot==TRUE) {
		outMessage=paste0(species,":  Skipping; plotting maps failed with: ", status$prog.plot)
		message(outMessage)
		exitSDMWorkflow(dirs=dirs,stats=stats, toSave=toSave,toDo=toDo,status=status)
		return(outMessage)
	}
	### 12.b. Make response curves
	#-------------------------------------------------------------------
	if(toDo$plotting$responseCurves){
		if(verbose>0) print(paste0(species,': plotting response curves ...'))
		prog.response.curves=try(responseCurves(dirs, 
																						env, 
																						stats, 
																						bg=bg,			
																						toDo=toDo,
																						best.var=best.var,
																						priors=NA,
																						pres=pres,
																						envMeans=env.means,
																						envSDs=env.sd,
																						extrapLimits=extrapLimits,
																						openFigs=FALSE))
		status$prog.response.response.curves=.statusTracker( status$prog.prog.response.curves,prog.response.curves, toSave=toSave,toDo=toDo, dirs=dirs,stats=stats,status=status)
		if(class(prog.response.curves)=='try-error') {
			outMessage=paste0(species,":  NOT Skipping the species because its not a super important issue, but plotting response curves failed with: ", status$prog.response.curves)
			message(outMessage)
		}																
	}
	
	if(verbose>0) print(paste0(species,': done plotting'))

	#-------------------------------------------------------------------
	## 13.  Clean up / Record model metadata
	#-------------------------------------------------------------------
	# TODO set this up for rmm objects	
	if(toDo$metadata$bienMetadata){	md=makeBIENMetadata(stats=stats,dirs=dirs,toDo=toDo)
		sp=statsPublic(stats, statPublicFile=paste0(dirs$algorithmDirs$PPM$statPublicDir, '/',species,'__StatsPub.csv'))
	}
	#rmm=rmmsBIEN(dirs=dirs,stats=stats,toDo=toDo)
 	if(verbose>2) print(paste0(species,': writing statistics and cleaning'))
 	status$time.all=stats$runTime=round((proc.time() - start.time)[3])
 	status$finished=TRUE 

 	ex=exitSDMWorkflow(dirs=dirs,stats=stats,status=status,toSave=toSave, toDo=toDo)
	
 	# clean up tmp files
	unlink(rasterOptions(overwrite=T)$tmpdir,recursive=T)
		
	if(verbose>0) print(paste0(species,': finished model in ', status$time.all))
	if(verbose>2) print(paste0(species,': -------------------------'))

	outMessage=paste0(species,': looks good')
	print(outMessage)
	return(outMessage)
  })
  return(out)
}

