
# Rscript "Final Iceberg Code/7_New Encounter GAMs.R" ####

library(raster); library(tidyverse); library(mgcv); library(cowplot)

NewIntersects <- readRDS("~/Albersnet/Iceberg Output Files/NewIntersects.rds")

PredReps <- c("Currents", paste0("Futures", 1:4))
PipelineReps <- LETTERS[1:4]

AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
AlberColours[length(AlberColours)+1:2] <- RColorBrewer::brewer.pal(11, AlberPalettes[[4]])[c(2,10)]

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

x = 2

# Resps <- c("OverlapSum", "SharingMean", "SharingSum")
Resps <- c("SomeBat", "NoBat")

NewIntersects %>% dplyr::select(-X,-Y, -contains("NoBat")) -
  NewIntersects %>% dplyr::select(contains("NoBat")) ->
  
  NewIntersects[,c(paste0("SomeBat", NewIntersects %>% dplyr::select(-X,-Y, -contains("NoBat")) %>% names))]

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

Elevation <- raster("Iceberg Input Files/altitude.tif")
Temperature <- raster("Iceberg Input Files/Temperature.tif")
Precipitation <- raster("Iceberg Input Files/Precipitation.tif")

NEPredictModels <- TestDFList <- FitList <- list()

CurrentsGridDF <- readRDS("~/Albersnet/Iceberg Output Files/CurrentsGridDF.rds") 

BlankDF <- CurrentsGridDF[,c("X", "Y", "Climate")] %>% rename(Richness = Climate)

BlankDF[,c("Elevation", "Temperature", "Precipitation")] <- cbind(
  
  Elevation %>% values,
  Temperature %>% values,  
  Precipitation %>% values
  
)[-Sea,]

NewIntersects[,c("Elevation", "Temperature", "Precipitation")] <- 
  BlankDF[,c("Elevation", "Temperature", "Precipitation")]

GridDF <- readRDS("~/Albersnet/Iceberg Output Files/CurrentsGridDF.rds") 

HabitatTypes <- c("Forest", "Nonforest", "Urban", "Cropland", "Rangeland")

for(Pipeline in PipelineReps[c(1,4)]){
  
  print(Pipeline)
  
  FitList[[Pipeline]] <- 
    NEPredictModels[[Pipeline]] <- 
    TestDFList[[Pipeline]] <- 
    list()
  
  for(x in (2:length(PredReps))[c(1,4)]){
    
    print(PredReps[x])
    
    NEPredictModels[[Pipeline]][[PredReps[x]]] <- FitList[[Pipeline]][[PredReps[x]]] <-  list()
    
    for(Resp in Resps){
      
      TestDF <- NewIntersects %>%  
        mutate(Overlaps = NewIntersects[,paste0(Resp,"Overlap.",PredReps[x],Pipeline)])
      
      if(Pipeline%in%c("A","C")){
        
        TestDF$Richness <- GridDF$ClimateLandUse
        
      }else{
        
        TestDF$Richness <- GridDF$Climate
        
      }
      
      LandUse <- brick(paste0("Iceberg Input Files/","ssp",((1:5)[-3])[x-1], "landusemax", ".tif"))
      
      TestDF$HabitatType <- HabitatTypes[values(LandUse)[-Sea]] %>% as.factor
      
      TestDF %>% na.omit() -> 
        TestDF ->
        TestDFList[[Pipeline]][[PredReps[x]]]
      
      Covar <- c(paste0("s(",c("Richness", 
                               "Elevation"),")"), "HabitatType")
      
      Formula <- Covar %>% paste(collapse = " + ") %>% paste0("Overlaps ~ ", .) %>% as.formula()
      
      NEPredictModels[[Pipeline]][[PredReps[x]]][[Resp]] <- bam(Formula,
                                                                data = TestDF, 
                                                                family = nb()
      )
      
      # Predicting lines ####
      
      Model <- NEPredictModels[[Pipeline]][[PredReps[x]]][[Resp]]
      
      # Effects ####
      
      RichnessRange <- seq(from = min(TestDFList[[Pipeline]][[PredReps[x]]]$Richness),
                           to = max(TestDFList[[Pipeline]][[PredReps[x]]]$Richness),
                           length = 51) %>% 
        c(mean(TestDFList[[Pipeline]][[PredReps[x]]]$Richness))
      
      ElevationRange <- seq(from = min(TestDFList[[Pipeline]][[PredReps[x]]]$Elevation),
                            to = max(TestDFList[[Pipeline]][[PredReps[x]]]$Elevation),
                            length = 51)  %>% 
        c(mean(TestDFList[[Pipeline]][[PredReps[x]]]$Elevation))
      
      FitList[[Pipeline]][[PredReps[x]]][[Resp]] <- expand.grid(
        Rep = paste0(PredReps[x],Pipeline),
        Richness = RichnessRange,
        Elevation = ElevationRange,
        HabitatType = levels(TestDF$HabitatType)
      )
      
      FitPredictions  <- predict.gam(Model, 
                                     newdata = FitList[[Pipeline]][[PredReps[x]]][[Resp]], 
                                     se.fit = T)
      
      FitList[[Pipeline]][[PredReps[x]]][[Resp]][,c("Fit","Lower", "Upper")] <- (with(FitPredictions, cbind(fit, fit - se.fit, fit + se.fit)))
      FitList[[Pipeline]][[PredReps[x]]][[Resp]][,c("Fit","Lower", "Upper")] <- exp(with(FitPredictions, cbind(fit, fit - se.fit, fit + se.fit)))
      
    }
    
  }
}

FullDataDF <- TestDFList %>% lapply(function(a){
  
  a %>% bind_rows(.id = "PredRep")
  
}) %>% bind_rows(.id = "Pipeline") %>% mutate(Rep = paste0(PredRep, Pipeline))

for(Pipeline in PipelineReps[c(1,4)]){
  print(Pipeline)
  for(x in (2:length(PredReps))[c(1,4)]){
    print(PredReps[x])
    for(Resp in Resps){
      Model <- NEPredictModels[[Pipeline]][[PredReps[x]]][[Resp]]
      
      TestDF <- FullDataDF %>% filter(Pipeline == Pipeline, PredRep == PredReps[x])
      
      RichnessRange <- seq(from = min(TestDF$Richness),
                           to = max(TestDF$Richness),
                           length = 51) %>%
        c(mean(TestDF$Richness))
      ElevationRange <- seq(from = min(TestDF$Elevation),
                            to = max(TestDF$Elevation),
                            length = 51)  %>%
        c(mean(TestDF$Elevation))
      FitList[[Pipeline]][[PredReps[x]]][[Resp]] <- expand.grid(
        Rep = paste0(PredReps[x],Pipeline),
        Richness = RichnessRange,
        Elevation = ElevationRange,
        HabitatType = levels(TestDF$HabitatType)
      )
      FitPredictions  <- predict.gam(Model,
                                     newdata = FitList[[Pipeline]][[PredReps[x]]][[Resp]],
                                     se.fit = T)
      FitList[[Pipeline]][[PredReps[x]]][[Resp]][,c("Fit","Lower", "Upper")] <- (with(FitPredictions, cbind(fit, fit - se.fit, fit + se.fit)))
      FitList[[Pipeline]][[PredReps[x]]][[Resp]][,c("Fit","Lower", "Upper")] <- exp(with(FitPredictions, cbind(fit, fit - se.fit, fit + se.fit)))
    }
  }
}


save(FitList, FullDataDF, NEPredictModels, file = "Iceberg Output Files/NEGams.Rdata")

FitList %>% unlist(recursive = F) %>% bind_rows() %>%
  filter(Elevation == last(unique(Elevation))) %>%
  ggplot(aes(Richness, Fit, group = HabitatType)) + 
  geom_point(data = FullDataDF, inherit.aes = F, aes(x = Richness, y = Overlaps), alpha = 0.1, colour = AlberColours[[2]]) + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, colour = NA) +
  geom_line() + #  scale_y_log10() +
  labs(y = "Overlaps", x = "Richness") +
  facet_wrap(~Rep, ncol = 4, scales = "free") +
  ggsave("NEGams_Richness.jpeg", units = "mm", dpi = 300, width = 400, height = 400)

FitList %>% unlist(recursive = F)  %>% unlist(recursive = F) %>% bind_rows() %>% 
  group_by(Rep) %>% summarise(Rich = last(Richness)) %>% 
  pull(Rich) %>% unique() -> Richnesses

FitList %>% unlist(recursive = F) %>% bind_rows() %>%
  filter(Richness %in% Richnesses) %>%
  ggplot(aes(Elevation, Fit, group = HabitatType)) + 
  geom_point(data = FullDataDF, inherit.aes = F, aes(x = Elevation, y = Overlaps), alpha = 0.1, colour = AlberColours[[2]]) + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, colour = NA) +
  geom_line() + #  scale_y_log10() +
  labs(y = "Overlaps", x = "Elevation") +
  facet_wrap(~Rep, ncol = 4, scales = "free") +
  ggsave("NEGams_Elevation.jpeg", units = "mm", dpi = 300, width = 400, height = 400)

plot_grid(FitList[[Pipeline]][[PredReps[x]]] %>% 
            filter(Elevation == last(unique(Elevation))) %>%
            ggplot(aes(Richness, Fit)) + 
            geom_point(data = TestDFList[[Pipeline]][[PredReps[x]]], inherit.aes = F, aes(x = Richness, y = OverlapSum), alpha = 0.1, colour = AlberColours[[2]]) + 
            geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, colour = NA) +
            geom_line() + #  scale_y_log10() +
            labs(y = "Overlaps", x = "Richness"),
          
          FitList[[Pipeline]][[PredReps[x]]] %>% 
            filter(Richness == last(unique(Richness))) %>%
            ggplot(aes(Elevation, Fit))+
            geom_point(data = TestDFList[[Pipeline]][[PredReps[x]]], inherit.aes = F, aes(x = Elevation, y = OverlapSum), alpha = 0.1, colour = AlberColours[[2]]) + 
            geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, colour = NA) +
            geom_line() + # scale_y_log10() +
            labs(y = "Overlaps", x = "Elevation")
)

plot_grid(FitList[[Pipeline]][[PredReps[x]]] %>% 
            filter(Elevation == last(unique(Elevation)),
                   Temperature == last(unique(Temperature)),
                   Precipitation == last(unique(Precipitation))) %>%
            ggplot(aes(Richness, Fit)) + 
            geom_point(data = TestDFList[[Pipeline]][[PredReps[x]]], inherit.aes = F, aes(x = Richness, y = OverlapSum), alpha = 0.1, colour = AlberColours[[2]]) + 
            geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, colour = NA) +
            geom_line() +
            labs(y = "Overlaps", x = "Richness"),
          
          FitList[[Pipeline]][[PredReps[x]]] %>% 
            filter(Richness == last(unique(Richness)),
                   Temperature == last(unique(Temperature)),
                   Precipitation == last(unique(Precipitation))) %>%
            ggplot(aes(Elevation, Fit))+
            geom_point(data = TestDFList[[Pipeline]][[PredReps[x]]], inherit.aes = F, aes(x = Elevation, y = OverlapSum), alpha = 0.1, colour = AlberColours[[2]]) + 
            geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, colour = NA) +
            geom_line() +
            labs(y = "Overlaps", x = "Elevation"), 
          
          FitList[[Pipeline]][[PredReps[x]]] %>% 
            filter(Richness == last(unique(Richness)),
                   Precipitation == last(unique(Precipitation)),
                   Elevation == last(unique(Elevation))) %>%
            ggplot(aes(Temperature, Fit)) + 
            geom_point(data = TestDFList[[Pipeline]][[PredReps[x]]], inherit.aes = F, aes(x = Temperature, y = OverlapSum), alpha = 0.1, colour = AlberColours[[2]]) + 
            geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, colour = NA) +
            geom_line() +
            labs(y = "Overlaps", x = "Temperature"),
          
          FitList[[Pipeline]][[PredReps[x]]] %>% 
            filter(Richness == last(unique(Richness)),
                   Temperature == last(unique(Temperature)),
                   Elevation == last(unique(Elevation))) %>%
            ggplot(aes(Precipitation, Fit))+
            geom_point(data = TestDFList[[Pipeline]][[PredReps[x]]], inherit.aes = F, aes(x = Precipitation, y = OverlapSum), alpha = 0.1, colour = AlberColours[[2]]) + 
            geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, colour = NA) +
            geom_line() +
            labs(y = "Overlaps", x = "Precipitation"), 
          
          ncol = 2
          
)

plot_grid(FitList[[Pipeline]][[PredReps[x]]] %>% 
            filter(Elevation == last(unique(Elevation)),
                   Temperature == last(unique(Temperature)),
                   Precipitation == last(unique(Precipitation))) %>%
            ggplot(aes(Richness, log(Fit+1))) + 
            geom_point(data = TestDFList[[Pipeline]][[PredReps[x]]], inherit.aes = F, aes(x = Richness, y = log(OverlapSum + 1)), alpha = 0.1, colour = AlberColours[[2]]) + 
            geom_ribbon(aes(ymin = log(Lower+1), ymax = log(Upper+1)), alpha = 0.2, colour = NA) +
            geom_line() +
            labs(y = "Overlaps", x = "Richness"),
          
          FitList[[Pipeline]][[PredReps[x]]] %>% 
            filter(Richness == last(unique(Richness)),
                   Temperature == last(unique(Temperature)),
                   Precipitation == last(unique(Precipitation))) %>%
            ggplot(aes(Elevation, log(Fit+1))) +
            geom_point(data = TestDFList[[Pipeline]][[PredReps[x]]], inherit.aes = F, aes(x = Elevation, y = log(OverlapSum + 1)), alpha = 0.1, colour = AlberColours[[2]]) + 
            geom_ribbon(aes(ymin = log(Lower+1), ymax = log(Upper+1)), alpha = 0.2, colour = NA) +
            geom_line() +
            labs(y = "Overlaps", x = "Elevation"), 
          
          FitList[[Pipeline]][[PredReps[x]]] %>% 
            filter(Richness == last(unique(Richness)),
                   Precipitation == last(unique(Precipitation)),
                   Elevation == last(unique(Elevation))) %>%
            ggplot(aes(Temperature, log(Fit+1))) + 
            geom_point(data = TestDFList[[Pipeline]][[PredReps[x]]], inherit.aes = F, aes(x = Temperature, y = log(OverlapSum + 1)), alpha = 0.1, colour = AlberColours[[2]]) + 
            geom_ribbon(aes(ymin = log(Lower+1), ymax = log(Upper+1)), alpha = 0.2, colour = NA) +
            geom_line() +
            labs(y = "Overlaps", x = "Temperature"),
          
          FitList[[Pipeline]][[PredReps[x]]] %>% 
            filter(Richness == last(unique(Richness)),
                   Temperature == last(unique(Temperature)),
                   Elevation == last(unique(Elevation))) %>%
            ggplot(aes(Precipitation, log(Fit+1)))+
            geom_point(data = TestDFList[[Pipeline]][[PredReps[x]]], inherit.aes = F, aes(x = Precipitation, y = log(OverlapSum + 1)), alpha = 0.1, colour = AlberColours[[2]]) + 
            geom_ribbon(aes(ymin = log(Lower+1), ymax = log(Upper+1)), alpha = 0.2, colour = NA) +
            geom_line() +
            labs(y = "Overlaps", x = "Precipitation"), 
          
          ncol = 2
          
)

PredDF <- TestDFList[[Pipeline]][[PredReps[x]]]
TestDFList[[Pipeline]][[PredReps[x]]]$Fit <- predict.gam(Model)

TestDFList[[Pipeline]][[PredReps[x]]] %>% ggplot(aes(X, Y, fill = exp(Fit))) + geom_tile() + scale_fill_continuous_sequential(palette = "Terrain")

# New version of GAM output figure ####

FitList %>% unlist(recursive = F) %>% unlist(recursive = F) %>% bind_rows() %>% 
  group_by(Rep) %>% summarise(Rich = last(Richness)) %>% 
  pull(Rich) %>% unique() -> Richnesses

lapply(Resps, function(a){
  
  FitList %>% unlist(recursive = F) %>%  map(a) %>% bind_rows(.id = "Pipeline") %>%
    filter(HabitatType == "Cropland", Richness %in% Richnesses)  %>%
    ggplot(aes(Elevation, Fit, fill = Pipeline, colour = Pipeline)) + 
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, colour = NA) + geom_line() +
    ggtitle(a) + scale_fill_discrete_sequential(palette = AlberPalettes[[1]]) + 
    labs(y = "Predicted Encounters") + theme(legend.position = c(0.5,0.8)) +
    scale_colour_discrete_sequential(palette = AlberPalettes[[1]])
  
}) %>% plot_grid(plotlist = ., ncol = 2)

lapply(Resps, function(a){
  
  FitList %>% unlist(recursive = F) %>%  map(a) %>% bind_rows(.id = "Pipeline") %>%
    filter(HabitatType == "Cropland", Elevation == last(unique(Elevation)))  %>%
    ggplot(aes(Richness, Fit, fill = Pipeline, colour = Pipeline)) + 
    geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, colour = NA) + geom_line() +
    ggtitle(a) + labs(y = "Predicted Encounters") + theme(legend.position = c(0.2,0.8)) + 
    scale_fill_discrete_sequential(palette = AlberPalettes[[2]]) + 
    scale_colour_discrete_sequential(palette = AlberPalettes[[2]])
  
}) %>% plot_grid(plotlist = ., ncol = 2)

# Deviance contributions #####

# Validating the model and getting deviance contributions 

load("Iceberg Output Files/NEGams.Rdata")

Iterations = 1

Resps <- paste0("A.", PredReps[c(2,5)])

RandomPredictionList <- DevianceList <- list()

RealPredictions <- InterceptPredictions <- list()

x <- 2
r <- "NoBat"

PredictCovar <- c("Richness", "Elevation", "Temperature", "Precipitation")[1:2]
CategoryCovar <- "HabitatType"

NEPredictModels <- NEPredictModels %>% unlist(recursive = F)

for(y in Resps[1]){
  
  print(y)
  
  if(str_detect(y, "1")){
    
    x <- 2
    r2 <- "SomeBatOverlap.Futures1A"
    r2 <- "NoBatOverlap.Futures1A"
    
  }else{
    
    x <- 5
    r2 <- "SomeBatOverlap.Futures4A"
    r2 <- "NoBatOverlap.Futures4A"
    
  }
  
  RandomPredictionList[[y]] <- DevianceList[[y]] <- list()
  
  RealOutcomes <- TestDFList[[Pipeline]][[PredReps[x]]][[r2]]
  
  RealPredictions[[y]] <- predict.bam(NEPredictModels[[y]][[r]], 
                                      newdata = TestDFList[[Pipeline]][[PredReps[x]]]) %>% exp
  
  InterceptPredictions[[y]] <- rep(mean(RealPredictions[[y]]), nrow(TestDFList[[Pipeline]][[PredReps[x]]]))
  
  for(z in PredictCovar){
    
    print(z)
    
    for(i in 1:Iterations){
      
      print(i)
      
      PredDF <- TestDFList[[Pipeline]][[PredReps[x]]]
      
      PredDF[,z] <- PredDF %>% slice(sample(1:n())) %>% pull(z)
      
      Predictions <- predict.bam(NEPredictModels[[y]][[r]], newdata = PredDF)
      
      RandomPredictionList[[z]][[i]] <- Predictions %>% exp
      
      ModelLikelihood = (RealOutcomes*log(RandomPredictionList[[z]][[i]]) - RandomPredictionList[[z]][[i]]) %>% sum
      
      Deviance = -2*ModelLikelihood
      
      DevianceList[[y]][[z]][[i]] <- Deviance
    }
  }
  
  for(z in CategoryCovar){
    
    print(z)
    
    for(i in 1:Iterations){
      
      print(i)
      
      PredDF <- TestDFList[[Pipeline]][[PredReps[x]]]
      
      PredDF[,z] <- PredDF %>% slice(sample(1:n())) %>% pull(z)
      
      Predictions <- predict.bam(NEPredictModels[[y]][["SomeBat"]], 
                                 newdata = PredDF)
      
      RandomPredictionList[[z]][[i]] <- Predictions %>% exp
      
      ModelLikelihood = (RealOutcomes*log(RandomPredictionList[[z]][[i]]) - RandomPredictionList[[z]][[i]]) %>% sum
      
      Deviance = -2*ModelLikelihood
      
      DevianceList[[y]][[z]][[i]] <- Deviance
    }
  }
}

RealDeviance <- InterceptDeviance <- list()

for(y in Resps[1]){
  
  print(y)
  
  if(str_detect(y, "1")){
    
    r <- "SomeBatOverlap.Futures1A"
    r <- "NoBatOverlap.Futures1A"
    
  }else{
    
    r <- "SomeBatOverlap.Futures1D"
    r <- "NoBatOverlap.Futures1D"
    
  }
  
  RealModelLikelihood = (TestDFList[[Pipeline]][[PredReps[x]]][[r]]*log(RealPredictions[[y]]) - RealPredictions[[y]]) %>% sum
  RealDeviance[[y]] = -2*RealModelLikelihood
  
  InterceptModelLikelihood = (TestDFList[[Pipeline]][[PredReps[x]]][[r]]*log(InterceptPredictions[[y]]) - InterceptPredictions[[y]]) %>% sum
  InterceptDeviance[[y]] = -2*InterceptModelLikelihood
  
}

DevianceList %>% lapply(function(a) sapply(a, mean)) %>% c(Real = RealDeviance, Intercept = InterceptDeviance)

lapply(Resps[1], function(a){
  DevianceDF <- data.frame(
    Var = names(DevianceList[[a]]),
    Model_Deviance = (((sapply(DevianceList[[a]], mean) - RealDeviance[[a]]) %>% prop.table())) %>% round(3)
  ) %>%
    mutate(Total_Deviance = Model_Deviance*(RealDeviance[[a]]/InterceptDeviance[[a]]))
  
}) -> DevianceDFList

names(DevianceDFList) <- Resps
