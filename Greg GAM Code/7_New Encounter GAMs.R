
# Rscript "Final Iceberg Code/07_New Encounter GAMs" ####

library(raster); library(tidyverse); library(mgcv)

NewIntersects <- readRDS("~/Albersnet/Iceberg Output Files/NewIntersects.rds")

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

x = 2

Resps <- c("OverlapSum", "SharingMean", "SharingSum")

Elevation <- raster("Iceberg Input Files/altitude.tif")
Temperature <- raster("Iceberg Input Files/Temperature.tif")
Precipitation <- raster("Iceberg Input Files/Precipitation.tif")

NEPredictModels <- TestDFList <- FitList <- list()

BlankDF <- CurrentsGridDF[,c("X", "Y", "Climate")] %>% rename(Richness = Climate)

BlankDF[,c("Elevation", "Temperature", "Precipitation")] <- cbind(
  
  Elevation %>% values,
  Temperature %>% values,  
  Precipitation %>% values
  
)[-Sea,]

NewIntersects[,c("Elevation", "Temperature", "Precipitation")] <- 
  BlankDF[,c("Elevation", "Temperature", "Precipitation")]

GridDF <- readRDS("~/Albersnet/Iceberg Output Files/CurrentsGridDF.rds") 

for(Pipeline in PipelineReps){
  
  print(Pipeline)
  
  FitList[[Pipeline]] <- 
    NEPredictModels[[Pipeline]] <- 
    TestDFList[[Pipeline]] <- 
    list()
  
  for(x in 2:length(PredReps)){
    
    print(PredReps[x])
    
    TestDF <- NewIntersects %>%  
      mutate(Overlaps = NewIntersects[,paste0("Overlap.",PredReps[x],Pipeline)])
    
    if(Pipeline%in%c("A","C")){
      
      TestDF$Richness <- GridDF$ClimateLandUse
      
    }else{
      
      TestDF$Richness <- GridDF$Climate
      
    }
    
    LandUse <- brick(paste0("Iceberg Input Files/","LandUse", PredReps[x], ".grd"))
    
    TestDF[,names(LandUse)] <- lapply(names(LandUse), function(a){
      
      (LandUse[[a]] %>% values)[-Sea]
      
    }) %>% bind_cols
    
    TestDF %>% na.omit() -> 
      TestDF ->
      TestDFList[[Pipeline]][[PredReps[x]]]
    
    Covar <- c(paste0("s(",c("Richness", 
                             "Elevation"),")"), names(LandUse))
    
    Formula <- Covar %>% paste(collapse = " + ") %>% paste0("Overlaps ~ ", .) %>% as.formula()
    
    NEPredictModels[[Pipeline]][[PredReps[x]]] <- list()
    
    NEPredictModels[[Pipeline]][[PredReps[x]]]$Overlaps <- bam(Formula,
                                                               data = TestDF, 
                                                               family = nb()
    )
    
    # Predicting lines ####
    
    Model <- NEPredictModels[[Pipeline]][[PredReps[x]]]$Overlaps
    
    # Effects ####
    
    RichnessRange <- seq(from = min(TestDFList[[Pipeline]][[PredReps[x]]]$Richness),
                         to = max(TestDFList[[Pipeline]][[PredReps[x]]]$Richness),
                         length = 21) %>% 
      c(mean(TestDFList[[Pipeline]][[PredReps[x]]]$Richness))
    
    ElevationRange <- seq(from = min(TestDFList[[Pipeline]][[PredReps[x]]]$Elevation),
                          to = max(TestDFList[[Pipeline]][[PredReps[x]]]$Elevation),
                          length = 21)  %>% 
      c(mean(TestDFList[[Pipeline]][[PredReps[x]]]$Elevation))
    
    FitList[[Pipeline]][[PredReps[x]]] <- expand.grid(
      Rep = paste0(PredReps[x],Pipeline),
      Richness = RichnessRange,
      Elevation = ElevationRange
    )
    
    FitList[[Pipeline]][[PredReps[x]]][,names(LandUse)] <- 0
    
    FitPredictions  <- predict.gam(Model, 
                                   newdata = FitList[[Pipeline]][[PredReps[x]]], 
                                   se.fit = T)
    
    FitList[[Pipeline]][[PredReps[x]]][,c("Fit","Lower", "Upper")] <- (with(FitPredictions, cbind(fit, fit - se.fit, fit + se.fit)))
    FitList[[Pipeline]][[PredReps[x]]][,c("Fit","Lower", "Upper")] <- exp(with(FitPredictions, cbind(fit, fit - se.fit, fit + se.fit)))
    
  }
}

FullDataDF <- TestDFList %>% lapply(function(a){
  
  a %>% bind_rows(.id = "PredRep")
  
}) %>% bind_rows(.id = "Pipeline") %>% mutate(Rep = paste0(PredRep, Pipeline))

save(FitList, FullDataDF, NEPredictModels, file = "Iceberg Output Files/NEGams.Rdata")

FitList %>% unlist(recursive = F) %>% bind_rows() %>%
  filter(Elevation == last(unique(Elevation))) %>%
  ggplot(aes(Richness, Fit)) + 
  geom_point(data = FullDataDF, inherit.aes = F, aes(x = Richness, y = Overlaps), alpha = 0.1, colour = AlberColours[[2]]) + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, colour = NA) +
  geom_line() + #  scale_y_log10() +
  labs(y = "Overlaps", x = "Richness") +
  facet_wrap(~Rep, ncol = 4, scales = "free") +
  ggsave("NEGams_Richness.jpeg", units = "mm", dpi = 300, width = 400, height = 400)

FitList %>% unlist(recursive = F) %>% bind_rows() %>% 
  group_by(Rep) %>% summarise(Rich = last(Richness)) %>% 
  pull(Rich) %>% unique() -> Richnesses

FitList %>% unlist(recursive = F) %>% bind_rows() %>%
  filter(Richness %in% Richnesses) %>%
  ggplot(aes(Elevation, Fit)) + 
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

# Deviance contributions #####

# Validating the model and getting deviance contributions 

Iterations = 10

Resps <- PredReps[2:5]

RandomPredictionList <- DevianceList <- list()

RealPredictions <- InterceptPredictions <- list()

y = "Futures1"
r = "OverlapSum"

PredictCovar <- c("Richness", "Elevation", "Temperature", "Precipitation")
CategoryCovar <- names(LandUse)

for(y in Resps){
  
  print(y)
  
  RandomPredictionList[[y]] <- DevianceList[[y]] <- list()
  
  RealOutcomes <- TestDFList[[Pipeline]][[PredReps[x]]][[r]]
  
  RealPredictions[[y]] <- predict.bam(NEPredictModels[[y]][[r]], 
                                      newdata = TestDFList[[Pipeline]][[PredReps[x]]]) %>% exp
  
  InterceptPredictions[[y]] <- rep(mean(RealPredictions[[y]]), nrow(TestDFList[[Pipeline]][[PredReps[x]]]))
  
  for(x in PredictCovar){
    
    print(x)
    
    for(i in 1:Iterations){
      
      print(i)
      
      PredDF <- TestDFList[[Pipeline]][[PredReps[x]]]
      
      PredDF[,x] <- PredDF %>% slice(sample(1:n())) %>% pull(x)
      
      Predictions <- predict.bam(NEPredictModels[[y]][[r]], newdata = PredDF)
      
      RandomPredictionList[[x]][[i]] <- Predictions %>% exp
      
      ModelLikelihood = (RealOutcomes*log(RandomPredictionList[[x]][[i]]) - RandomPredictionList[[x]][[i]]) %>% sum
      
      Deviance = -2*ModelLikelihood
      
      DevianceList[[y]][[x]][[i]] <- Deviance
    }
  }
  
  for(x in CategoryCovar){
    
    print(x)
    
    for(i in 1:Iterations){
      
      print(i)
      
      PredDF <- TestDFList[[Pipeline]][[PredReps[x]]]
      
      PredDF[,x] <- PredDF %>% slice(sample(1:n())) %>% pull(x)
      
      Predictions <- predict.bam(NEPredictModels[[y]][[r]], 
                                 newdata = PredDF)
      
      RandomPredictionList[[x]][[i]] <- Predictions %>% exp
      
      ModelLikelihood = (RealOutcomes*log(RandomPredictionList[[x]][[i]]) - RandomPredictionList[[x]][[i]]) %>% sum
      
      Deviance = -2*ModelLikelihood
      
      DevianceList[[y]][[x]][[i]] <- Deviance
    }
  }
}

RealDeviance <- InterceptDeviance <- list()

for(y in Resps){
  
  print(y)
  
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
