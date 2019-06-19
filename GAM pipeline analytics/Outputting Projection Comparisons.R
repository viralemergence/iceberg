
sp1 = AllMammaldf[2,"Sp"]
sp2 = AllMammaldf[2,"Sp2"]


MercatorAllMammaldf <- readRDS("~/Albersnet/Iceberg Output Files/MercatorAllMammaldf.rds")
load("~/Albersnet/Iceberg Output Files/MollweideAllMammaldf.Rdata")
MollweideAllMammaldf <- AllMammaldf
MercatorUncAllMammaldf <- readRDS("~/Albersnet/Iceberg Output Files/MercatorUncAllMammaldf.rds")


rbind(
  MercatorAllMammaldf %>% filter(Sp == sp1, Sp2 == sp2),
  
  MercatorUncAllMammaldf %>% filter(Sp == sp1, Sp2 == sp2),
  
  MollweideAllMammaldf %>% filter(Sp == sp1, Sp2 == sp2)
  
)


load("~/Albersnet/Iceberg Output Files/MollweideBAMList.Rdata")
MollweideBAMList <- BAMList

load("~/Albersnet/Iceberg Output Files/MercatorBAMList.Rdata")
MercatorBAMList <- BAMList

load("~/Albersnet/Iceberg Output Files/MercatorUncBAMList.Rdata")
MercatorUncBAMList <- BAMList

MercatorNewEncounters <- readRDS("~/Albersnet/Iceberg Output Files/MercatorNewEncounters.rds")
load("~/Albersnet/Iceberg Output Files/MollweideNewEncounters.Rdata")
MollweideNewEncounters <- NewEncountersList
MercatorUncNewEncounters <- readRDS("~/Albersnet/Iceberg Output Files/MercatorUncNewEncounters.rds")

sapply(MercatorNewEncounters, dim)
sapply(MercatorUncNewEncounters, dim)
sapply(MollweideNewEncounters, dim)

CompList <- list(MollweideBAMList, MercatorBAMList, MercatorUncBAMList)

FitList <- PostList <- DrawList <- list()

for(r in 1:length(CompList)){
  
  Model <- CompList[[r]][[1]]
  PostList[[r]] <- list()
  
  # Model Checking ####
  
  SpCoefNames <- names(Model$coef)[substr(names(Model$coef),1,5)=="SppSp"]
  SpCoef <- Model$coef[SpCoefNames]
  
  # Effects ####
  
  SpaceRange <- seq(from = 0,
                    to = 1,
                    length = 101) %>% 
    c(mean(DataList[[1]]$Space))
  
  PhyloRange <- seq(from = 0,
                    to = 1,
                    length = 101)  %>% 
    c(mean(DataList[[1]]$Phylo))
  
  FitList[[r]] <- expand.grid(Space = SpaceRange,
                              Phylo = PhyloRange,
                              MinCites = mean(DataList[[1]]$MinCites),
                              Domestic = 0
  ) %>%
    mutate(SpaceQuantile = ifelse(Space == last(unique(Space)), "1.5%",
                                  ifelse(Space == 0, "0%",
                                         ifelse(Space == 0.25, "25%",
                                                ifelse(Space == 0.5, "50%", NA)))),
           
           PhyloQuantile = ifelse(Phylo == last(unique(Phylo)), "0.1",
                                  ifelse(Phylo == 0, "0",
                                         ifelse(Phylo == 0.25, "0.25",
                                                ifelse(Phylo == 0.5, "0.5", NA)))),
           Gz = as.numeric(Space==0))
  
  FitList[[r]]$Spp <- matrix(0 , nrow = nrow(FitList[[r]]), ncol = length(SpCoef))
  
  FitPredictions  <- predict.gam(Model, 
                                 newdata = FitList[[r]], 
                                 se.fit = T)
  
  FitList[[r]][,c("Fit","Lower", "Upper")] <- logistic(with(FitPredictions, cbind(fit, fit - se.fit, fit + se.fit)))
  
  print("Getting posterior uncertainty!")
  
  # Posterior Uncertainty Simulation #### https://www.fromthebottomoftheheap.net/2014/06/16/simultaneous-confidence-intervals-for-derivatives/
  
  for(i in c("Space", "Phylo")){
    
    print(i)
    
    PredData <- FitList[[r]] 
    
    if(i == "Space") PredData <- PredData %>% filter(Phylo == last(unique(Phylo))) else{
      
      PredData <- PredData %>% filter(Space == 0)
      
    }
    
    lp <- predict(Model, newdata = PredData, 
                  type = "lpmatrix") %>% 
      as.data.frame()
    
    coefs <- coef(Model)
    vc <- vcov(Model)
    
    sim <- mvrnorm(100, mu = coefs, Sigma = vc)
    
    want <- lp %>% colnames
    
    lp <- lp %>% as.matrix #%>% logistic
    
    fits <- lp[, want] %*% t(sim[, want]) %>% as.data.frame() %>%
      mutate(i = PredData[,i])
    
    PostList[[r]][[i]] <- gather(fits, key = "Draw", value = "Fit", -i) %>%
      mutate(Fit = logistic(Fit))
    
  }
  
  Draw = F
  
  if(Draw){
    
    DrawList[[r]] <- list()
    
    for(i in c("Space", "Phylo")){
      
      print(i)
      
      lp = list()
      
      for(j in 1:100){
        
        print(j)
        
        PredData <- FitList[[r]]
        
        if(i == "Space"){
          
          PredData <- PredData %>% filter(Phylo == last(unique(Phylo))) %>%
            mutate(Phylo = sample(DataList[[r]]$Phylo, 1))
          
        } else {
          
          PredData <- PredData %>% filter(Space == last(unique(Space))) %>%
            mutate(Space = sample(DataList[[r]]$Space, 1))
          
        }
        
        lp[[j]] <- data.frame(Fit = predict(Model, newdata = PredData),
                              Iteration = as.factor(j),
                              i = PredData[,i])
        
      }
      
      DrawList[[r]][[i]] <- lp %>% bind_rows() %>% as.data.frame() %>%
        mutate(Fit = logistic(Fit))
      
    }
  }
  
}

save(FitList, PostList, DrawList, file = "Output Files/FitList.Rdata")

# Model Output Figure ####

for(x in 1:3){
  
  plot_grid(FitList[[x]] %>% 
              filter(!is.na(SpaceQuantile)) %>%
              ggplot(aes(Phylo, Fit, colour = SpaceQuantile)) + 
              geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = SpaceQuantile), alpha = 0.2, colour = NA) +
              geom_line(aes(group = as.factor(Space))) +
              labs(y = "Predicted Viral Sharing", x = "Phylogenetic Similarity", 
                   colour = "Overlap", fill = "Overlap") +
              lims(x = c(0,1), y = c(0,1)) +
              coord_fixed() +
              scale_color_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
              scale_fill_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
              theme(legend.position = c(0.1, 0.8), legend.background = element_rect(colour = "dark grey")) +
              geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Phylo), alpha = 0.01),
            
            FitList[[x]] %>% 
              filter(!is.na(PhyloQuantile)) %>%
              ggplot(aes(Space, Fit, colour = PhyloQuantile)) + 
              geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = PhyloQuantile), alpha = 0.2, colour = NA) +
              geom_line(aes(group = as.factor(Phylo))) +
              labs(y = "Predicted Viral Sharing", x = "Geographic Overlap", 
                   colour = "Relatedness", fill = "Relatedness") +
              lims(x = c(0,1), y = c(0,1)) +
              coord_fixed() +
              scale_color_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
              scale_fill_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
              theme(legend.position = c(0.1, 0.8), legend.background = element_rect(colour = "dark grey")) +
              geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Space), alpha = 0.01),
            
            FitList[[x]] %>% 
              filter(!Phylo == last(unique(Phylo)),
                     !Space == last(unique(Space))) %>%
              ggplot(aes(Space, Phylo)) + 
              geom_tile(aes(fill = Fit)) + 
              labs(x = "Geographic Overlap", 
                   y = "Phylogenetic Similarity",
                   fill = "Estimate") +
              #ggtitle("Tensor Field") +
              lims(x = c(0,1), y = c(0,1)) +
              coord_fixed() +
              theme(legend.position = "bottom") +
              scale_fill_continuous_sequential(palette = "Greens 2", cmax = 20, end = 1),
            
            DataList[[1]] %>%
              ggplot(aes(Space, Phylo)) + 
              labs(x = "Geographic Overlap", 
                   y = "Phylogenetic Similarity") +
              #ggtitle("Data Distribution") +
              scale_fill_continuous_sequential(palette = "purp", begin = 0.2) +
              lims(x = c(0,1), y = c(0,1)) +
              coord_fixed() +
              theme(legend.position = "bottom") +
              geom_hex(aes(fill = stat(log(count)))),
            
            nrow = 2, 
            rel_heights = c(1,1.23), 
            labels = "AUTO") %>% save_plot(filename = paste0(x,"Model Predictions.jpeg"), 
                                           #units = "mm", width = 200, height = 200,
                                           ncol = 2, # we're saving a grid plot of 2 columns
                                           nrow = 2, # and 2 rows
                                           # each individual subplot should have an aspect ratio of 1.3
                                           base_aspect_ratio = 1)
  
}


