
# Running Frequentist GAMS

# Rscript "Iceberg GAMs.R"

library(mgcv); library(tidyverse); library(ggregplot)

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

Covar <- c("s(Phylo, by = ordered(Gz))",
           "t2(Phylo, Space, by = ordered(!Gz))",
           "MinCites",
           "Domestic",
           "Spp")

BAMList <- DataList <- PPList <- list()

r = 1

for(r in 1:length(Resps[1])){
  
  print(Resps[r])
  
  DataList[[Resps[r]]] <- FinalHostMatrix[!NARows(FinalHostMatrix, Resps[r]),] %>% droplevels
  
  DataList[[Resps[r]]]$Sp <- factor(DataList[[Resps[r]]]$Sp, levels = sort(union(DataList[[Resps[r]]]$Sp,DataList[[Resps[r]]]$Sp2)))
  DataList[[Resps[r]]]$Sp2 <- factor(DataList[[Resps[r]]]$Sp2, levels = sort(union(DataList[[Resps[r]]]$Sp,DataList[[Resps[r]]]$Sp2)))
  
  DataList[[Resps[r]]] <- DataList[[Resps[r]]] %>% slice(order(Sp, Sp2))
  
  MZ1 <- model.matrix( ~ Sp - 1, data = DataList[[Resps[r]]]) %>% as.matrix
  MZ2 <- model.matrix( ~ Sp2 - 1, data = DataList[[Resps[r]]]) %>% as.matrix
  
  SppMatrix = MZ1 + MZ2
  
  DataList[[Resps[[r]]]]$Spp <- SppMatrix
  DataList[[Resps[[r]]]]$Cites <- rowSums(log(DataList[[Resps[r]]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1))
  DataList[[Resps[[r]]]]$MinCites <- apply(log(DataList[[Resps[r]]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1),1,min)
  DataList[[Resps[[r]]]]$Domestic <- ifelse(rowSums(cbind(2- DataList[[Resps[r]]]$hDom %>% as.factor %>% as.numeric,
                                                          2- DataList[[Resps[r]]]$hDom.Sp2 %>% as.factor %>% as.numeric))>0,1,0)
  
  PPList[[Resps[r]]] <- list(Spp = list(rank = nlevels(DataList[[Resps[r]]]$Sp), 
                                        diag(nlevels(DataList[[Resps[r]]]$Sp))))
  
  Formula = as.formula(paste0(Resps[r], 
                              " ~ ",
                              paste(Covar, collapse = " + ")
  ))
  
  BAMList[[Resps[r]]] <- bam(Formula,
                             data = DataList[[Resps[r]]], 
                             family = binomial(),
                             paraPen = PPList[[Resps[r]]], select = T
  )
  
}

save(DataList, PPList, BAMList, file = "Output Files/BAMList.Rdata")

FitList <- PostList <- DrawList <- list()

r = 1

for(r in 1:length(BAMList)){
  
  Model <- BAMList[[Resps[r]]]
  
  print(Resps[r])
  
  # Model Checking ####
  
  SpCoefNames <- names(Model$coef)[substr(names(Model$coef),1,5)=="SppSp"]
  SpCoef <- Model$coef[SpCoefNames]
  
  # Effects ####
  
  SpaceRange <- seq(from = 0,
                    to = 1,
                    length = 101) %>% 
    c(mean(DataList[[Resps[r]]]$Space))
  
  PhyloRange <- seq(from = 0,
                    to = 1,
                    length = 101)  %>% 
    c(mean(DataList[[Resps[r]]]$Phylo))
  
  FitList[[Resps[r]]] <- expand.grid(Space = SpaceRange,
                                     Phylo = PhyloRange,
                                     MinCites = mean(DataList[[Resps[r]]]$MinCites),
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
  
  FitList[[Resps[r]]]$Spp <- matrix(0 , nrow = nrow(FitList[[Resps[r]]]), ncol = length(SpCoef))
  
  FitPredictions  <- predict.gam(Model, 
                                 newdata = FitList[[Resps[r]]], 
                                 se.fit = T)
  
  FitList[[Resps[r]]][,c("Fit","Lower", "Upper")] <- logistic(with(FitPredictions, cbind(fit, fit - se.fit, fit + se.fit)))
  
  print("Getting posterior uncertainty!")
  
  # Posterior Uncertainty Simulation #### https://www.fromthebottomoftheheap.net/2014/06/16/simultaneous-confidence-intervals-for-derivatives/
  
  for(i in c("Space", "Phylo")){
    
    print(i)
    
    PredData <- FitList[[Resps[r]]] 
    
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
    
    PostList[[Resps[r]]][[i]] <- gather(fits, key = "Draw", value = "Fit", -i) %>%
      mutate(Fit = logistic(Fit))
    
  }
  
  Draw = F
  
  if(Draw){
    
    DrawList[[Resps[r]]] <- list()
    
    for(i in c("Space", "Phylo")){
      
      print(i)
      
      lp = list()
      
      for(j in 1:100){
        
        print(j)
        
        PredData <- FitList[[Resps[r]]]
        
        if(i == "Space"){
          
          PredData <- PredData %>% filter(Phylo == last(unique(Phylo))) %>%
            mutate(Phylo = sample(DataList[[Resps[r]]]$Phylo, 1))
          
        } else {
          
          PredData <- PredData %>% filter(Space == last(unique(Space))) %>%
            mutate(Space = sample(DataList[[Resps[r]]]$Space, 1))
          
        }
        
        lp[[j]] <- data.frame(Fit = predict(Model, newdata = PredData),
                              Iteration = as.factor(j),
                              i = PredData[,i])
        
      }
      
      DrawList[[Resps[r]]][[i]] <- lp %>% bind_rows() %>% as.data.frame() %>%
        mutate(Fit = logistic(Fit))
      
    }
  }
  
}

save(FitList, PostList, DrawList, file = "Output Files/FitList.Rdata")

# Model Output Figure ####

plot_grid(FitList[["VirusBinary"]] %>% 
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
          
          FitList[["VirusBinary"]] %>% 
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
          
          FitList[["VirusBinary"]] %>% 
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
          
          DataList$VirusBinary %>%
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
          labels = "AUTO") 
