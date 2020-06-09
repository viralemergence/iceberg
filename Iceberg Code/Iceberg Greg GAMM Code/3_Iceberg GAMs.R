# Rscript "3_Iceberg GAMs.R" ####

# Running Frequentist GAMS

# Rscript "3_Iceberg GAMs.R"

library(mgcv); library(tidyverse); library(ggregplot); library(MASS); library(cowplot); library(colorspace)
library(ggregplot); library(parallel); library(igraph); 
library(Matrix); library(ROCR)

# source("Iceberg Greg GAMM Code/2_Iceberg Data Import.R")

PipelineReps = LETTERS[1:4]

Resps <- c("VirusBinary","RNA","DNA","Vector","NVector")

Covar <- c("s(Phylo, by = ordered(Gz))",
           "t2(Phylo, Space, by = ordered(!Gz))",
           "MinCites",
           "Domestic",
           "Spp")

BAMList <- DataList <- PPList <- list()

r = 1

for(Pipeline in c("A", "B")){
  
  print(Pipeline)
  
  DataList[[Pipeline]] <- FinalHostMatrix[!NARows(FinalHostMatrix, "VirusBinary"),] %>% droplevels
  
  DataList[[Pipeline]]$Space <- DataList[[Pipeline]][,paste0("Space",Pipeline)]
  
  DataList[[Pipeline]]$Sp <- factor(DataList[[Pipeline]]$Sp, levels = sort(union(DataList[[Pipeline]]$Sp,DataList[[Pipeline]]$Sp2)))
  DataList[[Pipeline]]$Sp2 <- factor(DataList[[Pipeline]]$Sp2, levels = sort(union(DataList[[Pipeline]]$Sp,DataList[[Pipeline]]$Sp2)))
  
  DataList[[Pipeline]] <- DataList[[Pipeline]] %>% slice(order(Sp, Sp2))
  
  MZ1 <- model.matrix( ~ Sp - 1, data = DataList[[Pipeline]]) %>% as.matrix
  MZ2 <- model.matrix( ~ Sp2 - 1, data = DataList[[Pipeline]]) %>% as.matrix
  
  SppMatrix = MZ1 + MZ2
  
  DataList[[Pipeline]]$Spp <- SppMatrix
  DataList[[Pipeline]]$Cites <- rowSums(log(DataList[[Pipeline]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1))
  DataList[[Pipeline]]$MinCites <- apply(log(DataList[[Pipeline]][,c("hDiseaseZACites","hDiseaseZACites.Sp2")] + 1),1,min)
  DataList[[Pipeline]]$Domestic <- ifelse(rowSums(cbind(2- DataList[[Pipeline]]$hDom %>% as.factor %>% as.numeric,
                                                        2- DataList[[Pipeline]]$hDom.Sp2 %>% as.factor %>% as.numeric))>0,1,0)
  
  PPList[[Pipeline]] <- list(Spp = list(rank = nlevels(DataList[[Pipeline]]$Sp), 
                                        diag(nlevels(DataList[[Pipeline]]$Sp))))
  
  Formula = as.formula(paste0("VirusBinary", 
                              " ~ ",
                              paste(Covar, collapse = " + ")
  ))
  
  BAMList[[Pipeline]] <- bam(Formula,
                             data = DataList[[Pipeline]], 
                             family = binomial(),
                             paraPen = PPList[[Pipeline]], select = T
  )
  
}

save(DataList, PPList, BAMList, file = paste0("Iceberg Output Files/","BAMList.Rdata"))

FitList <- PostList <- DrawList <- list()

for(Pipeline in PipelineReps[1:length(BAMList)]){
  
  Model <- BAMList[[Pipeline]]
  
  print(Pipeline)
  
  # Model Checking ####
  
  SpCoefNames <- names(Model$coef)[substr(names(Model$coef),1,5)=="SppSp"]
  SpCoef <- Model$coef[SpCoefNames]
  
  # Effects ####
  
  SpaceRange <- seq(from = 0,
                    to = 1,
                    length = 101) %>% 
    c(mean(DataList[[Pipeline]]$Space))
  
  PhyloRange <- seq(from = 0,
                    to = 1,
                    length = 101)  %>% 
    c(mean(DataList[[Pipeline]]$Phylo))
  
  FitList[[Pipeline]] <- expand.grid(Space = SpaceRange,
                                     Phylo = PhyloRange,
                                     MinCites = mean(DataList[[Pipeline]]$MinCites),
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
  
  FitList[[Pipeline]]$Spp <- matrix(0 , nrow = nrow(FitList[[Pipeline]]), ncol = length(SpCoef))
  
  FitPredictions  <- predict.gam(Model, 
                                 newdata = FitList[[Pipeline]], 
                                 se.fit = T)
  
  FitList[[Pipeline]][,c("Fit","Lower", "Upper")] <- logistic(with(FitPredictions, cbind(fit, fit - se.fit, fit + se.fit)))
  
  print("Getting posterior uncertainty!")
  
  # Posterior Uncertainty Simulation #### https://www.fromthebottomoftheheap.net/2014/06/16/simultaneous-confidence-intervals-for-derivatives/
  
  for(i in c("Space", "Phylo")){
    
    print(i)
    
    PredData <- FitList[[Pipeline]] 
    
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
    
    PostList[[Pipeline]][[i]] <- gather(fits, key = "Draw", value = "Fit", -i) %>%
      mutate(Fit = logistic(Fit))
    
  }
}

save(FitList, PostList, DrawList, file = paste0("Iceberg Output Files/","FitList.Rdata"))

Output = F

if(Output){
  
  # Model Output Figure ####
  
  # pdf("GAMOutput.pdf", width = 9, height = 8)
  
  Pipeline <- "A"
  
  plot_grid(FitList[[Pipeline]] %>% 
              filter(!is.na(SpaceQuantile)) %>%
              ggplot(aes(Phylo, Fit, colour = SpaceQuantile)) + theme_cowplot() + 
              geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = SpaceQuantile), alpha = 0.2, colour = NA) +
              geom_line(aes(group = as.factor(Space))) +
              labs(y = "Viral sharing probability", x = "Phylogenetic similarity", 
                   colour = "Geographic overlap", fill = "Geographic overlap") +
              lims(x = c(0,1), y = c(0,1)) +
              coord_fixed() +
              scale_color_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
              scale_fill_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
              theme(legend.position = c(0.1, 0.8), 
                    legend.title = element_text(size = 10),
                    legend.background = element_rect(colour = "dark grey")) +
              geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Phylo), alpha = 0.01),
            
            FitList[[Pipeline]] %>% 
              filter(!is.na(PhyloQuantile)) %>%
              ggplot(aes(Space, Fit, colour = PhyloQuantile)) +  theme_cowplot() + 
              geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = PhyloQuantile), alpha = 0.2, colour = NA) +
              geom_line(aes(group = as.factor(Phylo))) +
              labs(y = "Viral sharing probability", x = "Geographic overlap", 
                   colour = "Phylogenetic similarity", fill = "Phylogenetic similarity") +
              lims(x = c(0,1), y = c(0,1)) +
              coord_fixed() +
              scale_color_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
              scale_fill_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
              theme(legend.position = c(0.1, 0.8), 
                    legend.title = element_text(size = 10),
                    legend.background = element_rect(colour = "dark grey")) +
              geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Space), alpha = 0.01),
            
            FitList[[Pipeline]] %>% 
              filter(!Phylo == last(unique(Phylo)),
                     !Space == last(unique(Space))) %>%
              ggplot(aes(Space, Phylo)) +  theme_cowplot() + 
              geom_tile(aes(fill = Fit)) + 
              labs(x = "Geographic overlap", 
                   y = "Phylogenetic similarity",
                   fill = "Viral sharing\nprobability") +
              #ggtitle("Tensor Field") +
              lims(x = c(0,1), y = c(0,1)) +
              coord_fixed() +
              theme(legend.position = "bottom",
                    legend.title = element_text(size = 10)) +
              geom_contour(aes(z = Fit), colour = "white", alpha = 0.8) + 
              metR::geom_text_contour(aes(z = Fit), colour = "white", size = 2.5, hjust = 0.5, vjust = 1.1, check_overlap = T) +
              scale_fill_continuous_sequential(palette = "ag_GrnYl",
                                               limits = c(0,1),
                                               breaks = c(0,0.5,1)),
            
            DataList[[Pipeline]] %>%
              ggplot(aes(Space, Phylo)) +  theme_cowplot() + 
              labs(x = "Geographic overlap", 
                   y = "Phylogenetic similarity") +
              #ggtitle("Data Distribution") +
              scale_fill_continuous_sequential(palette = "Heat 2", breaks = c(0:2*5)) +
              lims(x = c(0,1), y = c(0,1)) +
              coord_fixed() +
              theme(legend.position = "bottom") +
              geom_hex(aes(fill = stat(log(count)))),
            
            nrow = 2, 
            rel_heights = c(1,1.23), 
            labels = "AUTO") 
  
  dev.off()
  
  Resps <- c("VirusBinary")#,"RNA","DNA","Vector","NVector")
  
  # Simulating the predicted networks ####
  
  # Doing the simulating with random effects #####
  
  PredList1 <- list()
  
  Predictions1 <- predict.bam(BAMList[[1]], 
                              newdata = DataList[[1]],
                              type = "terms")
  
  Intercept1 <- attr(Predictions1, "constant")
  
  Predictions1 <- Predictions1 %>% as.data.frame()
  
  Predictions1$Intercept <- Intercept1
  
  N = nrow(DataList[[1]])
  
  PredList1 <- parallel::mclapply(1:1000, function(x){ # to do something non-specific
    
    BinPred <- rbinom(n = N,
                      prob = logistic(rowSums(Predictions1)),
                      size  = 1)
    
    BinPred
    
  }, mc.cores = 10)
  
  PredDF1 <- data.frame(PredList1)
  
  FinalHostMatrix$PredVirus1 <- apply(PredDF1, 1, mean)
  FinalHostMatrix$PredVirus1Q <- cut(FinalHostMatrix$PredVirus1,
                                     breaks = c(-1:10/10),
                                     labels = c(0:10/10))
  
  SimNets1 <- mclapply(1:length(PredList1), function(i){
    
    AssMat <- matrix(NA, 
                     nrow = nlevels(DataList[[1]]$Sp), 
                     ncol = nlevels(DataList[[1]]$Sp))
    
    AssMat[lower.tri(AssMat)] <- round(PredList1[[i]])
    AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
    diag(AssMat) <- 0
    dimnames(AssMat) <- list(levels(DataList[[1]]$Sp),
                             levels(DataList[[1]]$Sp))
    
    as(AssMat, "dgCMatrix")
    
  }, mc.cores = 10)
  
  SimGraphs1 <- mclapply(1:length(PredList1), function(i){
    
    graph.adjacency(SimNets1[[i]], mode = "undirected")
    
  }, mc.cores = 10)
  
  # Doing the simulating without random effects #####
  
  SpCoefNames <- names(BAMList[[1]]$coef)[substr(names(BAMList[[1]]$coef),1,5)=="SppSp"]
  SpCoef <- BAMList[[1]]$coef[SpCoefNames]
  
  Predictions1b <- predict.bam(BAMList[[1]], 
                               newdata = DataList[[1]],# %>% select(-Spp),
                               type = "terms",
                               exclude = "Spp")
  
  Intercept1b <- attr(Predictions1b, "constant")
  
  Predictions1b <- Predictions1b %>% as.data.frame
  
  N = nrow(DataList[[1]])
  
  PredList1b <- parallel::mclapply(1:1000, function(x){ # to do something non-specific
    
    Predictions1b[,"Spp"] <- sample(SpCoef, N, replace = T) + 
      sample(SpCoef, N, replace = T)
    
    Predictions1b[,"Intercept"] <- Intercept1b
    
    BinPred <- rbinom(n = N,
                      prob = logistic(rowSums(Predictions1b)),
                      size  = 1)
    
    BinPred
    
  }, mc.cores = 10)
  
  PredDF1b <- data.frame(PredList1b)
  FinalHostMatrix$PredVirus1b <- apply(PredDF1b, 1, mean)
  FinalHostMatrix$PredVirus1bQ <- cut(FinalHostMatrix$PredVirus1b,
                                      breaks = c(-1:10/10),
                                      labels = c(0:10/10))
  
  SimNets1b <- mclapply(1:length(PredList1b), function(i){
    
    AssMat <- matrix(NA, 
                     nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                     ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
    
    AssMat[lower.tri(AssMat)] <- round(PredList1b[[i]])
    AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
    diag(AssMat) <- 0
    dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                             union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
    
    as(AssMat, "dgCMatrix")
    
  }, mc.cores = 10)
  
  SimGraphs1b <- mclapply(1:length(PredList1b), function(i){
    
    graph.adjacency(SimNets1b[[i]], mode = "undirected")
    
  }, mc.cores = 10)
  
  # Doing the simulating with only random effects #####
  
  SpCoefNames <- names(BAMList[[1]]$coef)[substr(names(BAMList[[1]]$coef),1,5)=="SppSp"]
  SpCoef <- BAMList[[1]]$coef[SpCoefNames]
  
  Predictions1c <- predict.bam(BAMList[[1]], 
                               newdata = DataList[[1]] %>% 
                                 mutate_at(vars(Space, Phylo, MinCites), function(a) mean(a)) %>%
                                 mutate(Gz = 0),
                               type = "terms")
  
  Intercept1c <- attr(Predictions1c, "constant")
  
  Predictions1c <- Predictions1c %>% as.data.frame
  
  N = nrow(FinalHostMatrix)
  
  PredList1c <- parallel::mclapply(1:1000, function(x){ # to do something non-specific
    
    Predictions1c[,"Intercept"] <- Intercept1c
    
    BinPred <- rbinom(n = N,
                      prob = logistic(rowSums(Predictions1c)),
                      size  = 1)
    
    BinPred
    
  }, mc.cores = 10)
  
  PredDF1c <- data.frame(PredList1c)
  FinalHostMatrix$PredVirus1c <- apply(PredDF1c, 1, mean)
  FinalHostMatrix$PredVirus1cQ <- cut(FinalHostMatrix$PredVirus1c,
                                      breaks = c(-1:10/10),
                                      labels = c(0:10/10))
  
  SimNets1c <- mclapply(1:length(PredList1c), function(i){
    
    AssMat <- matrix(NA, 
                     nrow = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)), 
                     ncol = length(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2)))
    
    AssMat[lower.tri(AssMat)] <- round(PredList1c[[i]])
    AssMat[upper.tri(AssMat)] <- t(AssMat)[!is.na(t(AssMat))]
    diag(AssMat) <- 0
    dimnames(AssMat) <- list(union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2),
                             union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
    
    as(AssMat, "dgCMatrix")
    
  }, mc.cores = 10)
  
  SimGraphs1c <- mclapply(1:length(PredList1c), function(i){
    
    graph.adjacency(SimNets1c[[i]], mode = "undirected")
    
  }, mc.cores = 10)
  
  save(SimGraphs1, SimGraphs1b, SimGraphs1c, file = paste0("Iceberg Output Files/KnownSimGraphs.Rdata"))
  
  print(Sys.time())
  
  # Checking AUC
  
  AUCList <- list()
  
  for(x in c("PredVirus1", "PredVirus1b","PredVirus1c")){
    
    df <- FinalHostMatrix[,c("VirusBinary",x)]
    
    colnames(df) <- c('observed','predicted')
    
    df$observed[df$observed>1] = 1
    
    pred <- prediction(df$predicted, df$observed)
    
    AUCList[[x]] <- performance(pred,"auc")
      
  }
  
}

saveRDS(FinalHostMatrix, file = "Iceberg Output Files/KnownPredictedDegree.rds")
