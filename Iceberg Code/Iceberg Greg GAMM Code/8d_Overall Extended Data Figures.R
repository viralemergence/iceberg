
# Extended Data ####

library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr); library(SpRanger); library(cowplot)
library(fs); library(tidyverse); library(Matrix); library(parallel); library(mgcv); library(SpRanger); library(raster)
library(sf); library(fasterize);library(ggregplot); library(igraph);library(maptools); library(glue)
library(colorspace)

setwd("~/Albersnet/Iceberg Files/Summary")

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

dir_create("Iceberg Extended Data")

PredReps <- c("Currents", paste0("Futures", 1:4))

PipelineReps <- LETTERS[1:4]

SpaceVars <- paste0(paste("Space", PredReps, sep = "."),rep(PipelineReps, each = length(PredReps)))
SharingVars <- paste0(paste("Sharing",PredReps, sep = "."), rep(PipelineReps, each = length(PredReps)))

names(SpaceVars) <- names(SharingVars) <- paste0(PredReps,rep(PipelineReps, each = length(PredReps)))

Panth1 <- read.delim("~/Albersnet/data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order, hFamily = MSW05_Family)

Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

ClimateReps <- glue("Climate{1:9}")

"~/Albersnet/Iceberg Revision" %>% list.files(full.names = T) %>% extract2(1) %>%
  
  paste0("/PPM") %>% list.files(full.names = T) %>% extract2(2) %>%
  list.files %>%
  substr(1, 2) %>%
  unique ->
  
  CoryClimateReps

names(CoryClimateReps) <- ClimateReps

Unlist1 <- function(a) unlist(a, recursive = F)

# Overall trends ###

# AllMammaldf <- readRDS("~/Albersnet/Iceberg Output Files/AllMammaldf.rds")

glue("~/Albersnet/Iceberg Files/{ClimateReps}/Iceberg Output Files/NewEncounters.rds") %>% 
  map(readRDS) -> 
  NewEncountersListList

names(NewEncountersListList) <- ClimateReps

NewEncountersListList %<>% 
  map(Unlist1) %>% 
  Unlist1

NewEncountersListList %>% names 

# Import #####

CurrentsGridDF <- readRDS("~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/CurrentsGridDF.rds")
FuturesGridDF <- readRDS("~/Albersnet/Iceberg Files/Summary/Iceberg Output Files/FuturesGridDF.rds")

CurrentsGridDF %<>% arrange(X, Y)
FuturesGridDF %<>% arrange(X, Y)

FuturesGridDF$Y <- CurrentsGridDF$Y

PredReps <- c("Currents", paste0("Futures", 1:4))
PipelineReps <- LETTERS[1:4]

AreaRaster <- raster("Iceberg Input Files/LandArea.asc")
AreaValues <- raster::values(AreaRaster)

# IcebergAdjList <- readRDS(paste0("Iceberg Output Files/","IcebergAdjList.rds"))
# NewEncountersList <- readRDS(paste0("Iceberg Output Files/","NewEncounters.rds"))
#  AllMammaldf <- readRDS(paste0("Iceberg Output Files/","AllMammaldf.rds"))

SpaceVars <- paste0(paste("Space", PredReps, sep = "."),rep(PipelineReps, each = length(PredReps)))
SharingVars <- paste0(paste("Sharing",PredReps, sep = "."), rep(PipelineReps, each = length(PredReps)))

names(SpaceVars) <- names(SharingVars) <- paste0(PredReps,rep(PipelineReps, each = length(PredReps)))

# Setting up stuff ####

GridList <- NEGridList <- list()

CORES = 65

t1 <- Sys.time()

# Blanks
blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

Species <- 
  paste0("Iceberg Input Files/GretCDF/Currents") %>% list.files() %>% str_remove(".rds$") %>%
  sort

paste0("Iceberg Input Files/GretCDF/Currents") %>% 
  list.files(full.names = T) ->
  CurrentFiles

paste0("Iceberg Input Files/GretCDF/Currents") %>% 
  list.files() %>% str_remove(".rds$") ->
  names(CurrentFiles)

paste0("Iceberg Input Files/GretCDF/Futures") %>% 
  list.files(full.names = T) ->
  FutureFiles

paste0("Iceberg Input Files/GretCDF/Futures") %>% 
  list.files() %>% str_remove(".rds$") ->
  names(FutureFiles)

FutureCDFList <- list()

Species <- Species %>% sort %>% intersect(names(CurrentFiles)) #%>% intersect(names(FutureFiles))

AllNESpecies <- 
  
  NewEncountersList %>% unlist(recursive = F) %>% lapply(function(a){
    
    a %>% dplyr::select(Sp, Sp2) %>% unlist %>% unique
    
  }) %>% reduce(union) %>% sort

Species <- intersect(Species, AllNESpecies)

Panth1 %>% filter(Sp %in% Species, !hOrder == "Chiroptera") %>% 
  pull(hOrder) %>% unique -> MammalOrders

#MammalOrders %>% map(~readRDS(paste0("Iceberg Output Files/", .x, "NoBatNewIntersects.rds"))) -> 
#  OrderGridList 

PredReps[2:5] %>% rep(4) %>% 
  paste0(rep(LETTERS[1:4], each = 4)) -> NEReps

recode1 <- c(A = "Climate + Land + Dispersal",
             B = "Climate + Dispersal",
             C = "Climate + Land",
             D = "Climate Only")

recode2 <- c(Futures1 = 'RCP 2.6',
             Futures2 = 'RCP 4.5',
             Futures3 = 'RCP 7.0',
             Futures4 = 'RCP 8.5')

PipeLabels <- rep(recode1, each = 4) %>% paste(rep(recode2, 4))
names(PipeLabels) <- NEReps

# Overall maps ####

NewIntersects <- readRDS("~/Albersnet/Iceberg Files/Summary/Iceberg Output Files/NewIntersects.rds")

# New Encounters A and C ####

NewIntersects %>% 
  dplyr::select(X, Y, starts_with("Overlap.Futures")) %>%
  dplyr::select(-ends_with("B"), -ends_with("D")) %>%
  tidyr::gather("Key", "Encounters", -c(X,Y)) %>%
  mutate(Key = Key %>% str_remove("Overlap.")) %>%
  tidyr::separate(Key, sep = 8, into = c("RCP", "Pipeline")) %>%
  mutate(Pipeline = recode1[Pipeline],
         RCP = recode2[RCP]) %>%
  ggplot(., aes(X, Y, fill = log10(Encounters + 1))) + geom_tile() + 
  scale_fill_continuous_sequential(palette = AlberPalettes[[1]],
                                   labels = c(1, 11, 101, 1001, 10001)) +
  facet_grid(RCP~Pipeline) + 
  coord_fixed() +
  theme_void() + labs(fill = "Encounters") +
  theme(legend.position = "bottom", legend.text = element_text(angle = 45, hjust = 1)) +
  ggsave(filename = paste0("Iceberg Extended Data/", "A_C_Encounters", ".jpg"), width = 10, height = 8)


# New Sharings A and C ####

NewIntersects %>% dplyr::select(X, Y, starts_with("DeltaOverlapSharing.Futures")) %>%
  dplyr::select(-ends_with("B"), -ends_with("D")) %>%
  tidyr::gather("Key", "NewSharings", -c(X,Y)) %>%
  mutate(Key = Key %>% str_remove("DeltaOverlapSharing.")) %>%
  tidyr::separate(Key, sep = 8, into = c("RCP", "Pipeline")) %>%
  mutate(Pipeline = recode1[Pipeline],
         RCP = recode2[RCP]) %>%
  ggplot(., aes(X, Y, fill = log10(NewSharings + 1))) + geom_tile() + 
  scale_fill_continuous_sequential(palette = AlberPalettes[[1]],
                                   labels = c(1, 11, 101, 1001)) +
  facet_grid(RCP~Pipeline) + 
  coord_fixed() +
  theme_void() + labs(fill = "New sharings") +
  theme(legend.position = "bottom", legend.text = element_text(angle = 45, hjust = 1)) +
  ggsave(filename = paste0("Iceberg Extended Data/", "A_C_NewSharings", ".jpg"), width = 10, height = 8)

# Delta Sharings A and C ####

CurrentsGridDF %>% left_join(FuturesGridDF) %>% 
  dplyr::select(X, Y, 
                contains("ClimateLandUse"),
                Sharing.CurrentsA, starts_with("Sharing.Futures")) %>%
  dplyr::select(-ends_with("B"), 
                -ends_with("D")) %>%
  mutate_at(vars(Sharing.Futures1A:Sharing.Futures4C), 
            ~.x - Sharing.CurrentsA) %>%
  mutate(
    
    Sharing.Futures1A = Sharing.Futures1A/BufferClimateLandUse.Futures1,
    Sharing.Futures1C = Sharing.Futures1C/ClimateLandUse.Futures1,
    Sharing.Futures2A = Sharing.Futures2A/BufferClimateLandUse.Futures2,
    Sharing.Futures2C = Sharing.Futures2C/ClimateLandUse.Futures2,
    Sharing.Futures3A = Sharing.Futures3A/BufferClimateLandUse.Futures3,
    Sharing.Futures3C = Sharing.Futures3C/ClimateLandUse.Futures3,
    Sharing.Futures4A = Sharing.Futures4A/BufferClimateLandUse.Futures4,
    Sharing.Futures4C = Sharing.Futures4C/ClimateLandUse.Futures4
    
  ) %>% 
  dplyr::select(X, Y, contains("Sharing.Futures")) %>%
  tidyr::gather("Key", "DeltaSharing", -c(X,Y)) %>%
  mutate(Key = Key %>% str_remove("Sharing.")) %>%
  mutate(DeltaSharing = ifelse(DeltaSharing == "Inf", 0, DeltaSharing) %>% as.numeric) %>%
  tidyr::separate(Key, sep = 8, into = c("RCP", "Pipeline")) %>%
  mutate(Pipeline = recode1[Pipeline],
         RCP = recode2[RCP]) %>%
  filter(!RCP == "8.5") %>%
  ggplot(., aes(X, Y, fill = (DeltaSharing))) + geom_tile() + 
  scale_fill_continuous_sequential(palette = AlberPalettes[[1]]) +
  facet_grid(RCP~Pipeline) + coord_fixed() +
  theme(plot.title = element_text(size = 0.5))

CurrentsGridDF %>% left_join(FuturesGridDF) %>% 
  dplyr::select(X, Y, 
                contains("ClimateLandUse"),
                Sharing.CurrentsA, starts_with("Sharing.Futures")) %>%
  dplyr::select(-ends_with("B"), 
                -ends_with("D")) %>%
  mutate_at(vars(Sharing.Futures1A:Sharing.Futures4C), 
            ~.x - Sharing.CurrentsA) %>%
  # mutate(
  #    
  #   Sharing.Futures1A = Sharing.Futures1A/BufferClimateLandUse.Futures1,
  #   Sharing.Futures1C = Sharing.Futures1C/ClimateLandUse.Futures1,
  #    Sharing.Futures2A = Sharing.Futures2A/BufferClimateLandUse.Futures2,
  #    Sharing.Futures2C = Sharing.Futures2C/ClimateLandUse.Futures2,
  #    Sharing.Futures3A = Sharing.Futures3A/BufferClimateLandUse.Futures3,
  #    Sharing.Futures3C = Sharing.Futures3C/ClimateLandUse.Futures3,
  #    Sharing.Futures4A = Sharing.Futures4A/BufferClimateLandUse.Futures4,
  #    Sharing.Futures4C = Sharing.Futures4C/ClimateLandUse.Futures4
  #    
#  ) %>% 
dplyr::select(X, Y, contains("Sharing.Futures")) %>%
  tidyr::gather("Key", "DeltaSharing", -c(X,Y)) %>%
  mutate(Key = Key %>% str_remove("Sharing.")) %>%
  mutate(DeltaSharing = ifelse(DeltaSharing == "Inf", 0, DeltaSharing) %>% as.numeric) %>%
  mutate_at("DeltaSharing", ~kader:::cuberoot(.x)) %>%
  tidyr::separate(Key, sep = 8, into = c("RCP", "Pipeline")) %>%
  mutate(Pipeline = recode1[Pipeline],
         RCP = recode2[RCP]) %>%
  filter(RCP == "RCP 8.5") %>%
  ggplot(., aes(X, Y, fill = (DeltaSharing))) + geom_tile() + 
  scale_fill_continuous_diverging(palette = "Red-Green") +
  facet_grid(RCP~Pipeline) + coord_fixed() +
  theme(plot.title = element_text(size = 0.5)) +
  ggsave(filename = paste0("Iceberg Extended Data/", "A_C_MeanDeltaSharings", ".jpg"), width = 10, height = 8)


dplyr::select(X, Y, starts_with("Overlap.")) %>% 
  rename_all(~.x %>% str_remove("Overlap.")) %>%
  tidyr::gather("Pipeline", "Encounters", -c(X,Y)) %>%
  ggplot(., aes(X, Y, fill = log10(Encounters + 1))) + geom_tile() + 
  scale_fill_continuous_sequential(palette = AlberPalettes[[1]]) +
  facet_wrap(~Pipeline, labeller = as_labeller(PipeLabels)) + coord_fixed() +
  theme(plot.title = element_text(size = 0.5)) +
  ggsave(filename = paste0("Iceberg Extended Data/", "NE_", MammalOrders[i], ".jpg"), width = 10, height = 8)

# New encounters by order ####

AllMammals <- readRDS("~/Albersnet/Iceberg Files/Summary/AllMammals.rds")

AllMammalDF <- readRDS("~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/AllMammaldf.rds")

.x <- 1

Panth1 %>% 
  filter(Sp%in%(AllMammals)) %>%
  group_by(hOrder) %>% count -> 
  
  OrderN

NewEncountersListList[str_detect(names(NewEncountersListList), "A.Futures1")] ->
  NewNewEncountersListList

NullSlope <- AllMammals %>% length %>% 
  divide_by(sapply(NewNewEncountersListList, nrow)*2, .) %>% 
  mean

map(1:length(NewNewEncountersListList), function(.x){
  
  NewNewEncountersListList[[.x]] %>% 
    bind_rows(NewNewEncountersListList[[.x]] %>% 
                rename(Sp = Sp2, Sp2 = Sp,
                       hOrder.x = hOrder.y)) %>%
    group_by(hOrder.x) %>% 
    count() -> 
    NEOrderN
  
}) -> OrderList

OrderList %>% 
  bind_rows(.id = "ClimateRep") -> OrderListDF

OrderListDF %>% 
  #group_by(ClimateRep) %>% count()
  group_by(hOrder.x) %>% 
  summarise_at("n", mean) -> NEOrderN

OrderN %>% left_join(NEOrderN, by = c("hOrder" = "hOrder.x"), 
                     suffix = c(".Overall", ".NE")) %>% 
  ggplot(aes(n.Overall, n.NE)) +
  geom_abline(intercept = log10(NullSlope), lty = 2, alpha = 0.5) +
  # ggpubr::stat_cor() + 
  geom_point(alpha = 0.8, colour = AlberColours[["BuPu"]]) + 
  ggrepel::geom_label_repel(
    force = 50,
    data = OrderN %>% 
      left_join(NEOrderN, by = c("hOrder" = "hOrder.x"), 
                suffix = c(".Overall", ".NE")) %>%
      filter(n.Overall>100),
    aes(label = hOrder)) + 
  scale_y_log10() + scale_x_log10() +
  labs(x = "Order Richness", y = "New Encounters") +
  coord_equal() + 
  ggsave("Iceberg Extended Data/OrderNE.jpg", width = 7, height = 6)

# Range changes ####

if(file.exists("~/Albersnet/Iceberg Files/Summary/Iceberg Output Files/FullRangeChangeDF.rds")){
  
  RangeChangeDF <- readRDS("~/Albersnet/Iceberg Files/Summary/Iceberg Output Files/FullRangeChangeDF.rds")
  
}else{
  
  list.files("Iceberg Input Files/GretCDF/Currents") %>% str_remove(".rds$") -> 
    
    Species
  
  i <- 1
  
  map(1:length(Species), function(i){
    
    print(Species[i])
    
    FocalSp <- Species[i]
    
    GretCDF <- 
      
      readRDS(paste0("Iceberg Input Files/GretCDF/Currents/", FocalSp, ".rds")) %>% as.matrix %>% as.data.frame() %>%
      full_join(readRDS(paste0("Iceberg Input Files/GretCDF/Futures/", FocalSp, ".rds")) %>% 
                  as.matrix %>% as.data.frame() %>% dplyr::select(-c(Continent, BufferClimate)),
                by = c("X", "Y"), 
                suffix = c(".Currents", ".Futures"))
    
    GretCDF %>% dplyr::select(contains("LandUse"), -starts_with("LandUse"), -c(LandUse, ClimateLandUse)) %>%
      summarise_all(sum) %>% divide_by(sum(GretCDF$ClimateLandUse)) -> CLU
    
    GretCDF %>% dplyr::select(contains("Climate"), -contains("LandUse"), -c(X, Y, Climate, BufferClimate)) %>%
      summarise_all(sum) %>% divide_by(sum(GretCDF$Climate)) -> C
    
    bind_cols(CLU, C) %>% return
    
  }) -> RangeChange
  
  RangeChange %>% bind_rows %>% mutate(Species = Species) -> RangeChangeDF
  
  saveRDS(RangeChangeDF, file = "Iceberg Output Files/RangeChangeDF.rds")
  
  map(1:length(Species), function(i){
    
    print(Species[i])
    
    FocalSp <- Species[i]
    
    GretCDF <- 
      
      readRDS(paste0("Iceberg Input Files/GretCDF/Currents/", FocalSp, ".rds")) %>% as.matrix %>% as.data.frame()
    
    (sum(GretCDF$ClimateLandUse - GretCDF$Climate)/
        sum(GretCDF$Climate))
    
  }) -> LandUseLossCurrents
  
  
  map(1:length(Species), function(i){
    
    print(Species[i])
    
    FocalSp <- Species[i]
    
    GretCDF <- 
      
      readRDS(paste0("Iceberg Input Files/GretCDF/Currents/", FocalSp, ".rds")) %>% as.matrix %>% as.data.frame() %>%
      full_join(readRDS(paste0("Iceberg Input Files/GretCDF/Futures/", FocalSp, ".rds")) %>% 
                  as.matrix %>% as.data.frame() %>% dplyr::select(-c(Continent, BufferClimate)),
                by = c("X", "Y"), 
                suffix = c(".Currents", ".Futures"))
    
    GretCDF %>%
      as_tibble %>% 
      dplyr::select(-contains("LandUse"), 
                    -c(X,Y, BufferClimate, Climate, IUCN, Continent)) %>%
      names -> C
    
    C %>% str_replace("Climate", "ClimateLandUse") -> CLU
    
    (((GretCDF[,CLU] - GretCDF[,C])) %>% colSums(na.rm = T))/colSums(GretCDF[,C])->
      
      DF
    
    DF %>% return
    
  }) -> LandUseLoss
  
  LandUseLoss %>% bind_rows %>% mutate(Species = Species) -> LandUseLossDF
  
  LandUseLoss %>% map_dfc(~.x %>% as_tibble) %>% t %>% as.data.frame -> LandUseLossDF
  
  names(LandUseLossDF) <- names(LandUseLoss[[1]])
  
  LandUseLossDF$Species <- Species
  
}

LandUseLossDF$BufferClimateLandUse.Futures1 %>% mean(na.rm = T)

RangeChangeDF %>% 
  group_by(Rep, ClimateRep) %>% 
  dplyr::select(Rep, OverallChange, ClimateRep) %>% 
  summarise(OverallChange = median(OverallChange, na.rm = T)*100 - 100) %>% 
  ungroup %>% group_by(Rep) %>% summarise(OverallChangeSD = sd(OverallChange),
                                          OverallChange = mean(OverallChange)
  ) ->
  PercentChanges

PercentChanges %>% 
  mutate(Pipeline = rep(LETTERS[c(2, 1, 4, 3)], each = 4)) %>%
  mutate(RCP = rep(PredReps[2:5], 4)) ->
  PercentChanges

# Table making ####

glue("~/Albersnet/Iceberg Files/{ClimateReps}/Iceberg Output Files/NewEncounters.rds") %>% 
  map(readRDS) -> 
  NewEncountersListList

names(NewEncountersListList) <- ClimateReps

map(ClimateReps, function(.x){ 
  
  print(.x)
  
  if(.x != "Climate5"){
    
    lapply(PipelineReps, function(a){
      
      print(a)
      
      lapply(PredReps[2:5], function(b){
        
        print(b)
        
        NewEncountersListList[[.x]][[a]][[b]]$SharingVar <- NewEncountersListList[[.x]][[a]][[b]][,SharingVars[paste0(b,a)]]
        NewEncountersListList[[.x]][[a]][[b]]$DeltaSharingVar <- NewEncountersListList[[.x]][[a]][[b]][,paste0("Delta",SharingVars[paste0(b,a)])]
        
        NewEncountersListList[[.x]][[a]][[b]] %>% dplyr::summarise(
          Number = n(),
          Mean.Sharing = mean(SharingVar),
          Sum.Sharing = sum(SharingVar),
          New.Sharing = sum(DeltaSharingVar)
        )
        
      }) %>% bind_rows(.id = "PredRep")
      
    }) %>% bind_rows(.id = "Pipeline") %>%
      mutate(Pipeline = PipelineReps[as.numeric(Pipeline)],
             PredRep = PredReps[as.numeric(PredRep)+1])  %>%
      rename(RCP = PredRep) -> tablex
    
  }else{
    
    lapply(PipelineReps, function(a){
      
      print(a)
      
      lapply(PredReps[c(2, 4)], function(b){
        
        print(b)
        
        NewEncountersListList[[.x]][[a]][[b]]$SharingVar <- 
          NewEncountersListList[[.x]][[a]][[b]][,SharingVars[paste0(b,a)]]
        
        NewEncountersListList[[.x]][[a]][[b]]$DeltaSharingVar <- 
          NewEncountersListList[[.x]][[a]][[b]][,paste0("Delta",SharingVars[paste0(b,a)])]
        
        NewEncountersListList[[.x]][[a]][[b]] %>% 
          dplyr::summarise(
            Number = n(),
            Mean.Sharing = mean(SharingVar),
            Sum.Sharing = sum(SharingVar),
            New.Sharing = sum(DeltaSharingVar)
          )
        
      }) %>% bind_rows(.id = "PredRep")
      
    }) %>% bind_rows(.id = "Pipeline") %>%
      mutate(Pipeline = PipelineReps[as.numeric(Pipeline)],
             PredRep = PredReps[c(2, 4)][as.numeric(PredRep)])  %>%
      rename(RCP = PredRep) -> tablex
    
  }
  
}) -> TableXList

TableXList %>% bind_rows(.id = "ClimateRep") -> 
  
  TableXDF

TableXDF %>% group_by(Pipeline, RCP) %>% summarise_all(mean) -> tablex
TableXDF %>% group_by(Pipeline, RCP) %>% summarise_all(sd) %>% rename_all(~paste0(.x, ".sd")) -> tablexsd

tablex %<>% bind_cols(tablexsd)

tablex %<>% ungroup

FilllAlpha <- 0.7

PredReps[2:5] %>% rep(4) %>% 
  paste0(rep(LETTERS[1:4], each = 4)) -> NEReps

recode1 <- c(A = "Climate + Land + Dispersal",
             B = "Climate + Dispersal",
             C = "Climate + Land",
             D = "Climate Only")

recode2 <- c(Futures1 = 'RCP 2.6',
             Futures2 = 'RCP 4.5',
             Futures3 = 'RCP 7.0',
             Futures4 = 'RCP 8.5')

PipeLabels <- rep(recode1, each = 4) %>% paste(rep(recode2, 4))
names(PipeLabels) <- NEReps

tablex %>%
  mutate(
    Pipeline = recode1[Pipeline] %>% factor(levels = recode1[1:4]),
    RCP = recode2[RCP] %>% factor(levels = recode2[1:4]) %>% 
      str_remove_all("RCP ")
  ) -> tablex

PercentChanges %>%
  mutate(
    Pipeline = recode1[Pipeline] %>% factor(levels = recode1[1:4]),
    RCP = recode2[RCP] %>% factor(levels = recode2[1:4]) %>% 
      str_remove_all("RCP ")
  ) -> PercentChanges

glue("~/Albersnet/Iceberg Files/{ClimateReps}/Iceberg Output Files/AllMammaldf.rds") %>% 
  map(readRDS) -> 
  AllMammalDFList

names(AllMammalDFList) <- ClimateReps

AllMammalDFList %>% map(~.x %>% 
                          summarise_at(vars(starts_with("DeltaSharing")), ~mean(.x)) %>% 
                          reshape2::melt() %>% 
                          mutate(Variable = str_remove(variable, "DeltaSharing.")) %>%
                          rename(DeltaSharing = value) %>%
                          mutate(Pipeline = substr(Variable, nchar(Variable), nchar(Variable)),
                                 RCP = substr(Variable, 1, nchar(Variable)-1)) %>%
                          mutate(
                            RCP = recode2[RCP] %>% factor(levels = recode2[1:4]) %>% 
                              str_remove_all("RCP ")
                          )) %>% 
  bind_rows(.id = "ClimateRep") %>% 
  group_by(Variable, Pipeline, RCP, variable) %>% 
  summarise(Mean = mean(DeltaSharing), SD = sd(DeltaSharing)) -> 
  
  DeltasDF

DeltasDF %<>% mutate(Label = glue::glue("+ {(round(Mean, 4))}\n({round(SD, 4)})"))

# Heatmaps ####

# Range Change ####

PercentChanges %>% 
  #mutate_at("Pipeline")
  ggplot(aes(Pipeline, RCP)) + 
  geom_tile(colour = "black", alpha = FilllAlpha,
            aes(fill = OverallChange)) +
  geom_text(aes(label = paste0(round(OverallChange), "\n(", round(OverallChangeSD, 2), ")")),
            size = 3) + 
  labs(fill = "% Range\nchange") +
  scale_fill_continuous_sequential(palette = AlberPalettes[[2]], rev = T) +
  #theme(legend.position = "top") +
  scale_x_discrete(labels = recode1 %>% 
                     str_replace_all("[+] ", "+\n") %>% str_remove(" Only")) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45, size = 8)) +
  coord_fixed() +
  
  # 1E's
  
  tablex %>% ggplot(aes(Pipeline, RCP)) + 
  geom_tile(aes(fill = Number), colour = "black", alpha = FilllAlpha) + 
  geom_text(aes(label = paste0(round(Number), "\n(", round(Number.sd, 2), ")")),
            size = 3) + 
  labs(fill = "Encounters") +
  scale_fill_continuous_sequential(palette = AlberPalettes[[2]], rev = T) +
  #theme(legend.position = "top") +
  scale_x_discrete(labels = recode1 %>% 
                     str_replace_all("[+] ", "+\n") %>% str_remove(" Only")) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45, size = 8)) +
  coord_fixed() +
  
  # New sharings
  
  tablex %>% ggplot(aes(Pipeline, RCP)) + 
  geom_tile(aes(fill = New.Sharing), 
            colour = "black", alpha = FilllAlpha) + 
  geom_text(aes(label = paste0(round(New.Sharing), "\n(",round(New.Sharing.sd, 2), ")")),
            size = 3) + 
  labs(fill = "New\nsharings") +
  scale_fill_continuous_sequential(palette = AlberPalettes[[2]], rev = T) +
  #theme(legend.position = "top") +
  scale_x_discrete(labels = recode1 %>% 
                     str_replace_all("[+] ", "+\n") %>% str_remove(" Only")) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45, size = 8)) +
  coord_fixed() +
  
  # Overall deltas
  
  DeltasDF %>% 
  ggplot(aes(Pipeline, RCP, fill = Mean)) + 
  geom_tile(colour = "black", alpha = FilllAlpha) + 
  geom_text(aes(label = Label),
            size = 3) +
  labs(fill = "Change in\nsharing") +
  scale_fill_continuous_sequential(palette = AlberPalettes[[2]], rev = T) +
  #theme(legend.position = "top") +
  scale_x_discrete(labels = recode1 %>% 
                     str_replace_all("[+] ", "+\n") %>% str_remove(" Only")) +
  theme(axis.text.x = element_text(hjust = 1, angle = 45, size = 8)) +
  coord_fixed() + 
  plot_annotation(tag_levels = "A") +
  ggsave(filename = "Iceberg Extended Data/Heatmaps.jpg",
         units = "mm", 
         height = 250, width = 250, dpi = 600)

# Point and lines ####

# Process figure ####

FocalSp <- "Meles_meles"
FocalSp <- "Viverra_tangalunga"
FocalSp <- "Paguma_larvata"
FocalSp <- "Manis_pentadactyla"
FocalSp <- "Manis_javanica"
FocalSp <- "Arctonyx_collaris"
FocalSp <- "Taxidea_taxus"
FocalSp <- "Caracal_caracal"
FocalSp <- "Lycaon_pictus"
FocalSp <- "Felis_margarita"

GretCDF <- 
  readRDS(paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents/", 
                 FocalSp, ".rds")) %>% 
  as.matrix %>% as.data.frame() %>%
  full_join(readRDS(paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Futures/", FocalSp, ".rds")) %>% 
              as.matrix %>% as.data.frame() %>% dplyr::select(-c(Continent, BufferClimate)),
            by = c("X", "Y"), 
            suffix = c(".Currents", ".Futures"))

GretCDF %>% tidyr::gather("Step", "Value", -c(X, Y)) -> LongGretCDF

StepLevels <- c("IUCN", 
                "Climate", #"Continent", 
                "LandUse", "BufferClimate", #"LandUse.Futures1", 
                "Climate.Futures1", "BufferClimate.Futures1", 
                "ClimateLandUse.Futures1", "BufferClimateLandUse.Futures1")

names(StepLevels) <- c("Climate (IUCN Clipped)", "Land use", "Climate + dispersal buffer",
                       "Climate (RCP 2.6)", "Climate + dispersal",
                       "Climate + land use", "Climate + land use + dispersal")

FacetLabels <- c("Climate (IUCN Clipped)", "IUCN", "Land use", "Climate\n+ dispersal buffer",
                 "Climate (RCP 2.6)", "Climate + dispersal clipped",
                 "Climate + land use clipped", "Climate + land use\n+ dispersal clipped")

FacetLabels <- c("Climate (IUCN Clipped)", "IUCN", "Land use", "Climate + dispersal buffer",
                 "Climate (RCP 2.6)", "Climate + dispersal clipped",
                 "Climate + land use clipped", "Climate + land use + dispersal clipped")

names(FacetLabels) <- c("Climate", #"Continent", 
                        "IUCN", "LandUse", "BufferClimate", #"LandUse.Futures1", 
                        "Climate.Futures1", "BufferClimate.Futures1", 
                        "ClimateLandUse.Futures1", "BufferClimateLandUse.Futures1")

LongGretCDF %>% mutate(Step = factor(Step, levels = StepLevels)) %>% na.omit -> LongGretCDF

LongGretCDF %>% filter(!Step == "LandUse") %>% 
  filter(Value == 1) %>% 
  dplyr::select(X, Y) %>% ExtentGet() -> 
  Lims

LongGretCDF %>% 
  ggplot(aes(X, Y, fill = as.factor(Value))) + 
  geom_tile() + 
  facet_wrap(~Step, nrow = 5, labeller = as_labeller(FacetLabels)) + 
  scale_fill_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = c(2, 6)) +
  # coord_fixed() + 
  labs(fill = "Presence") + 
  #coord_cartesian(xlim = c(-40, 75), ylim = c(20, 80)) +
  coord_sf(xlim = c(Lims[[1]], Lims[[2]]), 
           ylim = c(Lims[[3]], Lims[[4]])) + 
  theme(strip.text = element_text(size = 10)) +
  theme_void() -> 
  
  BottomRow

Felis_margarita %>% data.frame %>% 
  mutate_at("PreClipClim", ~as.numeric(!is.na(.x))) %>% 
  ggplot(aes(X, Y, fill = as.factor(PreClipClim))) + 
  scale_fill_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = c(2, 6)) +
  coord_sf(xlim = c(Lims[[1]], Lims[[2]]), 
           ylim = c(Lims[[3]], Lims[[4]])) + 
  labs(fill = "Presence") +
  theme_void() +
  geom_tile() -> 
  
  TopRow

(TopRow/BottomRow) +
  plot_layout(heights = c(0.45, 1)) +
  ggsave(filename = paste0("Iceberg Extended Data/", 
                           "FocalSpecies_", FocalSp, ".jpg"), 
         width = 150, height = 250, 
         units = "mm", dpi = 600)

# GAMM panel ####

load(paste0("~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/","FitList.Rdata"))
load("~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/BAMList.Rdata")

Pipeline <- "A"

FitList[[Pipeline]] %>%
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
        legend.title = element_text(size = 10)) +
  geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Phylo), alpha = 0.01) +
  
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
        legend.title = element_text(size = 10)) +
  geom_rug(data = DataList[[1]], inherit.aes = F, aes(x = Space), alpha = 0.01) +
  
  FitList[[Pipeline]] %>%
  filter(!Phylo == last(unique(Phylo)),
         !Space == last(unique(Space))) %>%
  ggplot(aes(Space, Phylo)) +  theme_cowplot() +
  geom_tile(aes(fill = Fit)) +
  labs(x = "Geographic overlap",
       y = "Phylogenetic similarity",
       fill = "Viral sharing\nprobability") +
  lims(x = c(0,1), y = c(0,1)) +
  coord_fixed() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 10)) +
  geom_contour(aes(z = Fit), colour = "white", alpha = 0.8) +
  metR::geom_text_contour(aes(z = Fit), colour = "white", size = 2.5, hjust = 0.5, vjust = 1.1, check_overlap = T) +
  scale_fill_continuous_sequential(palette = "ag_GrnYl",
                                   limits = c(0,1),
                                   breaks = c(0,0.5,1)) +
  
  DataList[[Pipeline]] %>%
  ggplot(aes(Space, Phylo)) +  theme_cowplot() +
  labs(x = "Geographic overlap",
       y = "Phylogenetic similarity") +
  scale_fill_continuous_sequential(palette = "Heat 2", breaks = c(0:2*5)) +
  lims(x = c(0,1), y = c(0,1)) +
  coord_fixed() +
  theme(legend.position = "bottom") +
  geom_hex(aes(fill = stat(log(count)))) +
  
  ggsave(filename = paste0("Iceberg Extended Data/GAMMs.jpg"),
         units = "mm", height = 250, width = 250, dpi = 600)


# Degree predictions ####

FinalHostMatrix <- readRDS("Iceberg Output Files/KnownPredictedDegree.rds")

FinalHostMatrix %>% bind_rows(FinalHostMatrix %>% mutate(Sp = Sp2, Sp2 = Sp)) %>%
  group_by(Sp) %>% 
  summarise(Degree = sum(VirusBinary),
            PredDegree1 = sum(PredVirus1),
            PredDegree1b = sum(PredVirus1b),
            PredDegree1c = sum(PredVirus1c)) %>%
  tidyr::gather("Key", "Value", PredDegree1:PredDegree1c) ->
  
  HostsLong

HostsLong %>%
  ggplot(aes(Degree, Value, colour = Sp)) + 
  facet_wrap(~Key, 
             labeller = labeller(Key = c("PredDegree1" = "All effects",
                                         "PredDegree1b" = "Fixed effects",
                                         "PredDegree1c" = "Random effects"))) + 
  geom_abline(lty = 2, alpha = 0.3) +
  geom_point(alpha = 0.5) +
  coord_fixed() +
  theme(legend.position = "none", strip.background = element_rect(fill = "white")) +
  labs(x = "Observed degree", y = "Predicted degree") +
  #lims(x = c(0, MaxLim), y = c(0, MaxLim)) +
  scale_colour_discrete_sequential(palette = AlberPalettes[2])+
  ggsave("Iceberg Extended Data/DegreePredictions.jpg", units = "mm", width = 200, height = 100)


# New one ####

which(CoryClimateReps %in% c("ca", "mi")) %>% ClimateReps[.] %>% 
  map(~glue("~/Albersnet/Iceberg Files/{.x}/Iceberg Output Files/NewIntersects.rds") %>% readRDS) ->
  NIList

NIList %>% 
  bind_rows(.id = "ClimateReps") %>% 
  mutate_at("ClimateReps", ~c("CanESM5", "MIROC6")[as.numeric(.x)]) %>% 
  dplyr::select(X, Y, starts_with("Overlap.Futures"), ClimateReps) %>%
  dplyr::select(-ends_with("B"), -ends_with("D")) %>%
  tidyr::gather("Key", "Encounters", -c(X,Y, ClimateReps)) %>%
  mutate(Key = Key %>% str_remove("Overlap.")) %>%
  tidyr::separate(Key, sep = 8, into = c("RCP", "Pipeline")) %>%
  mutate(Pipeline = recode1[Pipeline],
         RCP = recode2[RCP]) %>%
  filter(str_detect(Pipeline, "Dispersal")) %>% 
  ggplot(., aes(X, Y, fill = log10(Encounters + 1))) + geom_tile() + 
  scale_fill_continuous_sequential(palette = AlberPalettes[[1]],
                                   labels = c(1, 11, 101, 1001, 10001)) +
  facet_grid(RCP~ClimateReps) + 
  coord_fixed() +
  theme_void() + labs(fill = "Encounters") +
  theme(legend.position = "bottom", legend.text = element_text(angle = 45, hjust = 1)) +
  ggsave(filename = paste0("Iceberg Extended Data/", "NE_TwoClimates", ".jpg"), width = 10, height = 8)


