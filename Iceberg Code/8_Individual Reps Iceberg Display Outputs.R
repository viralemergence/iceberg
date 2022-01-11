
# 8_Display Output #####

{
  
  library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr); library(SpRanger); library(cowplot);
  library(RColorBrewer); library(fs); library(colorspace)
  
  #devtools::install_github("gfalbery/ggregplot")
  library(ggregplot); library(ggplot2); library(RColorBrewer)
  
  ParasitePalettes<-c("PuRd","PuBu","BuGn","Purples","Oranges")
  ParasiteColours<-c("#DD1c77","#2B8CBE","#2CA25F",brewer.pal(5,"Purples")[4],brewer.pal(5,"Oranges")[4])
  
  AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
  AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
  AlberColours[length(AlberColours)+1:2] <- RColorBrewer::brewer.pal(11, AlberPalettes[[4]])[c(2,10)]
  
}

# setwd(paste0("~/Albersnet/Iceberg Files/", "Summary"))

dir_create("Iceberg Figures")
dir_create("Iceberg Figures/RasterDrafts")

LightSpectral <- rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))[10:100]

DarkSpectral <- rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))

PredReps <- c("Currents", paste0("Futures", 1:4))

PipelineReps <- LETTERS[1:4]

SpaceVars <- paste0(paste("Space", PredReps, sep = "."),rep(PipelineReps, each = length(PredReps)))
SharingVars <- paste0(paste("Sharing",PredReps, sep = "."), rep(PipelineReps, each = length(PredReps)))

names(SpaceVars) <- names(SharingVars) <- paste0(PredReps,rep(PipelineReps, each = length(PredReps)))

# Overall trends ###

AllMammaldf <- readRDS("Iceberg Output Files/AllMammaldf.rds")

# AllMammaldf %>% summarise_all(mean) %>% View()
# AllMammaldf %>% summarise_if(is.numeric, sum) %>% View()

# Import #####

CurrentsGridDF <- readRDS("~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/CurrentsGridDF.rds")
FuturesGridDF <- readRDS("Iceberg Output Files/FuturesGridDF.rds")

CurrentsGridDF %>%
  bind_cols(FuturesGridDF[,setdiff(names(FuturesGridDF), names(CurrentsGridDF))]) %>%
  rename(CurrentsAC = Climate, CurrentsBD = ClimateLandUse,
         Sharing.CurrentsAC = Sharing.CurrentsA, Sharing.CurrentsBD = Sharing.CurrentsB) %>%
  rename(Futures1D = Climate.Futures1, Futures2D = Climate.Futures2,
         Futures3D = Climate.Futures3, Futures4D = Climate.Futures4) %>%
  rename(Futures1C = ClimateLandUse.Futures1, Futures2C = ClimateLandUse.Futures2,
         Futures3C = ClimateLandUse.Futures3, Futures4C = ClimateLandUse.Futures4) %>%
  rename(Futures1B = BufferClimate.Futures1, Futures2B = BufferClimate.Futures2,
         Futures3B = BufferClimate.Futures3, Futures4B = BufferClimate.Futures4) %>%
  rename(Futures1A = BufferClimateLandUse.Futures1, Futures2A = BufferClimateLandUse.Futures2,
         Futures3A = BufferClimateLandUse.Futures3, Futures4A = BufferClimateLandUse.Futures4) %>%
  dplyr::select(-starts_with("LandUse"),X,Y, starts_with("Currents"),
                ends_with("A"), ends_with("B"), ends_with("C"), ends_with("D")) ->
  
  GretCDF

NewIntersects <- readRDS("Iceberg Output Files/NewIntersects.rds")
NewEncountersList <- readRDS(paste0("Iceberg Output Files/NewEncounters.rds"))

lapply(PipelineReps, function(a){
  
  lapply(PredReps[2:5], function(b){
    
    NewEncountersList[[a]][[b]]$SharingVar <- NewEncountersList[[a]][[b]][,SharingVars[paste0(b,a)]]
    NewEncountersList[[a]][[b]]$DeltaSharingVar <- NewEncountersList[[a]][[b]][,paste0("Delta",SharingVars[paste0(b,a)])]
    
    NewEncountersList[[a]][[b]] %>% dplyr::summarise(
      Number = n(),
      Mean.Sharing = mean(SharingVar),
      Sum.Sharing = sum(SharingVar),
      New.Sharing = sum(DeltaSharingVar)
    )
    
  }) %>% bind_rows(.id = "PredRep")
  
}) %>% bind_rows(.id = "Pipeline") %>%
  mutate(Pipeline = PipelineReps[as.numeric(Pipeline)],
         PredRep = PredReps[as.numeric(PredRep)+1]) -> tablex
tablex

recode1 <- c(A = "Climate + Land + Dispersal",
             B = "Climate + Land",
             C = "Climate + Dispersal",
             D = "Climate Only")

recode2 <- c(Futures1 = 'RCP 2.6',
             Futures2 = 'RCP 4.5',
             Futures3 = 'RCP 6.0',
             Futures4 = 'RCP 8.5')

tablex %>% as_tibble() %>% 
  mutate(Pipeline = recode(Pipeline, !!!recode1)) %>%
  mutate(PredRep = recode(PredRep, !!!recode2)) %>%
  group_by(PredRep) %>% 
  ggplot(aes(x = PredRep, y = Number, group = 1)) + 
  geom_line() + geom_point() + facet_wrap(~ Pipeline, scales = 'free') + xlab("Climate pathway") + 
  ylab("Number of new encounters")

####

nrow(NewEncountersList$A$Futures1[NewEncountersList$A$Futures1$hOrder.x == 'Chiroptera' |
                                    NewEncountersList$A$Futures1$hOrder.y == 'Chiroptera',])

nrow(NewEncountersList$A$Futures1)

#Ebola: ####
AssocsBase <- read_csv("https://raw.githubusercontent.com/ecohealthalliance/HP3/master/data/associations.csv") %>% data.frame()
HP3EbolaHosts <- AssocsBase %>% filter(vVirusNameCorrected == "Zaire_ebolavirus")
LauraEbolaHosts <- read.csv("Iceberg Input Files/ZEBOV hosts.csv")

EbolaHosts <- union(HP3EbolaHosts$Host, LauraEbolaHosts$bat_species) %>% 
  intersect(list.files("Iceberg Input Files/GretCDF/Currents") %>% str_remove(".rds$"))

EbolaHosts %>% setdiff(c("Miniopterus_schreibersii","Pipistrellus_pipistrellus", "Hipposideros_pomona",
                         "Cynopterus_sphinx","Acerodon_jubatus", "Rousettus_leschenaultii")) ->
  
  EbolaHosts

#  Current richness

EbolaGridDF <- readRDS("Iceberg Output Files/EbolaGridList.rds")

#New encounters (Ebola-nonEbola)
#Total new encounters & estimated sharing events A-D x 1-4

EbolaEncounters <- readRDS("Iceberg Output Files/EbolaEncounters.rds")

NonFigures <- F

if(NonFigures){
  
  names(EbolaEncounters) %>% lapply(function(a){
    
    list(
      Number = length(EbolaEncounters[[a]][,SharingVars[a]]),
      Mean.Sharing = mean(EbolaEncounters[[a]][,SharingVars[a]]),
      Sum.Sharing = sum(EbolaEncounters[[a]][,SharingVars[a]]),
      New.Sharing = sum(EbolaEncounters[[a]][,paste0("Delta",SharingVars[a])])
      
    ) %>% return
  }) %>% bind_rows()
  
  # Current probability of sharing matrix among _ONLY KNOWN and AFRICAN_ Ebola hosts
  
  AllMammaldf %>% filter(Sp%in%EbolaHosts&Sp2%in%EbolaHosts) %>% summarise_if(is.numeric, mean)
  AllMammaldf %>% filter(Sp%in%EbolaHosts&Sp2%in%EbolaHosts) %>% summarise_if(is.numeric, sum)
  
  # New encounters and sharing probability matrix between Ebola-nonEbola, including metadata on host order, in 2.6 A
  
  EbolaEncounters # Includes only one known African Ebola host and one other mammal, and possible data, for all 16 scenarios.
  
  # For example, all hosts in new encounters
  EbolaEncounters[[1]] %>% gather("Key", "Value", hOrder.x, hOrder.y) %>%
    rename(hOrder = Value) %>%
    group_by(hOrder) %>% summarise(N = n()) %>% arrange(N) %>%
    as.data.frame() %>% mutate(hOrder = factor(hOrder, levels = unique(hOrder))) %>%
    ggplot(aes(hOrder, N, fill = hOrder)) + geom_col(colour = "black")
  
  # Or if EXCLUDING KNOWN EBOLA HOSTS
  EbolaEncounters[[1]] %>% gather("Key", "Value", Sp, Sp2) %>% filter(!Value%in%EbolaHosts) %>%
    mutate(hOrder = ifelse(Key == "Sp",
                           as.character(hOrder.x),
                           as.character(hOrder.y))) %>%
    group_by(hOrder) %>% summarise(N = n()) %>% arrange(N) %>%
    as.data.frame() %>% mutate(hOrder = factor(hOrder, levels = unique(hOrder))) %>%
    ggplot(aes(hOrder, N, fill = hOrder)) + geom_col(colour = "black")
  
  # Then if looking at CHANGE IN SHARING of these hosts
  EbolaEncounters[[1]] %>% gather("Key", "Value", Sp, Sp2) %>% filter(!Value%in%EbolaHosts) %>%
    mutate(hOrder = ifelse(Key == "Sp",
                           as.character(hOrder.x),
                           as.character(hOrder.y))) %>%
    group_by(hOrder) %>% summarise(SharingChange = sum(DeltaSharing.Futures1A)) %>% arrange(SharingChange) %>%
    as.data.frame() %>% mutate(hOrder = factor(hOrder, levels = unique(hOrder))) %>%
    ggplot(aes(hOrder, SharingChange, fill = hOrder)) + geom_col(colour = "black")
  
  # Then if looking at CHANGE IN SHARING of these hosts
  EbolaEncounters[[4]] %>% gather("Key", "Value", Sp, Sp2) %>% filter(!Value%in%EbolaHosts) %>%
    mutate(hOrder = ifelse(Key == "Sp",
                           as.character(hOrder.x),
                           as.character(hOrder.y))) %>%
    group_by(hOrder) %>% summarise(SharingChange = sum(DeltaSharing.Futures4A)) %>% arrange(SharingChange) %>%
    as.data.frame() %>% mutate(hOrder = factor(hOrder, levels = unique(hOrder))) %>%
    ggplot(aes(hOrder, SharingChange, fill = hOrder)) + geom_col(colour = "black")
  
  #Tallies:
  #  New encounters by scenario A-D x 1-4
  
  EbolaEncounters %>% sapply(nrow) %>% 
    matrix(ncol = 4) %>% as.data.frame() ->
    EbolaNETable
  
  rownames(EbolaNETable) <- PredReps[2:length(PredReps)]
  colnames(EbolaNETable) <- PipelineReps
  
}

# Figure Prep ####

library(raster)

blank <- raster('~/Albersnet/Iceberg Files/CHELSA/UniversalBlank.tif')
Sea = which(is.na(raster::values(blank)))

rast <- function(x, Projection = NULL) {
  
  r <- blank
  
  values(r)[Sea] <- NA
  values(r)[-Sea] <- x
  
  if(!is.null(Projection)){
    
    r@crs <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
    
    r <- projectRaster(r, crs = Projection)
    
  }
  
  r
  
}

plot(rast(NewIntersects$Overlap.Futures1A))
plot(rast(NewIntersects$Overlap.Futures1D))

writeRaster(rast(NewIntersects$Overlap.Futures1A), paste0('Iceberg Output Files/RasterDrafts/',
                                                          'overlap1A.tif'),
            overwrite = T)

writeRaster(rast(NewIntersects$Overlap.Futures1D), paste0('Iceberg Output Files/RasterDrafts/',
                                                          'overlap1D.tif'),
            overwrite = T)

mean1D <- NewIntersects$OverlapSharing.Futures1D / NewIntersects$Overlap.Futures1D
mean1D[NewIntersects$Overlap.Futures1D<500] <- NA
plot(rast(mean1D))

# Grigures #### 

EbolaGridList <- readRDS('Iceberg Output Files/EbolaGridList.rds')

EboCurrents <- stack(lapply(EbolaGridList, function(x) {rast(x$ClimateLandUse)}))

conts <- raster('Iceberg Input Files/continents-madagascar.tif')
africa <- conts; africa[!(africa == 1)] <- NA; africa <- africa-1

EboFutures <- stack(lapply(EbolaGridList, function(x) {rast(x$BufferClimateLandUse.Futures1)}))

EbolaNew <- readRDS('Iceberg Output Files/EbolaNewIntersects.rds')

africashp <- rasterToPolygons(africa, dissolve = TRUE)
contsshp <- rasterToPolygons(conts, dissolve = TRUE)

# Figure 1 ####

Eckert <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

marker = c(color = colorRampPalette(brewer.pal(11,"Spectral"))(100))

pdf("Iceberg Figures/Figure1.pdf", width = 12, height = 10)

par(mfrow = c(2, 1), mar = c(0,0.5,0, 5))

plot(rast(NewIntersects$DeltaOverlapSharing.Futures1C),
     #legend.shrink = 0.8, 
     legend.width = 2, 
     box = FALSE,
     axes = FALSE,
     col = DarkSpectral,
     xlim = c(-12553230, 16774960), ylim = c(-8460601, 8373855),
     legend.args=list(text='Viral sharing events', side=2, line=0.4, cex=1.5))

plot(contsshp, lwd = 0.5, add = TRUE)

plot(rast(NewIntersects$DeltaOverlapSharing.Futures1A),
     #legend.shrink = 0.8, 
     legend.width = 2, 
     box = FALSE,
     axes = FALSE,
     col = DarkSpectral,
     xlim = c(-12553230, 16774960), ylim = c(-8460601, 8373855),
     legend.args=list(text='Viral sharing events', side=2, line=0.4, cex=1.5))

plot(contsshp, lwd = 0.5, add = TRUE)

dev.off()

# Figure 2 Prep ####

load("Iceberg Output Files/NEGams.Rdata")

unlist(NEPredictModels, recursive = F) %>% unlist(recursive = F) %>%
  
  lapply(function(a){
    
    a <- summary(a)
    N <- a$p.coeff %>% length
    
    data.frame(Mean = a$p.coeff, 
               SE = a$se[1:N],
               Var = names(a$p.coeff) %>% as.character() %>% str_remove("HabitatType")) %>%
      rbind(data.frame(Mean = 0, SE = 0, Var = "Cropland"), .) %>%
      filter(!Var == "(Intercept)")
    
  }) %>% bind_rows(.id = "Model") -> BAMFixed

FitList %>% unlist(recursive = F)  %>% unlist(recursive = F) %>% bind_rows() %>% 
  group_by(Rep) %>% summarise(Rich = last(Richness)) %>% 
  pull(Rich) %>% unique() -> Richnesses

FitList %>% unlist(recursive = F) %>%  map("NoBat") %>% bind_rows(.id = "Pipeline") %>%
  filter(HabitatType == "Cropland", Richness %in% Richnesses) %>% mutate(Bats = 'Non-bats') -> df1

FitList %>% unlist(recursive = F) %>%  map("SomeBat") %>% bind_rows(.id = "Pipeline") %>%
  filter(HabitatType == "Cropland", Richness %in% Richnesses)  %>% mutate(Bats = "Bat encounters") -> df2

named <- c(A.Futures1 = 'RCP 2.6 (CLD)',
           A.Futures4 = 'RCP 8.5 (CLD)',
           C.Futures1 = 'RCP 2.6 (CL)',
           C.Futures4 = 'RCP 8.5 (CL)') 

rbind(df1, df2) %>% mutate(Scenario = recode(Pipeline, !!!named)) %>% 
  ggplot(aes(Elevation, Fit, fill = Scenario, colour = Scenario)) + facet_wrap( ~ Bats) + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, colour = NA) + geom_line() +
  scale_fill_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = c(3,5,7,8)) + 
  labs(y = "Predicted Encounters",
       x = "Elevation",
       fill = 'Scenario', colour = 'Scenario') +
  theme_cowplot()  + 
  scale_colour_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = c(3,5,7,8)) + 
  #scale_y_log10(labels = scales::scientific) + #xlim(0,6100) + 
  scale_y_continuous(breaks = c((-1):6), labels = c(10^((-1):6))) +
  scale_x_continuous(breaks = c((-1):6), labels = c(10^((-1):6))) +
  theme(strip.background = element_blank()) -> g1; g1

FitList %>% unlist(recursive = F) %>%  map('NoBat') %>% bind_rows(.id = "Pipeline") %>%
  filter(HabitatType == "Cropland", Elevation == last(unique(Elevation))) %>% mutate(Bats = 'Non-bats')  -> df3

FitList %>% unlist(recursive = F) %>%  map('SomeBat') %>% bind_rows(.id = "Pipeline") %>%
  filter(HabitatType == "Cropland", Elevation == last(unique(Elevation))) %>% mutate(Bats = "Bat encounters") -> df4

named <- c(A.Futures1 = 'RCP 2.6 (CLD)',
           A.Futures4 = 'RCP 8.5 (CLD)',
           C.Futures1 = 'RCP 2.6 (CL)',
           C.Futures4 = 'RCP 8.5 (CL)') 

rbind(df3, df4) %>% mutate(Scenario = recode(Pipeline, !!!named)) %>% 
  ggplot(aes(Richness, Fit, fill = Scenario, colour = Scenario)) + facet_wrap(~ Bats) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, colour = NA) + geom_line() +
  labs(y = "Predicted Encounters",
       x = "Richness",
       fill = 'Scenario', colour = 'Scenario')+
  theme_cowplot()  +  
  scale_fill_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = c(3,5,7,8)) + 
  scale_colour_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = c(3,5,7,8)) + 
  #scale_y_log10(labels = scales::scientific) + 
  scale_y_continuous(breaks = c(seq(-1, 8, by = 2)), labels = c(10^(seq(-1, 8, by = 2)))) +
  scale_x_continuous(breaks = c(seq(-1, 8, by = 1)), labels = c(10^(seq(-1, 8, by = 1)))) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) -> g2; g2

plot_grid(g1, g2, nrow = 2)

#### BAM panels

named <- c(A.Futures1.NoBat = 'RCP 2.6 (CLD)',
           A.Futures4.NoBat = 'RCP 8.5 (CLD)',
           C.Futures1.NoBat = 'RCP 2.6 (CL)',
           C.Futures4.NoBat = 'RCP 8.5 (CL)',
           A.Futures1.SomeBat = 'RCP 2.6 (CLD)',
           A.Futures4.SomeBat = 'RCP 8.5 (CLD)',
           C.Futures1.SomeBat = 'RCP 2.6 (CL)',
           C.Futures4.SomeBat = 'RCP 8.5 (CL)') 

reclass <- c(Rangeland = 'Range',
             Nonforest = 'Other',
             Cropland = 'Crop', 
             Urban = 'Settled') 

BAMFixed %>% 
  filter(!Var %in% c("Elevation", "Richness")) %>%
  mutate(Bats = grepl("SomeBat",Model)) %>% 
  mutate(Bats = if_else(Bats == TRUE,'Bat encounters','Non-bat')) %>%
  mutate(Scenario = recode(Model, !!!named)) %>% 
  mutate(Var = recode(Var, !!!reclass)) %>%
  ggplot(aes(Var, Mean, colour = Scenario)) + theme_classic() + 
  facet_wrap(. ~ Bats) + 
  geom_hline(yintercept = 0, alpha = 0.7, lty = 2) +
  geom_point(position = position_dodge(w = 0.5)) + 
  geom_errorbar(aes(ymin = Mean - SE, 
                    ymax = Mean + SE), 
                position = position_dodge(w = 0.5), 
                width = 0.2) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(margin = margin(r = 5)),
        axis.title.y = element_text(margin = margin(r = 6, l = 8)),
        plot.margin = margin(0.5, 0.1, 0.5, 0, 'cm'),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.border = element_rect(colour = "black",
                                    fill = NA)) + 
  labs(y = "Effect size", x = "Land Use") + 
  scale_color_discrete_sequential(palette = AlberPalettes[[3]], nmax = 10, order = c(5,7,8,10)) +
  lims(x = rev(c("Crop", "Settled", "Range", "Forest", "Other"))) + 
  coord_flip()  -> g3; g3

TextSize = 14

TextTheme <- theme(axis.title.x = element_text(size = TextSize),
                   axis.title.y = element_text(size = TextSize),
                   axis.text.y = element_text(angle = 0, hjust = 1))

plot_grid(g1 + TextTheme,
          g2 + TextTheme,
          g3 + TextTheme,
          nrow = 3, rel_heights = c(1.15,1,1)) %>%
  save_plot(nrow = 3, ncol = 2, filename = "Iceberg Figures/Figure2CDEFG.jpeg",
            base_width = 4, base_height = 2.5)

# Adding bat and nonbat encounters ####

# Figure 2 Plotting ####


pdf("Iceberg Figures/Figure2AB.pdf", width = 12, height = 10)

par(mfrow = c(2, 1), mar = c(0,0.5,0, 5))

plot(rast(NewIntersects$Overlap.Futures1A - NewIntersects$NoBatOverlap.Futures1A), axes = FALSE, box = FALSE,
     legend.shrink = 0.8, 
     legend.width = 2, 
     axis.args = list(cex.axis = 1.2),
     col = DarkSpectral,
     legend.args=list(text='First encounters', side=2, line=0.4, cex=1.5))

plot(contsshp, lwd = 0.5, add = TRUE)

plot(rast(NewIntersects$NoBatOverlap.Futures1A), axes = FALSE, box = FALSE, legend.shrink = 0.8, 
     #legend.shrink = 0.8, 
     legend.width = 2, 
     axis.args = list(cex.axis = 1.2),
     col = DarkSpectral,
     legend.args=list(text='First encounters', side=2, line=0.4, cex=1.5))

plot(contsshp, lwd = 0.5, add = TRUE)

dev.off()

# Figure 3BCD ####

pdf("Iceberg Figures/Figure3BCD.pdf", width = 12, height = 5)

par(mfrow = c(1, 3), mar = c(0, 2, 0.2, 6), oma=c(1,1,1,2), cex=1.2)

plot(trim(sum(EboCurrents) + africa), xlim = c(-20, 55), ylim = c(-40,40), axes = FALSE, box = FALSE,
     col = DarkSpectral,
     legend.shrink = 0.5, 
     legend.width = 4,
     legend.args=list(text='Species', side=2, line=0.4, cex=1.5) ) # PANEL C

plot(africashp, add = TRUE, lwd = 0.5)

par(cex=1.2)
plot(trim(sum(EboFutures)-sum(EboCurrents) + africa), xlim = c(-20, 55), ylim = c(-40,40), axes = FALSE, box = FALSE,
     col = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)),
     legend.shrink = 0.5, 
     legend.width = 4, 
     legend.args=list(text=expression(paste(Delta, 'Species')), side=2, line=0.4, cex=1.5)) # PANEL C

plot(africashp, add = TRUE, lwd = 0.5)

par(cex=1.2)
plot(trim(rast(EbolaNew$Overlap.Futures1A)+africa), xlim = c(-20, 55), ylim = c(-40,40), axes = FALSE, box = FALSE,
     col = DarkSpectral,
     legend.shrink = 0.5, 
     legend.width = 4, 
     legend.args=list(text='First encounters', side=2, line=0.4, cex=1.5))

plot(africashp, add = TRUE, lwd = 0.5)

dev.off()

# Figure 3A ####

BatNew <- readRDS('Iceberg Output Files/BPNewIntersects.rds')

pdf("Iceberg Figures/Figure3A.pdf", width = 12, height = 10)

par(mfrow = c(1,1), mar=c(4,4,4,4.3))

plot(rast(BatNew$OverlapSharing.Futures1A), axes = FALSE, box = FALSE,
     col = DarkSpectral,
     legend.shrink = 0.4,
     legend.width = 1.5, 
     legend.args=list(text='Viral sharing events', side=2, line=0.3, cex=1.5))
plot(contsshp, lwd = 0.5,add = TRUE)

dev.off()

# Figure 3E ####

writeRaster(rast(NewIntersects$DeltaOverlapSharing.Futures1A), 
            file = "Iceberg Figures/Bivariate1A.tif", 
            overwrite = T)

if(0){
  
  # Extra stuff?
  
  par(fig = c(0,0.333,0,0.5), mar = c(0,0,0,4.5), new = TRUE)
  plot(trim(sum(EboCurrents) + africa), xlim = c(-20, 55), ylim = c(-40,40), axes = FALSE, box = FALSE,
       col = rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))[15:100]) # PANEL C
  par(fig = c(0,0.333,0,0.5), new = TRUE)
  plot(africashp, add = TRUE, lwd = 0.5)
  
  par(fig = c(0.333,0.666,0,0.5), mar = c(0,0,0,4.5), new = TRUE)
  plot(trim(sum(EboFutures)-sum(EboCurrents) + africa), xlim = c(-20, 55), ylim = c(-40,40), axes = FALSE, box = FALSE,
       col = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100))) # PANEL C
  par(fig = c(0.333,0.666,0,0.5), new = TRUE)
  plot(africashp, add = TRUE, lwd = 0.5)
  
  par(fig = c(0.666,1,0,0.5), mar = c(0,0,0,4.5), new = TRUE)
  plot(trim(rast(EbolaNew$Overlap.Futures1A)+africa), xlim = c(-20, 55), ylim = c(-40,40), axes = FALSE, box = FALSE,
       col = rev(colorRampPalette(brewer.pal(11,"Spectral"))(100)))
  par(fig = c(0.666,1,0,0.5), new = TRUE)
  plot(africashp, add = TRUE, lwd = 0.5)
  
  #### Sub-components
  
  
  plot.new()
  
  par(mfrow = c(1,2), mar = c(0,1,0,4))
  plot(rast(BatNew$Overlap.Futures1A), axes = FALSE, box = FALSE,
       col = rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))[10:100],
       legend.shrink = 0.25)
  plot(contsshp, lwd = 0.5, add = TRUE)
  
  plot(rast(BatNew$DeltaOverlapSharing.Futures1A), axes = FALSE, box = FALSE,
       col = rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))[10:100],
       legend.shrink = 0.25)
  plot(contsshp, lwd = 0.5,add = TRUE)
  
  plot.new()
  
  par(mfrow = c(1,3), mar = c(0,0.5,0,5.5))
  
  plot(trim(sum(EboCurrents) + africa), xlim = c(-20, 55), ylim = c(-40,40), axes = FALSE, box = FALSE,
       col = rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))[15:100],
       axis.args = list( cex.axis = 1.5), legend.shrink = 0.25) # PANEL C
  plot(africashp, add = TRUE, lwd = 0.5)
  
  plot(trim(sum(EboFutures)-sum(EboCurrents) + africa), xlim = c(-20, 55), ylim = c(-40,40), axes = FALSE, box = FALSE,
       col = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)),
       axis.args = list( cex.axis = 1.5), legend.shrink = 0.25) # PANEL C
  plot(africashp, add = TRUE, lwd = 0.5)
  
  plot(trim(rast(EbolaNew$Overlap.Futures1A)+africa), xlim = c(-20, 55), ylim = c(-40,40), axes = FALSE, box = FALSE,
       col = rev(colorRampPalette(brewer.pal(11,"Spectral"))(100)),
       axis.args = list( cex.axis = 1.5), legend.shrink = 0.25)
  plot(africashp, add = TRUE, lwd = 0.5)
  
}

ClimateReps %>% lapply(function(.x){
  
  Folder <- glue::glue("~/Albersnet/Iceberg Files/{.x}/Iceberg Extended Data")
  
  dir_copy(Folder, glue::glue("~/Albersnet/Iceberg Supplement/{.x}"))
  
})
