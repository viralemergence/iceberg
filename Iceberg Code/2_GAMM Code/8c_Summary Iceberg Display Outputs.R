
# 8_Display Output #####

{
  
  library(tidyverse); library(raster); library(parallel); library(sf); library(Matrix); library(magrittr); library(SpRanger); library(cowplot);
  library(RColorBrewer); library(fs); library(colorspace); library(glue); library(mgcv)
  library(patchwork)
  
  #devtools::install_github("gfalbery/ggregplot")
  library(ggregplot); library(ggplot2); library(RColorBrewer)
  
  ParasitePalettes<-c("PuRd","PuBu","BuGn","Purples","Oranges")
  ParasiteColours<-c("#DD1c77","#2B8CBE","#2CA25F",brewer.pal(5,"Purples")[4],brewer.pal(5,"Oranges")[4])
  
  AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
  AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
  AlberColours[length(AlberColours)+1:2] <- RColorBrewer::brewer.pal(11, AlberPalettes[[4]])[c(2,10)]
  
}

setwd(paste0("~/Albersnet/Iceberg Files/", "Summary"))

dir_create("Iceberg Figures")
dir_create("Iceberg Figures/RasterDrafts")

LightSpectral <- rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))[10:100]

DarkSpectral <- rev(colorRampPalette(brewer.pal(11,"Spectral"))(100))

PredReps <- c("Currents", paste0("Futures", 1:4))

PipelineReps <- LETTERS[1:4]

SpaceVars <- paste0(paste("Space", PredReps, sep = "."),rep(PipelineReps, each = length(PredReps)))
SharingVars <- paste0(paste("Sharing",PredReps, sep = "."), rep(PipelineReps, each = length(PredReps)))

names(SpaceVars) <- names(SharingVars) <- paste0(PredReps,rep(PipelineReps, each = length(PredReps)))

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

NewIntersects <- readRDS("Iceberg Output Files/NewIntersects.rds") %>% as.data.frame

NewIntersects %>% arrange(desc(Y), X) -> NewIntersects

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

EbolaGridDF <- readRDS("~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/EbolaGridList.rds")

# Setting up ####

Unlist1 <- function(a) unlist(a, recursive = F)

MeanSD <- function(a) list(Mean = mean(a, na.rm = T), 
                           SD = sd(a, na.rm = T))

glue::glue("Climate{1:9}") ->
  
  ClimateReps

# Importing New Encounters ####

glue("~/Albersnet/Iceberg Files/{ClimateReps}/Iceberg Output Files/NewEncounters.rds") %>% 
  map(readRDS) -> NewEncountersListList

names(NewEncountersListList) <- ClimateReps

NewEncountersListList %<>% 
  map(Unlist1) %>% 
  Unlist1

NewEncountersListList %>% names 

# Importing Ebola ####

glue("~/Albersnet/Iceberg Files/{ClimateReps}/Iceberg Output Files/EbolaEncounters.rds") %>% 
  map(readRDS) -> 
  EbolaEncountersListList

names(EbolaEncountersListList) <- ClimateReps

EbolaEncountersListList %<>% 
  #map(Unlist1) %>% 
  Unlist1

EbolaEncountersListList %>% names 

# Importing BP ####

glue("~/Albersnet/Iceberg Files/{ClimateReps}/Iceberg Output Files/BPEncounters.rds") %>% 
  map(readRDS) -> 
  BPEncountersListList

names(BPEncountersListList) <- ClimateReps

BPEncountersListList %<>% 
  #map(Unlist1) %>% 
  Unlist1

BPEncountersListList %>% names 

# Loading All All Mammal DF

#AllAllMammalDF <- 
#  readRDS("~/Albersnet/Iceberg Files/Summary/Iceberg Output Files/AllAllMammalDF.rds")

SummaryNetworkChanges <- 
  readRDS("~/Albersnet/Iceberg Files/Summary/Iceberg Output Files/SummaryNetworkChanges.rds")

# Figure Prep ####

blank <- matrix(0,360*2,720*2) # proper resolution
blank <- raster(blank)
extent(blank) <- c(-180, 180, -90, 90)
projection(blank) <- CRS("+proj=longlat +datum=WGS84")

UniversalBlank <- raster("Iceberg Input Files/UniversalBlank.tif")
Land = which(raster::values(UniversalBlank)==0)
Sea = which(is.na(raster::values(UniversalBlank)))

RasterCols <- 
  seq(from = extent(blank)[1], to = extent(blank)[2],
      length.out = ncol(blank)) %>%
  as.character()

RasterRows <- 
  seq(from = extent(blank)[3], to = extent(blank)[4],
      length.out = nrow(blank)) %>%
  as.character()

GregRast <- function(DF, Value, Projection = NULL){
  
  DF[,c("X", "Y", Value)] %>% 
    pivot_wider(names_from = "X", 
                values_from = Value) %>% 
    arrange(desc(Y)) %>% 
    as.data.frame -> 
    
    Matrix
  
  rownames(Matrix) <- Matrix$Y %>% as.character()
  
  Matrix[setdiff(RasterRows, rownames(Matrix)),] <- NA
  Matrix[,setdiff(RasterCols, colnames(Matrix))] <- NA
  
  Matrix %>% dplyr::select(-Y) %>% as.matrix -> Matrix2
  
  Matrix2[RasterRows, RasterCols] -> Matrix3
  
  Matrix3 %>% raster -> R1
  
  extent(R1) <- c(min(DF$X, na.rm = T), max(DF$X, na.rm = T),
                  min(DF$Y, na.rm = T), max(DF$Y, na.rm = T))
  
  R1 %>% resample(blank) -> R1
  
  return(R1)
  
}

GregRast <- function(DF, x, Projection = NULL) {
  
  r <- UniversalBlank
  
  values(r)[Sea] <- NA
  values(r)[-Sea] <- DF[,x]
  
  r
  
}

writeRaster(GregRast(NewIntersects, "Overlap.Futures1A"), paste0('Iceberg Output Files/RasterDrafts/',
                                                                 'overlap1A.tif'),
            overwrite = T)

writeRaster(GregRast(NewIntersects, "Overlap.Futures1D"), paste0('Iceberg Output Files/RasterDrafts/',
                                                                 'overlap1D.tif'),
            overwrite = T)

# Grigures #### 

EbolaGridList <- 
  readRDS('~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/EbolaGridList.rds')

EbolaGridList %>% 
  map(~GregRast(.x, "ClimateLandUse")) %>% stack ->
  
  EboCurrents

EbolaGridList %>% map(~GregRast(.x, "BufferClimateLandUse.Futures1")) %>% stack ->
  
  EboFutures

EbolaNew <- readRDS('Iceberg Output Files/EbolaNewIntersects.rds')

EbolaNew %>% arrange(desc(Y), X) %>% as.data.frame -> EbolaNew

conts <- raster('Iceberg Input Files/continents-madagascar.tif')
africa <- conts; africa[!(africa == 1)] <- NA; africa <- africa-1

africashp <- rasterToPolygons(africa, dissolve = TRUE)
contsshp <- rasterToPolygons(conts, dissolve = TRUE)

# Figure 1 ####

Eckert <- "+proj=eck4 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

marker = c(color = colorRampPalette(brewer.pal(11,"Spectral"))(100))

pdf("Iceberg Figures/Figure1.pdf", width = 12, height = 10)

par(mfrow = c(2, 1), mar = c(0,0.5,0, 5))

plot(GregRast(NewIntersects, "DeltaOverlapSharing.Futures1C"),
     #legend.shrink = 0.8, 
     legend.width = 2, 
     box = FALSE,
     axes = FALSE,
     col = DarkSpectral,
     xlim = c(-12553230, 16774960), ylim = c(-8460601, 8373855),
     legend.args=list(text='Viral sharing events', side=2, line=0.4, cex=1.5))

plot(contsshp, lwd = 0.5, add = TRUE)

plot(GregRast(NewIntersects, "DeltaOverlapSharing.Futures1A"),
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

i <- 1

FitListList <- ModelListList <- DataList <- list()

for(i in (i:length(ClimateReps))[-5]){
  
  print(i)
  
  glue("~/Albersnet/Iceberg Files/{ClimateReps[i]}/Iceberg Output Files/NEGams.rds") %>% 
    readRDS -> NEPredictModels
  
  glue("~/Albersnet/Iceberg Files/{ClimateReps[i]}/Iceberg Output Files/FitList.rds") %>% 
    readRDS -> FitList
  
  glue("~/Albersnet/Iceberg Files/{ClimateReps[i]}/Iceberg Output Files/FullDataDF.rds") %>% 
    readRDS -> FullDataDF
  
  ModelListList[[ClimateReps[i]]] <- NEPredictModels
  
  DataList[[ClimateReps[i]]] <- FullDataDF
  
  FitListList[[ClimateReps[i]]] <- FitList
  
}

names(ModelListList) <- ClimateReps[-5]

ModelListList %>% Unlist1 %>% Unlist1 %>% Unlist1 %>% 
  
  map_dfr(function(a){
    
    a <- summary(a)
    
    data.frame(Mean = a$p.coeff, 
               SE = a$se,
               Var = names(a$p.coeff) %>% as.character() %>% str_remove("HabitatType")) %>%
      rbind(data.frame(Mean = 0, SE = 0, Var = "Cropland"), .) %>%
      filter(!Var == "(Intercept)")
    
  }, .id = "ClimateRep") -> BAMFixed

FitListList %>% Unlist1 %>% Unlist1 %>% Unlist1 %>% bind_rows(.id = "ClimateRep") %>% 
  group_by(Rep, ClimateRep) %>% summarise(Rich = last(Richness)) %>% 
  pull(Rich) %>% unique() -> Richnesses

FitListList %>% Unlist1 %>% Unlist1 %>%
  map("NoBat", id = "ClimateRep") %>% 
  bind_rows(.id = "Pipeline") %>%
  filter(HabitatType == "Cropland", 
         Richness %in% Richnesses) %>% 
  mutate(Bats = 'Non-bats') -> df1

FitListList %>% Unlist1 %>% Unlist1 %>%
  map("SomeBat", id = "ClimateRep") %>% 
  bind_rows(.id = "Pipeline") %>%
  filter(HabitatType == "Cropland", 
         Richness %in% Richnesses) %>% 
  mutate(Bats = "Bat encounters") -> df2

named <- c(A.Futures1 = 'RCP 2.6 (CLD)',
           A.Futures4 = 'RCP 8.5 (CLD)',
           C.Futures1 = 'RCP 2.6 (CL)',
           C.Futures4 = 'RCP 8.5 (CL)') 

named <- c(Futures1A = 'RCP 2.6 (CLD)',
           Futures4A = 'RCP 8.5 (CLD)',
           Futures1C = 'RCP 2.6 (CL)',
           Futures4C = 'RCP 8.5 (CL)') 

rbind(df1, df2) %>% 
  mutate(ClimateRep = substr(Pipeline, 1, 8)) %>% 
  mutate(Scenario = recode(Rep, !!!named)) %>% 
  ggplot(aes(Elevation, Fit, 
             #fill = Scenario, 
             colour = Scenario,
             group = paste0(Scenario, ClimateRep))) + 
  facet_wrap(~ Bats) + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, colour = NA) + 
  geom_line(alpha = 0.4) + guides(alpha = F) +
  labs(y = "Predicted Encounters",
       x = "Elevation",
       fill = 'Scenario', 
       colour = 'Scenario') +
  theme_cowplot()  +  
  scale_fill_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = c(3,5,7,8)) + 
  scale_colour_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = c(3,5,7,8)) + 
  #scale_y_log10(labels = scales::scientific) + #xlim(0,6100) + 
  scale_y_continuous(breaks = c((-1):6), labels = c(10^((-1):6))) +
  scale_x_continuous(breaks = c((-1):6), labels = c(10^((-1):6))) +
  theme(strip.background = element_blank()) -> g1; g1

FitListList %>% Unlist1 %>% Unlist1 %>%
  map('NoBat') %>% bind_rows(.id = "Pipeline") %>%
  filter(HabitatType == "Cropland", Elevation == last(unique(Elevation))) %>% 
  mutate(Bats = 'Non-bats')  -> df3

FitListList %>% Unlist1 %>% Unlist1 %>%
  map('SomeBat') %>% bind_rows(.id = "Pipeline") %>%
  filter(HabitatType == "Cropland", 
         Elevation == last(unique(Elevation))) %>% 
  mutate(Bats = "Bat encounters") -> df4

rbind(df3, df4) %>% 
  mutate(ClimateRep = substr(Pipeline, 1, 8)) %>% 
  mutate(Scenario = recode(Rep, !!!named)) %>% 
  ggplot(aes(Richness, Fit, 
             #fill = Scenario, 
             colour = Scenario,
             group = paste0(Scenario, ClimateRep))) + 
  facet_wrap(~ Bats) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, colour = NA) + 
  geom_line(alpha = 0.4) +
  labs(y = "Predicted Encounters",
       x = "Richness",
       fill = NULL, 
       # fill = 'Scenario', 
       colour = NULL)+
  theme_cowplot()  +  
  scale_fill_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = c(3,5,7,8)) + 
  scale_colour_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = c(3,5,7,8)) + 
  scale_y_continuous(breaks = c(seq(-1, 8, by = 2)), labels = c(10^(seq(-1, 8, by = 2)))) +
  scale_x_continuous(breaks = c(seq(-1, 8, by = 1)), labels = c(10^(seq(-1, 8, by = 1)))) +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank()) -> g2; g2

#### BAM panels

named <- c(Futures1A = 'RCP 2.6 (CLD)',
           Futures4A = 'RCP 8.5 (CLD)',
           Futures1C = 'RCP 2.6 (CL)',
           Futures4C = 'RCP 8.5 (CL)')

reclass <- c(Rangeland = 'Range',
             Nonforest = 'Other',
             Cropland = 'Crop', 
             Urban = 'Settled') 

BAMFixed %>% 
  filter(!Var %in% c("Elevation", "Richness")) %>%
  separate(ClimateRep, "[.]", into = c("ClimateRep", "Pipeline", "PredRep", "Bats")) %>% 
  mutate_at("Bats", ~str_replace_all(.x, c("SomeBat" = "Bat encounters", 
                                           "NoBat" = "Non-bat encounters"))) %>% 
  mutate(Model = paste0(PredRep, Pipeline)) %>% 
  mutate(Scenario = recode(Model, !!!named)) %>% 
  mutate(Var = recode(Var, !!!reclass)) %>%
  ggplot(aes(Var, Mean, colour = Scenario, group = paste0(Scenario, ClimateRep))) + theme_classic() + 
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
  labs(y = "Effect size", x = "Land Use", colour = NULL, fill = NULL) + 
  scale_color_discrete_sequential(palette = AlberPalettes[[3]], nmax = 10, order = c(5,7,8,10)) +
  lims(x = rev(c("Crop", "Settled", "Range", "Forest", "Other"))) + 
  coord_flip()  -> g3; g3

TextSize = 14

TextTheme <- theme(axis.title.x = element_text(size = TextSize),
                   axis.title.y = element_text(size = TextSize),
                   axis.text.y = element_text(angle = 0, hjust = 1))

g1 + TextTheme +
  g2 + TextTheme +
  g3 + TextTheme +
  plot_layout(heights = c(1.15,1,1)) +
  ggsave(filename = "Iceberg Figures/Figure2CDEFG.jpeg", 
         units = "mm", height = 220, width = 180,
         dpi = 600)

# Adding bat and nonbat encounters ####

# Figure 2 Plotting ####

pdf("Iceberg Figures/Figure2AB.pdf", width = 12, height = 10)

par(mfrow = c(2, 1), mar = c(0, 0.5, 0, 5))

NewIntersects %>% 
  mutate(Diff = Overlap.Futures1A - NoBatOverlap.Futures1A) %>% 
  GregRast("Diff") %>% 
  plot(axes = FALSE, box = FALSE,
       legend.shrink = 0.8, 
       legend.width = 2, 
       axis.args = list(cex.axis = 1.2),
       col = DarkSpectral,
       legend.args=list(text='First encounters', side=2, line=0.4, cex=1.5))

plot(contsshp, lwd = 0.5, add = TRUE)

plot(GregRast(NewIntersects, "NoBatOverlap.Futures1A"), axes = FALSE, box = FALSE, legend.shrink = 0.8, 
     #legend.shrink = 0.8, 
     legend.width = 2, 
     axis.args = list(cex.axis = 1.2),
     col = DarkSpectral,
     legend.args=list(text='First encounters', side=2, line=0.4, cex=1.5))

plot(contsshp, lwd = 0.5, add = TRUE)

dev.off()

# Figure 3BCD ####

pdf("Iceberg Figures/Figure3BCD.pdf", width = 12, height = 5)

par(mfrow = c(1, 3), 
    mar = c(0, 2, 0.2, 6), 
    oma = c(1, 1, 1, 2), 
    cex = 1.2)

plot(trim(sum(EboCurrents) + africa), 
     xlim = c(-20, 55), 
     ylim = c(-40,40), axes = FALSE, box = FALSE,
     col = DarkSpectral,
     legend.shrink = 0.5, 
     legend.width = 4,
     legend.args=list(text='Species', side=2, line=0.4, cex=1.5) ) # PANEL C

plot(africashp, add = TRUE, lwd = 0.5)

par(cex = 1.2)

plot(trim(sum(EboFutures) - sum(EboCurrents) + africa), 
     xlim = c(-20, 55), ylim = c(-40,40), axes = FALSE, box = FALSE,
     col = rev(colorRampPalette(brewer.pal(11,"RdBu"))(100)),
     legend.shrink = 0.5, 
     legend.width = 4, 
     legend.args=list(text=expression(paste(Delta, 'Species')), side=2, line=0.4, cex=1.5)) # PANEL C

plot(africashp, add = TRUE, lwd = 0.5)

par(cex = 1.2)

plot(trim(GregRast(EbolaNew, "Overlap.Futures1A") + africa), 
     xlim = c(-20, 55), ylim = c(-40,40), axes = FALSE, box = FALSE,
     col = DarkSpectral,
     legend.shrink = 0.5, 
     legend.width = 4, 
     legend.args=list(text='First encounters', side=2, line=0.4, cex=1.5))

plot(africashp, add = TRUE, lwd = 0.5)

dev.off()

# Figure 3A ####

BatNew <- readRDS('Iceberg Output Files/BPNewIntersects.rds')

BatNew %>% arrange(desc(Y), X) %>% as.data.frame -> BatNew

pdf("Iceberg Figures/Figure3A.pdf", width = 12, height = 10)

par(mfrow = c(1,1), mar=c(4,4,4,4.3))

plot(GregRast(BatNew, "OverlapSharing.Futures1A"), axes = FALSE, box = FALSE,
     col = DarkSpectral,
     legend.shrink = 0.4,
     legend.width = 1.5, 
     legend.args=list(text='Viral sharing events', side=2, line=0.3, cex=1.5))

plot(contsshp, lwd = 0.5,add = TRUE)

dev.off()

# Figure 3E ####

writeRaster(GregRast(NewIntersects, "DeltaOverlapSharing.Futures1A"), 
            file = "Iceberg Figures/Bivariate1A.tif", 
            overwrite = T)
