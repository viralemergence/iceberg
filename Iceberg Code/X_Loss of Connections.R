
# Loss of Connections ####

library(tidyverse); library(cowplot)

theme_set(theme_cowplot() + theme(strip.background = element_rect(fill = "white")))

CLC_database_hosts <- read.csv("~/Albersnet/CLC_database_hosts.csv")

CLC_Parasites <- read.csv("~/Albersnet/CLC_database_lifehistory.csv")

CLC_Parasites %<>% 
  
  mutate_at(c("Parasite.species"), ~.x %>% str_trim %>% str_replace_all(" ", "_"))

CLC_database_hosts %>% 
  
  mutate_at(c("Parasite.species", "Host.species"), ~.x %>% str_trim %>% str_replace_all(" ", "_")) ->
  
  CLC

CLC %<>% filter(Host.species %in% Panth1$Sp)

CLC %<>% rename_all(~.x %>% str_remove(".species$"))

# CLC %<>% dplyr::select(Parasite, Host, Stage)

CLC %>% group_by(Parasite) %>% summarise(N = nunique(Stage)) %>% filter(N>1) %>%
  pull(Parasite) -> N2Parasites

CLC %>% filter(Parasite %in% N2Parasites) ->
  
  CLC2

# Convert CLC to network ####

CLC2$Parasite %>% unique %>% 
  sort %>% 
  lapply(function(a){
    
    CLC2 %>% filter(Parasite == a) %>% 
      dplyr::select(Host, Stage) -> SubDF
    
    SubDF %>% table() %>% graph.incidence(weighted = NULL) %>%
      bipartite.projection %>% extract2("proj1") %>% get.adjacency %>%
      as.matrix %>% subtract(1, .) %>% reshape2::melt() %>%
      filter(value == 1, !(Var1 == Var2)) %>%
      rename(Sp = Var1, Sp2 = Var2) %>%
      mutate(Parasite.species = a)
    
  }) %>% bind_rows() %>% unique -> StedgeList

StedgeList %<>% mutate(Pair = paste0(Sp, ".", Sp2))

# OldEncounters <- readRDS("~/Albersnet/Iceberg Files/Climate1/Iceberg Output Files/OldEncounters.rds")

ClimateReps <- glue::glue("Climate{1:9}")

OldEncounters <- 
  ClimateReps %>% map(~readRDS(glue::glue("~/Albersnet/Iceberg Files/{.x}/Iceberg Output Files/OldEncounters.rds")))

names(OldEncounters) <- ClimateReps

OldEncounters %>% 
  unlist(recursive = F) %>% #unlist(recursive = F) %>%
  unlist(recursive = F) %>% #unlist(recursive = F) %>%
  map(~.x %>% 
        mutate(Pair = paste0(Sp, ".", Sp2)) %>% 
        filter(Pair %in% StedgeList$Pair) %>% dplyr::select(Sp, Sp2, hOrder.x, hOrder.y) %>%
        left_join(StedgeList, by = c("Sp", "Sp2"))) %>%
  bind_rows(.id = "Scenario") %>% 
  left_join(CLC %>% dplyr::select(Parasite.species = Parasite, Parasite.group) %>% unique) ->
  LostConnections

LostConnections$Parasite.species %>% setdiff(N2Parasites, .) %>% length %>% 
  divide_by(nunique(LostConnections$Parasite.species), .) %>% 
  multiply_by(100)

GroupCounts <- 
  CLC2 %>% dplyr::select(Parasite.species = Parasite, Parasite.group) %>% 
  unique %>% group_by(Parasite.group) %>% 
  count() %>% rename(N.Total = n)

LostSummary <- 
  LostConnections %>% 
  dplyr::select(Scenario, Parasite.species, Parasite.group) %>% unique %>% 
  group_by(Scenario, Parasite.group) %>% count() %>% 
  full_join(GroupCounts) %>% 
  mutate(N.Prop = n/N.Total) 

LostSummary %>% 
  ggplot(aes(Scenario, n, fill = Parasite.group)) + geom_col(position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

LostSummary %>% 
  ggplot(aes(Scenario, N.Prop, fill = Parasite.group)) + geom_col(position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

LostSummary %>% 
  ggplot(aes(n, fill = Parasite.group)) + 
  geom_histogram(colour = "black") + 
  labs(y = "Number of scenarios", x = "Number of lost connections") +
  facet_wrap(~Parasite.group, scales = "free") +
  lims(x = c(0, NA)) + theme(legend.position = "none")

LostSummary %>% 
  ggplot(aes(N.Prop, colour = Parasite.group)) + 
  geom_density()

LostSummary %<>% 
  mutate_at("Scenario", ~str_remove(.x, "Climate")) %>% 
  separate(Scenario, sep = "[.]", into = c("GCM", "Pipeline", "RCP")) %>% 
  mutate_at(c("GCM", "Pipeline", "RCP"), as.character)

Recode <- c("2.6", "4.5", "7.0", "8.5")
names(Recode) <- unique(LostSummary$RCP)

LostSummary %>% ungroup %>% 
  group_by(GCM, RCP, Pipeline) %>% 
  summarise(N.Prop = sum(n)/sum(N.Total)) %>% 
  ungroup %>% 
  group_by(RCP, Pipeline) %>% 
  summarise(Mean = mean(N.Prop),
            Lower = min(N.Prop),
            Upper = max(N.Prop)) %>% 
  mutate_at("RCP", ~str_replace_all(.x, Recode)) -> 
  LostSummarySummary

LostSummarySummary %<>% mutate_at(c("Mean", "Lower", "Upper"), ~.x*100)

LostSummarySummary %<>% mutate(Label = paste0(round(Mean, 1), "%\n(", round(Lower, 1), ", ", round(Upper, 1), ")"))

LostSummarySummary %>% 
  ggplot(aes(RCP, Pipeline)) + 
  geom_tile(aes(fill = Mean)) + 
  geom_text(aes(label = Label)) + 
  scale_y_discrete(labels = c("CLD", "CD", "CL", "C")) +
  scale_fill_continuous_sequential(palette = AlberPalettes[[2]], limits = c(0, NA)) +
  coord_fixed() + labs(fill = "% lost")

StedgeList %>% 
  inner_join(AllMammaldf %>% mutate(Pair = paste0(Sp, ".", Sp2)), 
             by = c("Pair")) -> 
  
  OverallStedgeDF

OverallStedgeDF$DeltaSpace.Futures1D %>% mean(na.rm = T)

# Stats for new paper ####

MeanSD <- function(x){
  
  list(Mean = mean(x), SD = sd(x))
  
}

# Scale of biotic network collapse due to disjunct host ranges
# # of species which lose an overlap
# Breakdown by RCP

LostSummary %>% group_by(GCM, RCP, Pipeline) %>% 
  summarise(N.Prop = sum(n)/sum(N.Total),
            n = sum(n)) %>% 
  ungroup %>% group_by(RCP) %>% 
  summarise(MeanSD = list(MeanSD(n))) %>% 
  pull(MeanSD)

# Breakdown by scenario

LostSummary %>% group_by(GCM, RCP, Pipeline) %>% 
  summarise(N.Prop = sum(n)/sum(N.Total),
            n = sum(n)) %>% 
  ungroup %>% group_by(Pipeline) %>% 
  summarise(MeanSD = list(MeanSD(n))) %>% 
  pull(MeanSD)

# Breakdown by parasite group (hypothesis about specificity)

LostSummary %>% 
  group_by(Parasite.group) %>% 
  summarise(MeanSD = list(MeanSD(N.Prop))) %>% 
  pull(MeanSD)

# Total extinction
# # of species which lose every overlap for at least one pair of life stages

CurrentsRangeAdjA <- readRDS("~/Albersnet/Iceberg Files/Climate8/Iceberg Output Files/CurrentsRangeAdjA.rds")
CurrentsRangeAdjB <- readRDS("~/Albersnet/Iceberg Files/Climate8/Iceberg Output Files/CurrentsRangeAdjB.rds")

CurrentsRangeAdjA %>% reshape2::melt() %>% 
  left_join(CurrentsRangeAdjB %>% reshape2::melt(), by = c("Var1", "Var2")) %>% 
  rename(Sp = Var1, Sp2 = Var2) -> AllOverlaps

AllOverlaps %<>% rename(Space.CurrentsA = value.x, Space.CurrentsB = value.y)

AllOverlaps %>% 
  dplyr::select(Sp, Sp2, A = Space.CurrentsA, B = Space.CurrentsB) %>% 
  tidyr::gather("Pipeline", "Space", A:B) -> 
  AllOverlapDF

AllOverlapDF %<>% filter(!is.na(Space))

LostConnections %<>%
  mutate_at("Scenario", ~str_remove(.x, "Climate")) %>%
  separate(Scenario, sep = "[.]", into = c("GCM", "Pipeline", "RCP")) %>%
  mutate_at(c("GCM", "Pipeline", "RCP"), as.character)

Recode <- c("2.6", "4.5", "7.0", "8.5")
names(Recode) <- unique(LostSummary$RCP)

ToCompare <- StedgeList %>% inner_join(AllOverlapDF, ., by = c("Sp", "Sp2")) %>% 
  arrange(Parasite.species, Sp, Sp2)

GoCompare <- ToCompare %>% filter(Space>0)

StageCounts <- GoCompare %>% group_by(Pipeline, Parasite.species) %>% count() %>% 
  rename(NPairs = n)

StageCounts <- 
  StageCounts %>% 
  bind_rows(StageCounts %>% mutate_at("Pipeline", ~str_replace_all(.x, c("A" = "C", "B" = "D"))))

LostConnections %>% 
  filter(Sp %in% colnames(CurrentsRangeAdjA), Sp2 %in% colnames(CurrentsRangeAdjA)) %>% 
  group_by(GCM, RCP, Pipeline, Parasite.species) %>%
  count %>%
  # mutate_at("Scenario", ~str_remove(.x, "Climate")) %>%
  # separate(Scenario, sep = "[.]", into = c("GCM", "Pipeline", "RCP")) %>%
  # mutate_at(c("GCM", "Pipeline", "RCP"), as.character) %>% 
  left_join(StageCounts, by = c("Parasite.species", "Pipeline")) %>% 
  mutate(PairsLost.Prop = n/NPairs) -> 
  ComparisonDF

# ComparisonDF$PairsLost.Prop %>% mean(na.rm = T)
# 
# ComparisonDF$PairsLost.Prop %>% qplot
# 
# ComparisonDF %>% filter(PairsLost.Prop > 0.2)
# 
# Failed <- ComparisonDF %>% filter(is.na(NPairs)) %>% pull(Parasite.species) %>% unique
# 
# ToCompare %>% filter(Parasite.species %in% Failed) %>% pull(Space) %>% qplot
# 
# GoCompare %>% filter(Parasite.species %in% Failed)# %>% pull(Space) %>% qplot
# 
# StageCounts %>% filter(Parasite.species %in% Failed) %>% # %>% pull(Space) %>% qplot
# # bind_rows(StageCounts)
#   pull(Pipeline)
#   
# LostConnections %>% 
#   filter(Sp %in% colnames(CurrentsRangeAdjA), Sp2 %in% colnames(CurrentsRangeAdjA)) %>% 
#   group_by(GCM, RCP, Pipeline, Parasite.species) %>%
#   count %>% 
#   filter(Parasite.species %in% Failed) %>% 
#   full_join(StageCounts %>% filter(Parasite.species %in% Failed), 
#             by = c("Pipeline", "Parasite.species")) %>% 
#   filter(is.na(NPairs))


# This will be small enough they can be individually identified and discussed
# As a %: what % of species lose all overlap range (compared to #/% in Carlson 2017 Sci Adv)

# Looking at the overlap between parasite pairs ####

i <- 1

print(i)

FocalGCM <- LostConnections[i, "GCM"]

FocalSpp <- LostConnections[i,c("Sp", "Sp2")] %>% unlist 

FocalSpp %>% 
  map(~glue::glue("~/Albersnet/Iceberg Files/Climate{FocalGCM}/Iceberg Input Files/GretCDF/Currents/{.x}.rds"))

UniqueLost <- LostConnections[,c("Sp", "Sp2", "Pair")] %>% unique

1:nrow(UniqueLost) %>% 
  
  lapply(function(i){
    
    print(i)
    
    FocalSpp <- LostConnections[i,c("Sp", "Sp2")] %>% unlist 
    
    FocalSpp %>% 
      map(~glue::glue("~/Albersnet/Iceberg Files/Climate1/Iceberg Input Files/GretCDF/Currents/{.x}.rds") %>% 
            readRDS)
    
  }) -> LostOverlaps

names(LostOverlaps) <- UniqueLost$Pair

RCPs <- glue::glue("Futures{1:4}")
GCMs <- 1:9
Pipelines <- LETTERS[1:4]

i <- 1
j <- RCPs[1]
k <- Pipelines[1]

LostTileList <- list()

for(i in GCMs){
  
  print(i)
  
  for(j in RCPs){
    
    print(j)
    
    for(k in Pipelines){
      
      print(k)
      
      LostConnections %>% filter(GCM == i, RCP == j, Pipeline == k) %>% 
        pull(Pair) %>% 
        map(function(b){
          
          if(length(b)>0){
            
            LostOverlaps[[b]] %>% 
              lapply(function(a){
                
                a %>% as.matrix %>% data.frame %>% mutate(A = ClimateLandUse,
                                                          B = Climate,
                                                          C = ClimateLandUse,
                                                          D = Climate)
                
              }) %>% reduce(~full_join(.x, .y, by = c("X", "Y"))) -> WideDF
            
            WideDF %<>% 
              dplyr::select(X, Y, matches(paste0("^", k, "..$"))) %>% 
              filter_all(~.x != 0)
            
            # colnames(WideDF)[3:4] <- b %>% str_split("[.]") %>% unlist
            
          }
        }) -> ListOne
      
      if(length(ListOne) > 0){
        
        ListOne %>% bind_rows %>% group_by(X, Y) %>% count -> LostTileDF
        
      }else{
        
        LostTileDF <- NULL
        
      }
      
      LostTileList[[paste(i, j, k, sep = ".")]] <- LostTileDF
      
    }
    
    # LostTileDF %>% ggplot(aes(X, Y, fill = AsBinary(n))) + geom_tile() + coord_sf()
    
  }
}

LostTileList %>% saveRDS("LostTileDF.rds")

# Hypothetical extinction from overlap loss
#   OPTIONAL! Would only work if it uses (all species, not just those that lose a total link) but (subset to those where 100% of their hosts at every stage are included)
#   Get area, before and after (every scenario), of maximum possible extent for at least one pathway through all life stages
#   Run range losses through Thomas et al. 2004 estimator
#   Compare to 8-10% statistic (from RCP 8.5?) from Sci Adv paper


