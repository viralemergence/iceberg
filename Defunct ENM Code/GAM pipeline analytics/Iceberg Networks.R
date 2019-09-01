
# Iceberg Networks ####

# Getting taxonomic patterns of predictions ####

Patterndf <- data.frame(Sp = AllMammals2) %>% 
  left_join(Panth1)

Patterndf[,SpaceVars] <-
  lapply(IcebergAdjList, function(a) rowSums(a[AllMammals2,AllMammals2])) %>% bind_cols

Patterndf[,SharingVars] <-
  lapply(PredNetworkList, function(a) rowSums(a[AllMammals2,AllMammals2])) %>% bind_cols

Patterndf[,paste0("Delta",SpaceVars[2:5])] <- Patterndf[,SpaceVars[2:5]] - Patterndf$Space.Currents
Patterndf[,paste0("Delta",SharingVars[2:5])] <- Patterndf[,SharingVars[2:5]] - Patterndf$Sharing.Currents

BarGraph(Patterndf, "hOrder", "DeltaSharing.Futures1", Order = T, Just = T, Text = "N") + scale_fill_discrete_sequential(palette = AlberPalettes[[1]])
BarGraph(Patterndf, "hOrder", "DeltaSpace.Futures1", Order = T, Just = T, Text = "N") + scale_fill_discrete_sequential(palette = AlberPalettes[[1]])
BarGraph(Patterndf, "hOrder", "CurrentSharing", Order = T, Just = T, Text = "N")
BarGraph(Patterndf, "hOrder", "FutureSharing", Order = T, Just = T, Text = "N")

colnames(Patterndf)[2:7] %>% lapply(function(a){
  
  BarGraph(Patterndf, "hOrder", a, Order = T, Just = T, Text = "N") + 
    theme(legend.position = "none") +
    scale_fill_discrete_sequential(palette = AlberPalettes[[1]])
  
}) %>% arrange_ggplot2(ncol = 3)

# Just new encounters

NewEncountersGraph <- graph_from_edgelist(NewEncounters[,c("Sp", "Sp2")] %>% as.matrix, directed = F)
E(NewEncountersGraph)$weights <- NewEncounters$CurrentSharing
E(NewEncountersGraph)$weights2 <- NewEncounters$FutureSharing

NewEncountersAdj <- get.adjacency(NewEncountersGraph, 
                                  attr = "weights"
) %>% 
  as.matrix #%>% as.data.frame

NewEncountersAdj2 <- get.adjacency(NewEncountersGraph, 
                                   attr = "weights2"
) %>% 
  as.matrix #%>% as.data.frame

NewEncountersPatterns <- NewEncountersAdj %>% 
  reshape2::melt(value.name = "CurrentSharing") %>% 
  rename(Sp = Var1, Sp2 = Var2) %>%
  inner_join(Panth1[,c("Sp", "hFamily", "hOrder")], by = c("Sp" = "Sp")) %>%
  inner_join(Panth1[,c("Sp", "hFamily", "hOrder")], by = c("Sp2" = "Sp")) %>% 
  right_join(NewEncountersAdj2 %>% 
               reshape2::melt(value.name = "FutureSharing") %>%
               rename(Sp = Var1, Sp2 = Var2)) %>%
  filter(!FutureSharing == 0) %>%
  mutate(#DeltaOverlap = 
    DeltaSharing = FutureSharing - CurrentSharing)

NewEncountersPatterns %>% group_by(hOrder.x, hOrder.y) %>% 
  summarise(FutureSharing = mean(FutureSharing), N = n()) %>% 
  ggplot(aes(hOrder.x, hOrder.y, fill = log(FutureSharing))) + 
  geom_tile() + 
  coord_fixed() + 
  geom_text(aes(label = N)) + 
  scale_fill_continuous_sequential(palette = "reds") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

colnames(Patterndf)[2:7] %>% lapply(function(a){
  
  BarGraph(Patterndf, "hOrder", a, Order = T, Just = T, Text = "N") + 
    theme(legend.position = "none") +
    scale_fill_discrete_sequential(palette = AlberPalettes[[1]])
  
}) %>% arrange_ggplot2(ncol = 3)

# Species level network ####

SpeciesGraph <- graph_from_edgelist(AllMammalMatrix[,c("Sp","Sp2")] %>% as.matrix(), directed = F)

SpeciesAdjList <- SpeciesGraphList <- list()

for(i in (which(colnames(AllMammalMatrix)=="Space.Currents")):ncol(AllMammalMatrix)){
  
  E(SpeciesGraph)$weight <- AllMammalMatrix[,i]
  
  SpeciesAdjList[[colnames(AllMammalMatrix)[i]]] <- get.adjacency(SpeciesGraph, attr = "weight")
  
  SpeciesGraphList[[colnames(AllMammalMatrix)[i]]] <- SpeciesGraph
  
}

# Order level network #### 

OrderSharingMeans <- AllMammalMatrix %>% 
  group_by(hOrder.x, hOrder.y) %>%
  summarise_at(vars(starts_with("Space"),
                    starts_with("Sharing"), 
                    starts_with("Delta")), 
               mean) %>% as.data.frame()

OrderSharingSums <- AllMammalMatrix %>% 
  group_by(hOrder.x, hOrder.y) %>%
  summarise_at(vars(starts_with("Space"),
                    starts_with("Sharing"), 
                    starts_with("Delta")), 
               sum) %>% as.data.frame()

OrderSharing <- full_join(OrderSharingMeans, OrderSharingSums, by = c("hOrder.x", "hOrder.y"), suffix = c(".Mean", ".Sum"))

Lower <- which(lower.tri(matrix(NA, nlevels(OrderSharing$hOrder.x), nlevels(OrderSharing$hOrder.x))))

OrderSharing <- OrderSharing[Lower,]

OrderGraph <- graph_from_edgelist(OrderSharing[,1:2] %>% as.matrix(), directed = F)

OrderAdjList <- OrderGraphList <- list()

for(i in (which(colnames(OrderSharing)=="Space.Currents.Mean")):ncol(OrderSharing)){
  
  E(OrderGraph)$weight <- OrderSharing[,i]
  
  OrderAdjList[[colnames(OrderSharing)[i]]] <- get.adjacency(OrderGraph, attr = "weight")
  
  OrderGraphList[[colnames(OrderSharing)[i]]] <- OrderGraph
  
}

# Misc plots for now #####
AllMammaldf2 %>% dplyr::select(Sharing.Currents:Sharing.Futures2) %>% ggpairs(lower = list(continuous = "smooth"))

