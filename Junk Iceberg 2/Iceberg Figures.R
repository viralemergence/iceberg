
# Overall SI Iceberg Figures ####

# Map: First contacts, sum probabilities, mean probabilities (three panel) x 4 scenarios

jpeg("Iceberg Figures/NewOverlaps.jpeg", units = "mm", width = 150, height = 250, res = 300)

list(
  NEGridList %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
    geom_tile(aes(fill = OverlapSum)) + 
    coord_fixed() + 
    theme_void() +
    facet_grid(~PredRep),
  
  NEGridList %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
    geom_tile(aes(fill = SharingSum)) + 
    coord_fixed() +
    theme_void() +
    facet_grid(~PredRep),
  
  NEGridList %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
    geom_tile(aes(fill = SharingMean)) + 
    coord_fixed() +
    theme_void() +
    facet_grid(~PredRep)
  
) %>% arrange_ggplot2(nrow = 3)

dev.off()

FillVars <- c("OverlapSum", "SharingSum", "SharingMean")
VarLabs = c("Richness", "SharingSum", "SharingMean")
names(VarLabs) <- FillVars

MaxFills <- VarLabs[FillVars] %>% sapply(function(a) (GridList %>% bind_rows)[,a] %>% max)
NEMaxFills <- FillVars %>% sapply(function(a) (NEGridList %>% bind_rows)[,a] %>% max)

for(a in PredReps){
  print(a)
  
  for(b in FillVars){   
    print(b)  
    
    print(MaxFills[b])
    
    GridList[[a]][,"Fill"] <- GridList[[a]][,VarLabs[b]]
    
    GridList[[a]] %>% 
      mutate(Long = X, Lat = Y) %>%
      ggplot(aes(Long, Lat)) + 
      geom_tile(aes(fill = Fill)) +
      theme_void() + 
      scale_fill_continuous_sequential(palette = AlberPalettes[[1]], limits = c(0, MaxFills[b])) +
      labs(fill = VarLabs[b]) + coord_fixed() +
      ggsave(file = paste0("Iceberg Figures/",a,".",b,".jpeg"), units = "mm", width = 200, height = 100, dpi = 300)
    
    if(a != "Currents"){
      
      NEGridList[[a]][,"Fill"] <- NEGridList[[a]][,b]
      
      print(NEMaxFills[b])
      
      NEGridList[[a]] %>% 
        mutate(Long = X, Lat = Y) %>%
        ggplot(aes(Long, Lat)) + 
        geom_tile(aes(fill = Fill)) +
        theme_void() + 
        scale_fill_continuous_sequential(palette = AlberPalettes[[1]], limits = c(0, NEMaxFills[b])) +
        labs(fill = b) + coord_fixed() +
        ggsave(file = paste0("Iceberg Figures/NE.",a,".",b,".jpeg"), units = "mm", width = 200, height = 100, dpi = 300)
      
    }
  }
}

# Map: Species richness, sum probabilities, mean probabilities (three panel) x 4 scenarios

jpeg("Iceberg Figures/OverallMaps.jpeg", units = "mm", width = 150, height = 250, res = 300)

list(
  GridList %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
    geom_tile(aes(fill = Richness)) + 
    coord_fixed() + 
    theme_void() +
    facet_grid(~PredRep),
  
  GridList %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
    geom_tile(aes(fill = SharingSum)) + 
    coord_fixed() +
    theme_void() +
    facet_grid(~PredRep),
  
  GridList %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
    geom_tile(aes(fill = SharingMean)) + 
    coord_fixed() +
    theme_void() +
    facet_grid(~PredRep)
  
) %>% arrange_ggplot2(nrow = 3)

dev.off()

# Map: bat-primate new contacts weighted sum for four scenarios

jpeg("Iceberg Figures/BPNewEncounters.jpeg", units = "mm", width = 350, height = 200, res = 300)

list(
  BPNEGridList %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
    geom_tile(aes(fill = OverlapSum)) + 
    coord_fixed() +
    theme_void() + 
    theme(legend.position = "top") +
    scale_fill_continuous_sequential(palette = AlberPalettes[[1]]) +
    facet_grid(PredRep~.),
  
  BPNEGridList %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
    geom_tile(aes(fill = SharingSum)) + 
    coord_fixed() +
    theme_void() +
    theme(legend.position = "top") +
    scale_fill_continuous_sequential(palette = AlberPalettes[[1]]) +
    facet_grid(PredRep~.),
  
  BPNEGridList %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
    geom_tile(aes(fill = SharingMean)) + 
    coord_fixed() +
    theme_void() +
    theme(legend.position = "top") +
    scale_fill_continuous_sequential(palette = AlberPalettes[[1]]) +
    facet_grid(PredRep~.)
  
) %>% arrange_ggplot2(nrow = 1)

dev.off()

# Map: new ebola host-mammal meetings ####

jpeg("Iceberg Figures/EbolaNewEncounters.jpeg", units = "mm", width = 350, height = 200, res = 300)

list(
  EbolaNEGridList[c(1,4)] %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
    geom_tile(aes(fill = OverlapSum)) + 
    coord_fixed() +
    theme_void() + 
    theme(legend.position = "top") +
    scale_fill_continuous_sequential(palette = AlberPalettes[[1]]) +
    facet_grid(PredRep~.),
  
  EbolaNEGridList %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
    geom_tile(aes(fill = SharingSum)) + 
    coord_fixed() +
    theme_void() +
    theme(legend.position = "top") +
    scale_fill_continuous_sequential(palette = AlberPalettes[[1]]) +
    facet_grid(PredRep~.),
  
  EbolaNEGridList %>% bind_rows() %>%  ggplot(aes(X, Y)) + 
    geom_tile(aes(fill = SharingMean)) + 
    coord_fixed() +
    theme_void() +
    theme(legend.position = "top") +
    scale_fill_continuous_sequential(palette = AlberPalettes[[1]]) +
    facet_grid(PredRep~.)
  
) %>% arrange_ggplot2(nrow = 1)

dev.off()


# Map: two-panel x 4 scenario within vs between order link sharing



# Map: human risk red-blue against sum probabilities



# Misc plots that colin didn't ask for

# Number of encounters per species ####

PredReps[2:5] %>% lapply(function(a){
  
  table(c(NewEncountersList[[a]]$Sp,NewEncountersList[[a]]$Sp2)) %>% qplot
  
}) %>% arrange_ggplot2




