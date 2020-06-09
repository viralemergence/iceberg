
# Loss of Connections ####

CLC_database_hosts <- read.csv("~/Albersnet/CLC_database_hosts.csv")

CLC_Parasites <- read.csv("~/Albersnet/CLC_database_lifehistory.csv")

CLC_Parasites %<>% 
  
  mutate_at(c("Parasite.species"), ~.x %>% str_trim %>% str_replace_all(" ", "_"))

CLC_database_hosts %>% 
  
  mutate_at(c("Parasite.species", "Host.species"), ~.x %>% str_trim %>% str_replace_all(" ", "_")) ->
  
  CLC

CLC %<>% filter(Host.species %in% Panth1$Sp)

CLC %<>% rename_all(~.x %>% str_remove(".species$"))

CLC %<>% dplyr::select(Parasite, Host, Stage)

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

OldEncounters %>% 
  unlist(recursive = F) %>% #unlist(recursive = F) %>%
  map(~.x %>% 
        mutate(Pair = paste0(Sp, ".", Sp2)) %>% 
        filter(Pair %in% StedgeList$Pair) %>% dplyr::select(Sp, Sp2, hOrder.x, hOrder.y) %>%
        left_join(StedgeList, by = c("Sp", "Sp2"))) %>%
  bind_rows(.id = "Scenario") %>% 
  left_join(CLC %>% dplyr::select(Parasite.species, Parasite.group) %>% unique) ->
  LostConnections

LostConnections %>% group_by(Scenario, Parasite.group) %>% count() %>% 
  ggplot(aes(Scenario, n, fill = Parasite.group)) + geom_col(position = "dodge") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

StedgeList %>% 
  inner_join(AllMammaldf %>% mutate(Pair = paste0(Sp, ".", Sp2)), 
            by = c("Pair")) -> 
  
  OverallStedgeDF

OverallStedgeDF$DeltaSpace.Futures1D %>% mean(na.rm = T)







