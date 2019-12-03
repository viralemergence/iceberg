

Currents <- IcebergAdjList[[1]]$Currents[lower.tri(IcebergAdjList[[1]]$Currents)] %>% c

Futures1 <- IcebergAdjList[[1]]$Futures1 %>% c

ggregplot::Prev(Currents==0&Futures1>0)

Futures2 <- IcebergAdjList[[1]]$Futures2 %>% c

lapply(IcebergAdjList, function(a){
  
  Currents <- (a$Currents[lower.tri(a$Currents)] %>% c)>0
  
  Futures <- lapply(a[2:5], function(b){
    
    length(which(Currents==0&c(b[lower.tri(b)])>0))
    
  })
  
}) -> YNList

lapply(IcebergAdjList, function(a){
  
  sapply(a, function(b) mean(c(b[lower.tri(b)])))  
  
}) -> MeanList

YNList %>% bind_rows(.id = "Pipeline") %>% 
  gather("RCP", "New Encounters", Futures1:Futures4) -> NETable

NETable %>% ggplot(aes(RCP, Pipeline, fill = `New Encounters`)) + 
  geom_tile() + coord_fixed() + 
  geom_text(aes(label = `New Encounters`)) +
  scale_y_discrete(labels = c("Full", "Dispersal", "Land Use", "Climate")) +
  scale_x_discrete(labels = c(RCPs)) +
  scale_fill_continuous_sequential(palette = AlberPalettes[[2]], 
                                   limits = c(0, 40000))



lapply()


