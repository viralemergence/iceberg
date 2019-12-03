
head(ChangeDF)

ChangeDF %>% filter(Rep == 'Climate.Futures1') %>% 
  pull(OverallChange) %>% na.omit() %>% mean()


ChangeDF %>% filter(Rep == 'BufferClimateLandUse.Futures1') %>% 
  pull(OverallChange) %>% na.omit() %>% mean()
