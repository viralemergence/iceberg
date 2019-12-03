
Sp <- "Nyctalus_noctula"

Currents <- readRDS(paste0("Iceberg Input Files/GretCDF/Currents/", Sp, ".rds")) %>% as.matrix() %>% as.data.frame()

Futures <- readRDS(paste0("Iceberg Input Files/GretCDF/Futures/", Sp, ".rds")) %>% as.matrix() %>% as.data.frame()

Currents %>% #%>% #bind_cols(Futures) %>% 
  ggplot(aes(X,Y,fill = ClimateLandUse)) + geom_tile()

Currents %>% bind_cols(Futures) %>% 
  ggplot(aes(X,Y,fill = BufferClimateLandUse.Futures1)) + geom_tile()
