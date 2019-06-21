cont <- readOGR(dsn="C:/Users/cjcar/Desktop/FinalData", layer='contnew')

blank <- matrix(0,360*2,720*2)
blank <- raster(blank)
extent(blank) <- c(-180,180,-90,90)

cont <- clgeo_Clean(cont)

africa <- cont.r <- fasterize(st_as_sf(cont[1,]), blank)
aus <- cont.r <- fasterize(st_as_sf(cont[3,]), blank)
euras <- cont.r <- fasterize(st_as_sf(cont[4,]), blank)
n.am <- cont.r <- fasterize(st_as_sf(cont[5,]), blank)
s.am <- cont.r <- fasterize(st_as_sf(cont[6,]), blank)

