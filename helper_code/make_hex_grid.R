hex_grid_spsample <- function(spdata, cellsize, buffer) {
  spdata@bbox <- bbox(spdata) + matrix(c(rep(-buffer,2), rep(buffer,2)), 2, 2)
  p4s = proj4string(spdata)
  
  hex_cent <- spsample(spdata, cellsize=cellsize, type="hexagonal", offset=c(0.5,0.5))
  hex_poly <- HexPoints2SpatialPolygons(hex_cent)
  hex_nb = dnearneigh(hex_cent, cellsize-0.01*cellsize, cellsize+0.01*cellsize)
  
  visit = unique(over(spdata, hex_poly))
  neighbor_df = data.frame(
    hex_to = unlist(hex_nb),
    hex_from = rep(1:length(hex_cent), sapply(hex_nb, length))
  ) %>% filter(hex_from %in% visit)
  
  idx = sort(unique(c(neighbor_df$hex_from, neighbor_df$hex_to)))
  
  hex_poly = hex_poly[idx]
  row.names(hex_poly) = as.character(1:length(hex_poly))
  
  hex_cent = hex_cent[idx]
  hex_nb = dnearneigh(hex_cent, cellsize-0.01*cellsize, cellsize+0.01*cellsize)
  
  visit = unique(over(spdata, hex_poly))
  neighbor_df = data.frame(
    hex_from = rep(1:length(hex_cent), sapply(hex_nb, length)),
    hex_to = unlist(hex_nb)
  ) 
  neighbor_df$visit = neighbor_df$hex_from %in% visit
  
  to_next =  coordinates(hex_cent[neighbor_df$hex_to,]) - coordinates(hex_cent[neighbor_df$hex_from,])
  neighbor_df$bearing_to_next = (atan2(to_next[,1], to_next[,2])*180/pi) %>% ifelse(.<0, 360+., .) 
  unit_xy = cbind(cos(pi*(90-neighbor_df$bearing_to_next)/180), sin(pi*(90-neighbor_df$bearing_to_next)/180))
  neighbor_df$north = as.vector(round(unit_xy %*% c(0,1),3))
  neighbor_df$east = as.vector(round(unit_xy %*% c(1,0),3))
  
  sp::proj4string(hex_poly) = p4s
  sp::proj4string(hex_cent) = p4s
  
  output <- list(
    poly=hex_poly,
    neighbor_df = neighbor_df,
    hex_centroids = hex_cent
  )
}