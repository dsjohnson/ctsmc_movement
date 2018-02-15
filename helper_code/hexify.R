hexify = function(x, hex_list, grid,idx){
  require(raster)
  require(sp)
  require(stringr)
  if(!missing(idx)){
    x = x[[idx]]
    grid = grid[idx]
  }
  proj_h = proj4string(hex_list$poly)
  if(is.na(proj_h)) stop("hex_list poly element must have proj4string")
  proj_x = proj4string(x)
  if(str_count(proj_x, "\\+proj=longlat")==1){
    if(extent(x)@xmax > 180){
      proj_x = paste0(proj_x, " +lon_wrap=180")
    }
  }
  poly = sp::spTransform(hex_list$poly, proj_x)
  x = crop(x, extend(extent(poly), res(x)))
  out = raster::extract(x, poly, weights=T)
  foo = function(m){
    if(is.null(m)) stop("error")
    if(nrow(m)>1){
      return(apply(m[,-ncol(m)], 2, weighted.mean, w=m[,ncol(m)]))
    } else return(m[,-ncol(m)])
  }   
  out = t(sapply(out, foo))
  out = data.frame(
    hex=as.numeric(sapply(slot(poly, "polygons"), function(x) slot(x, "ID"))), 
    out
  ) 
  colnames(out)[-1] = format(grid,"%Y-%m-%d %H:%M:%S")
  out %>% gather(key=Time, value=val, -hex) %>% arrange(hex, Time) %>% 
    mutate(Time = lubridate::ymd_hms(Time)) -> out
  return(out)
}

hexify_df = function(raster_list, hex_list, quad_list, bound){
  hex_cov = vector(mode="list", length(raster_list))
  for(j in 1:length(raster_list)){
    if(length(quad_list)==1) grid=quad_list[[1]]
    else grid=quad_list[[j]]
    idx = seq(min(which(grid>=bound[[1]]))-1, max(which(grid<=bound[[2]]))+1, 1)
    hex_cov[[j]] =  hexify(raster_list[[j]], hex_list, grid, idx)
    
    message(names(raster_list)[j], " hexified ...")
  }
  hex_cov_df = hex_cov[[1]] 
  if(length(raster_list)>1) for(j in 2:length(hex_cov)) hex_cov_df %>% full_join(hex_cov[[j]], by=c("hex", "Time")) -> hex_cov_df
  names(hex_cov_df)[-c(1:2)] = names(raster_list)
  hex_cov_df %>% arrange(hex, Time) -> hex_cov_df
  hex_cov_df = lapply(hex_cov_df, zoo::na.locf, na.rm=FALSE) %>% as.data.frame() -> hex_cov_df
  return(hex_cov_df)
}
