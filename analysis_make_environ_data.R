
### Load necessary packages
source("load_packages.R")

### Download data from figshare
if(!file.exists("raw_data/environmental_data")){
  dir.create("raw_data/environmental_data")
  path = "raw_data/environmental_data/"
  # 2005 wind u
  download.file("https://ndownloader.figshare.com/files/10364952", paste0(path, "uwnd.10m.gauss.2005.nc"))
  # 2005 wind v
  download.file("https://ndownloader.figshare.com/files/10364955", paste0(path, "vwnd.10m.gauss.2005.nc"))
  # 2006 wind u
  download.file("https://ndownloader.figshare.com/files/10364964", paste0(path, "uwnd.10m.gauss.2006.nc"))
  # 2005 wind v
  download.file("https://ndownloader.figshare.com/files/10364967", paste0(path, "vwnd.10m.gauss.2006.nc"))
  #Geo current 1
  download.file("https://ndownloader.figshare.com/files/10365015", paste0(path, "geo_current1.nc"))
  #Geo current 2
  download.file("https://ndownloader.figshare.com/files/10365018", paste0(path, "geo_current2.nc"))
  #Geo current 3
  download.file("https://ndownloader.figshare.com/files/10365021", paste0(path, "geo_current3.nc"))
  #SST 1
  download.file("https://ndownloader.figshare.com/files/10365024", paste0(path, "sst1.nc"))
  #SST 2
  download.file("https://ndownloader.figshare.com/files/10365027", paste0(path, "sst2.nc"))
}

### Read netcdf data for wind and sst into a list of raster bricks there is one for each year,
### So, they need to be stiched together temporally
cov_list = list(
  surface_wind_u  = addLayer(
    brick("raw_data/environmental_data/uwnd.10m.gauss.2005.nc"), 
    brick("raw_data/environmental_data/uwnd.10m.gauss.2006.nc")
  ),
  surface_wind_v = addLayer(
    brick("raw_data/environmental_data/vwnd.10m.gauss.2005.nc"), 
    brick("raw_data/environmental_data/vwnd.10m.gauss.2006.nc")
  ),
  sst = addLayer(
    brick("raw_data/environmental_data/sst1.nc"), 
    brick("raw_data/environmental_data/sst2.nc")
  )
)

### Read in netcdf data for currents
### The geo. currents data resides in 2 netcdf files that need to be combined spatially

geo1_u = brick("raw_data/environmental_data/geo_current1.nc", varname="u")
geo1_v = brick("raw_data/environmental_data/geo_current1.nc", varname="v")
conn = nc_open("raw_data/environmental_data/geo_current1.nc")
curr_grid = ncdf4.helpers::nc.get.time.series(conn) %>% crawl::intToPOSIX()
nc_close(conn)
idx = which(curr_grid >= ymd("2005-09-01") & curr_grid <= ymd("2006-09-02"))
geo1_u = geo1_u[[idx]]
geo1_v = geo1_v[[idx]]

geo2_u = brick("raw_data/environmental_data/geo_current2.nc", varname="u")
geo2_v = brick("raw_data/environmental_data/geo_current2.nc", varname="v")
conn = nc_open("raw_data/environmental_data/geo_current2.nc")
curr_grid = ncdf4.helpers::nc.get.time.series(conn) %>% crawl::intToPOSIX()
nc_close(conn)
idx = which(curr_grid >= ymd("2005-09-01") & curr_grid <= ymd("2006-09-02"))
geo2_u = geo2_u[[idx]]
geo2_v = geo2_v[[idx]]

### Add geo currents to cov_list after merging together for a full N. Pacific raster stack
cov_list$geo_curr_u = raster::merge(geo1_u, geo2_u)
cov_list$geo_curr_v = raster::merge(geo1_v, geo2_v)


### This wasn't use in the fur seal analysis, but I left it here for future reference
### compute gradiants for SST
cov_list$sst_u = cov_list$sst
cov_list$sst_v = cov_list$sst
for(l in 1:nlayers(cov_list$sst)){
  tmp_grad = terrain(cov_list$sst[[l]], opt=c('slope','aspect'))
  cov_list$sst_u[[l]] = tmp_grad$slope*sin(tmp_grad$aspect)
  cov_list$sst_v[[l]] = tmp_grad$slope*cos(tmp_grad$aspect)
}

### obtain quadrature grid time for all covariates
quad_list = vector(mode="list", length=length(cov_list))
names(quad_list) = names(cov_list)

# surface winds u,v
conn = nc_open("raw_data/environmental_data/uwnd.10m.gauss.2005.nc")
quad_list$surface_wind_u = ncdf4.helpers::nc.get.time.series(conn)
nc_close(conn)

conn = nc_open("raw_data/environmental_data/uwnd.10m.gauss.2006.nc")
quad_list$surface_wind_u = c(quad_list$surface_wind_u, ncdf4.helpers::nc.get.time.series(conn))
nc_close(conn)
quad_list$surface_wind_u %>% intToPOSIX() -> quad_list$surface_wind_u

quad_list$surface_wind_v = quad_list$surface_wind_u

# sst, sst_u, sst_v
conn = nc_open("raw_data/environmental_data/sst1.nc")
quad_list$sst = ncdf4.helpers::nc.get.time.series(conn)
nc_close(conn)

conn = nc_open("raw_data/environmental_data/sst2.nc")
quad_list$sst = c(quad_list$sst, ncdf4.helpers::nc.get.time.series(conn))
nc_close(conn)
quad_list$sst %>% intToPOSIX() -> quad_list$sst
quad_list$sst_u <- quad_list$sst_v <- quad_list$sst

# geo curr
quad_list$geo_curr_u <- quad_list$geo_curr_v <- curr_grid[idx]

if(!dir.exists("data_products")) dir.create("data_products")
saveRDS(cov_list, file = "data_products/cov_list.rds")
saveRDS(quad_list, file = "data_products/quad_list.rds")

### Make some movies
# The code below can make some animations of the 3 different spatio-temporal environmental data sets
# if(!dir.exists("plots")) dir.create("plots")
# library(rasterVis)
# library(animation)
# load("data_products/pup_frame.RData")
# source("helper_code/make_hex_grid.R")
# sims = pup_frame %>% mutate(sims=purrr::map(sims,~as(.x,"data.frame"))) %>% 
#   dplyr::select(dbid, data, sims, hex_grid) %>% dplyr::select(dbid, sims) %>% unnest() %>% 
#   mutate(dbid=factor(dbid), group=paste0(dbid, reps)) %>% 
#   group_by(dbid) %>% mutate(mins = minute(Time) %>% factor() %>% as.numeric()) %>% ungroup() %>% 
#   filter(mins==4)
# coordinates(sims) = ~mu.x+mu.y
# proj4string(sims) = st_crs(npac_poly)[[2]]
# sims = spTransform(sims, paste0(proj4string(cov_list$sst), " +lon_wrap=180"))
# 
# 
# cov_list_sub = cov_list %>% map(crop, y=extent(sims))
# npac_sp = as(npac(resolution="i"), "Spatial") %>% spTransform(proj4string(sims)) %>% 
#   ms_clip(target=., bbox=sims@bbox %>% as.vector())
# 
# 
# ani.options(
#   interval = 0.12, 
#   ani.dev=function(...){png(res=72*4,...)},
#   ani.width=480*4,
#   ani.height=480*4
# )
# 
# ### SST
# myColorkey = list(
#   at=seq(-2,25,length=100),
#   labels=list(at=seq(-2,25,length=10))
# )
# saveVideo({
#   for(t in 1:dim(cov_list_sub$sst)[3]){
#   # for(t in 1:30){
#     tmp=levelplot(
#       cov_list_sub$sst[[t]], xlab="Longitude", ylab="Latitude", par.settings = viridisTheme, 
#       main=paste0("SST: ", as.character(quad_list$sst[t])), 
#       at=seq(-2,25,length=100),
#       colorkey=myColorkey,
#       margin=F) + layer(sp.polygons(npac_sp, fill='gray70'))
#     print(tmp)
#     cat(t, "  ")
#   }
# }, video.name = "plots/sst_video.mp4")
# 
# 
# ### Surface wind
# myColorkey = list(
#   at=seq(0,32,length=100),
#   labels=list(at=round(seq(0,32,3)))
# )
# saveVideo({
# for(t in 1:dim(cov_list_sub$surface_wind_u)[3]){
#   # for(t in 1:33){
#   w = brick(cov_list_sub$surface_wind_u[[t]], cov_list_sub$surface_wind_v[[t]])
#   tmp=vectorplot(
#     w, isField ="dXY", xlab="Longitude", ylab="Latitude", #par.settings = viridisTheme,
#     main=paste0("Surface wind: ", date(quad_list$surface_wind_u[t])), 
#     par.settings = viridisTheme, 
#     at=seq(0,32,length=100),
#     colorkey=myColorkey,
#     lwd.arrows=0.5,
#     col.arrows='white',
#     margin=F) + layer(sp.polygons(npac_sp, fill='gray70', col='gray70'))
#   print(tmp)
#   cat(t, "  ")
# }
# }, video.name="plots/winds_video.mp4")
# 
# 
# ### Ocean currents
# curr_vel = as.vector(sqrt(values(cov_list_sub$geo_curr_u)^2 + values(cov_list_sub$geo_curr_v)^2)) %>% 
#   .[!is.na(.)]
# 
# myColorkey = list(
#   at=seq(0,1.05, length=100),
#   labels=list(at=seq(0,0.7, 0.1))
# )
# saveVideo({
#   for(t in 1:dim(cov_list_sub$geo_curr_u)[3]){
#     w = brick(cov_list_sub$geo_curr_u[[t]], cov_list_sub$geo_curr_v[[t]])
#     tmp=vectorplot(
#       w, isField ="dXY", xlab="Longitude", ylab="Latitude", #par.settings = viridisTheme,
#       main=paste0("Geostrophic current: ", date(quad_list$geo_curr_u[t])), 
#       par.settings = viridisTheme, 
#       at=seq(0,1.05, length=100),
#       # colorkey=myColorkey,
#       lwd.arrows=0.5,
#       col.arrows='white',
#       margin=F) + layer(sp.polygons(npac_sp, fill='gray70', col='gray70'))
#     print(tmp)
#     cat(t, "  ")
#   }
# }, video.name="plots/curr_video.mp4")




