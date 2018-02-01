
### Load necessary packages
source("load_packages.R")

### Load North Pacific map and projection
npac_poly = npac()
proj = st_crs(npac_poly)$proj4string %>% CRS()

### Source helper code
source("helper_code/make_hex_grid.R")
source("helper_code/crawl_funcs.R")

### Read in pup telemetry data
pup_frame = read_csv("raw_data/Pups_2005.csv", 
                     col_types=cols(
                       id = col_character(),
                       dbid = col_integer(),
                       instrumenttype = col_character(),
                       ptt = col_integer(),
                       site = col_character(),
                       GMT = col_datetime(format = "%m/%d/%Y %H:%M:%S"),
                       depcode = col_integer(),
                       Habitat = col_character(),
                       DateDeploy = col_date(format = "%m/%d/%y"),
                       loc_class = col_factor(levels = c('3', '2', '1', '0', 'A', 'B','Z')),
                       lat = col_double(),
                       long = col_double(),
                       long2 = col_double(),
                       sex = col_factor(levels = c('F','M'))
                     )
) %>% 
  dplyr::filter(site!="SMP", loc_class!='Z', loc_class!="B") %>% 
  droplevels() %>% 
  dplyr::select(dbid, site, sex, GMT, long, lat, loc_class) 

### Extract long deployment females animals
pup_frame %>% group_by(dbid, sex) %>% summarise(max_date=max(GMT)) %>% 
  dplyr::filter(month(max_date)>=6 & year(max_date)==2006 & sex=='F') %>% dplyr::select(dbid) %>% 
  left_join(pup_frame) -> pup_frame


### Project pup telemetry data
coordinates(pup_frame) = ~ long+lat
proj4string(pup_frame) = CRS("+init=epsg:4326")
pup_frame = spTransform(pup_frame, proj) %>% as("data.frame") %>% 
  rename(x=long, y=lat) 

### Create some deployment summaries for futher filtering
pup_frame %>% group_by(dbid, site, sex) %>% nest() %>% 
  mutate(
    duration = map_dbl(data, ~max(.x$GMT)-min(.x$GMT)),
    first_obs = map_dbl(data, ~min(.x$GMT)) %>% intToPOSIX(),
    last_obs = map_dbl(data, ~max(.x$GMT)) %>% intToPOSIX()
  ) -> pup_frame

### Set seed for duplication
set.seed(123)

### Fit crawl model
pup_frame %>% 
  mutate(crw_fit = map(data, crawl_fit, 
                       err.model=list(x=~loc_class-1), 
                       drift=TRUE,
                       Time.name="GMT",
                       fixPar=c(log(250), log(500), log(1500), NA, NA, NA,4, NA, NA),
                       theta=c(log(2000), log(2000), 5, 0, 0),
                       constr=list(
                         lower=c(rep(log(1500),2), rep(-Inf,3)), upper=rep(Inf,5)
                       ),
                       method='L-BFGS-B'
  )) -> pup_frame


### Create path posterior simulations...
pup_frame %>% mutate(
  simulator = map(crw_fit, crwSimulator, predTime="15 mins", parIS=0),
  sims = map(simulator, sim_tracks, reps=20, CRSobj=proj)
) -> pup_frame


### Form hex grid around simulated paths
pup_frame %>% 
  mutate(hex_grid = map(sims, hex_grid_spsample, cellsize=40000, buffer=150000) ) -> pup_frame

if(!dir.exists("data_products")) dir.create("data_products")
save(pup_frame, npac_poly, proj, file="data_products/pup_frame.RData")


### Make data plots
### map and data plot

if(!dir.exists("plots")) dir.create("plots")

for(i in 1:length(unique(pup_frame$dbid))){
  locs = pup_frame$data[[i]] 
  hexes = pup_frame$hex_grid[[i]]$poly %>% st_as_sf()
  bb = st_bbox(hexes)
  
  p_data = ggplot() + 
    geom_sf(data=npac_poly, fill='black', color='black') +
    geom_sf(data=hexes, fill=NA) +
    coord_sf(xlim = c(bb["xmin"]-100000, bb["xmax"]+100000), ylim = c(bb["ymin"]-100000, bb["ymax"]+100000)) +
    xlab("Longitude") + ylab("Latitude") +
    theme(panel.grid.major = element_line(colour = "white"), legend.position="none") +
    scale_x_continuous(breaks = seq(-180,180,5)) +
    scale_y_continuous(breaks = seq(0,90,5))
  
  sims=pup_frame$sims[[i]] %>% as(.,"data.frame")
  for(r in 1:max(sims$reps)){
    p_data = p_data + geom_path(data=sims[sims$reps==r,], aes(x=mu.x, y=mu.y), color='red', alpha=0.2)
  }
  
  p_data = p_data + geom_point(data=locs, aes(x=x, y=y), cex=0.5, col='blue')
  nm = paste0("plots/p_",unique(pup_frame$dbid)[i],"_data.jpg")
  ggsave(nm, p_data, dpi="retina", width=7, height=7, units = 'in')
  
  rm('locs', 'hexes', 'bb', 'p_data', 'sims')
}

sims = pup_frame %>% mutate(sims=purrr::map(sims,~as(.x,"data.frame"))) %>% 
  select(dbid, data, sims, hex_grid) %>% select(dbid, sims) %>% unnest() %>% 
  mutate(dbid=factor(dbid), group=paste0(dbid, reps)) %>% 
  group_by(dbid) %>% mutate(mins = minute(Time) %>% factor() %>% as.numeric()) %>% ungroup()

p_sims = ggplot() + 
  geom_sf(data=npac_poly, fill='black', color='black') +
  coord_sf(xlim = c(min(sims$mu.x), max(sims$mu.x)), ylim = c(min(sims$mu.y), max(sims$mu.y))) +
  # scale_x_continuous(breaks = seq(0, 360, by=10)) +
  # scale_y_continuous(breaks = seq(0, 180, by=5)) +
  geom_path(aes(y=mu.y, x=mu.x, color=dbid, group=group), alpha=0.2, data=sims %>% filter(mins==1)) +
  xlab("Longitude") + ylab("Latitude")+
  theme(panel.grid.major = element_line(colour = "white"), legend.position="none")
ggsave("plots/sims_plot.jpg", p_sims, dpi="retina", width=7, height=4, units = "in")
