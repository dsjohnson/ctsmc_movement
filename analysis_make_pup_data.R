
### Load necessary packages
source("load_packages.R")

### Set up parallel evaluation
registerDoFuture()
plan(multisession)

### Load North Pacific map and projection
npac_poly = npac()
proj = st_crs(npac_poly)$proj4string %>% CRS()

### Source helper code
source("helper_code/make_hex_grid.R")
source("helper_code/crawl_funcs.R")

### Download data from figshare
if(!file.exists("raw_data/nfs_pups_2005.csv")){
  if(!file.exists("raw_data")) dir.create("raw_data")
  path = "raw_data/nfs_pups_2005.csv"
  download.file("https://ndownloader.figshare.com/files/10364766", "raw_data/nfs_pups_2005.csv")
}



### Read in pup telemetry data
pup_frame = read_csv("raw_data/nfs_pups_2005.csv",
                     col_types=cols(
                       dbid = col_integer(),
                       site = col_character(),
                       GMT = col_datetime(format = "%m/%d/%y %H:%M"),
                       Habitat = col_character(),
                       DateDeploy = col_date(format = "%m/%d/%y"),
                       loc_class = col_factor(levels = c('3', '2', '1', '0', 'A', 'B','Z')),
                       lat = col_double(),
                       long = col_double(),
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
  ) %>% ungroup() -> pup_frame

### Set seed for duplication
set.seed(123)

### Fit crawl model
# pup_frame %>%
#   mutate(crw_fit = future_map(data, crawl_fit,
#                               err.model=list(x=~loc_class-1),
#                               drift=TRUE,
#                               Time.name="GMT",
#                               fixPar=c(log(250), log(500), log(1500), NA, NA, NA,4, NA, NA),
#                               theta=c(log(2000), log(2000), 5, 2, 6),
#                               constr=list(
#                                 lower=c(rep(log(1500),2), rep(-Inf,3)), upper=rep(Inf,5)
#                               ),
#                               method='L-BFGS-B'
#   )
#   
#   ) -> pup_frame

sim_list <- foreach(i = 1:nrow(pup_frame)) %dopar% {
  ### Fit CTCRW model for MI
  crw_fit <- 
    crawl_fit(pup_frame$data[[i]],
              err.model=list(x=~loc_class-1),
              drift=TRUE,
              Time.name="GMT",
              fixPar=c(log(250), log(500), log(1500), NA, NA, NA,4, NA, NA),
              theta=c(log(2000), log(2000), 5, 2, 6),
              constr=list(
                lower=c(rep(log(1500),2), rep(-Inf,3)), upper=rep(Inf,5)
              ),
              method='L-BFGS-B')
  ### Create path posterior simulations...
  simulator <- crwSimulator(crw_fit, predTime="15 mins", parIS=0)
  sims <- sim_tracks(simulator,reps=20, CRSobj=proj)
  ### Form hex grid around simulated paths
  hex_grid <- hex_grid_spsample(sims, cellsize=40000, buffer=150000)
  list(crw_fit=crw_fit, simulator=simulator, sims=sims, hex_grid=hex_grid)
}

pup_frame$crw_fit <- map(sim_list, ~{.x[[1]]})
pup_frame$simulator <- map(sim_list, ~{.x[[2]]})
pup_frame$sims <- map(sim_list, ~{.x[[3]]})
pup_frame$hex_grid <- map(sim_list, ~{.x[[4]]})

rm(sims_list); gc()

if(!dir.exists("data_products")) dir.create("data_products")
save(pup_frame, npac_poly, proj, file="data_products/pup_frame.RData")


### Make data plots
### map and data plot

if(!dir.exists("plots")) dir.create("plots")

for(i in 1:length(unique(pup_frame$dbid))){
  message("plotting ", unique(pup_frame$dbid)[i], " ...")
  locs = pup_frame$data[[i]]
  hexes = pup_frame$hex_grid[[i]]$poly %>% st_as_sf()
  bb = st_bbox(hexes)
  
  p_data = ggplot() +
    geom_sf(data=npac_poly, fill='black', color='black') +
    geom_sf(data=hexes, fill=NA, color='darkgray', lwd=0.5) +
    coord_sf(xlim = c(bb["xmin"]-100000, bb["xmax"]+100000), ylim = c(bb["ymin"]-100000, bb["ymax"]+100000)) +
    xlab("Longitude") + ylab("Latitude") +
    # theme(panel.grid.major = element_line(colour = "white"), legend.position="none")
    cowplot::theme_cowplot() + theme(legend.position="none", panel.border = element_rect(colour = "black", fill=NA))
  
  sims=pup_frame$sims[[i]] %>% as(.,"data.frame")
  for(r in 1:max(sims$reps)){
    p_data = p_data + geom_path(data=sims[sims$reps==r,], aes(x=mu.x, y=mu.y), color='darkred', alpha=1, lwd=0.15)
  }
  
  p_data = p_data + geom_point(data=locs, aes(x=x, y=y), cex=0.5, col='blue')
  nm = paste0("plots/p_",unique(pup_frame$dbid)[i],"_data.pdf")
  ggsave(nm, p_data, dpi="retina", width=7, height=7, units = 'in')
  
}

sims <- pup_frame %>% mutate(sims=purrr::map(sims,~as(.x,"data.frame"))) %>%
  dplyr::select(dbid, data, sims, hex_grid) %>% dplyr::select(dbid, sims) %>% unnest() %>%
  mutate(dbid=factor(dbid), group=paste0(dbid, reps)) %>%
  group_by(dbid) %>% mutate(mins = minute(Time) %>% factor() %>% as.numeric()) %>% ungroup()

p_sims <- ggplot() +
  geom_sf(data=npac_poly, fill='black', color='black') +
  coord_sf(xlim = c(min(sims$mu.x), max(sims$mu.x)), ylim = c(min(sims$mu.y), max(sims$mu.y))) +
  geom_path(aes(y=mu.y, x=mu.x, color=dbid, group=group), alpha=1, lwd=0.15, data=sims %>% filter(mins==1)) +
  xlab("Longitude") + ylab("Latitude") +
  # theme(panel.grid.major = element_line(colour = "white"), legend.position="none")
  theme_cowplot() + theme(legend.position="none", panel.border = element_rect(colour = "black", fill=NA))
ggsave("plots/sims_plot.pdf", p_sims, dpi="retina", width=7, height=4, units = "in")
