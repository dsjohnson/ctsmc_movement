
### Load necessary packages
source("load_packages.R")

### Setup parallel processing
registerDoFuture()
plan(multisession)
options(future.progress=TRUE)

### Source 'helper' code
source("helper_code/fit_ctsmc.R")
source("helper_code/hexify.R")
source("helper_code/vector2scalar.R")
source("helper_code/make_disc_path.R")
source("helper_code/make_glm_data.R")
source("helper_code/pmp.R")

### Load pup data
load("data_products/pup_frame.RData")

### Load environmental data
# load("data_products/environ_cov_data.rda")

### create output directories and lists
if(!dir.exists("results")) dir.create("results")

### Model set
base = "z ~offset(log(delta)) + log(tsm)*prev_move + s(elapsed_time,by=north,k=4) + s(elapsed_time,by=east,k=4)"
model_forms = c(
  base,
  paste0(base, " + s(elapsed_time,by=wind, k=4)"),
  paste0(base, " + sst"),
  paste0(base, " + s(elapsed_time,by=geo_curr, k=4)"),
  paste0(base, " + s(elapsed_time,by=wind, k=4) + sst"),
  paste0(base, " + s(elapsed_time,by=wind, k=4) + s(elapsed_time,by=geo_curr, k=4)"),
  paste0(base, " + sst + s(elapsed_time,by=geo_curr, k=4)"),
  paste0(base, " + s(elapsed_time,by=wind, k=4) + s(elapsed_time,by=geo_curr, k=4) + sst")
)


### Convert environmental data to hex grids and save the file names to the pup data
pup_frame %>% mutate(
  hex_df = pmap(.,~future(save_hex_df(...))) %>% values
) -> pup_frame


### Fit models
pup_frame %>% mutate(
  fit = pmap(.,~fit_ctsmc(...,model_forms = model_forms))
) -> pup_frame







#################################################################
for(i in 1:nrow(pup_frame)){
  pmp_store[[i]] = matrix(NA, nrow=max(pup_frame$sims[[i]]$reps), ncol=length(model_forms))
  fit_frame = expand.grid(
    dbid=pup_frame$dbid[i], 
    model=1:length(model_forms), 
    rep=1:max(pup_frame$sims[[i]]$reps), 
    fit=vector("list",1)) %>% 
    as.tibble()
  
  # initial processing of covariates to glm form
  hex_list = pup_frame$hex_grid[[i]]
  bound = list(pup_frame$first_obs[i], pup_frame$last_obs[i])
  hex_cov_df = hexify_df(cov_list, hex_list, quad_list, bound=bound) %>% 
    mutate(Time = with_tz(Time, "GMT"))
  saveRDS(hex_cov_df, paste0("data_products/",pup_frame$dbid[i], "_hex_cov_df.rds"))
  
  for(k in 1:max(pup_frame$sims[[i]]$reps)){
    message("**** ",i, "-", k, " ****")
    # for(k in 1:10){
    # create discrete path
    quad_times = unique(hex_cov_df$Time)
    disc_path = make_disc_path(pup_frame$sims[[i]], k, hex_list$poly, quad_times)
    
    # create glm data
    glm_data = make_glm_data(disc_path, hex_list$neighbor_df, hex_cov_df) %>% 
      mutate(
        wind = vector2scalar(surface_wind_u, surface_wind_v, bearing_to_next),
        geo_curr = vector2scalar(geo_curr_u, geo_curr_v, bearing_to_next)
      )
    
    # Fit models 
    for(p in 1:length(model_forms)){
      fit_frame %>% mutate(
        fit = replace(fit, rep==k & model==p, list(gam(as.formula(model_forms[p]), family = "poisson", data=glm_data, control=list(scalePenalty=F), method="REML")))
      ) -> fit_frame
      message("       p = ",p)
    }
    # Calculate post. mod. prob.
    pmp_store[[i]][k,] = fit_frame %>% filter(rep==k) %>% pull(fit) %>% imputation_pmp_gam()
  }
  saveRDS(fit_frame, paste0("results/",pup_frame$dbid[i], "_fit_frame.rds"))
}

saveRDS(pmp_store, "results/pmp.rds")
