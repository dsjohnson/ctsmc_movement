save_hex_df = function(dbid, data, sims, hex_grid, ...){
  cov_list = readRDS("data_products/cov_list.rds")
  quad_list = readRDS("data_products/quad_list.rds")
  bound = data %>% pull(GMT) %>% range() %>% as.numeric() %>% crawl::intToPOSIX()
  hex_cov_df = hexify_df(cov_list, hex_grid, quad_list, bound=bound) %>% mutate(Time = with_tz(Time, "GMT"))
  saveRDS(hex_cov_df, paste0("data_products/",dbid, "_hex_df.rds"))
  return(paste0("data_products/",dbid, "_hex_df.rds"))
}


fit_ctsmc = function(dbid, sims, hex_df, hex_grid,...,model_forms){
  message("Starting ",dbid,"...")
  hex_cov_df = readRDS(hex_df)
  num_reps = max(sims$reps)
  fit_frame = expand.grid(
    dbid=dbid, 
    rep=1:num_reps) %>% as.tibble()
  quad_times = unique(hex_cov_df$Time)
  fit_frame = fit_frame %>% mutate(
    fit = foreach(k=1:num_reps) %dopar% {
      disc_path = make_disc_path(sims, k, hex_grid$poly, quad_times)
      glm_data = make_glm_data(disc_path, hex_grid$neighbor_df, hex_cov_df) %>% 
        mutate(
          wind = vector2scalar(surface_wind_u, surface_wind_v, bearing_to_next),
          geo_curr = vector2scalar(geo_curr_u, geo_curr_v, bearing_to_next)
        )
      fit_rep = vector("list", length(model_forms))
      for(p in 1:length(model_forms)){
        fit_rep[[p]] = gam(as.formula(model_forms[p]), family = "poisson", data=glm_data, select=T, control=list(scalePenalty=F), method="REML")
      }
      tibble(model=1:length(model_forms), fit = fit_rep, pmp=imputation_pmp_gam(fit_rep))
    }
  ) %>% unnest()
  
  if(!file.exists("results")) dir.create("results")
  saveRDS(fit_frame, paste0("results/",dbid, "_fit_frame.rds"))
  return(paste0("results/",dbid, "_fit_frame.rds"))
}