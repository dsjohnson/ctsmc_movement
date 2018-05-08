save_hex_df = function(dbid, data, sims, hex_grid, ...){
  cov_list = readRDS("data_products/cov_list.rds")
  quad_list = readRDS("data_products/quad_list.rds")
  bound = data %>% pull(GMT) %>% range() %>% as.numeric() %>% crawl::intToPOSIX()
  hex_cov_df = hexify_df(cov_list, hex_grid, quad_list, bound=bound) %>% mutate(Time = with_tz(Time, "GMT"))
  saveRDS(hex_cov_df, paste0("data_products/",dbid, "_hex_df.rds"))
  return(paste0("data_products/",dbid, "_hex_df.rds"))
}


fit_ctsmc = function(dbid, sims, hex_df, hex_grid,...){
  message("Starting ",dbid,"...")
  hex_cov_df = readRDS(hex_df)
  num_reps = max(sims$reps)
  # create storage data frame 
  fit_frame = expand.grid(
    dbid=dbid, 
    rep=1:num_reps) %>% as.tibble()
  quad_times = unique(hex_cov_df$Time)
  
  # perform fitting in parallel over imputations
  fit_frame = fit_frame %>% mutate(
    fit = foreach(k=1:num_reps) %dopar% {
      
      # spatially discritize path
      disc_path = make_disc_path(sims, k, hex_grid$poly, quad_times)
      
      # create GLM model data
      glm_data = make_glm_data(disc_path, hex_grid$neighbor_df, hex_cov_df) %>% 
        mutate(
          wind = vector2scalar(surface_wind_u, surface_wind_v, bearing_to_next),
          geo_curr = vector2scalar(geo_curr_u, geo_curr_v, bearing_to_next)
        ) %>% filter(!is.na(tsm))
      
      # Create knots and basis for RBF temporally varying coefs
      df=3
      kts = seq(min(glm_data$elapsed_time), max(glm_data$elapsed_time), length=df) %>% 
        c(.[1]-(.[2]-.[1]), ., .[df]+(.[df]-.[df-1])) 
      r=max(diff(kts))
      R = rbf_matrix(glm_data$elapsed_time, knots=kts, r)
      
      # Define model set
      base = " ~ log(tsm)*prev_move + R:north + R:east"
      model_forms= c("", " + R:wind", " + R:geo_curr", " + R:sst", 
                     " + R:wind + R:geo_curr", " + R:wind + R:sst", 
                     " + R:geo_curr + R:sst",
                     " + R:wind + R:geo_curr + R:sst") %>% paste(base,.) 
      
      # Fit each model and output to tibble data frame 
      tibble(model=1:length(model_forms)) %>% mutate(
        fit = foreach(p = 1:length(model_forms)) %do% {
          X = model.matrix(as.formula(model_forms[p]), data=glm_data)
          S = crossprod(X)
          gam(z ~ 0 + offset(log(delta)) + X, family = "poisson", paraPen=list(X=list(S)), data=glm_data, method="REML")
        }
      ) 
      
    }
  ) %>% unnest()
  
  # Save complete fitting data for an individual
  if(!file.exists("results")) dir.create("results")
  saveRDS(fit_frame, paste0("results/",dbid, "_fit_frame.rds"))
  return(paste0("results/",dbid, "_fit_frame.rds"))
}


rbf_matrix = function(x, knots, r){
  # X = x
  X=NULL
  for(k in 1:length(knots)){
    X = cbind(X, exp(-((x-knots[k])/r)^2))
  }
  return(X)
}

