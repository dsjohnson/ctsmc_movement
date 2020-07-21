
### Load necessary packages
source("load_packages.R")

### Setup parallel processing
registerDoFuture()
plan(multisession, workers=6)

### Source 'helper' code
source("helper_code/fit_ctsmc.R")
source("helper_code/hexify.R")
source("helper_code/vector2scalar.R")
source("helper_code/make_disc_path.R")
source("helper_code/make_glm_data.R")
source("helper_code/pmp.R")

### Load pup data
pup_355 <- readRDS("results/pup_frame_fitted.rds") %>% filter(dbid=="355")

### create output directories and lists
if(!dir.exists("results")) dir.create("results")

hex_cov_df <- readRDS("data_products/355_hex_df.rds")
sims <- pup_355$sims[[1]]
num_reps <- max(sims$reps)
# create storage data frame 
fit_frame <- expand.grid(
  dbid="355", 
  rep=1:num_reps) %>% as_tibble()
quad_times <- unique(hex_cov_df$Time)
hex_grid <- pup_355$hex_grid[[1]]


# perform fitting in parallel over imputations
fit_frame <- fit_frame %>% mutate(
  fit = foreach(k=1:num_reps) %dopar% {
    
    # spatially discrete path
    disc_path <- make_disc_path(sims, k, hex_grid$poly, quad_times)
    
    # create GLM model data
    glm_data <- make_glm_data(disc_path, hex_grid$neighbor_df, hex_cov_df) %>% 
      mutate(
        wind = vector2scalar(surface_wind_u, surface_wind_v, bearing_to_next),
        geo_curr = vector2scalar(geo_curr_u, geo_curr_v, bearing_to_next)
      ) %>% filter(!is.na(tsm))
    
    glm_data <- glm_data %>% mutate(
      full_move = zoo::na.locf(move, na.rm=F, fromLast=T)
    )
    
    # Create constant covariates between moves for fitting CTMC model version
    ctmc_data <- glm_data %>% 
      dplyr::select(tsm, hex_to, bearing_to_next, full_move, wind, geo_curr, sst, elapsed_time) %>% 
      group_by(full_move) %>% nest() %>% mutate(
        data = map(data, ~{filter(.x, tsm==min(tsm))})
      ) %>% unnest(cols=c(data)) %>% 
      dplyr::rename(wind_ctmc=wind, geo_curr_ctmc=geo_curr, sst_ctmc=sst, elapsed_time_ctmc=elapsed_time) %>% 
      dplyr::select(-tsm)
    
    glm_data <- glm_data %>% left_join(ctmc_data)
    
    
    
    # Create knots and basis for RBF temporally varying coefs
    df <- 3
    kts <- seq(min(glm_data$elapsed_time_ctmc), max(glm_data$elapsed_time_ctmc), length=df) %>% 
      c(.[1]-(.[2]-.[1]), ., .[df]+(.[df]-.[df-1])) 
    r <- max(diff(kts))
    R_ctmc <- rbf_matrix(glm_data$elapsed_time_ctmc, knots=kts, r)
    
    # Define model set
    base_ctmc <- " ~ prev_move + R_ctmc:north + R_ctmc:east"
    model_forms_ctmc <- c("", " + R_ctmc:wind_ctmc", " + R_ctmc:geo_curr_ctmc", " + R_ctmc:sst_ctmc", 
                          " + R_ctmc:wind_ctmc + R_ctmc:geo_curr_ctmc", " + R_ctmc:wind_ctmc + R_ctmc:sst_ctmc", 
                          " + R_ctmc:geo_curr_ctmc + R_ctmc:sst_ctmc",
                          " + R_ctmc:wind_ctmc + R_ctmc:geo_curr_ctmc + R_ctmc:sst_ctmc") %>% paste(base_ctmc,.) 
    
    
    kts <- seq(min(glm_data$elapsed_time), max(glm_data$elapsed_time), length=df) %>% 
      c(.[1]-(.[2]-.[1]), ., .[df]+(.[df]-.[df-1])) 
    r <- max(diff(kts))
    R <- rbf_matrix(glm_data$elapsed_time, knots=kts, r)
    
    base = " ~ log(tsm)*prev_move + R:north + R:east"
    model_forms= c("", " + R:wind", " + R:geo_curr", " + R:sst", 
                   " + R:wind + R:geo_curr", " + R:wind + R:sst", 
                   " + R:geo_curr + R:sst",
                   " + R:wind + R:geo_curr + R:sst") %>% paste(base,.) 
    model_forms <- c(model_forms_ctmc, model_forms)
    
    # Fit each model and output to tibble data frame 
    tibble(model=1:length(model_forms)) %>% mutate(
      fit = foreach(p = 1:length(model_forms)) %do% {
        X <- model.matrix(as.formula(model_forms[p]), data=glm_data)
        S <- crossprod(X)
        gam(z ~ 0 + offset(log(delta)) + X, family = "poisson", paraPen=list(X=list(S)), data=glm_data, method="REML")
      }
    ) 
  }
) %>% unnest()

if(!file.exists("results")) dir.create("results")
saveRDS(fit_frame, paste0("results/355_ctmc_fit_frame.rds"))
plan('sequential')

fit_frame <- fit_frame %>% group_by(dbid, rep) %>% nest() %>% 
  dplyr::mutate(
    pmp = purrr::map(data, ~calc_pmp(.x$fit,"laplace")),
    pmp_bic = purrr::map(data, ~calc_pmp(.x$fit,"bic"))
  ) %>% tidyr::unnest(cols = c('data', 'pmp', 'pmp_bic'))

pmp_355 <- fit_frame %>% group_by(dbid, model) %>% 
  summarise(pmp_est=mean(pmp), pmp_est_bic=mean(pmp_bic)) %>% ungroup()

model_forms <- c("{\\tt base} only", " {\\tt wind$_{ijt}$}", " {\\tt curr$_{ijt}$}", " {\\tt sst$_{it}$}", 
               " {\\tt wind$_{ijt}$ + curr$_{ijt}$}", "{\\tt wind$_{ijt}$ + sst$_{it}$}", 
               " {\\tt curr$_{ijt}$ + sst$_{it}$}",
               " {\\tt wind$_{ijt}$ + curr$_{ijt}$ + sst$_{it}$}") 
pmp_355$model <- rep(model_forms,2)
pmp_355$class <- c(rep("Markov", 8), rep("semi-Markov", 8))

saveRDS(pmp_355, "results/pmp_355.rds")



