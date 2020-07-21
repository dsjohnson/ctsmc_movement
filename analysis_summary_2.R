source("load_packages.R")
source("helper_code/pmp.R")
source("helper_code/fit_ctsmc.R")

### Set up parallel evaluation
# registerDoFuture()
# plan(multisession)


pup_frame <- readRDS("results/pup_frame_fitted.rds") %>% 
  dplyr::select(dbid, fit)
pup_frame$tilde_beta <- vector("list", nrow(pup_frame))
for(i in 1:nrow(pup_frame)){
  # i <- 5 # for animal 355
  fit_frame <- readRDS(pup_frame$fit[[i]]) %>% filter(model==8)
  gc()
  ps <- fit_frame %>% mutate(
    b = map(fit, coef),
    V = map(fit, vcov),
    ps = map2(b, V, ~rmvnorm(5000, .x, .y)) %>% map(., as.data.frame)
  ) %>% dplyr::select(ps) %>% unnest(cols=c(ps)) %>% as.matrix() %>% {.[,c(2,4)]}
  out_df <- data.frame(previous = c(-1, -0.5, 0.5, 1))
  for(j in 1:5){
    tilde_beta <- mcmc(ps[,1] + out_df$previous[j]*ps[,2])
    ci <- round(HPDinterval(tilde_beta),2)
    out_df$Estimate[j] <- round(mean(tilde_beta),2)
    out_df$CI[j] <- paste0("[",ci[1],", ", ci[2], "]")
  }
  pup_frame$tilde_beta[[i]] <- out_df
  cat(i, "  ")
}

tilde_beta <- pup_frame %>% dplyr::select(dbid, tilde_beta) %>% 
  unnest(cols = c(tilde_beta))

saveRDS(tilde_beta,file="results/tilde_beta_data.rds")

