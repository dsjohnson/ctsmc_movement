library(tidyverse)
library(mgcv)
# library(cowplot)
library(sf)
library(coda)

# load("pup_frame.RData")

# load("reps20.RData")

unc_ci = function(eff_mat, se_mat, prob=0.95, sim=50000){
  m = colMeans(eff_mat)
  sd = sqrt(colMeans(se_mat^2) + apply(eff_mat, 2, var))
  out_df = data.frame(est=m, se=sd, lowCI=NA, upCI=NA)
  for(p in 1:ncol(eff_mat)){
    out_df[p,3:4] = mapply(function(eff, se){rnorm(sim, eff, se)}, eff=eff_mat[,p], se=se_mat[,p]) %>% 
      as.vector() %>% mcmc() %>% HPDinterval(prob=prob)
  }
  return(out_df)
}



### Covariate effect plots

fits = fit_frame %>% filter(model==8)

wind_eff <- north_eff <- east_eff <- sst_eff <- curr_eff <- matrix(NA, nrow(fits), 100)
wind_se <- north_se <- east_se <- sst_se <- curr_se <- matrix(NA, nrow(fits), 100)

for(k in 1:nrow(fits)){
  plot_data <- {
    pdf(NULL)
    res <- plot(fits$fit[[k]], pages=1, scale=0)
    invisible(dev.off())
    res
  }
  north_eff[k,]= plot_data[[1]]$fit
  north_se[k,]= plot_data[[1]]$se
  east_eff[k,]= plot_data[[2]]$fit
  east_se[k,]= plot_data[[2]]$se
  wind_eff[k,] = plot_data[[3]]$fit
  wind_se[k,] = plot_data[[3]]$se
  sst_eff[k,] = plot_data[[5]]$fit
  sst_se[k,] = plot_data[[5]]$se
  curr_eff[k,] = plot_data[[4]]$fit
  curr_se[k,] = plot_data[[4]]$se
  
}

north_ci = unc_ci(north_eff, north_se)
east_ci = unc_ci(east_eff, east_se)
wind_ci = unc_ci(wind_eff, wind_se)
curr_ci = unc_ci(curr_eff, curr_se)
sst_ci = unc_ci(sst_eff, sst_se)

p_north = ggplot(data=north_ci)+geom_path(aes(x=plot_data[[1]]$x, y=est), lwd=1) + 
  geom_ribbon(aes(ymin=lowCI, ymax=upCI, x=plot_data[[1]]$x), alpha=0.3) + geom_abline(slope=0, intercept = 0, color='grey60') + 
  xlab('Time (h)') + ylab('North')
p_east = ggplot(data=east_ci)+geom_path(aes(x=plot_data[[2]]$x, y=est), lwd=1) + 
  geom_ribbon(aes(ymin=lowCI, ymax=upCI, x=plot_data[[2]]$x), alpha=0.3) + geom_abline(slope=0, intercept = 0, color='grey60') + 
  xlab('Time (h)') + ylab('East')
p_wind = ggplot(data=wind_ci)+geom_path(aes(x=plot_data[[3]]$x, y=est), lwd=1) + 
  geom_ribbon(aes(ymin=lowCI, ymax=upCI, x=plot_data[[3]]$x), alpha=0.3) + geom_abline(slope=0, intercept = 0, color='grey60') + 
  xlab('Time (h)') + ylab('Wind')
p_curr = ggplot(data=curr_ci)+geom_path(aes(x=plot_data[[3]]$x, y=est), lwd=1) + 
  geom_ribbon(aes(ymin=lowCI, ymax=upCI, x=plot_data[[3]]$x), alpha=0.3) + geom_abline(slope=0, intercept = 0, color='grey60') + 
  xlab('Time (h)') + ylab('Current')
p_sst = ggplot(data=sst_ci)+geom_path(aes(x=plot_data[[3]]$x, y=est), lwd=1) + 
  geom_ribbon(aes(ymin=lowCI, ymax=upCI, x=plot_data[[3]]$x), alpha=0.3) + geom_abline(slope=0, intercept = 0, color='grey60') + 
  xlab('Time (h)') + ylab('SST')

p_eff = plot_grid(p_north, p_east, p_wind, p_curr, p_sst,ncol = 2, labels="AUTO")
ggsave("p_eff.pdf", p_eff, width=7, height=10, units = 'in')


fits = fit_frame %>% filter(model==9)

wind_eff <- sst_eff <- curr_eff <- matrix(NA, nrow(fits), 100)
wind_se <- sst_se <- curr_se <- matrix(NA, nrow(fits), 100)

for(k in 1:nrow(fits)){
  plot_data <- {
    pdf(NULL)
    res <- plot(fits$fit[[k]], pages=1, scale=0)
    invisible(dev.off())
    res
  }
  wind_eff[k,]= plot_data[[1]]$fit
  wind_se[k,]= plot_data[[1]]$se
  curr_eff[k,]= plot_data[[2]]$fit
  curr_se[k,]= plot_data[[2]]$se
  sst_eff[k,] = plot_data[[3]]$fit
  sst_se[k,] = plot_data[[3]]$se
}

wind_ci = unc_ci(wind_eff, wind_se)
curr_ci = unc_ci(curr_eff, curr_se)
sst_ci = unc_ci(sst_eff, sst_se)

p_wind = ggplot(data=wind_ci)+geom_path(aes(x=plot_data[[1]]$x, y=est), lwd=1) + 
  geom_ribbon(aes(ymin=lowCI, ymax=upCI, x=plot_data[[1]]$x), alpha=0.3) + geom_abline(slope=0, intercept = 0, color='grey60') + 
  xlab('Time (h)') + ylab('WIND')
p_curr = ggplot(data=curr_ci)+geom_path(aes(x=plot_data[[2]]$x, y=est), lwd=1) + 
  geom_ribbon(aes(ymin=lowCI, ymax=upCI, x=plot_data[[2]]$x), alpha=0.3) + geom_abline(slope=0, intercept = 0, color='grey60') + 
  xlab('Time (h)') + ylab('CURR')
p_sst = ggplot(data=sst_ci)+geom_path(aes(x=plot_data[[3]]$x, y=est), lwd=1) + 
  geom_ribbon(aes(ymin=lowCI, ymax=upCI, x=plot_data[[3]]$x), alpha=0.3) + geom_abline(slope=0, intercept = 0, color='grey60') + 
  xlab('Time (h)') + ylab('SST')

p_eff_maxPMP = plot_grid(p_wind, p_curr, p_sst,ncol = 2, labels="AUTO")
ggsave("p_eff_maxPMP.pdf", p_eff_maxPMP, width=7, height=10, units = 'in')



### Baseline intensity function
lambda0 <- lambda0_se <- matrix(NA, nrow(fits), 1)
for(k in 1:nrow(fits)){
  lambda0[k] = coef(fits$fit[[k]])[2]
  lambda0_se[k] = vcov(fits$fit[[k]])[2,2]
}

lambda_df = unc_ci(lambda0, lambda0_se)

prev <- prev_se <- matrix(NA, nrow(fits), 1)
for(k in 1:nrow(fits)){
  prev[k] = coef(fits$fit[[k]])[3]
  prev_se[k] = vcov(fits$fit[[k]])[3,3]
}

prev_df = unc_ci(prev, prev_se)




