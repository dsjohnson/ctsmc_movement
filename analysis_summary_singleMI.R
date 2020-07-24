
source("load_packages.R")
source("helper_code/pmp.R")
source("helper_code/fit_ctsmc.R")

pup_frame <- readRDS("results/pup_frame_fitted.rds") %>% ungroup() %>% 
  filter(dbid=='355')

dt = difftime(pup_frame$last_obs[i],pup_frame$first_obs[i], units = "hours")
time = seq(0,dt, length=100)
df = 3
kts = seq(0, dt, length=df)
kts = c(kts[1]-(kts[2]-kts[1]), kts, kts[df]+(kts[df]-kts[df-1]))
R = rbf_matrix(time, kts, min(diff(kts)))
x_nms=c('north','east','wind','geo_curr','sst')
eff_nms = c('North','East','Wind','Geostrophic current','SST')

fit_frame <- readRDS(pup_frame$fit[[1]]) %>% filter(model==8) %>%
  mutate(
    b = map(fit, coef),
    V = map(fit, vcov, unconditional=T),
    ps = map2(b, V, ~rmvnorm(5000, .x, .y)) %>% map(., as.data.frame)
  )

fit_frame <- fit_frame %>% 
  mutate(
    plot_data = map2(b, V, ~{
      plot_data = NULL
      for(e in 1:length(x_nms)){
        idx =  names(.x) %>% grep(x_nms[e],.)         
        b = .x[idx]
        est = as.vector(R%*%b) 
        V = R%*%.y[idx,idx]%*%t(R)
        out <- tibble(effect=eff_nms[e], time=time, est=est, v=diag(V), 
               lower=as.vector(est-1.96*sqrt(diag(V))), 
               upper=as.vector(est+1.96*sqrt(diag(V)))
        )
        plot_data <- bind_rows(plot_data, out)
      }
      plot_data
    })
  )

plot_data <- fit_frame %>% dplyr::select(rep, plot_data) %>% 
  unnest(cols="plot_data") %>% ungroup()
avg_plot_data <- plot_data %>% group_by(effect, time) %>% 
  summarise(
    est_mi = mean(est),
    vb = var(est),
    ev = mean(v)
  ) %>% 
  mutate(
    sd = sqrt(vb+ev),
    lower_mi = est_mi-1.96*sd,
    upper_mi = est_mi+1.96*sd
  ) %>% dplyr::select(-vb, -ev)



ggplot(data=plot_data)+
  geom_abline(slope=0, intercept = 0, color="gray70", lwd=1) + 
  geom_path(aes(x=time, y=est, group=rep)) +
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper, group=rep), alpha=0.07)+
  xlab('Time (h)') + ylab('Effect coefficient') + 
  facet_wrap(~effect, nrow=3, scales="free") +
  geom_path(aes(x=time, y=est_mi), color='red', data=avg_plot_data) +
  geom_path(aes(x=time, y=upper_mi), color='red', data=avg_plot_data) +
  geom_path(aes(x=time, y=lower_mi), color='red', data=avg_plot_data) +
  theme_classic()



ggsave(paste0('plots/',pup_frame$dbid[[1]], '_eff_MI_plot.pdf'), width=6.5, height=7.7, units='in', dpi = "retina")

paste0(pup_frame$dbid[[i]], ' complete')
