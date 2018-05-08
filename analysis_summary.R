
source("load_packages.R")
source("helper_code/pmp.R")

pup_frame = readRDS("results/pup_frame_fitted.rds") 

pmp_frame = pup_frame %>% dplyr::select(dbid, fit) 
pmp_frame = pmp_frame %>% mutate(
  pmp = foreach(i = 1:n()) %do% { 
    message(dbid[i], "...")
    readRDS(fit[[i]]) %>% 
      dplyr::select(dbid, rep, model, fit) %>% 
      group_by(dbid, rep) %>% tidyr::nest() %>% 
      dplyr::mutate(
        pmp = purrr::map(data, ~calc_pmp(.x$fit,"laplace")),
        pmp_bic = purrr::map(data, ~calc_pmp(.x$fit,"bic"))
        ) %>%
      tidyr::unnest() %>%
      dplyr::select(-fit)
  }
) %>% dplyr::select(-fit) %>% tidyr::unnest()

## Summarize for plotting
pmp_frame = pmp_frame %>% mutate(dbid=factor(dbid), model=factor(model))
pmp_summ = pmp_frame %>% group_by(dbid, model) %>% summarise(pmp_est=mean(pmp), pmp_est_bic=mean(pmp_bic)) %>% ungroup()  
# pos_mods = pmp_summ %>% group_by(model) %>% summarise(pmp=mean(pmp_est)) %>% filter(pmp>0) %>% pull(model)
# pmp_summ = pmp_summ %>% filter(model%in%pos_mods)
pmp_pop = pmp_summ %>% group_by(model) %>% 
  summarise(
    pmp_pop=mean(pmp_est), range=list(range(pmp_est)),
    pmp_pop_bic=mean(pmp_est_bic), range_bic=list(range(pmp_est_bic))
    ) %>% ungroup()
pmp_frame %>% right_join(pmp_summ) -> pmp_frame

saveRDS(pmp_frame, "results/pmp_data.rds")

ggplot(data=pmp_frame) + geom_point(aes(x=model, y=pmp, color=dbid), alpha=0.1, cex=5, position = position_dodge(width=0.75)) +
  geom_point(data=pmp_summ, aes(x=model, y=pmp_est, color=dbid), alpha=1, cex=5, position = position_dodge(width=0.75)) +
  geom_point(data=pmp_pop, aes(x=model, y=pmp_pop), cex=5) + #theme_light() +
  xlab("\nModel") + ylab("Posterior model probability (PMP)\n") +
  theme(
    axis.text = element_text(size=rel(2)),
    axis.title = element_text(size=rel(2))
  )
ggsave("plots/pmp_fig.jpg", dpi="retina", width=12, height=8, units = "in")

### Make table all animals
library(xtable)

pmp_frame %>% dplyr::select(dbid, model, pmp_est) %>% distinct() %>% 
  spread(model, pmp_est) %>% xtable()

pmp_frame %>% dplyr::select(dbid, model, pmp_est_bic) %>% distinct() %>% 
  spread(model, pmp_est_bic) %>% xtable()

### table for animal 355
pmp_frame %>% filter(dbid=='355') %>% dplyr::select(model, rep, pmp) -> pmp_355
pmp_355 %>% group_by(model) %>% summarise(pmp_est=mean(pmp), min=min(pmp), max=max(pmp))


foreach(i = 1:nrow(pup_frame)) %do% {
  fit_frame = readRDS(pup_frame$fit[[i]]) %>% filter(model==8)
  fit_frame %>% mutate(
    b = map(fit, coef),
    V = map(fit, vcov),
    ps = map2(b, V, ~rmvnorm(5000, .x, .y)) %>% map(., as.data.frame)
  ) %>% dplyr::select(ps) %>% unnest() -> ps
  dt = difftime(pup_frame$last_obs[i],pup_frame$first_obs[i], units = "hours")
  time = seq(0,dt, length=100)
  df = 3
  kts = seq(0, dt, length=df)
  kts = c(kts[1]-(kts[2]-kts[1]), kts, kts[df]+(kts[df]-kts[df-1])) 
  R = rbf_matrix(time, kts, min(diff(kts)))
  plot_data = NULL
  x_nms=c('north','east','wind','geo_curr','sst')
  eff_nms = c('North','East','Wind','Geostrophic current','SST')
  for(e in 1:5){
    idx =  colnames(ps) %>% grep(x_nms[e],.)         
    b = ps[,idx]
    ps_eff =R%*%t(b) %>% t()
    est = tibble(est=colMeans(ps_eff))
    est = ps_eff %>% mcmc() %>% HPDinterval() %>% as.tibble() %>% bind_cols(est,.) %>% 
      mutate(effect=eff_nms[e], time=time)
    plot_data = bind_rows(plot_data, est)
  }

  ggplot(data=plot_data)+
    geom_abline(slope=0, intercept = 0, color="gray70", lwd=1) + 
    geom_path(aes(x=time, y=est)) +
    geom_ribbon(aes(x=time, ymin=lower, ymax=upper), alpha=0.2)+ 
    xlab('Time (h)') + ylab('Effect coef.') + 
    facet_wrap(~effect, nrow=3, scales="free")

  ggsave(paste0('plots/',pup_frame$dbid[[i]], '_eff_plot.pdf'), width=6.5, units='in', dpi = "retina")
  
  paste0(pup_frame$dbid[[i]], ' complete')
}
