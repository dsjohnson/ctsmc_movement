
source("load_packages.R")
source("helper_code/pmp.R")

pup_frame = readRDS("results/pup_frame_fitted.rds")

pmp_frame = pup_frame %>% dplyr::select(dbid, fit)
pmp_frame = pmp_frame %>% mutate(
  pmp = foreach(i = 1:n()) %do% { 
    out = readRDS(fit[[i]]) %>% 
      dplyr::select(dbid, rep, model, fit) %>% 
      group_by(dbid, rep) %>% tidyr::nest() %>% 
      dplyr::mutate(pmp = purrr::map(data, ~imputation_pmp_gam(.x$fit,1))) %>%
      tidyr::unnest() %>% dplyr::select(-fit)
    out
  }
) %>% dplyr::select(-fit) %>% tidyr::unnest()

pmp_frame = pmp_frame %>% mutate(dbid=factor(dbid), model=factor(model))
pmp_summ = pmp_frame %>% group_by(dbid, model) %>% summarise(pmp_est=mean(pmp)) %>% ungroup()  
pos_mods = pmp_summ %>% group_by(model) %>% summarise(pmp=mean(pmp_est)) %>% filter(pmp>0) %>% pull(model)
pmp_summ = pmp_summ %>% filter(model%in%pos_mods)
pmp_pop = pmp_summ %>% group_by(model) %>% summarise(pmp_pop=mean(pmp_est), range=list(range(pmp_est))) %>% ungroup()
pmp_frame %>% right_join(pmp_summ) -> pmp_frame
ggplot(data=pmp_frame) + geom_point(aes(x=model, y=pmp, color=dbid), alpha=0.1, cex=5, position = position_dodge(width=0.75)) +
  geom_point(data=pmp_summ, aes(x=model, y=pmp_est, color=dbid), alpha=1, cex=5, position = position_dodge(width=0.75)) +
  geom_point(data=pmp_pop, aes(x=model, y=pmp_pop), cex=5) + theme_light() +
  xlab("\nModel") + ylab("Posterior model probability (PMP)\n") +
  theme(
    axis.text = element_text(size=rel(2)),
    axis.title = element_text(size=rel(2))
  )
ggsave("plots/pmp_fig.jpg", dpi="retina", width=6.5, units = "in")




### Summarize model 8 effects
get_plot_data = function(object){
  pdf(NULL)
  res <- plot(object, pages=1, scale=0)
  invisible(dev.off())
  res
}

fit_list = pull(pup_frame, fit)
dbid = pull(pup_frame, dbid)

foreach(i = 1:length(fit_list)) %do% {
  fit_frame = readRDS(fit_list[[i]]) %>% filter(model==8)
  fit_frame %>% mutate(
    plot_data = purrr::map(fit, get_plot_data),
    north_pd = purrr::map(plot_data, ~.x[[1]]),
    east_pd = purrr::map(plot_data, ~.x[[2]]),
    wind_pd = purrr::map(plot_data, ~.x[[3]]),
    curr_pd = purrr::map(plot_data, ~.x[[4]]),
    sst_pd = purrr::map(plot_data, ~.x[[5]])
  ) %>% dplyr::select(-plot_data) -> fit_frame
  
  north_eff = map(fit_frame$north_pd, ~.x$fit) %>% unlist() %>% matrix(nrow(fit_frame), byrow = T)
  north_se = map(fit_frame$north_pd, ~.x$se) %>% unlist() %>% matrix(nrow(fit_frame), byrow = T)
  north_ci = unc_ci(north_eff, north_se, sim=50000) %>% mutate(effect="North") %>% 
    mutate(time=fit_frame$north_pd[[1]]$x)
 
  east_eff = map(fit_frame$east_pd, ~.x$fit) %>% unlist() %>% matrix(nrow(fit_frame), byrow = T)
  east_se = map(fit_frame$east_pd, ~.x$se) %>% unlist() %>% matrix(nrow(fit_frame), byrow = T)
  east_ci = unc_ci(east_eff, east_se, sim=50000) %>% mutate(effect="East")%>% 
    mutate(time=fit_frame$east_pd[[1]]$x)
  
  wind_eff = map(fit_frame$wind_pd, ~.x$fit) %>% unlist() %>% matrix(nrow(fit_frame), byrow = T)
  wind_se = map(fit_frame$wind_pd, ~.x$se) %>% unlist() %>% matrix(nrow(fit_frame), byrow = T)
  wind_ci = unc_ci(wind_eff, wind_se, sim=50000) %>% mutate(effect="Surface wind")%>% 
    mutate(time=fit_frame$wind_pd[[1]]$x)
  
  curr_eff = map(fit_frame$curr_pd, ~.x$fit) %>% unlist() %>% matrix(nrow(fit_frame), byrow = T)
  curr_se = map(fit_frame$curr_pd, ~.x$se) %>% unlist() %>% matrix(nrow(fit_frame), byrow = T)
  curr_ci = unc_ci(curr_eff, curr_se, sim=50000) %>% mutate(effect="Geostrophic current")%>% 
    mutate(time=fit_frame$curr_pd[[1]]$x)
  
  sst_eff = map(fit_frame$sst_pd, ~.x$fit) %>% unlist() %>% matrix(nrow(fit_frame), byrow = T)
  sst_se = map(fit_frame$sst_pd, ~.x$se) %>% unlist() %>% matrix(nrow(fit_frame), byrow = T)
  sst_ci = unc_ci(sst_eff, sst_se, sim=50000) %>% mutate(effect="Sea surface temperature")%>% 
    mutate(time=fit_frame$sst_pd[[1]]$x)
  
  eff_ci = bind_rows(north_ci, east_ci, wind_ci, curr_ci, sst_ci) %>% 
    mutate(dbid = dbid[[i]])
  
  ggplot(data=eff_ci) + 
    geom_path(aes(x=time, y=est), lwd=2) +
    geom_ribbon(aes(ymin=lowCI, ymax=upCI, x=time), alpha=0.2) + 
    geom_abline(slope=0, intercept = 0, color="gray70", lwd=2) + 
    xlab('Time (h)') + ylab('Effect coef.') + facet_wrap(~effect, nrow=3, scales="free")
  
  ggsave(paste0('plots/',dbid[[i]], '_eff_plot.pdf'), width=6.5, units='in', dpi = "retina")
  
  paste0(dbid[[i]], ' complete')
}
