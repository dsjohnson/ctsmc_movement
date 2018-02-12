
source("load_packages.R")
load("data_products/pup_frame.RData")
source("helper_code/pmp.R")

### Look at the PMP
pmp = readRDS("results/pmp.rds")
nms = c("rep",paste0("M",1:8))
pmp = pmp %>%  tibble(dbid = names(.), pmp=.) %>% 
  mutate(
    pmp=map(pmp,~cbind(1:20,.x)) %>% map(function(x){colnames(x)=nms; return(x)}) %>% 
      map(as_tibble)
  ) %>% unnest() %>% gather(key=model, pmp, -dbid, -rep)
pmp_summ = pmp %>% group_by(dbid, model) %>% summarise(pmp_est=mean(pmp)) %>% ungroup()  
pos_mods = pmp_summ %>% group_by(model) %>% summarise(pmp=mean(pmp_est)) %>% filter(pmp>0) %>% pull(model)
pmp_summ = pmp_summ %>% filter(model%in%pos_mods)
pmp_pop = pmp_summ %>% group_by(model) %>% summarise(pmp_pop=mean(pmp_est), range=list(range(pmp_est))) %>% ungroup()
pmp %>% right_join(pmp_summ) -> pmp
ggplot(data=pmp) + geom_point(aes(x=model, y=pmp, color=dbid), alpha=0.1, cex=5, position = position_dodge(width=0.75)) +
  geom_point(data=pmp_summ, aes(x=model, y=pmp_est, color=dbid), alpha=1, cex=5, position = position_dodge(width=0.75)) +
  geom_point(data=pmp_pop, aes(x=model, y=pmp_pop), cex=5) + theme_light() +
  xlab("\nModel") + ylab("Posterior model probability (PMP)\n") +
  theme(
    axis.text = element_text(size=rel(2)),
    axis.title = element_text(size=rel(2))
  )
ggsave("plots/pmp_fig.jpg", dpi="retina", width=7, height=7, units = "in")




### Summarize model 6 effects

fit_list = list.files("results") %>% .[grepl("fit_frame",.)]
# fit_list = fit_list[1:5]
wind_eff <- north_eff <- east_eff <- curr_eff <- matrix(NA, nrow(fit_frame), 100)
wind_se <- north_se <- east_se <- curr_se <- matrix(NA, nrow(fit_frame), 100)

north_ci <- east_ci <- wind_ci <- curr_ci <- NULL
for(i in 1:length(fit_list)){
  fit_frame = readRDS(paste0("results/",fit_list[[i]])) %>% filter(model==6)
  for(k in 1:nrow(fit_frame)){
    plot_data <- {
      pdf(NULL)
      res <- plot(fit_frame$fit[[k]], pages=1, scale=0)
      invisible(dev.off())
      res
    }
    north_eff[k,]= plot_data[[1]]$fit
    north_se[k,]= plot_data[[1]]$se
    east_eff[k,]= plot_data[[2]]$fit
    east_se[k,]= plot_data[[2]]$se
    wind_eff[k,] = plot_data[[3]]$fit
    wind_se[k,] = plot_data[[3]]$se
    curr_eff[k,] = plot_data[[4]]$fit
    curr_se[k,] = plot_data[[4]]$se
    
  }
  unc_ci(north_eff, north_se, sim=10000) %>% cbind(dbid=fit_frame$dbid[[1]],time=plot_data[[1]]$x,.) %>% rbind(north_ci,.) -> north_ci
  unc_ci(east_eff, east_se, sim=10000) %>% cbind(dbid=fit_frame$dbid[[1]],time=plot_data[[1]]$x,.) %>% rbind(east_ci,.) -> east_ci
  unc_ci(wind_eff, wind_se, sim=10000) %>% cbind(dbid=fit_frame$dbid[[1]],time=plot_data[[1]]$x,.) %>% rbind(wind_ci,.) -> wind_ci
  unc_ci(curr_eff, curr_se, sim=10000) %>% cbind(dbid=fit_frame$dbid[[1]],time=plot_data[[1]]$x,.) %>% rbind(curr_ci,.) -> curr_ci
  
  cat(fit_frame$dbid[[1]], "  ")
}
north_ci %>% mutate(dbid=factor(dbid)) -> north_ci
east_ci %>% mutate(dbid=factor(dbid)) -> east_ci
wind_ci %>% mutate(dbid=factor(dbid)) -> wind_ci
curr_ci %>% mutate(dbid=factor(dbid)) -> curr_ci

save(list=c("north_ci","east_ci","wind_ci","curr_ci"), file="results/ci_data.RData")

library(cowplot)
for(i in 1:length(fit_list)){
  p_north = ggplot(data=north_ci %>% filter(dbid==levels(dbid)[i])) + 
    geom_path(aes(x=time, y=est), lwd=2) +
    geom_ribbon(aes(ymin=lowCI, ymax=upCI, x=time), alpha=0.2) + 
    geom_abline(slope=0, intercept = 0, color="gray70", lwd=2) + 
    xlab('Time (h)') + ylab(expression(gamma*" North"))
  
  p_east = ggplot(data=east_ci %>% filter(dbid==levels(dbid)[i])) + 
    geom_path(aes(x=time, y=est), lwd=2) + 
    geom_ribbon(aes(ymin=lowCI, ymax=upCI, x=time), alpha=0.2) + 
    geom_abline(slope=0, intercept = 0, color="gray70", lwd=2) + 
    xlab('Time (h)') + ylab(expression(gamma*' East'))
  
  p_wind = ggplot(data=wind_ci %>% filter(dbid==levels(dbid)[i])) + 
    geom_path(aes(x=time, y=est), lwd=2) + 
    geom_ribbon(aes(ymin=lowCI, ymax=upCI, x=time), alpha=0.2) + 
    geom_abline(slope=0, intercept = 0, color="gray70", lwd=2) + 
    xlab('Time (h)') + ylab(expression(gamma*' Wind'))
  
  p_curr = ggplot(data=curr_ci %>% filter(dbid==levels(dbid)[i])) + 
    geom_path(aes(x=time, y=est), lwd=2) + 
    geom_ribbon(aes(ymin=lowCI, ymax=upCI, x=time), alpha=0.2) + 
    geom_abline(slope=0, intercept = 0, color="gray70", lwd=2) + 
    xlab('Time (h)') + ylab(expression(gamma*' Current'))
  
  p_eff = plot_grid(p_north, p_east, p_wind, p_curr,ncol = 2, labels="AUTO")
  ggsave(paste0("plots/",dbid[i],"_p_eff.pdf"), p_eff, width=7, height=5, units = 'in')
}

p_north = ggplot(data=north_ci) + 
  geom_abline(slope=0, intercept = 0, color=1, lty="dashed") + 
  geom_path(aes(x=time, y=est, color=dbid)) +
  xlab('Time (h)') + ylab(expression(gamma*" North")) + theme(legend.position="none")

p_east = ggplot(data=east_ci) + 
  geom_abline(slope=0, intercept = 0, color=1, lty="dashed") + 
  geom_path(aes(x=time, y=est, color=dbid)) +
  xlab('Time (h)') + ylab(expression(gamma*" East"))+ theme(legend.position="none")

p_wind = ggplot(data=wind_ci) + 
  geom_abline(slope=0, intercept = 0, color=1, lty="dashed") + 
  geom_path(aes(x=time, y=est, color=dbid)) +
  xlab('Time (h)') + ylab(expression(gamma*" Wind"))+ theme(legend.position="none")

p_curr = ggplot(data=curr_ci) + 
  geom_abline(slope=0, intercept = 0, color=1, lty="dashed") + 
  geom_path(aes(x=time, y=est, color=dbid)) +
  xlab('Time (h)') + ylab(expression(gamma*" Current"))+ theme(legend.position="none")

p_eff = plot_grid(p_north, p_east, p_wind, p_curr,ncol = 2, labels="AUTO")
ggsave(paste0("plots/pop_p_eff.pdf"), p_eff, width=7, height=7, units = 'in')



