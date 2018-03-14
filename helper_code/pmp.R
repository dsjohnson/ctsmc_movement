imputation_pmp_gam = function(model_list, sp.scale=10){
  out = rep(NA, length(model_list))
  for(i in 1:length(model_list)){
    lambda = (model_list[[i]]$sp * unlist(lapply(model_list[[i]]$smooth, "[[", "S.scale"))) 
    out[i] = 
      as.numeric(
        logLik(model_list[[i]]) + ln_prior_gam(model_list[[i]]) +
          0.5*determinant(vcov(model_list[[i]]), unconditional=T, logarithm = T)$modulus +
          0.5*determinant(sp.vcov(model_list[[i]]), logarithm = T)$modulus +
          0.5*(length(coef(model_list[[i]])) + length(lambda))*log(2*pi) + 
          sum(dnorm(1/sqrt(lambda), 0, sp.scale, log=T) - (3/2)*log(lambda) -log(2))
      )
  }
  out = exp(out-max(out))
  out = round(out/sum(out),2)
  return(out)
}

ln_prior_gam = function(object){
  lambda = (object$sp * unlist(lapply(object$smooth, "[[", "S.scale"))) 
  if(length(lambda)==0) return(0)
  S = lapply(object$smooth, "[[", "S")
  if(length(lambda)!=2*length(S)) stop("In gam fitting: must have 'select=TRUE' for evaluating PMP!")
  lambda = lambda  %>% split(ceiling(seq_along(.)/2))
  S = S %>% map2(lambda, ~.x[[1]]*.y[[1]] + .x[[2]]*.y[[2]]) %>% .bdiag(.)
  E = eigen(S, only.values = T)$values
  M = sum(E!=0)
  ln_det_S = sum(log(E[E!=0]))
  b_sm = tail(coef(object), nrow(S))
  return(as.numeric(-(t(b_sm) %*% S %*%  b_sm)/2 + ln_det_S/2 - M*log(2*pi)/2))
}

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
