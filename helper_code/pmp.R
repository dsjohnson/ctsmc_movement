imputation_pmp_gam = function(model_list){
  out = rep(NA, length(model_list))
  for(i in 1:length(model_list)){
    out[i] = 
      as.numeric(
        -2*logLik(model_list[[i]]) -
          determinant(vcov(model_list[[i]]), logarithm = T)$modulus -
          length(coef(model_list[[i]]))*log(2*pi) -
          2*ln_prior_gam(model_list[[i]])
        )
  }
  out = exp(-0.5*(out-min(out)))
  out = round(out/sum(out),2)
  return(out)
}

ln_prior_gam = function(object){
  # return(0)
  lambda = object$sp * unlist(lapply(object$smooth, "[[", "S.scale"))
  if(length(lambda)==0){
    return(0)
  } else{
    S = lapply(lapply(object$smooth, "[[", "S"), "[[", 1) %>% 
      mapply('*', lambda, ., SIMPLIFY = F) %>% .bdiag(.)
    E = eigen(S, only.values = T)$values
    M = sum(E!=0)
    ln_det_S = sum(log(E[E!=0]))
    b_sm = tail(coef(object), nrow(S))
    return(as.numeric(-(t(b_sm) %*% S %*%  b_sm)/2 + ln_det_S/2 - M*log(2*pi)/2))
  }
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
