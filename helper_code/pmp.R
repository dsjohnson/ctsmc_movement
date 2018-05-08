calc_pmp = function(model_list, method="laplace"){
  out = rep(NA, length(model_list))
  for(p in 1:length(model_list)){
    if(method=="laplace"){
    S = model_list[[p]]$sp * model_list[[p]]$paraPen$S[[1]]
    out[p] =
      as.numeric(
        logLik(model_list[[p]]) +
        mvtnorm::dmvnorm(coef(model_list[[p]]), sigma=solve(S), log=T) +
        0.5*determinant(vcov.gam(model_list[[p]], unconditional = T), logarithm = T)$modulus +
        0.5*length(coef(model_list[[p]]))*log(2*pi)
      )
    } else if(method=='bic') {
      out[p] = logLik(model_list[[p]]) - 0.5*log(sum(model_list[[p]]$y))*length(coef(model_list[[p]]))
    } else{ stop("Method not recognized")}
  }
  out = exp(out-max(out))
  out = out/sum(out) %>% round(2)
  return(out)
}

# ln_prior_gam = function(object){
#   lambda = (object$sp * unlist(lapply(object$smooth, "[[", "S.scale")))
#   if(length(lambda)==0) return(0)
#   S = lapply(object$smooth, "[[", "S")
#   if(length(lambda)!=2*length(S)) stop("In gam fitting: must have 'select=TRUE' for evaluating PMP!")
#   lambda = lambda  %>% split(ceiling(seq_along(.)/2))
#   S = S %>% map2(lambda, ~.x[[1]]*.y[[1]] + .x[[2]]*.y[[2]]) %>% .bdiag(.)
#   D = S %>% solve() %>% diag() %>% sqrt()
#   lambda_tilde = D/sqrt(unlist(lambda))
#   E = eigen(S, only.values = T)$values
#   M = sum(E!=0)
#   ln_det_S = sum(log(E[E!=0]))
#   b_sm = tail(coef(object), nrow(S))
#   return(
#     list(
#       lambda_tilde=lambda_tilde,
#       R = S %>% solve() %>% cov2cor(),
#       dens = as.numeric(-(t(b_sm) %*% S %*%  b_sm)/2 + ln_det_S/2 - M*log(2*pi)/2)
#     )
#   )
# }

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
