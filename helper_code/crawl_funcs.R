# fit crawl and draw random imputations from a fitted track
crawl_fit = function(data, mov.model = ~1, err.model = NULL, activity = NULL, 
                     drift = FALSE, coord = c("x", "y"), Time.name,
                     theta, fixPar, method = "Nelder-Mead", control = NULL, 
                     constr = list(lower = -Inf, upper = Inf), prior = NULL, need.hess = TRUE, 
                     initialSANN = list(maxit = 200), 
                     attempts = 1){
  
  # load("~/Desktop/NMML_Phase_2/R_scripts/NAP_Hex_Grid_and_Covariate_Operations/an355.RData")
  
  if(drift){
    initial = list(
      a=c(data[1,coord[1]][[1]],0,0,data[1,coord[2]][[1]],0,0),
      P=diag(c(10000^2,5400^2,5400^2,10000^2,5400^2,5400^2))
    )
  } else {
    initial = list(
      a=c(data[1,coord[1]][[1]],0,data[1,coord[2]][[1]],0),
      P=diag(c(10000^2,5400^2,10000^2,5400^2))
    )
  }
  
  fit1 <- crwMLE(
    mov.model, err.model, activity, 
    drift, data, coord, Time.name, initial.state=initial, 
    theta, fixPar, method, control, constr, prior, need.hess, initialSANN, 
    attempts
  )
  
return(fit1)
}



sim_tracks = function(simObj, reps, CRSobj){
  require(lubridate)
  require(sp)
  samp = vector("list",reps)
  predTime = crawl::intToPOSIX(simObj$Time[simObj$locType=="p"])
  for(i in 1:reps){
    samp[[i]] <- crwPostIS(simObj, fullPost = FALSE)$alpha.sim[simObj$locType=="p",c("mu.x","mu.y")]
  }
  samp=do.call(rbind, samp)
  samp =data.frame(reps=rep(1:reps, each=nrow(samp)/reps), Time=rep(predTime, reps), samp)
  coordinates(samp) = ~mu.x+mu.y
  if(!missing(CRSobj)) proj4string(samp) = CRSobj
  return(samp)
}


