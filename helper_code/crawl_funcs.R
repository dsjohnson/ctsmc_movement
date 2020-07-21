# fit crawl and draw random imputations from a fitted track
crawl_fit = function(data, mov.model = ~1, err.model = NULL, activity = NULL, 
                     drift = FALSE, coord = c("x", "y"), Time.name, time.scale="hours",
                     theta, fixPar, method = "Nelder-Mead", control = NULL, 
                     constr = list(lower = -Inf, upper = Inf), prior = NULL, need.hess = TRUE, 
                     initialSANN = list(maxit = 200), 
                     attempts = 1){
  
  
  fit1 <- crwMLE(
    mov.model=mov.model, err.model=err.model, activity=activity, 
    drift=drift, data=data, coord=coord, Time.name=Time.name, time.scale = time.scale,
    theta=theta, fixPar=fixPar, method=method, control=control, constr=constr, 
    prior=prior, need.hess=need.hess, initialSANN=initialSANN, 
    attempts
  )
  
return(fit1)
}



sim_tracks <- function(simObj, reps, CRSobj){
  require(lubridate)
  require(sp)
  samp <- vector("list",reps)
  if(packageVersion("crawl")=='2.2.2'){
    predTime <- crawl::intToPOSIX(simObj$TimeNum[simObj$locType=="p"]*3600)
  } else{
    predTime <- crawl::intToPOSIX(simObj$TimeNum[simObj$locType=="p"])
  }
  for(i in 1:reps){
    samp[[i]] <- crwPostIS(simObj, fullPost = FALSE)$alpha.sim[simObj$locType=="p",c("mu.x","mu.y")]
  }
  samp <- do.call(rbind, samp)
  samp <- data.frame(reps=rep(1:reps, each=nrow(samp)/reps), Time=rep(predTime, reps), samp)
  coordinates(samp) <- ~mu.x+mu.y
  if(!missing(CRSobj)) proj4string(samp) <- CRSobj
  return(samp)
}


