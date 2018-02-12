make_disc_path = function(sims, rep, hex_poly, quad_times){
  bound = list(min(sims$Time) , max(sims$Time))
  cont_path_time = slot(sims,"data") %>% filter(reps==1) %>% dplyr::select(Time)
  disc_path = data.frame(
    cont_path_time, 
    move_hex_to=over(sims[sims@data$reps==k,], hex_poly)
  ) %>% 
    mutate(move_hex_from = c(NA,move_hex_to[-n()])) %>% 
    full_join(data.frame(Time=quad_times, cov_trans=1), by="Time") %>% dplyr::arrange(Time) %>% 
    mutate(grid_time = ifelse(cov_trans==1, Time, NA) %>% zoo::na.locf() %>% crawl::intToPOSIX()) %>%
    filter(Time>=bound[[1]]) %>% 
    mutate(
      elapsed_time = c(0,cumsum(diff(as.numeric(Time))))/3600,
      move = ifelse(move_hex_from!=move_hex_to, 1, 0) %>% ifelse(is.na(.), 0, .) %>% 
        cumsum() %>% ifelse(move_hex_from==move_hex_to, NA, .),
      tempA = zoo::na.locf(move, fromLast=TRUE, na.rm=FALSE),
      tempB = c(0, diff(elapsed_time))
    ) %>% group_by(tempA) %>% 
    mutate(
      tempC = c(NA, diff(elapsed_time)) %>% ifelse(is.na(.), tempB, .),
      tsm = cumsum(tempB) %>% ifelse(.==0, NA, .)
    ) %>% ungroup() %>% dplyr::select(-contains("temp"))
  disc_path %>% 
    mutate(
      move_hex_to = zoo::na.locf(move_hex_to), 
      move_hex_from = ifelse(is.na(move_hex_from), move_hex_to, move_hex_from),
      delta = c(NA, diff(elapsed_time))
    ) -> disc_path
  
  return(disc_path)
  }