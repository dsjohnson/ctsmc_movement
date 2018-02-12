make_glm_data = function(disc_path, neighbor_df, hex_cov_df){
  
  glm_data = left_join(disc_path, neighbor_df, by=c("move_hex_from"="hex_from")) %>% 
    mutate(z = ifelse(!is.na(move) & hex_to == move_hex_to, 1, 0)) %>% 
    left_join(hex_cov_df, by=c("move_hex_from"="hex", "grid_time"="Time")) 
  
  zorig = glm_data$z[nrow(glm_data)]
  if(!zorig){
    glm_data$z[nrow(glm_data)] = 1
    glm_data$move[(nrow(glm_data)-6+1):nrow(glm_data)] = max(glm_data$move, na.rm=T)+1
  }
  prev_move_data = glm_data %>% filter(z==1) %>% dplyr::select(move, bearing_to_next) %>% 
    mutate(prev_bearing = c(NA, bearing_to_next[-n()])) %>% dplyr::select(-bearing_to_next)
  
  glm_data %>% left_join(prev_move_data, by=c("move")) %>% 
    mutate(
      prev_bearing = zoo::na.locf(prev_bearing, fromLast=TRUE, na.rm=FALSE)
      )-> glm_data
  if(!zorig){
    glm_data$z[nrow(glm_data)] = 0
    glm_data$move[(nrow(glm_data)-6+1):nrow(glm_data)] = NA
  }
  t1 = glm_data$Time[glm_data$z==1 & glm_data$move==1]
  glm_data %>% mutate(
    prev_bearing=ifelse(Time<=t1, NA, prev_bearing),
    prev_move_u = cos(pi*(90-prev_bearing)/180),
    prev_move_v = sin(pi*(90-prev_bearing)/180),
    prev_move = vector2scalar(prev_move_u, prev_move_v, bearing_to_next)
    ) -> glm_data
  glm_data$prev_move = ifelse(is.na(glm_data$prev_move), 0, glm_data$prev_move)
  return(glm_data)
}