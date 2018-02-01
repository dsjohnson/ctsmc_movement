vector2scalar = function(u, v, bearing, rotation_deg=NULL){
  unit_xy = cbind(cos(pi*(90-bearing)/180), sin(pi*(90-bearing)/180))
  if(is.null(rotation_deg)) uv = cbind(u,v)
  else uv = rotate_vector(u, v, rotation_deg)
  scalar = rowSums(unit_xy * uv)
  return(scalar)
}

rotate_vector = function(u,v, deg){
  theta = deg*pi/180
  R = matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  out = t(t(R) %*% t(cbind(u,v)))
  return(out)
}