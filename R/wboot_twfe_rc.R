# Bootrstapped TWFE Difference-in-Differences with Repeated Cross Section Data
# 2 periods and 2 groups

wboot_twfe_rc <- function(nn, n, y, dd, post, x, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  #Compute the TWFE Regression
  if(!is.null(x)){
    reg.b <- stats::lm(y ~  dd:post + post + dd + x, weights = b.weights)
    #reg.b <- stats::glm(y ~  dd:post + post + dd + x + offset(n), family = poisson(link = "log))
  } else{
    reg.b <- stats::lm(y ~  dd:post + post + dd, weights = b.weights)
    #reg.b <- stats::glm(y ~  dd:post + post + dd + offset(n), family = poisson(link = "log))
  }
  twfe.att.b <- reg.b$coefficients["dd:post"]
  #-----------------------------------------------------------------------------
  return(twfe.att.b)
}
