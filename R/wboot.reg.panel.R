# Bootrstapped Regression-based Robust Difference-in-Differences with Panel Data
# 2 periods and 2 groups

wboot.reg.panel <- function(nn, n, deltaY, D, offset, int.cov, i.weights){
  #-----------------------------------------------------------------------------
  v <- stats::rexp(n)
  #v <- v / mean(v)
  #weights for the bootstrap
  b.weights <- as.vector(i.weights * v)
  #Compute the Outcome regression for the control group
  # reg.coeff.b <- stats::coef(stats::lm(deltaY ~ -1 + int.cov,
  #                                      subset = D==0,
  #                                      weights = b.weights))
  control_filter <- (D == 0)
  reg.coeff.b <- stats::coef(fastglm::fastglm(
                              x = int.cov[control_filter, , drop = FALSE],
                              y = deltaY[control_filter],
                              #weights = b.weights[control_filter],
                              offset = n, 
                                  family = poisson(link = "log")
  ))
  out.reg.b <- as.vector(tcrossprod(reg.coeff.b, int.cov))
  # Compute OR estimator
  att.b <- mean(b.weights * D * (deltaY - out.reg.b))/mean(b.weights * D)
  #-----------------------------------------------------------------------------
  return(att.b)
}
