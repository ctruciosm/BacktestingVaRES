################################
##           DGPs             ##
################################

library(quantreg)
VaR_VQR = function(r,VaR, alpha){
  
  fit1 = suppressWarnings(summary(rq(r ~ VaR, tau = alpha, method = "fn"), method="fn" , se="nid" , cov=TRUE))
  
  a1 = fit1$coefficients[1]
  a2 = fit1$coefficients[2]
  
  M = matrix(nrow = 2 , ncol=1)
  M[1,1] = a1
  M[2,1] = (a2-1)
  
  icov = matrix(nrow = 2 , ncol = 2)
  aa = fit1$cov[1,1]
  bb = fit1$cov[1,2]
  cc = fit1$cov[2,1]
  dd = fit1$cov[2,2]
  icov[2,1] = 1/(bb-aa*dd/cc)
  icov[2,2] = 1/(dd-cc*bb/aa)
  icov[1,1] = -icov[2,1]*dd/cc
  icov[1,2] = -icov[2,2]*bb/aa
  
  statistic = (t(M)) %*% icov %*% M 
  # Added by myself, only in case of computational problems
  if(is.na(statistic)){
    fit1 = suppressWarnings(summary(rq(r ~ VaR, tau = alpha, method = "fn"), method="fn" , se="boot" , cov=TRUE))
    
    a1 = fit1$coefficients[1]
    a2 = fit1$coefficients[2]
    
    M = matrix(nrow = 2 , ncol=1)
    M[1,1] = a1
    M[2,1] = (a2-1)
    
    icov = matrix(nrow = 2 , ncol = 2)
    aa = fit1$cov[1,1]
    bb = fit1$cov[1,2]
    cc = fit1$cov[2,1]
    dd = fit1$cov[2,2]
    icov[2,1] = 1/(bb-aa*dd/cc)
    icov[2,2] = 1/(dd-cc*bb/aa)
    icov[1,1] = -icov[2,1]*dd/cc
    icov[1,2] = -icov[2,2]*bb/aa
    
    statistic = (t(M)) %*% icov %*% M 
  }
  
  p.value = 1-pchisq(statistic[1,1], df=2)
  
  return(p.value)
}

simulate_garch <- function(n, omega, alpha, beta, nu, risk_level) {
  nbur <- 500
  ntot <- n + nbur
  s <- sqrt((nu - 2) / nu)
  e <-  rt(ntot, nu) * s
  r <- rep(0, ntot)
  h <- rep(0, ntot)
  var <- rep(0, ntot)
  es <- rep(0, ntot)
  r[1] <- e[1]
  h[1] <- 1
  f = function(x, sigma_) x*rugarch::ddist(distribution = "std", y = x, mu = 0, sigma = sigma_, shape = nu)
  
  for (i in 2:ntot) {
    h[i] <- omega + alpha * r[i - 1]^2 + beta*h[i - 1]
    r[i] <- sqrt(h[i])*e[i]
    var[i] <- qt(risk_level, nu) * sqrt(h[i]) * s
    es[i] <- integrate(f, -Inf, var[i], sigma_ = sqrt(h[i]))$value/risk_level
  }
  out <- data.frame(r = r[-c(1:nbur)], h = h[-c(1:nbur)], var = var[-c(1:nbur)], es = es[-c(1:nbur)])
  return(out)
}



