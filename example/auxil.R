## (C) Michael Pokojovy (2023)

mu.hat <- function(x, df, estimator = "KWA") {
  est = if (estimator == "usual")
    mean(x)
  else if (estimator == "Qn")
    robustbase::s_Qn(x, mu.too = TRUE)[1]
  else if (estimator == "MCD.50") {
    robustbase::covMcd(x = x, raw.only = TRUE, use.correction = FALSE, alpha = 0.50)$raw.center
  } else if (estimator == "MCD.75") {
    robustbase::covMcd(x = x, raw.only = TRUE, use.correction = FALSE, alpha = 0.75)$raw.center
  } else if (estimator == "KWA")
    KWA1D::KWA1D.loc(x, df)
  else if (estimator == "HL")
    0.5*median(outer(x, x, FUN = "+"))
  
  return(as.numeric(est))
}

sigma.hat <- function(x, df, estimator = "KWA") {
  c.IQR = qt(0.75, df) - qt(0.25, df)
  
  p = 1; bdp = 0.50
  I = integrate(function(u) qf(u, p, df)*p, 0, 1 - bdp)$value
  c.MCD.50 = 1/(I/(1 - bdp)/p)

  p = 1; bdp = 0.25
  I = integrate(function(u) qf(u, p, df)*p, 0, 1 - bdp)$value
  c.MCD.75 = 1/(I/(1 - bdp)/p)
  
  Qn.constant = 1/(2.0^(0.5 + 0.5/df)*qt(5/8, df = df))
  
  est = if (estimator == "usual")
    if (df > 2) sd(x)/sqrt((df - 2)/df) else NaN
  else if (estimator == "Qn")
    robustbase::s_Qn(x, constant = Qn.constant, finite.corr = FALSE, mu.too = FALSE)
  else if (estimator == "MCD.50") {
    mcd = robustbase::covMcd(x = x, raw.only = TRUE, use.correction = FALSE, alpha = 0.50) 
    mcd$raw.cov/mcd$raw.cnp2[1]*c.MCD.50
  } else if (estimator == "MCD.75") {
    mcd = robustbase::covMcd(x = x, raw.only = TRUE, use.correction = FALSE, alpha = 0.75) 
    mcd$raw.cov/mcd$raw.cnp2[1]*c.MCD.75
  } else if (estimator == "KWA")
    KWA1D::KWA1D.scale(x, df)
  else if (estimator == "IQR")
    IQR(x)/c.IQR
  else if (estimator == "HL") {
    log.x = log(abs(x))
    exp(0.5*median(outer(log.x, log.x, FUN = "+")))
  }
  
  return(as.numeric(est))
}