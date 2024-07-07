# (C) Michael Pokojovy (2024)

do.simulation <- function(n, df = 1, nrep = 500, type = "loc", estimator = "KWA", 
                          logarithm.flag = TRUE, parallel.flag = TRUE) {
  
  c.IQR = qcauchy(0.75) - qcauchy(0.25)
  
  c.MCD.50 = {alpha = 0.50; a = qcauchy(0.5*(1 - alpha)); b = qcauchy(1 - 0.5*(1 - alpha));
              (b - a)/pi/alpha - 1}
  
  c.MCD.75 = {alpha = 0.75; a = qcauchy(0.5*(1 - alpha)); b = qcauchy(1 - 0.5*(1 - alpha));
              (b - a)/pi/alpha - 1}
  
  iter <- function(n, type = "loc", estimator = "KWA", logarithm.flag = TRUE) {
    x = rt(n, df)
    
    if (type == "loc") {
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
    } else {
      est = if (estimator == "usual")
        sd(x)
      else if (estimator == "Qn")
        robustbase::s_Qn(x, constant = 1.2071, finite.corr = FALSE, mu.too = FALSE)
      else if (estimator == "MCD.50") {
        mcd = robustbase::covMcd(x = x, raw.only = TRUE, use.correction = FALSE, alpha = 0.50) 
        mcd$raw.cov/mcd$raw.cnp2[1]/c.MCD.50
      } else if (estimator == "MCD.75") {
        mcd = robustbase::covMcd(x = x, raw.only = TRUE, use.correction = FALSE, alpha = 0.75) 
        mcd$raw.cov/mcd$raw.cnp2[1]/c.MCD.75
      } else if (estimator == "KWA")
        KWA1D::KWA1D.scale(x, df)
      else if (estimator == "IQR")
        IQR(x)/c.IQR
      else if (estimator == "HL") {
        log.x = log(abs(x))
        exp(0.5*median(outer(log.x, log.x, FUN = "+")))
      }
    }
    
    est = abs(est)
    if (abs(est) == 0) {
      est = 1E-9
    }
    
    if ((type == "scale") && (logarithm.flag)) {
      est = log(est)
    }
      
    return(c(est, est^2, est^4)/nrep)
  }
  
  res = if (parallel.flag) {
    foreach(i = 1:nrep, .combine = "+", .inorder = FALSE) %dopar% {
      err = iter(n, type, estimator, logarithm.flag)
    }
  } else {
    foreach(i = 1:nrep, .combine = "+", .inorder = FALSE) %do% {
      err = iter(n, type, estimator, logarithm.flag)
    }
  }
  
  #res = c(res[2], max(0.0, res[2] - nrep/(nrep - 1)*res[1]^2), sqrt(max(0.0, res[3] - nrep/(nrep - 1)*res[2]^2)))
  res = c(res[2], max(0.0, res[2] - nrep/(nrep - 1)*res[1]^2), sqrt(max(0.0, res[3] - nrep/(nrep - 1)*res[2]^2)))
      
  return(res)
}