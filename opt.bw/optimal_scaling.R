# (C) Michael Pokojovy (2024)

Rcpp::sourceCpp("KWA_1D.cpp")
source("auxil.R")
    
set.seed(0)

nrep = 100000L

n.grid = c(30, 50, 100, 200)
h.loc.opt   = rep(0.0, length(n.grid))
h.scale.opt = rep(0.0, length(n.grid))
  
h.grid.loc   = seq(0.10, 1.0, length.out = 100L)
h.grid.scale = seq(0.05, 0.5, length.out = 100L)

for (i.n in 1:length(n.grid)) {
  ptm = proc.time()
  
  n = n.grid[i.n]
  cat("n = ", n, ": computing...\n", sep = "")
  
  par(mfrow = c(1, 2))
  
  # mu
  MSE.loc = rep(0.0, length(h.grid.loc))
  
  for (i.h in 1:length(h.grid.loc)) {
    MSE.loc[i.h] = mean(replicate(nrep, {x = rcauchy(n); c = Qn_scale(x); KWA_loc_(x, c*h.grid.loc[i.h])^2}))
  }
  
  plot(h.grid.loc, MSE.loc, type = "p", xlab = "Bandwidth h", ylab = "MSE", main = "KWA location")
  
  fit = cobs::conreg(h.grid.loc, MSE.loc, convex = TRUE)
  h.out = seq(min(h.grid.loc), max(h.grid.loc), length.out = 10000L)
  fit = smooth.spline(h.grid.loc[fit$iKnots], fit$yf[fit$iKnots], spar = 0.6)
  fit = predict(fit, h.out)
  lines(fit$x, fit$y, col = "red")
  
  h.loc.opt[i.n] = fit$x[which.min(fit$y)]
  
  cat("KWA location: h_opt = ", h.loc.opt[i.n], "\n", sep = "")
  cat("KWA location: objective = ", min(fit$y), "\n", sep = "")
  
  # sigma
  MSE.scale = rep(0.0, length(h.grid.scale))
  
  for (i.h in 1:length(h.grid.scale)) {
    MSE.scale[i.h] = mean(replicate(nrep, {x = rcauchy(n); c = Qn_scale(x); log(KWA_scale_(x, c*h.grid.scale[i.h]))^2}))
  }
  
  plot(h.grid.scale, MSE.scale, type = "p", xlab = "Bandwidth h", ylab = "MSE", main = "KWA scale")
  
  fit = cobs::conreg(h.grid.scale, MSE.scale, convex = TRUE)
  h.out = seq(min(h.grid.scale), max(h.grid.scale), length.out = 10000L)
  fit = smooth.spline(h.grid.scale[fit$iKnots], fit$yf[fit$iKnots], spar = 0.4)
  fit = predict(fit, h.out)
  lines(fit$x, fit$y, col = "red")
  
  h.scale.opt[i.n] = fit$x[which.min(fit$y)]
  
  cat("KWA scale: h_opt = ", h.scale.opt[i.n], "\n", sep = "")
  cat("KWA scle: objective = ", min(fit$y), "\n", sep = "")
  
  print(proc.time() - ptm)
  cat("\n")
}

cat("\n\n")

# Optimal scaling
par(mfrow = c(1, 2))

par = lm(log(h.loc.opt) ~ log(n.grid))$coefficients
plot(log(n.grid), log(h.loc.opt),   type = "p", xlab = "log(n)", ylab = "log(h_opt)", 
     main = "KWA location", sub = paste("C_hat = ", zapsmall(exp(par[1]), 4), ", alpha_hat = ", zapsmall(par[2], 4), sep = ""))
abline(a = par[1], b = par[2], col = "red")
cat("Optimal bandwidth for KWA location: h = ", zapsmall(exp(par[1]), 4), "*n^(", zapsmall(par[2], 4), ")*sigma_hat_Qn\n", sep = "")

par = lm(log(h.scale.opt) ~ log(n.grid))$coefficients
plot(log(n.grid), log(h.scale.opt), type = "p", xlab = "log(n)", ylab = "log(h_opt)",
     main = "KWA scale", sub = paste("C_hat = ", zapsmall(exp(par[1]), 4), ", alpha_hat = ", zapsmall(par[2], 4), sep = ""))
abline(a = par[1], b = par[2], col = "red")
cat("Optimal bandwidth for KWA scale: h = ", zapsmall(exp(par[1]), 4), "*n^(", zapsmall(par[2], 4), ")*sigma_hat_Qn\n", sep = "")