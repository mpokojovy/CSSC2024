# (C) Michael Pokojovy (2024)

source("auxil.R")

load("opt.bw.loc.RData")

h.loc.opt.mat = h.opt.mat

set.seed(1)

parallel.flag = TRUE
save.flag = TRUE
plot.flag = TRUE

library(doParallel)
library(robustbase)

if (parallel.flag) {
  nclust = max(detectCores() - 4L, 1L)
  HPC = makeCluster(nclust)
  registerDoParallel(HPC)
}

df.grid = c(1:5, 10, 20, 30)
n.grid = c(30, 50, 100, 150, 200, 300)

h.opt.mat = matrix(0.0, nrow = length(df.grid), ncol = length(n.grid))
rownames(h.opt.mat) = sapply(df.grid, function(df) paste("df=", df, sep = ""))
colnames(h.opt.mat) = sapply(n.grid,  function(n)  paste("n=",   n, sep = ""))

batch = 500000L
#batch = 500L

h.grid  = signif(c(exp(seq(from = log(0.05), to = log(1.00), length.out = 120L)), 
                   seq(from = 1.05, to = 15.0, by = 0.05)), digits = 4)

ptm = proc.time()

for (i.df in 1:length(df.grid)) {
  df = df.grid[i.df]
  
  model = lm(log(h.loc.opt.mat[i.df, ]) ~ log(n.grid))
  p.val = summary(model)$coefficients[2, 4]
  
  if (p.val < 0.01) {
    coeff.opt.loc = model$coefficients
    coeff.opt.loc[1] = exp(coeff.opt.loc[1])
    cat("df = ", df, ": p-val = ", p.val, ", h_opt_loc = ", coeff.opt.loc[1], "*n^(", coeff.opt.loc[2], ")\n", sep = "")
  } else {
    coeff.opt.loc = c(0.0, 0.0)
    coeff.opt.loc[2] = 0.0
    coeff.opt.loc[1] = exp(mean(log(h.loc.opt.mat[i.df, ])))
    cat("df = ", df, ": p-val = ", p.val, ", h_opt_loc = ", coeff.opt.loc[1], "\n", sep = "")
  }
  
  for (i.n  in 1:length(n.grid)) {
    n  = n.grid[i.n]
    
    bw.loc = coeff.opt.loc[1]*n^coeff.opt.loc[2]
    
    err <- function(h, batch = 100L, eps = 1E-6) {
      xbar  = 0.0
      x2bar = 0.0
      for (rep in 1:batch) {
        x = rt(n, df); c = Qn_scale(x, df)
        obj = log(KWA::KWA.1D.scale(x, df, bw.loc = c*bw.loc, bw.scale = c*h))
        xbar  = xbar  + obj/batch
        x2bar = x2bar + obj^2/(batch - 1)
      }
      err = x2bar - batch/(batch - 1)*xbar^2
      return(err)
    }
    
    if (parallel.flag) {
      err.grid = foreach(h = h.grid, .inorder = TRUE, .export = c("err", "Qn_scale", "bw.loc"), .packages = c("KWA"), .combine = "c") %dopar% {
                         err(h, batch = batch)}
    } else {
      err.grid = foreach(h = h.grid, .inorder = TRUE, .export = c("err", "Qn_scale", "bw.loc"), .packages = c("KWA"), .combine = "c") %do% {
                         err(h, batch = batch)}
    }
    
    err.grid = log(err.grid)
    
    ## optimal scaling
    fit = smooth.spline(h.grid, err.grid) #, spar = 0.9)
    h.opt   = fit$x[which.min(fit$y)]
    err.opt = fit$y[which.min(fit$y)]
    
    se = sqrt(fit$cv.crit/fit$df)
    ind = min(which(fit$y <= err.opt + 2.0*se))
    h.opt.2se   = fit$x[ind]
    err.opt.2se = fit$y[ind]
    
    h.opt.mat[i.df, i.n] = h.opt.2se
    ##
    
    if (save.flag) {
      save(h.grid, err.grid, h.opt, err.opt, h.opt.2se, err.opt.2se,
           file = paste("output/scale.bw.df=", df, ".n=", n, ".RData", sep = ""))
    }
    
    if (plot.flag) {
      file.name = paste("fig/scale.bw.df=", df, ".n=", n, ".pdf", sep = "")
      grDevices::pdf(file.name, width = 10, height = 6)
      
      par(mfrow = c(1L, 1L))
      par(mar = c(4.0, 4.0, 2.0, 0.5))
      
      dy = max(err.grid) - min(err.grid)
      
      plot(h.grid, err.grid, type = "p", xlab = bquote("Bandwidth"~italic(h)), ylab = "log(Var)", 
           ylim = c(min(err.grid), min(err.grid) + 1.15*dy),
           main = bquote("KWA scale ("*italic(df)*" = "*.(df)*", "*italic(n)*" = "*.(n)*")"))
      lines(fit$x, fit$y, col = "blue")
      points(h.opt, err.opt, col = "red", pch = 8)
      points(h.opt.2se, err.opt.2se, col = "green", pch = 8)
      
      legend("topright",  legend = c(bquote(italic(h)["opt"]~"="~.(h.opt)), bquote(italic(h)["opt 2se"]~"="~.(h.opt.2se))),
             col = c("red", "green"), pch = c(8, 8))
      
      grDevices::dev.off()
    }
  }
}

print(proc.time() - ptm)

cat("optimal bandwidth for KWA scale:\n")
print(h.opt.mat)

if (save.flag) {
  save(h.opt.mat, file = "opt.bw.scale.RData")
}

try(dev.off(dev.list()["RStudioGD"]), silent = TRUE)
try(dev.off(), silent = TRUE)

if (parallel.flag)
  stopCluster(HPC)