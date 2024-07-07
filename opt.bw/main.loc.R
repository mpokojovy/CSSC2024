# (C) Michael Pokojovy (2024)

source("auxil.R")

set.seed(0)

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

# h.grid  = signif(c(exp(seq(from = log(0.05), to = log(1.00), length.out = 120L)), 
#                    seq(from = 1.05, to = 15.0, by = 0.05)), digits = 4)

h.grid  = signif(c(exp(seq(from = log(0.05), to = log(1.00), length.out = 120L)), 
                   seq(from = 1.025, to = 8.0, by = 0.025)), digits = 4)

ptm = proc.time()

for (i.df in 1:length(df.grid))
for (i.n  in 1:length(n.grid)) {
  df = df.grid[i.df]
  n  = n.grid[i.n]
  
  MSE <- function(h, batch = 100L, eps = 1E-6) {
    MSE = 0.0
    for (rep in 1:batch) {
      x = rt(n, df); c = Qn_scale(x, df)
      obj = KWA::KWA.1D.loc(x, df, bw = c*h)^2
      MSE = MSE + obj/batch
    }
    return(MSE)
  }
  
  if (parallel.flag) {
    MSE.grid = foreach(h = h.grid, .inorder = TRUE, .export = c("MSE", "Qn_scale"), .packages = c("KWA"), .combine = "c") %dopar% {
                       MSE(h, batch = batch)}
  } else {
    MSE.grid = foreach(h = h.grid, .inorder = TRUE, .export = c("MSE", "Qn_scale"), .packages = c("KWA"), .combine = "c") %do% {
                       MSE(h, batch = batch)}
  }

  MSE.grid = log(MSE.grid)
  
  ## optimal scaling
  fit = smooth.spline(h.grid, MSE.grid) #, spar = 0.9)
  h.opt   = fit$x[which.min(fit$y)]
  MSE.opt = fit$y[which.min(fit$y)]
  
  se = sqrt(fit$cv.crit/fit$df)
  ind = min(which(fit$y <= MSE.opt + 2.0*se))
  h.opt.2se   = fit$x[ind]
  MSE.opt.2se = fit$y[ind]

  h.opt.mat[i.df, i.n] = h.opt.2se
  ##
  
  if (save.flag) {
    save(h.grid, MSE.grid, h.opt, MSE.opt, h.opt.2se, MSE.opt.2se,
         file = paste("output/loc.bw.df=", df, ".n=", n, ".RData", sep = ""))
  }
  
  if (plot.flag) {
    file.name = paste("fig/loc.bw.df=", df, ".n=", n, ".pdf", sep = "")
    grDevices::pdf(file.name, width = 10, height = 6)
    
    par(mfrow = c(1L, 1L))
    par(mar = c(4.0, 4.0, 2.0, 0.5))
    
    dy = max(MSE.grid) - min(MSE.grid)
    
    plot(h.grid, MSE.grid, type = "p", xlab = bquote("Bandwidth"~italic(h)), ylab = "log(MSE)", 
         ylim = c(min(MSE.grid), min(MSE.grid) + 1.15*dy),
         main = bquote("KWA location ("*italic(df)*" = "*.(df)*", "*italic(n)*" = "*.(n)*")"))
    lines(fit$x, fit$y, col = "blue")
    points(h.opt, MSE.opt, col = "red", pch = 8)
    points(h.opt.2se, MSE.opt.2se, col = "green", pch = 8)
    
    legend("topright",  legend = c(bquote(italic(h)["opt"]~"="~.(h.opt)), bquote(italic(h)["opt 2se"]~"="~.(h.opt.2se))),
           col = c("red", "green"), pch = c(8, 8))
    
    grDevices::dev.off()
  }
}

print(proc.time() - ptm)

cat("optimal bandwidth for KWA location:\n")
print(h.opt.mat)

if (save.flag) {
  save(h.opt.mat, file = "opt.bw.loc.RData")
}

try(dev.off(dev.list()["RStudioGD"]), silent = TRUE)
try(dev.off(), silent = TRUE)

if (parallel.flag)
  stopCluster(HPC)