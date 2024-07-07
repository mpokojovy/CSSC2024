# (C) Michael Pokojovy (2024)

source("auxil.R")

load("opt.bw.loc.RData")
h.loc.opt.mat = h.opt.mat
rm(h.opt.mat)

load("opt.bw.scale.RData")
h.scale.opt.mat = h.opt.mat
rm(h.opt.mat)

set.seed(2)

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

bias.mat = matrix(0.0, nrow = length(df.grid), ncol = length(n.grid))
rownames(bias.mat) = sapply(df.grid, function(df) paste("df=", df, sep = ""))
colnames(bias.mat) = sapply(n.grid,  function(n)  paste("n=",   n, sep = ""))

nrep = 1E6

ptm = proc.time()

for (i.df in 1:length(df.grid)) {
  df = df.grid[i.df]

  ## location
  model = lm(log(h.loc.opt.mat[i.df, ]) ~ log(n.grid))
  p.val = summary(model)$coefficients[2, 4]

  if (p.val < 0.01) {
    coeff.opt.loc = model$coefficients
    coeff.opt.loc[1] = exp(coeff.opt.loc[1])
  } else {
    coeff.opt.loc[2] = 0.0
    coeff.opt.loc[1] = exp(mean(log(h.loc.opt.mat[i.df, ])))
  }

  ## scale
  model = lm(log(h.scale.opt.mat[i.df, ]) ~ log(n.grid))
  p.val = summary(model)$coefficients[2, 4]

  if (p.val < 0.01) {
    coeff.opt.scale = model$coefficients
    coeff.opt.scale[1] = exp(coeff.opt.scale[1])
  } else {
    coeff.opt.scale = c(0.0, 0.0)
    coeff.opt.scale[2] = 0.0
    coeff.opt.scale[1] = exp(mean(log(h.scale.opt.mat[i.df, ])))
  }

  for (i.n in 1:length(n.grid)) {
    n = n.grid[i.n]

    bw.loc   = coeff.opt.loc[1]*n^coeff.opt.loc[2]
    bw.scale = coeff.opt.scale[1]*n^coeff.opt.scale[2]

    sim <- function(n, df) {
      x = rt(n, df); c = Qn_scale(x, df)
      obj = log(KWA::KWA.1D.scale(x, df, bw.loc = c*bw.loc, bw.scale = c*bw.scale))
      return(obj)
    }

    if (parallel.flag) {
      bias.mat[i.df, i.n] = foreach(rep = 1:nrep, .inorder = FALSE, .export = c("sim", "Qn_scale", "bw.loc", "bw.scale"), .packages = c("KWA"), .combine = "+") %dopar% {
                                    sim(n, df)/nrep}
    } else {
      bias.mat[i.df, i.n] = foreach(rep = 1:nrep, .inorder = FALSE, .export = c("sim", "Qn_scale", "bw.loc", "bw.scale"), .packages = c("KWA"), .combine = "+") %do% {
                                    sim(n, df)/nrep}
    }
  }
}

print(proc.time() - ptm)

if (parallel.flag)
  stopCluster(HPC)

cat("bias for KWA scale:\n")
print(bias.mat)

if (save.flag) {
  save(bias.mat, file = "bias.RData")
}

for (i.df in 1:length(df.grid)) {
  df = df.grid[i.df]
  
  ## location 
  model = lm(log(h.loc.opt.mat[i.df, ]) ~ log(n.grid))
  p.val = summary(model)$coefficients[2, 4]
  
  if (p.val < 0.01) {
    coeff.opt.loc = model$coefficients
    coeff.opt.loc[1] = exp(coeff.opt.loc[1])
    cat("df = ", df, ": p-val = ", p.val, ", h_opt_loc = ", coeff.opt.loc[1], "*n^(", coeff.opt.loc[2], ")\n", sep = "")
  } else {
    coeff.opt.loc[2] = 0.0
    coeff.opt.loc[1] = exp(mean(log(h.loc.opt.mat[i.df, ])))
    cat("df = ", df, ": p-val = ", p.val, ", h_opt_loc = ", coeff.opt.loc[1], "\n", sep = "")
  }
  
  ## scale
  model = lm(log(h.scale.opt.mat[i.df, ]) ~ log(n.grid))
  p.val = summary(model)$coefficients[2, 4]
  
  if (p.val < 0.01) {
    coeff.opt.scale = model$coefficients
    coeff.opt.scale[1] = exp(coeff.opt.scale[1])
    cat("df = ", df, ": p-val = ", p.val, ", h_opt_scale = ", coeff.opt.scale[1], "*n^(", coeff.opt.scale[2], ")\n", sep = "")
  } else {
    coeff.opt.scale = c(0.0, 0.0)
    coeff.opt.scale[2] = 0.0
    coeff.opt.scale[1] = exp(mean(log(h.scale.opt.mat[i.df, ])))
    cat("df = ", df, ": p-val = ", p.val, ", h_opt_scale = ", coeff.opt.scale[1], "\n", sep = "")
  }
  
  ## bias
  if (plot.flag) {
    .df = data.frame(n = n.grid, scaling = bias.mat[i.df, ])
    model = nls(scaling ~ a + b/(n + c) + d/(n^2 + e), data = .df, algorithm = "port",
                control = nls.control(maxiter = 500L, printEval = TRUE),
                start = c(a = 1.0, b = 1.0, c = 1.0, d = 1.0, e = 50.0^2), lower = c(-Inf, -Inf, 1.0, -Inf, 1.0))
    theta = summary(model)$coefficients
    
    cat("df = ", df, ": bias = ", theta[1], " + (", theta[2], ")/(n + ", theta[3], ")",
                                  " + (", theta[4], ")/(n^2 + ", sqrt(theta[5]), "^2)\n", sep = "")
    
    file.name = paste("fig/bias.df=", df, ".pdf", sep = "")
    grDevices::pdf(file.name, width = 10, height = 6)
    
    par(mfrow = c(1L, 1L))
    par(mar = c(4.0, 4.0, 2.0, 0.5))
    
    dy = max(bias.mat[i.df, ]) - min(bias.mat[i.df, ])
    
    plot(n.grid, bias.mat[i.df, ], type = "p", xlab = bquote("Sample size"~italic(n)), ylab = "Bias", 
         ylim = c(min(bias.mat[i.df, ]), min(bias.mat[i.df, ]) + 1.15*dy),
         main = bquote("Bias of KWA scale ("*italic(df)*" = "*.(df)*")"))
    points(n.grid, bias.mat[i.df, ])
    
    n.range = seq(20, max(n.grid), by = 1)
    scale = predict(model, newdata = data.frame(n = n.range))
    lines(n.range, scale, col = "blue")
    
    legend("topright", legend = c(bquote("simulated"),
                                  bquote("fitted: "~.(theta[1])*" + ("*.(theta[2])*")/("*italic(n)*" + "*.(theta[3])*")"*
                                       " + ("*.(theta[4])*")/("*italic(n)^2*" + "*.(sqrt(theta[5]))^2*")")),
           lty = c(NA, 1), pch = c(1, NA), col = c("black", "blue"))
    
    grDevices::dev.off()
  }
}

try(dev.off(dev.list()["RStudioGD"]), silent = TRUE)
try(dev.off(), silent = TRUE)
