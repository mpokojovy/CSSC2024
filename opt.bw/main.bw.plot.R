# (C) Michael Pokojovy (2024)

PATH = NULL # set path

setwd(PATH)

load("opt.bw.loc.RData")
h.loc.opt.mat = h.opt.mat
rm(h.opt.mat)

load("opt.bw.scale.RData")
h.scale.opt.mat = h.opt.mat
rm(h.opt.mat)

df.grid = c(1:5, 10, 20, 30)
n.grid = c(30, 50, 100, 150, 200, 300)

for (i.df in 1:length(df.grid)) {
  df = df.grid[i.df]

  file.name = paste("fig/opt.bw.df=", df, ".pdf", sep = "")
  grDevices::pdf(file.name, width = 10, height = 6)
  
  par(mfrow = c(1L, 2L))
  par(mar = c(4.0, 5.0, 2.0, 0.5))
  
  ## location
  plot(log(n.grid), log(h.loc.opt.mat[i.df, ]), 
       xlab = bquote("log("*italic(n)*")"), ylab = bquote("log("*italic(h)["opt"]*")"),
       main = bquote("KWA location ("*italic(df)*" = "*.(df)*")"))
  
  model = lm(log(h.loc.opt.mat[i.df, ]) ~ log(n.grid))
  p.val = summary(model)$coefficients[2, 4]

  if (p.val < 0.01) {
    coeff.opt.loc = model$coefficients
    coeff.opt.loc[1] = exp(coeff.opt.loc[1])
    
    abline(model$coefficients[1], model$coefficients[2], col = "red")
  } else {
    coeff.opt.loc[2] = 0.0
    coeff.opt.loc[1] = exp(mean(log(h.loc.opt.mat[i.df, ])))
    
    abline(h = log(coeff.opt.loc[1]), col = "red")
  }

  ## scale
  plot(log(n.grid), log(h.scale.opt.mat[i.df, ]),
       xlab = bquote("log("*italic(n)*")"), ylab = bquote("log("*italic(h)["opt"]*")"),
       main = bquote("KWA scale ("*italic(df)*" = "*.(df)*")"))
  
  model = lm(log(h.scale.opt.mat[i.df, ]) ~ log(n.grid))
  p.val = summary(model)$coefficients[2, 4]

  if ((p.val < 0.01) && (df > 1)) {
    coeff.opt.scale = model$coefficients
    coeff.opt.scale[1] = exp(coeff.opt.scale[1])
    
    abline(model$coefficients[1], model$coefficients[2], col = "red")
  } else {
    coeff.opt.scale = c(0.0, 0.0)
    coeff.opt.scale[2] = 0.0
    coeff.opt.scale[1] = exp(mean(log(h.scale.opt.mat[i.df, ])))
    
    abline(h = log(coeff.opt.scale[1]), col = "red")
  }
  
  grDevices::dev.off()
}

try(dev.off(dev.list()["RStudioGD"]), silent = TRUE)
try(dev.off(), silent = TRUE)