## (C) Michael Pokojovy (2023)

setwd("???")

source("auxil.R")

# Select the dataset
dataset.name = "AMC"
#dataset.name = "GME"
#dataset.name = "FB"


if (dataset.name == "FB") {
  dataset = read.csv("FB.csv", sep = ",")
} else if (dataset.name == "AMC") {
  dataset = read.csv("AMC.csv", sep = ",")
} else {
  dataset = read.csv("GME.csv", sep = ",")
}

dataset$Date = as.Date(dataset$Date, tryFormats = "%Y-%m-%d")
dataset = dataset[order(dataset$Date), ]

## Prepare log-diffs

x = diff(log(dataset$Close), differences = 1)
n = length(x)

## Do search

df.search = seq(from = 1.00, to = 15.0, by = 0.05)

pval.max  = 0
pvals = rep(0.0, length(df.search))
df    = NA
loc   = NA
scale = NA

for (i.df in 1:length(df.search)) {
  .df = df.search[i.df]
  
  est = KWA1D::KWA1D.loc.scale(x, .df)
  
  z = (x - est$loc)/est$scale
  
  #pval = ks.test(z, "pt", .df, exact = TRUE)$p.value
  pval = goftest::cvm.test(z, "pt", df = .df)$p.value
  
  pvals[i.df] = pval
  
  if (pval >= pval.max) {
    pval.max  = pval
    df    = .df
    loc   = est$loc
    scale = est$scale
  }
}

## df selection
file.name = paste(dataset.name, ".df.selection.pdf", sep = "")
grDevices::pdf(file.name, width = 6, height = 5)
par(mfcol = c(1, 1))

plot(df.search, pvals, type = "l", ylim = c(0.0, 1.0), 
     main = "df estimation using Cramer-von Mises test",
     xlab = "df", ylab = "p-value")
points(df, pval.max)

grDevices::dev.off()

## density plot
file.name = paste(dataset.name, ".density.pdf", sep = "")
grDevices::pdf(file.name, width = 6, height = 5)
par(mfcol = c(1, 1))

x.grid = seq(from = min(x), to = max(x), length.out = 500)
y.grid = dt((x.grid - loc)/scale, df = df)/scale
plot(density(x, adjust = 1.0), main = bquote("Density plot ("*hat(df)~"="~.(df)*")"), lwd = 2)
lines(x.grid, y.grid, col = "red", lty = 2, lwd = 2)

grDevices::dev.off()

##
try(dev.off(dev.list()["RStudioGD"]), silent = TRUE)
try(dev.off(), silent = TRUE)

## Table

cat("Stock: ", dataset.name, "\n\n", sep = "")
cat("Location estimators:\n")

est.list = c("KWA", "Qn", "usual", "MCD.50", "MCD.75", "HL")
for (est in est.list) {
  cat(est, " = ", mu.hat(x, df = df, estimator = est), "\n")
}

cat("\n")
cat("Scale estimators:\n")

est.list = c("KWA", "Qn", "IQR", "usual", "MCD.50", "MCD.75", "HL")
for (est in est.list) {
  cat(est, " = ", sigma.hat(x, df = df, estimator = est), "\n")
}

cat("\n")
cat("AIC:\n")

est.list = c("KWA", "Qn", "usual", "MCD.50", "MCD.75", "HL")
for (est in est.list) {
  mu    = mu.hat(x, df = df, estimator = est)
  sigma = sigma.hat(x, df = df, estimator = est)
  
  npar = 3L
  z = (x - mu)/sigma
  AIC = 2.0*npar - 2.0*sum(dt(z, df = df, log = TRUE))
  
  cat(est, " = ", AIC, "\n")
}

cat("\n")
cat("BIC:\n")

est.list = c("KWA", "Qn", "usual", "MCD.50", "MCD.75", "HL")
for (est in est.list) {
  mu    = mu.hat(x, df = df, estimator = est)
  sigma = sigma.hat(x, df = df, estimator = est)
  
  npar = if (est == "HL") 2L else 3L
  z = (x - mu)/sigma
  BIC = npar*log(n) - 2.0*sum(dt(z, df = df, log = TRUE))
  
  cat(est, " = ", BIC, "\n")
}

cat("\n")
cat("CvM p-value:\n")

est.list = c("KWA", "Qn", "usual", "MCD.50", "MCD.75", "HL")
for (est in est.list) {
  mu    = mu.hat(x, df = df, estimator = est)
  sigma = sigma.hat(x, df = df, estimator = est)
  
  if (is.na(sigma))
    next
  
  npar = if (est == "HL") 2L else 3L
  z = (x - mu)/sigma
  pval = goftest::cvm.test(z, "pt", df = .df, estimated = TRUE)$p.value
  
  cat(est, " = ", pval, "\n")
}