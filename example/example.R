PATH = NULL # set path

setwd(PATH)

source("auxil.R")

## Choose the dataset
dataset.name = "AMC"
#dataset.name = "GME"
#dataset.name = "FB"

if (dataset.name == "FB") {
  df = read.csv("FB.csv", sep = ",")
  df$Date = as.Date(df$Date, tryFormats = "%Y-%m-%d")
  df = df[order(df$Date), ]
} else if (dataset.name == "AMC") {
  df = read.csv("AMC.csv", sep = ",")
  df$Date = as.Date(df$Date, tryFormats = "%Y-%m-%d")
  df = df[order(df$Date), ]
} else {
  df = read.csv("GME.csv", sep = ",")
  df$Date = as.Date(df$Date, tryFormats = "%Y-%m-%d")
  df = df[order(df$Date), ]
}

## Prepare log-diffs

#price.out = approx(df$Date - min(df$Date), df$Close)$y
#x = diff(log(price.out), differences = 1)

x = diff(log(df$Close), differences = 1)

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
  pval = (goftest::cvm.test(z, "pt", df = .df)$p.value)
  
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
plot(density(x, adjust = 1.0), main = bquote("Density plot ("*hat(df)~"="~.(df)*")"))
lines(x.grid, y.grid, col = "red")

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