## (C) Michael Pokojovy (2023)

library("KWA1D")
library(foreach)

setwd("???")

source("auxil.R")

set.seed(123)

nrep = 100000L # 1E5L

logarithm.flag = TRUE
parallel.flag  = TRUE

if (parallel.flag) {
  library(doParallel)
  HPC <- makeCluster(detectCores())
  registerDoParallel(HPC)
}

##
MSE.file         = file("MSE.txt")
err.MSE.file     = file("err.MSE.txt")
paired.test.file = file("paired.test.txt")

cat("", file = MSE.file,         append = FALSE)
cat("", file = err.MSE.file,     append = FALSE)
cat("", file = paired.test.file, append = FALSE)

n.array  = c(30, 50, 100, 200, 500)
df.array = c(1, 2, 3, 4, 5, 10, 20, 30, Inf)

tic = proc.time()

res.df = data.frame()

for (type in c("loc", "scale"))
for (estimator in if (type == "loc") c("KWA", "MCD.50", "MCD.75", "Qn", "HL", "usual") else 
                                     c("KWA", "MCD.50", "MCD.75", "Qn", "HL", "IQR", "usual")) {
  MSE.matrix     = matrix(0, nrow = length(df.array), ncol = length(n.array))
  err.MSE.matrix = matrix(0, nrow = length(df.array), ncol = length(n.array))
  
  rownames(MSE.matrix) = sapply(df.array, function(df) paste("df=", df, sep = ""))
  colnames(MSE.matrix) = sapply(n.array,  function(n)  paste("n=",  n,  sep = ""))
  
  rownames(err.MSE.matrix) = sapply(df.array, function(df) paste("df=", df, sep = ""))
  colnames(err.MSE.matrix) = sapply(n.array,  function(n)  paste("n=",  n,  sep = ""))
  
  for (i.df in 1:length(df.array))
  for (i.n  in 1:length(n.array)) {
    df = df.array[i.df]
    n  = n.array[i.n]
    
    res = do.simulation(n = n, df = df, nrep = nrep, type = type, estimator = estimator, 
                        logarithm.flag = logarithm.flag, parallel.flag = parallel.flag)
    
    if (type == "loc") {
      MSE.matrix[i.df, i.n]     = res[1]
      err.MSE.matrix[i.df, i.n] = res[3]/sqrt(nrep)
    } else {
      MSE.matrix[i.df, i.n]     = res[2]
      err.MSE.matrix[i.df, i.n] = res[3]/sqrt(nrep - 1)
    }
    
    res = data.frame(param = type, est = estimator, n = n, df = df, err = MSE.matrix[i.df, i.n], err.ste = err.MSE.matrix[i.df, i.n])
    
    res.df = rbind(res.df, res)
  }

  cat("Estimator: ", estimator, " ", type, "\n", sep = "", file = "MSE.txt", append = TRUE)
  cat(capture.output(print(MSE.matrix)), "", sep = "\n", file = "MSE.txt", append = TRUE)
  
  cat("Estimator: ", estimator, " ", type, "\n", sep = "", file = "err.MSE.txt", append = TRUE)
  cat(capture.output(print(err.MSE.matrix)), "", sep = "\n", file = "err.MSE.txt", append = TRUE)
  
  if (estimator == "KWA") {
    MSE.matrix.KWA     = MSE.matrix
    err.MSE.matrix.KWA = err.MSE.matrix
  } else {
    MSE.diff         = MSE.matrix.KWA - MSE.matrix
    std.err.MSE.diff = sqrt((err.MSE.matrix.KWA^2 + err.MSE.matrix^2)/nrep)
    
    test.matrix = ifelse(MSE.diff - qnorm(0.975)*std.err.MSE.diff > 0, "competitor", 
                         ifelse(MSE.diff + qnorm(0.975)*std.err.MSE.diff < 0, "ours", "tie"))
    
    cat("Estimator: ", estimator, " ", type, "\n", sep = "", file = "paired.test.txt", append = TRUE)
    cat(capture.output(print(test.matrix)), "", sep = "\n", file = "paired.test.txt", append = TRUE)
  }
}

save(res.df, file = "results.RData")

print(proc.time() - tic)

close(MSE.file)
close(err.MSE.file)
close(paired.test.file)

if (parallel.flag)
  stopCluster(HPC)