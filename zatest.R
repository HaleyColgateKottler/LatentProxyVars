library(lavaan)
library(cubature)
library(mvtnorm)
library(dplyr)
library(reshape2)
library(ggplot2)
source("datagen2.R")
# source("computeTrue.R")
source("comparisonMethods2.R")
source("latentMethod2.R")

dataImport <- function(tag, n) {
  A <- read.csv(file.path("Data", "Samples", paste("A_", tag, "_", n, ".csv", sep = "")), header = FALSE)
  colnames(A) <- "A"
  Z <- read.csv(file.path("Data", "Samples", paste("Z_", tag, "_", n, ".csv", sep = "")), header = FALSE)
  colnames(Z) <- paste0("V", 1:ncol(Z))
  Y <- read.csv(file.path("Data", "Samples", paste("Y_", tag, "_", n, ".csv", sep = "")), header = FALSE)
  colnames(Y) <- "Y"
  H <- read.csv(file.path("Data", "Samples", paste("H_", tag, "_", n, ".csv", sep = "")), header = FALSE)
  colnames(H) <- paste0("H", 1:ncol(H))
  df <- cbind(Z, H, A, Y)
  return(df)
}

tag <- "1DwithAZ"

k <- 1
p <- 7
savemarker <- 10
dataGenNCE(k, p, 3, 1, 1000, .3, c(.5, .7, .9, -.6, -.8, .2, -.4), 2, matrix(c(1), nrow = 1,
                                                                 byrow = TRUE),
        .75, c(.5), c(-.2, .9, .7, .5, -.4, -.8, .1), tag)
sampleSize <- 1000
load(file.path("Data", "Parameters", paste("param_", tag, ".RData", sep = "")))
# lambda
est.df <- data.frame(matrix(nrow = 0, ncol = 33))
colnames(est.df) <- c(
  paste0("Latent", .1 * (0:10)), paste0("Linear", .1 * (0:10)),
  paste0("IPW", .1 * (0:10))
)
set.seed <- 17

latent <- matrix(NaN, nrow = savemarker, ncol = 11)
IPW <- matrix(NaN, nrow = savemarker, ncol = 11)
linear <- matrix(NaN, nrow = savemarker, ncol = 11)
for (j in 0:9) {
  i <- j %% savemarker + 1
  
  rawData <- dataImport(tag, sampleSize)
  ZH <- rawData[, c(paste0("V", 1:p), paste0("H", 1:k))]
  
  rawData <- rawData[, which(!(colnames(rawData) %in% paste0("H", 1:k)))]
  
  z_weights = c(.8, -.6, .7, .9, -.8, .4, -.9)
  
  for (h.weight.base in 0:10) {
    h.weight <- h.weight.base * .1
    
    for (entry.row in 1:1000) {
      zhrow <- ZH[entry.row, ]
      epsilon_A <- rnorm(1, mean = 0, sd = 1)
      epsilon_Y <- rnorm(1, mean = 0, sd = sqrt(var_Y))
      rawData$A[entry.row] <- (1 - h.weight) * (z_weights %*% t(zhrow[1:p])) +
        h.weight * (lambda[p+1, 1:k] %*% t(zhrow[(p + 1):(p + k)])) + 
        A_intercept + epsilon_A
      rawData$Y[entry.row] <- (C1 %*% t(zhrow[(p + 1):(p + k)])) + 
        C0 * rawData$A[entry.row] + Y_intercept + epsilon_Y
    }
    
    linear.ATE <- linearEst(rawData)
    linear[i, h.weight.base + 1] <- linear.ATE
    IPW.ATE <- IPWest(rawData)
    IPW[i, h.weight.base + 1] <- IPW.ATE
    ATEest <- latent.ATE.NCE(rawData, k, p)
    latent[i, h.weight.base + 1] <- ATEest
  }
  
  ATE.true <- C0
  
  if ((j + 1) %% savemarker == 0) {
    print(j)
    est.df <- rbind(est.df, cbind(latent, linear, IPW))
    write.table(est.df, file.path("Data", "Estimates", paste("ests_", tag, "_HZweighted.csv", sep = "")), sep = ",", row.names = FALSE)
    latent <- matrix(NaN, nrow = savemarker, ncol = 11)
    IPW <- matrix(NaN, nrow = savemarker, ncol = 11)
    linear <- matrix(NaN, nrow = savemarker, ncol = 11)
  }
}

temp.df <- read.csv(file.path("Data", "Estimates", paste("ests_", tag, "_HZweighted.csv", sep = "")),
                    colClasses = "numeric"
)

mean.ests <- data.frame(matrix(nrow = 0, ncol = 5))
sampSize <- 1000
for (j in 0:10) {
  hWeight <- j * .1
  ciMult <- qt(.975, sampSize - 1)
  mean.ests[3 * j + 3, ] <- c(
    hWeight, "linear", mean(as.numeric(temp.df[, 11 + j + 1])),
    quantile(temp.df[, 11 + j + 1], probs = c(.05, .95), na.rm = TRUE)
  )
  mean.ests[3 * j + 2, ] <- c(
    hWeight, "IPW", mean(as.numeric(temp.df[, 22 + j + 1])),
    quantile(temp.df[, 22 + j + 1], probs = c(.05, .95), na.rm = TRUE)
  )
  mean.ests[3 * j + 1, ] <- c(
    hWeight, "latent", mean(as.numeric(temp.df[, j + 1])),
    quantile(temp.df[, j + 1], probs = c(.05, .95), na.rm = TRUE)
  )
}

colnames(mean.ests) <- c("HWeight", "Type", "Mean", "Q.05", "Q.95")
mean.ests$HWeight <- as.numeric(mean.ests$HWeight)
mean.ests$Type <- factor(mean.ests$Type)
mean.ests$Mean <- as.numeric(mean.ests$Mean)
mean.ests$mins <- as.numeric(mean.ests$Q.05)
mean.ests$maxs <- as.numeric(mean.ests$Q.95)

mean.ests <- mean.ests[order(mean.ests$HWeight), ]

ggplot(mean.ests) +
  geom_ribbon(aes(x = HWeight, ymin = mins, ymax = maxs, fill = Type, alpha = .1)) +
  geom_line(aes(x = HWeight, y = Mean, group = Type, color = Type), linewidth = 1) +
  geom_point(aes(x = HWeight, y = Mean, color = Type), size = 3) +
  geom_hline(yintercept = 0.5) +
  annotate("text", x = .5, y = .51, label = "True ATE") +
  xlab("H Weight") +
  ylab("Average ATE Estimate")
ggsave(file.path("Data", "Figures", paste("Errs_by_ZHweight_", tag, ".png", sep = "")))
