library(lavaan)
library(cubature)
library(mvtnorm)
library(dplyr)
library(reshape2)
library(ggplot2)

source("computeTrue.R")
source("comparisonMethods2.R")
source("latentMethod2.R")

savemarker <- 50
k <- 1
p <- 3
tag <- "IV_multiplicative"

for (sample.size in 0:4) {
  sample.size <- 200 + 200 * sample.size
  est.df <- data.frame(matrix(nrow = 0, ncol = 3))
  colnames(est.df) <- c("latent", "linear", "IPW")
  set.seed <- 17
  latent <- c()
  ipw <- c()
  linear <- c()
  iv <- c()

  for (j in 1:100) {
    i <- j %% savemarker + 1

    U <- rnorm(sample.size)
    Z <- rmvnorm(sample.size, mean = c(.5, .7, -.8))
    A <- -.9 * U * (c(.1, -.7, .4) %*% t(Z)) - 1 + .4 * rnorm(sample.size)
    Y <- -.8 * U + .5 * A + .2 + rnorm(sample.size)
    raw.data <- data.frame("Z" = Z, "A" = A, "Y" = Y)
    colnames(raw.data) <- c("V1", "V2", "V3", "A", "Y")

    ate.true <- .5

    linear.ate <- linearEst(raw.data)
    linear[i] <- linear.ate
    ipw.ate <- IPWest(raw.data)
    ipw[i] <- ipw.ate
    iv.ate <- IVest(raw.data)
    iv[i] <- iv.ate

    ate.est <- latent.ATE(raw.data, k, p)
    latent[i] <- ate.est

    if (j %% savemarker == 0) {
      print(sample.size)
      print(j)
      est.df <- rbind(est.df, cbind(latent, linear, ipw, iv))
      write.table(est.df,
        file.path(
          "Data", "Estimates",
          paste("ests_", tag, as.character(sample.size),
            "_IVtrial.csv",
            sep = ""
          )
        ),
        sep = ",",
        row.names = FALSE
      )
      latent <- c()
      ipw <- c()
      linear <- c()
      iv <- c()
    }
  }
}


mean.ests <- data.frame(matrix(nrow = 0, ncol = 5))

for (j in 0:4) {
  samp.size <- j * 200 + 200
  temp.df <- read.csv(
    file.path(
      "Data", "Estimates",
      paste("ests_", tag, as.character(samp.size),
        "_IVtrial.csv",
        sep = ""
      )
    ),
    colClasses = "numeric"
  )
  mean.ests[4 * j + 1, ] <- c(
    samp.size, "linear", mean(temp.df$linear, na.rm = TRUE),
    quantile(temp.df$linear, probs = c(.05, .95), na.rm = TRUE)
  )
  mean.ests[4 * j + 2, ] <- c(
    samp.size, "IPW", mean(temp.df$ipw, na.rm = TRUE),
    quantile(temp.df$ipw, probs = c(.05, .95), na.rm = TRUE)
  )
  mean.ests[4 * j + 3, ] <- c(
    samp.size, "latent", mean(temp.df$latent, na.rm = TRUE),
    quantile(temp.df$latent, probs = c(.05, .95), na.rm = TRUE)
  )
  mean.ests[4 * j + 4, ] <- c(
    samp.size, "IV", mean(temp.df$iv, na.rm = TRUE),
    quantile(temp.df$iv, probs = c(.05, .95), na.rm = TRUE)
  )
}

colnames(mean.ests) <- c("SampleSize", "Type", "Mean", "Q.05", "Q.95")
mean.ests$Type <- factor(mean.ests$Type)
mean.ests$SampleSize <- as.numeric(mean.ests$SampleSize)
mean.ests$Mean <- as.numeric(mean.ests$Mean)
mean.ests$Q.05 <- as.numeric(mean.ests$Q.05)
mean.ests$Q.95 <- as.numeric(mean.ests$Q.95)

mean.ests <- mean.ests[order(mean.ests$SampleSize), ]

ggplot(mean.ests) +
  geom_ribbon(aes(
    x = SampleSize, ymin = Q.05, ymax = Q.95, fill = Type,
    alpha = .025
  )) +
  geom_line(aes(x = SampleSize, y = Mean, group = Type, color = Type),
    linewidth = 2
  ) +
  geom_point(aes(x = SampleSize, y = Mean, color = Type), size = 3) +
  geom_hline(yintercept = 0.5) +
  annotate("text", x = 215, y = .51, label = "True ATE") +
  xlab("Sample Size") +
  ylab("Average ATE Estimate")
ggsave(file.path("Data", "Figures", paste("Errs_by_sample.size_",
  tag, "_IVtrial.png",
  sep = ""
)))
