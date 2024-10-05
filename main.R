library(lavaan)
library(cubature)
library(mvtnorm)
library(dplyr)
library(reshape2)
library(ggplot2)
source("datagen.R")
source("computeTrue.R")
source("comparisonMethods.R")
source("latentMethod.R")

data.import <- function(tag, n) {
  A <- read.csv(
    file.path(
      "Data", "Samples",
      paste("A_", tag, "_", n, ".csv", sep = "")
    ),
    header = FALSE
  )
  Z <- read.csv(
    file.path(
      "Data", "Samples",
      paste("Z_", tag, "_", n, ".csv", sep = "")
    ),
    header = FALSE
  )
  Y <- read.csv(
    file.path(
      "Data", "Samples",
      paste("Y_", tag, "_", n, ".csv", sep = "")
    ),
    header = FALSE
  )
  df <- Z
  df$A <- A$V1
  df$Y <- Y$V1
  return(df)
}

tag <- "1D"
savemarker <- 100
k <- 1
p <- 9

dataGen(
  k, p, 1, .25, 1000, .3, c(.5, .25, .3, -.1, -.2, .4, .7, .8, -.6),
  1, c(.8, -.7),
  c(.5, .6), tag
)

for (sample.size in 0:4) {
  sample.size <- 200 + 200 * sample.size
  est.df <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(est.df) <- c("latent", "linear", "IPW", "IV")
  set.seed <- 17
  latent <- c()
  ipw <- c()
  linear <- c()
  iv <- c()

  for (j in 1:100) {
    i <- j %% savemarker + 1

    azGen(tag, sample.size)

    raw.data <- data.import(tag, sample.size)

    ate.true <- .5

    linear.ate <- linearEst(raw.data)
    linear[i] <- linear.ate
    ipw.ate <- continuousIPWest(raw.data)
    ipw[i] <- ipw.ate
    iv.ate <- IVest(raw.data)
    iv[i] <- iv.ate

    ate.est <- latent.ATE(raw.data, p)
    latent[i] <- ate.est

    if (j %% savemarker == 0) {
      print(sample.size)
      print(j)
      est.df <- rbind(est.df, cbind(latent, linear, ipw, iv))
      write.table(est.df,
        file.path(
          "Data", "Estimates",
          paste("ests_", tag, as.character(sample.size),
            ".csv",
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
        ".csv",
        sep = ""
      )
    ),
    colClasses = "numeric"
  )
  mean.ests[4 * j + 1, ] <- c(
    samp.size, "linear", mean(temp.df$linear),
    quantile(temp.df$linear, probs = c(.05, .95))
  )
  mean.ests[4 * j + 2, ] <- c(
    samp.size, "IPW", mean(temp.df$ipw),
    quantile(temp.df$ipw, probs = c(.05, .95))
  )
  mean.ests[4 * j + 3, ] <- c(
    samp.size, "latent", mean(temp.df$latent),
    quantile(temp.df$latent, probs = c(.05, .95))
  )
  mean.ests[4 * j + 4, ] <- c(
    samp.size, "IV", mean(temp.df$iv),
    quantile(temp.df$iv, probs = c(.05, .95))
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
    alpha = .05
  )) +
  geom_line(aes(x = SampleSize, y = Mean, group = Type, color = Type),
    linewidth = 2
  ) +
  geom_point(aes(x = SampleSize, y = Mean, color = Type), size = 3) +
  geom_hline(yintercept = 0.5) +
  annotate("text", x = 215, y = .51, label = "True ATE") +
  xlab("Sample Size") +
  ylab("Average ATE Estimate") +
  facet_wrap(Type ~ ., scales = "free")
ggsave(file.path("Data", "Figures", paste("Errs_by_sample_size_",
  tag, ".png",
  sep = ""
)))
