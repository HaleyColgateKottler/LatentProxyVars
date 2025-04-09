library(boot)
source("datagen.R")
source("latentMethod.R")
library(ggplot2)
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
# saveParams(2, 6, 3, .1, 0, rep(0,6), diag(2),
#            c(.7, -.5, -.8), c(.5, .7, .2), 'boot')

latent.func <- function(data.df, indices){
  df <- data.df[indices,]
  latent.ATE(df, 1:6, 6)$estimate
}

save.name <- 'coverage.csv'
if (file.exists(save.name)){
  df <- read.csv(save.name)
} else {
  df <- data.frame(matrix(0, nrow = 0, ncol = 7))
}
colnames(df) <- c('est', 'norm_low', 'norm_high',
            'perc_low', 'perc_high',
            'base_low', 'base_high')

for (i in 11:100){
  print(i)
  set.seed(i)
  azGen('boot', 500)
  raw.data <- data.import('boot', 500)
  boot.out <- boot(data = raw.data, statistic = latent.func, R = 2000)
  cis <- boot.ci(boot.out, type = c("norm", "perc", "basic"))
  df[nrow(df) + 1, ] <- c(boot.out$t0[[1]], cis$normal[2:3], cis$percent[4:5], cis$basic[4:5])
  if (i %% 10 == 0){
    write.csv(df, save.name, row.names = FALSE)
  }
}

write.csv(df, save.name, row.names = FALSE)

norm.cov <- c()
perc.cov <- c()
base.cov <- c()
for (i in 1:nrow(df)){
  if (df$norm_low[i] < .5 && df$norm_high[i] > .5){
    norm.cov[i] <- 1
  } else {
    norm.cov[i] <- 0
  }
  if (df$perc_low[i] < .5 && df$perc_high[i] > .5){
    perc.cov[i] <- 1
  } else {
    perc.cov[i] <- 0
  }
  if (df$base_low[i] < .5 && df$base_high[i] > .5){
    base.cov[i] <- 1
  } else {
    base.cov[i] <- 0
  }
}
df$norm.cov <- factor(norm.cov)
df$perc.cov <- factor(perc.cov)
df$base.cov <- factor(base.cov)
tag <- 'norm'
ggplot(df) + geom_hline(yintercept = .5) +
    geom_point(aes(x=1:nrow(df), y = est, color = norm.cov)) +
    geom_errorbar(aes(x = 1:nrow(df), ymin = norm_low, ymax = norm_high,
                      color = norm.cov))
ggsave(paste("coverage_plot_", tag, "_norm.png", sep = ""))
df$norm.cov

tag <- 'perc'
ggplot(df) + geom_hline(yintercept = .5) +
  geom_point(aes(x=1:nrow(df), y = est, color = perc.cov)) +
  geom_errorbar(aes(x = 1:nrow(df), ymin = perc_low, ymax = perc_high,
                    color = perc.cov))
ggsave(paste("coverage_plot_", tag, "_perc.png", sep = ""))

tag <- 'base'
ggplot(df) + geom_hline(yintercept = .5) +
  geom_point(aes(x=1:nrow(df), y = est, color = base.cov)) +
  geom_errorbar(aes(x = 1:nrow(df), ymin = base_low, ymax = base_high,
                    color = base.cov))
ggsave(paste("coverage_plot_", tag, "_base.png", sep = ""))
