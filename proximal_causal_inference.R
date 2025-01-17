# proximal causal inference
library(lavaan)
library(cubature)
library(mvtnorm)
library(dplyr)
library(pracma)
library(ggplot2)
source("latentMethod.R")
source("datagen.R")

sample.size <- 500

output <- sigmaSim(1, 6, 2)
lambda <- output[[1]]
psi <- output[[2]]

ests <- data.frame(matrix(0, nrow = 0, ncol = 2))

for (sample.size in seq(from = 500, to = 5000, by = 500)){
  temp.ests <- c()
  for (rep in 1:1000){
    U <- rnorm(sample.size)
    W <- t(lambda[1:3, 1] %*% t(U)) + rmvnorm(sample.size, mean = rep(0, 3), sigma = psi[1:3, 1:3])
    V <- t(lambda[4:6, 1] %*% t(U)) + rmvnorm(sample.size, mean = rep(0, 3), sigma = psi[4:6, 4:6])
    A <- t(lambda[7, 1] %*% t(U) + c(.9, .3, .7) %*% t(V) + .1*rnorm(sample.size))
    
    Y <- 2 + .5 * A + .75 * U + t(c(-.7, 1, .6) %*% t(W)) + .1*rnorm(sample.size)
    
    # calc h(w)
    
    df <- data.frame(cbind(A, V, Y))
    colnames(df) <- c("A", "V1", "V2", "V3", "Y")
    
    m.hW <- lm(W ~ A + V1 + V2 + V3 + A * V1 + A * V2 + A * V3, df)
    
    wav <- predict(m.hW, df)
    
    df <- cbind(df, wav)
    colnames(df) <- c("A", "V1", "V2", "V3", "Y", "wav1", "wav2", "wav3")
    
    # estimate tau_a
    m1 <- lm(Y ~ A + wav1 + wav2 + wav3 + A * wav1 + A * wav2 + A * wav3, df)
    alpha <- m1$coefficients[c("(Intercept)", "wav1", "wav2", "wav3")]
    gamma <- m1$coefficients[c("A", "A:wav1", "A:wav2", "A:wav3")]
    
    # calc CATE
    m.WV <- lm(W ~ V1 + V2 + V3, df)
    WV <- predict(m.WV, df)
    CATE <- gamma[1] + gamma[2:4] %*% t(WV)
    ATE <- mean(CATE)
    temp.ests <- c(temp.ests, ATE)
  }
  ests[nrow(ests) + 1,] <- c(sample.size, mean(temp.ests))
}

colnames(ests) <- c("SampleSize", "Estimate")
ggplot(ests) + geom_line(aes(x = SampleSize, y = Estimate))
