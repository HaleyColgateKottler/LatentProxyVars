uzuauy_az_gen <- function(tag, sample.size) {
  load(file.path(
    "Data", "Parameters",
    paste("param_", tag, ".RData", sep = "")
  ))

  H <- rmvnorm(sample.size, mean = rep(0, k))
  Z <- t(lambda[1:p, ] %*% t(H)) + rmvnorm(sample.size, sigma = psi[1:p, 1:p])
  colnames(Z) <- paste0("V", 1:p)
  A <- t(b %*% t(H)) + rnorm(sample.size)
  Y <- rowSums((alpha + A %*% gamma) * cbind(1, H)) + rnorm(sample.size, sd = sqrt(.1))

  write.table(H, file = file.path("Data", "Samples", paste("H_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Z, file = file.path("Data", "Samples", paste("Z_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(A, file = file.path("Data", "Samples", paste("A_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Y, file = file.path("Data", "Samples", paste("Y_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
}

uzuauy_test <- function(kvals, p, sample.size, trials, tag, savemarker = 100) {
  indexing <- 0
  for (z.level in c("low", "mid", "high")) {
    for (a.level in c("low", "mid", "high")) {
      for (y.level in c(20, 40, 60)) {
        est.df <- data.frame(matrix(nrow = 0, ncol = 4))
        colnames(est.df) <- c("latent", "linear", "IPW", "IV")
        latent <- c()
        ipw <- c()
        linear <- c()
        iv <- c()
        proximal <- c()

        for (j in 1:trials) {
          i <- j %% savemarker + 1
          temp.tag <- paste("z", z.level, "_a", a.level, "_y", y.level, sep = "")

          uzuauy_az_gen(temp.tag, sample.size)

          raw.data <- data.import(temp.tag, sample.size)

          linear.ate <- linearEst(raw.data)
          linear[i] <- linear.ate
          ipw.ate <- continuousIPWest(raw.data)
          ipw[i] <- ipw.ate
          iv.ate <- IVest(raw.data)
          iv[i] <- iv.ate
          proximal.ate <- proximal_causal(raw.data, paste("V", 1:floor(p/2), sep = ""),
                                          paste("V", (floor(p/2)+1):p, sep = ""))
          proximal[i] <- proximal.ate

          ate.est <- latent.ATE(raw.data, kvals, p)
          latent[i] <- ate.est$estimate

          if (j %% savemarker == 0 | j == trials) {
            print(c(z.level, a.level, y.level))
            print(j)
            est.df <- rbind(est.df, cbind(latent, linear, ipw, iv, proximal))
            write.table(est.df,
              file.path(
                "Data", "Estimates",
                paste("ests_", temp.tag, as.character(sample.size),
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
            proximal <- c()
          }
        }
      }
    }
  }
}

graph.uzuauy <- function(tag, sample.sizes) {
  mean.ests <- data.frame(matrix(nrow = 0, ncol = 7))
  j <- 0
  for (z.level in c("low", "mid", "high")) {
    for (a.level in c("low", "mid", "high")) {
      for (y.level in c(20, 40, 60)) {
        j <- j + 1
        temp.tag <- paste("z", z.level, "_a", a.level, "_y", y.level, sep = "")
        temp.df <- read.csv(
          file.path(
            "Data", "Estimates",
            paste("ests_", temp.tag, as.character(sample.size),
              ".csv",
              sep = ""
            )
          ),
          colClasses = "numeric"
        )
        mean.ests[5 * j - 4, ] <- c(
          z.level, a.level, y.level, "linear", mean(temp.df$linear),
          quantile(temp.df$linear, probs = c(.05, .95))
        )
        mean.ests[5 * j - 3, ] <- c(
          z.level, a.level, y.level, "IPW", mean(temp.df$ipw),
          quantile(temp.df$ipw, probs = c(.05, .95))
        )
        mean.ests[5 * j - 2, ] <- c(
          z.level, a.level, y.level, "latent", mean(temp.df$latent),
          quantile(temp.df$latent, probs = c(.05, .95))
        )
        mean.ests[5 * j - 1, ] <- c(
          z.level, a.level, y.level, "IV", mean(temp.df$iv),
          quantile(temp.df$iv, probs = c(.05, .95))
        )
        mean.ests[5 * j, ] <- c(
          z.level, a.level, y.level, "proximal", mean(temp.df$iv),
          quantile(temp.df$iv, probs = c(.05, .95))
        )
      }
    }
  }
  colnames(mean.ests) <- c("z.level", "a.level", "y.level", "Type", "Mean", "Q.05", "Q.95")
  mean.ests$Type <- factor(mean.ests$Type)
  mean.ests$Mean <- as.numeric(mean.ests$Mean)
  mean.ests$Q.05 <- as.numeric(mean.ests$Q.05)
  mean.ests$Q.95 <- as.numeric(mean.ests$Q.95)

  mean.ests$marker <- ifelse(mean.ests$z.level == "low", "1",
    ifelse(mean.ests$z.level == "mid", "2", "3")
  )
  mean.ests$marker <- ifelse(mean.ests$a.level == "low",
    paste(mean.ests$marker, "1", sep = ""),
    ifelse(mean.ests$a.level == "mid",
      paste(mean.ests$marker, "2", sep = ""),
      paste(mean.ests$marker, "3", sep = "")
    )
  )
  mean.ests$marker <- paste(mean.ests$marker, as.numeric(mean.ests$y.level) / 20,
    sep = ""
  )
  mean.ests$marker <- factor(mean.ests$marker)
  mean.ests$z.level <- paste(mean.ests$z.level, "UZ")
  mean.ests$a.level <- paste(mean.ests$a.level, "UA")
  mean.ests$y.level <- ifelse(mean.ests$y.level == 20, "low UY",
    ifelse(mean.ests$y.level == 40,
      "mid UY", "high UY"
    )
  )
  mean.ests$z.level <- factor(mean.ests$z.level, c("low UZ", "mid UZ", "high UZ"),
    ordered = TRUE
  )
  mean.ests$a.level <- factor(mean.ests$a.level, c("low UA", "mid UA", "high UA"),
    ordered = TRUE
  )
  mean.ests$y.level <- factor(mean.ests$y.level, c("low UY", "mid UY", "high UY"),
    ordered = TRUE
  )
  load(file.path(
    "Data", "Parameters",
    paste("param_", temp.tag, ".RData", sep = "")
  ))
  true_ate <- gamma[1]
  p1 <- ggplot(mean.ests) +
    geom_boxplot(aes(x = y.level, y = Mean)) +
    facet_wrap(z.level ~ a.level) +
    geom_hline(yintercept = true_ate) +
    ylab("Average ATE Estimate")
  ggsave(file.path("Data", "Figures", paste("UZUAUY_byY",
    tag, ".png",
    sep = ""
  )), p1)

  p2 <- ggplot(mean.ests) +
    geom_boxplot(aes(x = a.level, y = Mean)) +
    facet_wrap(y.level ~ z.level) +
    geom_hline(yintercept = true_ate) +
    ylab("Average ATE Estimate")
  ggsave(file.path("Data", "Figures", paste("UZUAUY_byA",
    tag, ".png",
    sep = ""
  )), p2)

  p3 <- ggplot(mean.ests) +
    geom_boxplot(aes(x = z.level, y = Mean)) +
    facet_wrap(a.level ~ y.level) +
    geom_hline(yintercept = true_ate) +
    ylab("Average ATE Estimate")
  ggsave(file.path("Data", "Figures", paste("UZUAUY_byZ",
    tag, ".png",
    sep = ""
  )), p3)
  return(list(p1,p2,p3))
}




# source('datagen.R')
# set.seed(127)
# k = 1
# p = 9
# 
# uz.lows <- list()
# uz.mids <- list()
# uz.highs <- list()
# infos <- data.frame(matrix(0, nrow = 0, ncol = 3))
# 
# for (i in 1:100){
#   uz.low <- sigmaSim(k, p-1, 1)
#   uz.mid <- sigmaSim(k, p-1, 2)
#   uz.high <- sigmaSim(k, p-1, 3)
#   
#   uz.lows[[i]] <- uz.low
#   uz.mids[[i]] <- uz.mid
#   uz.highs[[i]] <- uz.high
#   
#   info.uz <- c(det(uz.low[[1]] %*% t(uz.low[[1]]) + uz.low[[2]]) / det(uz.low[[2]]),
#                det(uz.mid[[1]] %*% t(uz.mid[[1]]) + uz.mid[[2]]) / det(uz.mid[[2]]),
#                det(uz.high[[1]] %*% t(uz.high[[1]]) + uz.high[[2]]) / det(uz.high[[2]])
#                )
#   info.uz <- .5*log(info.uz)
#   infos[i,] <- info.uz
# }
# 
# min.info <- min(infos[,1])
# max.info <- max(infos[,3])
# mid.info <- (min.info + max.info)/2
# 
# min.loc <- which(infos[,1] == min.info)
# max.loc <- which(infos[,3] == max.info)
# mid.difs <- abs(infos[,2] - mid.info)
# mid.loc <- which(mid.difs == min(mid.difs))
# 
# uz.low <- uz.lows[[min.loc]]
# uz.mid <- uz.mids[[mid.loc]]
# uz.high <- uz.highs[[max.loc]]
# 
# 
# ua.lows <- list()
# ua.mids <- list()
# ua.highs <- list()
# infos <- data.frame(matrix(0, nrow = 0, ncol = 3))
# 
# for (i in 1:100){
#   ua.low <- sigmaSim(k, 0, 1)[[1]]
#   ua.mid <- sigmaSim(k, 0, 2)[[1]]
#   ua.high <- sigmaSim(k, 0, 3)[[1]]
#   
#   ua.lows[[i]] <- ua.low
#   ua.mids[[i]] <- ua.mid
#   ua.highs[[i]] <- ua.high
#   
#   info.az <- c(ua.low %*% t(ua.low) + 1,
#                ua.mid %*% t(ua.mid) + 1,
#                ua.high %*% t(ua.high) + 1
#   )
#   info.az <- .5*log(info.az)
#   infos[i,] <- info.az
# }
# 
# min.info <- min(infos[,1])
# max.info <- max(infos[,3])
# mid.info <- (min.info + max.info)/2
# 
# min.loc <- which(infos[,1] == min.info)
# max.loc <- which(infos[,3] == max.info)
# mid.difs <- abs(infos[,2] - mid.info)
# mid.loc <- which(mid.difs == min(mid.difs))
# 
# ua.low <- ua.lows[[min.loc[[1]]]]
# ua.mid <- ua.mids[[mid.loc[[1]]]]
# ua.high <- ua.highs[[max.loc[[1]]]]
# 
# alpha <- c(.5, -.8)
# get.gamma <- function(alpha, gamma, b, target){
#    numerat <- alpha[2]^2 + 2*gamma*alpha[2]*b + gamma^2*b^2# - .1 * (target - 1)
#    const <- 2*(target-1)*b^2
#    sqrt(numerat/const)
# }
# 
# lambdas <- list('low' = uz.low[[1]],
#                 'mid' = uz.mid[[1]],
#                 'high' = uz.high[[1]])
# psis <- list('low' = uz.low[[2]],
#                 'mid' = uz.mid[[2]],
#                 'high' = uz.high[[2]])
# 
# betas <- list('low' = ua.low,
#              'mid' = ua.mid,
#              'high' = ua.high)
# gammas <- data.frame(matrix(0, nrow = 0, ncol = 3))
# i=1
# for (alev in c('low', 'mid', 'high')){
#   for (ylev in c(10, 20, 30)){
#     gammas[i,] <- c(alev, ylev, get.gamma(alpha, .7, betas[[alev]], ylev))
#     i = i+1
#   }
# }
# 
# for (zlev in c('low', 'mid', 'high')){
#   for (alev in c('low', 'mid', 'high')){
#     for (ylev in c(10, 20, 30)){
#       lambda <- lambdas[[zlev]]
#       psi <- psis[[zlev]]
#       b <- betas[[alev]]
#       gamma <- c(.7, get.gamma(alpha, .7, b, ylev))
# 
#       save(k, p, lambda, psi, alpha, gamma, b,
#           file = file.path("Data", "Parameters",
#           paste("param_z", zlev, "_a", alev, "_y", ylev, ".RData", sep = "")))
#     }
#   }
# }
