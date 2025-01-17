ratio_param_gen <- function(k, pvals, alpha, gamma, tag) {
  max.p <- max(pvals)
  factor.loadings <- sigmaSim(k, max.p, 1)
  lambda.full <- factor.loadings[[1]]
  psi.full <- factor.loadings[[2]]
  var_Y <- .1
  H_covar <- diag(k)
  A_intercept <- 0
  for (p in pvals) {
    lambda <- lambda.full[1:(p + 1), ]
    lambda[p + 1, ] <- lambda.full[max.p + 1, ]
    psi <- psi.full[1:(p + 1), 1:(p + 1)]
    psi[p + 1, p + 1] <- psi.full[max.p, max.p]
    Z_intercept <- rep(1, p)
    save(k, p, var_Y, lambda, psi, alpha, gamma, H_covar,
      Z_intercept, A_intercept,
      file = file.path(
        "Data", "Parameters",
        paste("param_", tag, p,
          ".RData",
          sep = ""
        )
      )
    )
  }
}

ratio_test <- function(kvals, pvals, sample.size, reps, tag, savemarker = 100) {
  for (p in pvals) {
    tag.new <- paste(tag, p, sep = "")
    file.name <- file.path(
      "Data", "Estimates",
      paste("ests_", tag.new, as.character(sample.size),
            ".csv",
            sep = ""
      )
    )
    if (file.exists(file.name)){
      est.df <- read.csv(file.name)
    } else {
      est.df <- data.frame(matrix(nrow = 0, ncol = 5))
      colnames(est.df) <- c("latent", "linear", "IPW", "IV", "proximal")
    }
    latent <- c()
    ipw <- c()
    linear <- c()
    iv <- c()
    proximal <- c()

    for (j in 1:reps) {
      i <- j %% savemarker + 1

      azGen(tag.new, sample.size)

      raw.data <- data.import(tag.new, sample.size)

      linear.ate <- linearEst(raw.data)
      linear[i] <- linear.ate
      ipw.ate <- continuousIPWest(raw.data)
      ipw[i] <- ipw.ate
      iv.ate <- IVest(raw.data)
      iv[i] <- iv.ate
      proximal.ate <- proximal_causal(raw.data, paste("V", 1:floor(p/2), sep = ""),
                                      paste("V", (floor(p/2)+1):p, sep = ""))
      proximal[i] <- proximal.ate

      ate.est <- latent.ATE(raw.data, kvals[kvals <= p], p)
      latent[i] <- ate.est$estimate

      if (j %% savemarker == 0 | j == reps) {
        print(p)
        print(j)
        est.df <- rbind(est.df, cbind(latent, linear, ipw, iv, proximal))
        write.table(est.df,
          file.path(
            "Data", "Estimates",
            paste("ests_", tag.new, as.character(sample.size),
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

graph.ratio <- function(tag, pvals) {
  mean.ests <- data.frame(matrix(nrow = 0, ncol = 5))
  j <- 0
  for (p in pvals) {
    j <- j + 1
    temp.df <- read.csv(
      file.path(
        "Data", "Estimates",
        paste("ests_", tag, p, sample.size,
          ".csv",
          sep = ""
        )
      ),
      colClasses = "numeric"
    )
    
    mean.ests[5 * j - 4, ] <- c(
      p/k, "linear", mean(temp.df$linear, na.rm = TRUE),
      quantile(temp.df$linear, probs = c(.05, .95), na.rm = TRUE)
    )
    mean.ests[5 * j - 3, ] <- c(
      p/k, "IPW", mean(temp.df$ipw, na.rm = TRUE),
      quantile(temp.df$ipw, probs = c(.05, .95), na.rm = TRUE)
    )
    mean.ests[5 * j, ] <- c(
      p/k, "latent", mean(temp.df$latent, na.rm = TRUE),
      quantile(temp.df$latent, probs = c(.05, .95), na.rm = TRUE)
    )
    mean.ests[5 * j - 1, ] <- c(
      p/k, "IV", mean(temp.df$iv, na.rm = TRUE),
      quantile(temp.df$iv, probs = c(.05, .95), na.rm = TRUE)
    )
    mean.ests[5 * j - 2, ] <- c(
      p/k, "proximal", mean(temp.df$proximal, na.rm = TRUE),
      quantile(temp.df$proximal, probs = c(.05, .95), na.rm = TRUE)
    )
  }

  colnames(mean.ests) <- c("pk", "Method", "Mean", "Q.05", "Q.95")
  mean.ests$Method <- factor(mean.ests$Method)
  mean.ests$pk <- as.numeric(mean.ests$pk)
  mean.ests$Mean <- as.numeric(mean.ests$Mean)
  mean.ests$Q.05 <- as.numeric(mean.ests$Q.05)
  mean.ests$Q.95 <- as.numeric(mean.ests$Q.95)

  mean.ests <- mean.ests[order(mean.ests$pk), ]
  load(file.path(
    "Data", "Parameters",
    paste("param_", tag, p, ".RData", sep = "")
  ))
  true_ate <- gamma[1]
  p1 <- ggplot(mean.ests) +
    geom_hline(yintercept = true_ate) +
    geom_ribbon(aes(
      x = pk, ymin = Q.05, ymax = Q.95, fill = Method,
      alpha = .05
    )) +
    geom_line(aes(x = pk, y = Mean, group = Method, color = Method),
      linewidth = 2
    ) +
    geom_point(aes(x = pk, y = Mean, color = Method), size = 3) +
    xlab("p/k") +
    ylab("ATE Estimate") + guides(alpha = "none") +
    ggtitle("Ratio of p to k") +
    theme(legend.position = c(.87, .5)) +
    coord_cartesian(ylim = c(-2, 2)) +
    scale_x_continuous(breaks = 1:10)
  ggsave(file.path("Data", "Figures", paste(tag, ".png",
    sep = ""
  )), p1)
  
  mean.ests <- mean.ests[mean.ests$Method != 'IV', ]
  p2 <- ggplot(mean.ests) +
    geom_hline(yintercept = true_ate) +
    geom_ribbon(aes(
      x = pk, ymin = Q.05, ymax = Q.95, fill = Method,
      alpha = .05
    )) +
    geom_line(aes(x = pk, y = Mean, group = Method, color = Method),
              linewidth = 2
    ) +
    geom_point(aes(x = pk, y = Mean, color = Method), size = 3) +
    xlab("p/k") +
    ylab("ATE Estimate") + guides(alpha = "none") +
    ggtitle("Ratio of p to k") +
    theme(legend.position = c(.87, .5)) +
    coord_cartesian(ylim = c(-2, 2)) +
    scale_x_continuous(breaks = 1:10)
  ggsave(file.path("Data", "Figures", paste(tag, "_noIV.png",
                                            sep = ""
  )), p2)
  return(list(p1,p2))
}
