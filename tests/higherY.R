squared_param_gen <- function(k, p, alpha, gamma, tag) {
  factor.loadings <- sigmaSim(k, p, 3)
  lambda <- factor.loadings[[1]]
  psi <- factor.loadings[[2]]

  var_Y <- .1

  save(k, p, var_Y, lambda, psi, alpha, gamma,
    file = file.path(
      "Data", "Parameters",
      paste("param_", tag,
        ".RData",
        sep = ""
      )
    )
  )
}

# squared_param_gen(k, p, alpha, c(gamma[1:2], -.8), tag)
squared_az_gen <- function(tag, sample.size) {
  load(file.path(
    "Data", "Parameters",
    paste("param_", tag, ".RData", sep = "")
  ))
  df <- data.frame(matrix(0, nrow = 0, ncol = k + p + 2))
  for (samp in 1:sample.size) {
    epsilon_AZ <- MASS::mvrnorm(n = 1, mu = rep(0, p + 1), Sigma = psi)
    epsilon_Z <- epsilon_AZ[1:p]
    epsilon_A <- epsilon_AZ[p + 1]
    epsilon_Y <- rnorm(1, mean = 0, sd = sqrt(var_Y))

    H <- MASS::mvrnorm(n = 1, mu = rep(0, k), Sigma = diag(k))
    ZA <- lambda %*% H
    Z <- ZA[1:p, 1] + epsilon_Z
    A <- ZA[p + 1, 1] + epsilon_A
    Y <- (alpha + A * gamma) %*% c(1, H, H^2) + epsilon_Y

    df[samp, ] <- c(H, Z, A, Y)
  }
  colnames(df) <- c(paste0("H", 1:k), paste0("V", 1:p), "A", "Y")
  write.table(df[, paste0("H", 1:k)], file = file.path("Data", "Samples", paste("H_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(df[, paste0("V", 1:p)], file = file.path("Data", "Samples", paste("Z_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(df$A, file = file.path("Data", "Samples", paste("A_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(df$Y, file = file.path("Data", "Samples", paste("Y_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
}

squared_test <- function(kvals, p, sample.size, reps, tag, savemarker = 100) {
  file.name <- file.path(
    "Data", "Estimates",
    paste("ests_", tag, as.character(sample.size),
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
    squared_az_gen(tag, sample.size)
    raw.data <- data.import(tag, sample.size)

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

    if (j %% savemarker == 0 | j == reps) {
      print(sample.size)
      print(j)
      est.df <- rbind(est.df, cbind(latent, linear, ipw, iv, proximal))
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
      proximal <- c()
    }
  }
}
