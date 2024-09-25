no_latent_param_gen <- function(p, alpha, gamma, tag){
  factor.loadings <- sigmaSim(p, 0, 3)
  lambda <- sigmaSim(p, 0, 3)[[1]]
  psi <- sigmaSim(2, p, 3)[[2]][1:p, 1:p]
  
  var_Y <- .1
  
  save(p, var_Y, lambda, psi, alpha, gamma,
       file = file.path("Data", "Parameters",
                        paste("param_", tag,
                              ".RData", sep = "")))
}

no_latent_az_gen <- function(tag, sample.size){
  load(file.path("Data", "Parameters",
                 paste("param_", tag, ".RData", sep = "")))
  
  Z <- rmvnorm(sample.size, sigma = psi[1:p, 1:p])
  colnames(Z) <- paste0("V", 1:p)
  A <- t(lambda %*% t(Z)) + rnorm(sample.size)
  Y <- rowSums((alpha + A%*%gamma) * cbind(1, Z)) + rnorm(sample.size, sd = var_Y)
  
  write.table(Z, file = file.path("Data", "Samples", paste("Z_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(A, file = file.path("Data", "Samples", paste("A_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Y, file = file.path("Data", "Samples", paste("Y_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
}

no_latent_test <- function(kvals, sample.size, reps, tag, savemarker = 100){
  est.df <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(est.df) <- c("latent", "linear", "IPW", "IV")
  latent <- c()
  ipw <- c()
  linear <- c()
  iv <- c()
  reps = trials
  for (j in 1:reps) {
    i <- j %% savemarker

    no_latent_az_gen(tag, sample.size)
    raw.data <- data.import(tag, sample.size)

    linear.ate <- linearEst(raw.data)
    linear[i] <- linear.ate
    ipw.ate <- continuousIPWest(raw.data)
    ipw[i] <- ipw.ate
    iv.ate <- IVest(raw.data)
    iv[i] <- iv.ate
    
    ate.est <- latent.ATE(raw.data, kvals, p)
    latent[i] <- ate.est$estimate
    
    if (j %% savemarker == 0 | j == reps) {
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
