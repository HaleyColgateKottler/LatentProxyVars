binary_param_gen <- function(k, p, alpha, gamma, tag){
  factor.loadings <- sigmaSim(k, p, 2)
  lambda <- factor.loadings[[1]]
  psi <- factor.loadings[[2]]
  psi[p+1, p+1] <- 1
  
  var_Y <- .1
  
  save(k, p, var_Y, lambda, psi, alpha, gamma,
       file = file.path("Data", "Parameters",
                        paste("param_", tag,
                              ".RData", sep = "")))
}

binary_az_gen <- function(tag, sample.size){
  load(file.path("Data", "Parameters",
                 paste("param_", tag, ".RData", sep = "")))
  
  H <- rmvnorm(sample.size, mean = rep(0, k))
  Z <- t(lambda[1:p,] %*% t(H)) + rmvnorm(sample.size, sigma = psi[1:p, 1:p])
  colnames(Z) <- paste0("V", 1:p)
  A.temp <- t(lambda[p+1,] %*% t(H)) + rnorm(sample.size)
  A <- ifelse(A.temp > 0, 1, 0)
  Y <- rowSums((alpha + A%*%gamma) * cbind(1, H)) + rnorm(sample.size, sd = var_Y)
  
  write.table(H, file = file.path("Data", "Samples", paste("H_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Z, file = file.path("Data", "Samples", paste("Z_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(A, file = file.path("Data", "Samples", paste("A_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Y, file = file.path("Data", "Samples", paste("Y_", tag, "_", sample.size, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
}

binary_test <- function(kvals, sample.size, reps, tag, savemarker = 100){
  est.df <- data.frame(matrix(nrow = 0, ncol = 4))
  colnames(est.df) <- c("latent", "linear", "IPW", "IV")
  latent <- c()
  ipw <- c()
  linear <- c()
  iv <- c()
  
  for (j in 1:trials) {
    i <- j %% savemarker
    binary_az_gen(tag, sample.size)
    raw.data <- data.import(tag, sample.size)
    
    linear.ate <- linearEst(raw.data)
    linear[i] <- linear.ate
    ipw.ate <- binaryIPWest(raw.data)
    ipw[i] <- ipw.ate
    iv.ate <- IVest(raw.data)
    iv[i] <- iv.ate
    
    ate.est <- latent.ATE(raw.data, kvals, p)
    latent[i] <- ate.est$estimate
    
    if (j %% savemarker == 0 | j == trials) {
      print(sample.size)
      print(j)
      est.df <- rbind(est.df, cbind(latent, linear, ipw, iv))
      write.table(est.df,
                  file.path(
                    "Data", "Estimates",
                    paste('ests_', tag, as.character(sample.size),
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
