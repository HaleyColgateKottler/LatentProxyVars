sigmaSim <- function(k, p, communality) {
  # p number of variables
  # k number of factors
  # communality 1 low, 2 wide, 3 high
  
  p <- p + 1
  
  A <- matrix(0, p, k)
  
  # need row entries to sum to k-1
  A[, 1] <- sample(0:(k - 1), p, replace = TRUE)
  
  for (i in 1:p) {
    if (k - 2 < 1) break
    for (j in 2:(k - 2)) {
      current_sum <- sum(A[i, ])
      if (current_sum == k - 1) {
        A[i, j] <- 0
      } else {
        A[i, j] <- sample(0:(k - 1 - current_sum), 1, replace = TRUE)
      }
    }
  }
  A[, k] <- (k - 1) * matrix(1, 1, p) - rowSums(A)
  
  # Add normal deviation
  c <- .1 * sample(7:9, k)
  x <- rnorm(p * k)
  x1 <- matrix(x, nrow = p)
  d <- c(matrix(1, ncol = p, nrow = 1) / wordspace::rowNorms(x1))
  
  Y <- A * c + d * x1 * sqrt(1 - c**2)
  
  # Apply skewing function
  Y2 <- Y + abs(Y) + 0.2
  Y3 <- abs(Y) + 0.2
  Z <- (1.2 / 2.2) * (Y * Y2) / Y3
  
  g <- rep(1, p) / wordspace::rowNorms(Z)
  
  # Scale to set communality
  if (communality == 1) {
    B1 <- diag(.1 * sample(2:4, p, replace = TRUE), p)
  } else if (communality == 2) {
    B1 <- diag(.1 * sample(2:8, p, replace = TRUE), p)
  } else if (communality == 3) {
    B1 <- diag(.1 * sample(6:8, p, replace = TRUE), p)
  } else {
    B1 <- matrix(0, p, p)
  }
  
  B2 <- diag(p) - B1
  
  # Final factor loading matrix
  lambda <- c(sqrt(B1) %*% g) * Z
  psi <- sqrt(B2)
  
  return(list(lambda, psi))
}

saveParams <- function(k, p, communality, var_Y, A_intercept, Z_intercept,
                       H_covar, alpha, gamma, tag) {
  output <- sigmaSim(k, p, communality)
  lambda <- output[[1]]
  psi <- output[[2]]
  
  save(k, p, communality, var_Y, lambda, psi, H_covar, alpha, gamma,
       A_intercept, Z_intercept,
       file = file.path("Data", "Parameters",
                        paste("param_", tag, ".RData", sep = "")))
}

azGen <- function(tag, n) {
  load(file.path("Data", "Parameters",
                 paste("param_", tag, ".RData", sep = "")))
  H_set <- matrix(, nrow = n, ncol = k)
  Z_set <- matrix(, nrow = n, ncol = p)
  A_set <- c()
  Y_set <- c()
  
  for (i in 1:n) {
    epsilon_AZ <- MASS::mvrnorm(n = 1, mu = rep(0, p + 1), Sigma = psi)
    epsilon_Z <- epsilon_AZ[1:p]
    epsilon_A <- epsilon_AZ[p + 1]
    epsilon_Y <- rnorm(1, mean = 0, sd = sqrt(var_Y))
    
    H <- MASS::mvrnorm(n = 1, mu = rep(0, k), Sigma = H_covar)
    ZA <- lambda %*% H
    Z <- ZA[1:p, 1] + Z_intercept + epsilon_Z
    A <- ZA[p + 1, 1] + A_intercept + epsilon_A
    Y <- (alpha + A*gamma) %*% c(1, H) + epsilon_Y
    
    H_set[i, ] <- t(H)
    Z_set[i, ] <- Z
    A_set[i] <- A
    Y_set[i] <- Y
  }
  
  write.table(H_set, file = file.path("Data", "Samples", paste("H_", tag, "_", n, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Z_set, file = file.path("Data", "Samples", paste("Z_", tag, "_", n, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(A_set, file = file.path("Data", "Samples", paste("A_", tag, "_", n, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Y_set, file = file.path("Data", "Samples", paste("Y_", tag, "_", n, ".csv", sep = "")), row.names = FALSE, col.names = FALSE, sep = ",")
}

# dataGen(1, 3, 3, .5,1000, -.3, c(.2, -.1, .6), 0,
# c(1), c(.9), .6, c(.8), c(-.2), tag)
dataGen <- function(k, p, communality, var_Y, n, A_intercept, Z_intercept,
                    H_covar, alpha, gamma, tag) {
  saveParams(k, p, communality, var_Y, A_intercept, Z_intercept,
             H_covar, alpha, gamma, tag)
  azGen(tag, n)
}

# dataGen(2, 6, 2, .5, 1000, .3, c(.5,.25, .3, -.1, -.2, .4), 2,
# matrix(c(1,.5,.5,1), nrow = 2, byrow = TRUE), 2, c(4,1), "2D")

# dataGen(1, 3, 1, .25, 500, -.2, c(.2, -.1, .4), 2,
# c(1), .5, c(-.8, .5, .3), "1Da")
dataGenNCE <- function(k, p, communality, var_Y, n, A_intercept, Z_intercept,
                       Y_intercept, H_covar, C0, C1, C2, tag) {
  saveParams(k, p, communality, var_Y, A_intercept, Z_intercept, Y_intercept,
             H_covar, C0, C1, C2, tag)
  azGenNCE(tag, n)
}

azGenNCE <- function(tag, n) {
  load(file.path("Data", "Parameters",
                 paste("param_", tag, ".RData", sep = "")))
  H_set <- matrix(, nrow = n, ncol = k)
  Z_set <- matrix(, nrow = n, ncol = p)
  A_set <- c()
  Y_set <- c()
  
  for (i in 1:n) {
    epsilon_AZ <- MASS::mvrnorm(n = 1, mu = rep(0, p + 1), Sigma = psi)
    epsilon_Z <- epsilon_AZ[1:p]
    epsilon_A <- epsilon_AZ[p + 1]
    epsilon_Y <- rnorm(1, mean = 0, sd = sqrt(var_Y))
    
    H <- MASS::mvrnorm(n = 1, mu = rep(0, k), Sigma = H_covar)
    ZA <- lambda %*% H
    Z <- ZA[1:p, 1] + Z_intercept + epsilon_Z
    A <- ZA[p + 1, 1] + (C2 %*% Z) + A_intercept + epsilon_A
    Y <- (C1 %*% H) + C0 * A + Y_intercept + epsilon_Y
    
    H_set[i, ] <- t(H)
    Z_set[i, ] <- Z
    A_set[i] <- A
    Y_set[i] <- Y
  }
  
  file.tag = paste(tag, "_", n, ".csv", sep = "")
  write.table(H_set, file = file.path("Data", "Samples",
                                      paste("H_", file.tag, sep = "")),
              row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Z_set, file = file.path("Data", "Samples",
                                      paste("Z_", file.tag, sep = "")),
              row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(A_set, file = file.path("Data", "Samples",
                                      paste("A_", file.tag, sep = "")),
              row.names = FALSE, col.names = FALSE, sep = ",")
  write.table(Y_set, file = file.path("Data", "Samples",
                                      paste("Y_", file.tag, sep = "")),
              row.names = FALSE, col.names = FALSE, sep = ",")
}
