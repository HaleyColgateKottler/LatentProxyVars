test.consistency <- function(tag, sample.sizes, kvals, savemarker = 100, reps = 100){
  indexing <- 0
  for (sample.size in sample.sizes) {
    indexing <- indexing + 1
    est.df <- data.frame(matrix(nrow = 0, ncol = 4))
    colnames(est.df) <- c("latent", "linear", "IPW", "IV")
    latent <- c()
    ipw <- c()
    linear <- c()
    iv <- c()
    
    for (j in 1:reps) {
      i <- j %% savemarker
      
      azGen(tag, sample.size)
      
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
}

graph.consistency <- function(tag, sample.sizes){
  mean.ests <- data.frame(matrix(nrow = 0, ncol = 5))
  j <- 0
  for (samp.size in sample.sizes) {
    j <- j + 1
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
    mean.ests[4 * j - 3, ] <- c(
      samp.size, "linear", mean(temp.df$linear),
      quantile(temp.df$linear, probs = c(.05, .95))
    )
    mean.ests[4 * j - 2, ] <- c(
      samp.size, "IPW", mean(temp.df$ipw),
      quantile(temp.df$ipw, probs = c(.05, .95))
    )
    mean.ests[4 * j - 1, ] <- c(
      samp.size, "latent", mean(temp.df$latent),
      quantile(temp.df$latent, probs = c(.05, .95))
    )
    mean.ests[4 * j, ] <- c(
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
  load(file.path("Data", "Parameters",
                 paste("param_", tag, ".RData", sep = "")))
  true_ate <- gamma[1]
  ggplot(mean.ests) +
    geom_ribbon(aes(
      x = SampleSize, ymin = Q.05, ymax = Q.95, fill = Type,
      alpha = .05
    )) +
    geom_line(aes(x = SampleSize, y = Mean, group = Type, color = Type),
              linewidth = 2
    ) +
    geom_point(aes(x = SampleSize, y = Mean, color = Type), size = 3) +
    geom_hline(yintercept = true_ate) +
    annotate("text", x = 215, y = .51, label = "True ATE") +
    xlab("Sample Size") +
    ylab("Average ATE Estimate")
  ggsave(file.path("Data", "Figures", paste("Errs_by_sample_size_",
                                            tag, ".png",
                                            sep = ""
  )))
}
