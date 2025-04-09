test.consistency <- function(tag, p, sample.sizes, kvals, savemarker = 100, reps = 100) {
  indexing <- 0
  for (sample.size in sample.sizes) {
    set.seed(127+sample.size/100)
    indexing <- indexing + 1
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
      print(j)
      i <- j %% savemarker + 1

      azGen(tag, sample.size)

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

      ate.est <- latent.ATE(raw.data, kvals, p)$estimate
      latent[i] <- ate.est

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
}

graph.consistency <- function(tag, p, sample.sizes) {
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
    mean.ests[5 * j - 4, ] <- c(
      samp.size, "linear", median(temp.df$linear, na.rm = TRUE),
      quantile(temp.df$linear, probs = c(.25, .75), na.rm = TRUE)
    )
    mean.ests[5 * j - 3, ] <- c(
      samp.size, "IPW", median(temp.df$ipw, na.rm = TRUE),
      quantile(temp.df$ipw, probs = c(.25, .75), na.rm = TRUE)
    )
    mean.ests[5 * j, ] <- c(
      samp.size, "latent", median(temp.df$latent, na.rm = TRUE),
      quantile(temp.df$latent, probs = c(.25, .75), na.rm = TRUE)
    )
    mean.ests[5 * j - 1, ] <- c(
      samp.size, "IV", median(temp.df$iv, na.rm = TRUE),
      quantile(temp.df$iv, probs = c(.25, .75), na.rm = TRUE)
    )
    mean.ests[5 * j - 2, ] <- c(
      samp.size, "proximal", median(temp.df$proximal, na.rm = TRUE),
      quantile(temp.df$proximal, probs = c(.25, .75), na.rm = TRUE)
    )
  }

  colnames(mean.ests) <- c("SampleSize", "Method", "Mean", "Q.25", "Q.75")
  mean.ests$Method <- factor(mean.ests$Method)
  mean.ests$SampleSize <- as.numeric(mean.ests$SampleSize)
  mean.ests$Mean <- as.numeric(mean.ests$Mean)
  mean.ests$Q.25 <- as.numeric(mean.ests$Q.25)
  mean.ests$Q.75 <- as.numeric(mean.ests$Q.75)

  mean.ests <- mean.ests[order(mean.ests$SampleSize), ]
  load(file.path(
    "Data", "Parameters",
    paste("param_", tag, ".RData", sep = "")
  ))
  true_ate <- gamma[1]
  main.plot1 <- ggplot(mean.ests) +
    geom_ribbon(aes(
      x = SampleSize, ymin = Q.25, ymax = Q.75, fill = Method),
      alpha = .3
    ) +
    geom_hline(yintercept = true_ate) +
    geom_line(aes(x = SampleSize, y = Mean, group = Method, color = Method),
              linewidth = 2
    ) +
    geom_point(aes(x = SampleSize, y = Mean, color = Method), size = 3) +
    xlab("Sample Size") +
    ylab("ATE Estimate") + guides(alpha = "none") +
    theme(legend.position = c(.87, .5))
  ggsave(file.path("Data", "Figures", paste("Errs_by_sample_size_",
    tag, ".png",
    sep = ""
  )), plot = main.plot1)
  
  mean.ests <- mean.ests[mean.ests$Method != 'IV', ]
  main.plot2 <- ggplot(mean.ests) +
    geom_ribbon(aes(
      x = SampleSize, ymin = Q.25, ymax = Q.75, fill = Method),
      alpha = .3
    ) +
    geom_hline(yintercept = true_ate) +
    geom_line(aes(x = SampleSize, y = Mean, group = Method, color = Method),
                  linewidth = 2
    ) +
    geom_point(aes(x = SampleSize, y = Mean, color = Method), size = 3) +
    xlab("Sample Size") +
    ylab("ATE Estimate") + guides(alpha = "none") +
    theme(legend.position = c(.87, .5))
  ggsave(file.path("Data", "Figures", paste("Errs_by_sample_size_",
                                            tag, "_noIV.png",
                                            sep = ""
  )), plot = main.plot2)
  return(list(main.plot1, main.plot2))
}
