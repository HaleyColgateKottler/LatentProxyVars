library(dplyr)
library(plotly)

latent = c()
latent.conf = c()
nieve = c()
nieve.conf = c()
IPW = c()
IPW.conf = c()
match = c()
match.conf = c()
linear = c()
linear.conf = c()
sampleSizes = c(1000,2000,3000,4000,5000)
num.trials = c()

for (sampleSize in sampleSizes){
  tag = paste("1DsampSize",sampleSize,sep="")
  df = read.csv(paste("errors_", tag, ".csv", sep = ""), header=FALSE, col.names = c("latent", "IPW", "linear", "match", "nieve"))
  num.trials = append(num.trials, nrow(df))

  sds = sapply(df, sd)
  means = colMeans(df)
  latent = append(latent, means[1])
  latent.conf = append(latent.conf, sds[1]*1.96/sqrt(nrow(df)))
  IPW = append(IPW, means[2])
  IPW.conf = append(IPW.conf, sds[2]*1.96/sqrt(nrow(df)))
  linear = append(linear, means[3])
  linear.conf = append(linear.conf, sds[3]*1.96/sqrt(nrow(df)))
  match = append(match, means[4])
  match.conf = append(match.conf, sds[4]*1.96/sqrt(nrow(df)))
  nieve = append(nieve, means[5])
  nieve.conf = append(nieve.conf, sds[5]*1.96/sqrt(nrow(df)))
}

if (length(unique(num.trials)) > 1){
  print("unequal number of samples")
  print(num.trials)
}
plot(sampleSizes, linear, type = "o", col = "red", ylim=c(.5,.9), ann=FALSE)
lines(sampleSizes, linear-linear.conf, col="red",lty=2)
lines(sampleSizes, linear+linear.conf, col="red",lty=2)
lines(sampleSizes, nieve, type = "o", col = "blue")
lines(sampleSizes, nieve-nieve.conf, col="blue",lty=2)
lines(sampleSizes, nieve+nieve.conf, col="blue",lty=2)
lines(sampleSizes, match, type = "o", col = "orange")#, ylim=c(.07,.12), ann=FALSE)
lines(sampleSizes, match-match.conf, col="orange",lty=2)
lines(sampleSizes, match+match.conf, col="orange",lty=2)
lines(sampleSizes, IPW, type = "o", col = "green")
lines(sampleSizes, IPW-IPW.conf, col="green",lty=2)
lines(sampleSizes, IPW+IPW.conf, col="green",lty=2)
lines(sampleSizes, latent, type = "o", col = "magenta")
lines(sampleSizes, latent-latent.conf, col="magenta",lty=2)
lines(sampleSizes, latent+latent.conf, col="magenta",lty=2)
lines(c(1000,5000), c(.6,.6))
text(1250, .61, "True ATE")
title(main = "Average ATE Estimate (500 trials)", sub="Dotted lines give 95% confidence intervals", xlab = "Sample Size", ylab = "Average ATE Estimate")
legend("topright",c("linear", "nieve", "matching", "IPW", "latent"), title = "Estimator", lty = 1, col = c("red", "blue", "orange", "green", "magenta"))

# box plots and error bars
hist(df$IPW)
hist(df$latent)

plot(sampleSizes, match, type = "o", col = "orange", ylim=c(.5,.725),ann=FALSE)
lines(sampleSizes, match-match.conf, col="orange",lty=2)
lines(sampleSizes, match+match.conf, col="orange",lty=2)
lines(sampleSizes, IPW, type = "o", col = "green")
lines(sampleSizes, IPW-IPW.conf, col="green",lty=2)
lines(sampleSizes, IPW+IPW.conf, col="green",lty=2)
lines(sampleSizes, latent, type = "o", col = "magenta")
lines(sampleSizes, latent-latent.conf, col="magenta",lty=2)
lines(sampleSizes, latent+latent.conf, col="magenta",lty=2)
lines(c(1000,5000), c(.6,.6))
text(1250, .61, "True ATE")
title(main = "Average ATE Estimate (500 trials)", sub="Dotted lines give 95% confidence intervals", xlab = "Sample Size", ylab = "Average ATE Estimate")
legend("topright",c("matching", "IPW", "latent"), title = "Estimator", lty = 1, col = c("orange", "green", "magenta"))


hist(df$latent, col="magenta",xlim=c(min(df$latent), max(df$IPW)))
hist(df$IPW, add=TRUE, col="green")
lines(c(.6,.6),c(0,150), lwd=5)

sampleSize = 1000
tag = paste("1DsampSize",sampleSize,sep="")
df = read.csv(paste("errors_", tag, ".csv", sep = ""), header=FALSE, col.names = c("latent", "IPW", "linear", "match", "nieve"))
fig1 <- plot_ly() %>%
                 add_histogram(x = df$latent, name = "Latent", nbinsx = 20, marker = list(color = "magenta",
                                                                                          line = list(color = "magenta",
                                                                                                      width = 2))) %>%
                 add_histogram(x = df$IPW, name = "IPW", nbinsx = 20,marker = list(color = "lightgreen",
                                                                                   line = list(color = "lightgreen",
                                                                                               width = 2))) %>%
                 add_trace(x = c(.6,.6), y = c(0,220), name = 'True ATE',mode = 'lines', line=list(color="black", width=4)) %>%
                 layout(title = paste("ATE Estimates with 1000 Samples"),
                        xaxis = list(title = "ATE Estimate"),
                        yaxis = list(title = "Frequency"))

sampleSize = 2000
tag = paste("1DsampSize",sampleSize,sep="")
df = read.csv(paste("errors_", tag, ".csv", sep = ""), header=FALSE, col.names = c("latent", "IPW", "linear", "match", "nieve"))
fig2 <- plot_ly() %>%
  add_histogram(x = df$latent, name = "Latent", nbinsx = 20, marker = list(color = "magenta",
                                                                           line = list(color = "magenta",
                                                                                       width = 2))) %>%
  add_histogram(x = df$IPW, name = "IPW", nbinsx = 20,marker = list(color = "lightgreen",
                                                                    line = list(color = "lightgreen",
                                                                                width = 2))) %>%
  add_trace(x = c(.6,.6), y = c(0,220), name = 'True ATE',mode = 'lines', line=list(color="black", width=4)) %>%
  layout(title = paste("ATE Estimates with 2000 Samples"),
         xaxis = list(title = "ATE Estimate"),
         yaxis = list(title = "Frequency"))

sampleSize = 3000
tag = paste("1DsampSize",sampleSize,sep="")
df = read.csv(paste("errors_", tag, ".csv", sep = ""), header=FALSE, col.names = c("latent", "IPW", "linear", "match", "nieve"))
fig3 <- plot_ly() %>%
  add_histogram(x = df$latent, name = "Latent", nbinsx = 20, marker = list(color = "magenta",
                                                                           line = list(color = "magenta",
                                                                                       width = 2))) %>%
  add_histogram(x = df$IPW, name = "IPW", nbinsx = 20,marker = list(color = "lightgreen",
                                                                    line = list(color = "lightgreen",
                                                                                width = 2))) %>%
  add_trace(x = c(.6,.6), y = c(0,220), name = 'True ATE',mode = 'lines', line=list(color="black", width=4)) %>%
  layout(title = paste("ATE Estimates with 3000 Samples"),
         xaxis = list(title = "ATE Estimate"),
         yaxis = list(title = "Frequency"))

sampleSize = 4000
tag = paste("1DsampSize",sampleSize,sep="")
df = read.csv(paste("errors_", tag, ".csv", sep = ""), header=FALSE, col.names = c("latent", "IPW", "linear", "match", "nieve"))
fig4 <- plot_ly() %>%
  add_histogram(x = df$latent, name = "Latent", nbinsx = 20, marker = list(color = "magenta",
                                                                           line = list(color = "magenta",
                                                                                       width = 2))) %>%
  add_histogram(x = df$IPW, name = "IPW", nbinsx = 20,marker = list(color = "lightgreen",
                                                                    line = list(color = "lightgreen",
                                                                                width = 2))) %>%
  add_trace(x = c(.6,.6), y = c(0,220), name = 'True ATE',mode = 'lines', line=list(color="black", width=4)) %>%
  layout(title = paste("ATE Estimates with 4000 Samples"),
         xaxis = list(title = "ATE Estimate"),
         yaxis = list(title = "Frequency"))


sampleSize = 5000
tag = paste("1DsampSize",sampleSize,sep="")
df = read.csv(paste("errors_", tag, ".csv", sep = ""), header=FALSE, col.names = c("latent", "IPW", "linear", "match", "nieve"))
fig5 <- plot_ly() %>%
  add_histogram(x = df$latent, name = "Latent", nbinsx = 20, marker = list(color = "magenta",
                                                                           line = list(color = "magenta",
                                                                                       width = 2))) %>%
  add_histogram(x = df$IPW, name = "IPW", nbinsx = 20,marker = list(color = "lightgreen",
                                                                    line = list(color = "lightgreen",
                                                                                width = 2))) %>%
  add_trace(x = c(.6,.6), y = c(0,220), name = 'True ATE',mode = 'lines', line=list(color="black", width=4)) %>%
  layout(title = paste("ATE Estimates with 5000 Samples"),
         xaxis = list(title = "ATE Estimate"),
         yaxis = list(title = "Frequency"))
fig <- subplot(fig1,fig5) %>% 
  
  layout(title = "ATE Estimates by Sample Size (500 Trials)")

fig
