library(dplyr)

latent.errors = c()
latent.est = c()
nieve = c()
IPW.est = c()
IPW.errors = c()
match = c()
linear = c()
sampleSizes = c(1000,2000,3000,4000,5000)

for (sampleSize in sampleSizes){
  tag = paste("1DsampSize",sampleSize,sep="")
  df = read.csv(paste("errors_", tag, ".csv", sep = ""))
  
  stats = colMeans(df)
  latent.est = append(latent.est, stats[1])
  latent.errors = append(latent.errors, stats[2])
  IPW.est = append(IPW.est, stats[3])
  IPW.errors = append(IPW.errors, stats[4])
  linear = append(linear, stats[5])
  match = append(match, stats[6])
  nieve = append(nieve, stats[7])
}
plot(sampleSizes, linear, type = "o", col = "red", ylim=c(.07,.3), ann=FALSE)
lines(sampleSizes, nieve, type = "o", col = "blue")
lines(sampleSizes, match, type = "o", col = "orange")#, ylim=c(.07,.12), ann=FALSE)
lines(sampleSizes, IPW.errors, type = "o", col = "green")
lines(sampleSizes, latent.errors, type = "o", col = "magenta")
title(main = "Average Absolute ATE Error (500 trials)", xlab = "Sample Size", ylab = "Average Error")
legend("topright",c("linear", "nieve", "matching", "IPW", "latent"), title = "Estimator", lty = 1, col = c("red", "blue", "orange", "green", "magenta"))

# box plots and error bars
hist(df$IPW.est)
hist(df$latent.est)

plot(sampleSizes, match, type = "o", col = "orange", ylim=c(.07,.11), ann=FALSE)
lines(sampleSizes, IPW.errors, type = "o", col = "green")
lines(sampleSizes, latent.errors, type = "o", col = "magenta")
title(main = "Average Absolute ATE Error (500 trials)", xlab = "Sample Size", ylab = "Average Error")
legend("topright",c("matching", "IPW", "latent"), title = "Estimator", lty = 1, col = c("orange", "green", "magenta"))

