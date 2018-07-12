# BayesianVariableSelection2.R
# R version 3.4.3 (2017-11-30)
# July 12, 2018. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Multiple linear regression using shrinkage priors. Code adapted
# from: https://www4.stat.ncsu.edu/~reich/ST590/code/

#-----------------------------------------------------------------------
require(rjags)
library(data.table)
#-----------------------------------------------------------------------

# Set working directory.
#setwd("/Users/mblai/Documents/Thesis/Hyperspectral Imaging")

# Import .csv files. 
chl <- fread('Chlorophyll.csv', skip = 2, stringsAsFactors = T)

# Subset 80% of data for training. 
rand <- sample(round(dim(chl)[1]*0.8))
trainChl <- chl[rand, ]
testChl <- chl[-rand, ]

# Subset and standardize response variable. 
chlorophyll <- trainChl$Measured_Chl 
chlorophyll <- (chlorophyll - mean(chlorophyll))/sd(chlorophyll)

# Subset and scale predictors. 
wavelengths     <- trainChl[, 7:1907]
wavelengths     <- scale(wavelengths)

# Define n as the number of rows, or observations. 
n     <- dim(trainChl)[1]
# Define p as the number of predictors, or wavelengths. 
p     <- ncol(wavelengths)


model_string <- "model{

  # Likelihood
  for(i in 1:n){

  chlorophyll[i]   ~ dnorm(mu[i], inv.var)
  mu[i] <- alpha + inprod(wavelengths[i,], beta[])

  }

  # Prior for beta
  for(j in 1:p){

  beta[j] ~ ddexp(0, inv.var.b)

  }

  # Prior for the inverse variance
  inv.var   ~ dgamma(0.01, 0.01)
  inv.var.b ~ dgamma(0.01, 0.01)
  alpha     ~ dnorm(0, 0.01)

}"


model <- jags.model(textConnection(model_string), 
                     data = list(chlorophyll = chlorophyll,
                                 n = n, p = p, 
                                 wavelengths = wavelengths))

update(model, 1000)

samp <- coda.samples(model, 
                      variable.names = c("beta"), 
                      n.iter = 2000)

summary(samp)
saveRDS(samp, file = "samp.rds")

# Extract the MCMC samples from each fit:
s <- samp[[1]]
saveRDS(s, file = "s.rds")

# Plot the posterior for each beta:
par(mfrow = c(5, 3))

for(index in 1:p){
  
  d <- density(s[,index])
  mx <- max(d$y)
  
  plot(d, ylim = c(0, mx),
       xlab = "Beta", 
       ylab = "Posterior density", 
       main = colnames(wavelengths)[index])
  
}
