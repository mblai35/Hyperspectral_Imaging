# JeffreysPrior.R
# R version 3.4.3 (2017-11-30)
# July 12, 2018. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Multiple linear regression using adaptive shrinkage with Jeffrey's 
# prior. Code adapted from: 
# http://evolvedmicrobe.com/Literature/2009_Review%20of%20Bayesian%20Variable%20selection%20methods.pdf

#-----------------------------------------------------------------------
require(rjags)
library(data.table)
#-----------------------------------------------------------------------

# Set working directory.
#setwd("/Users/mblai/Documents/Thesis/Hyperspectral Imaging")

# Import .csv files. 
chl <- fread('Chlorophyll.csv', skip = 2, stringsAsFactors = T)

# Subset 80% of data for training. 
rand <- sample(round(dim(chl)[1] * 0.8))
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

chlorophyll[i]   ~ dnorm(mu[i], tau)
mu[i] <- alpha + inprod(wavelengths[i,], beta[])

}

# Prior for beta
for(j in 1:p){

InTau[j] ~ dunif(-50, 50)
TauM[j] <- exp(InTau[j])
beta[j] ~ ddexp(0, TauM[j])
Ind[j] <- step(abs(beta[j]) - 0.05)

}

# Prior for the inverse variance
tau   ~ dgamma(0.0001, 0.0001)
alpha     ~ dnorm(0, 0.000001)


}"


model <- jags.model(textConnection(model_string), 
                    data = list(chlorophyll = chlorophyll,
                                n = n, p = p, 
                                wavelengths = wavelengths))

update(model, 1000)

samp <- coda.samples(model, 
                     variable.names = c("beta", "Ind"), 
                     n.iter = 2000)

summary(samp)
saveRDS(samp, file = "samp.rds")

# Extract the MCMC samples from each fit:
s <- samp[[1]]
saveRDS(s, file = "s.rds")

