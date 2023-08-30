
# Single-species, single-season occupancy model
# This species distribution model accounts for imperfect detection, hence produces a true distribution
# Ugyen Penjor, 2018 ugyenpenjor.bt@gmail.com
# Link to relevant document: 10.13140/RG.2.2.31667.91689 
# Nature Conservation Division, Department of Forests and Park Services, Bhutan
# Distribution and winter habitat use of Bhutan Takin, Burdocas taxicolor whitei, in Bhutan

rm(list=ls())
ls()

# Load packages
library(jagsUI)
library(wiqid)

# Load data
load("single_species_single_season_occupancy_model_data.RData")
ls()

# Bundle data for JAGS
jagsDataL <- list(y=y1, n=n, nSites=length(n),
                  mea=as.vector(scale(meadow)), con=as.vector(scale(conifer)), cti=as.vector(scale(cti)), roa=as.vector(scale(road)), 
                  roug=as.vector(scale(roughness)), slo=as.vector(scale(slope)), sno=as.vector(scale(snow)), sno2=as.vector(scale(snow))^2,
                  wea=wea, veg=veg)
str(jagsDataL)

# Occupany model in BUGS language

# JAGS model
modelText <- "model{
  # Likelihood
  for(i in 1:nSites) {

    # Ecological model
    z[i] ~ dbern(psi[i])
      logit(psi[i]) <- b0 + bMea * mea[i] + bCon * con[i] + bCti * cti[i] + bRoa * roa[i] + bRou * roug[i] + bSlo * slo[i] +
                       bSno * sno[i] + bSno2 * sno2[i]

    # Observation model
    y[i] ~ dbin(z[i] * p[i], n[i]) 
      logit(p[i]) <- a0 + aWea[wea[i]] + aVeg[veg[i]]
      logLik[i] <- log(psi[i] * dbin(y[i], p[i], n[i]) + step(-y[i]) * (1-psi[i]))
  } # i

  # Priors
  # Priors on detection intercept and covariates (broad uniform priors)
  a0 ~ dunif(-10, 10)   
  aWea[1]  <- 0 # Reference category - sunny weather (observations mostly made on sunnay day)
  for (k in 2:4) {
    aWea[k] ~ dunif(-10, 10) 
  } 
  aVeg[1]  <- 0 # Reference category - broadleaf forest (BL): likewise most observation plots in broadleaved forest
  for (k in 2:6) { 
    aVeg[k] ~ dunif(-10, 10) 
  } 

  # Priors on occupancy intercept and covariates
  b0 ~ dunif(-10, 10)
  bMea ~ dunif(-5, 5)
  bCon ~ dunif(-5, 5)
  bCti ~ dunif(-5, 5)
  bRoa ~ dunif(-5, 5)
  bRou ~ dunif(-5, 5)
  bSlo ~ dunif(-5, 5)
  bSno ~ dunif(-5, 5)
  bSno2 ~ dunif(-5, 5)

  # Derived variable
  Nsite.hat <- sum(z[]) # estimated no. of sites being occupied
  
}"
writeLines(modelText, "Single_species_single_season_model.jags")

# Initial values to select optimum step size during MCMC
inits <- function() list(z=rep(1, length(y1)))    # initial observations given 1: all observed

# Parameters wanted
wanted <- c("b0", "bMea", "bCon", "bCti", "bRoa", "bRou", "bSlo", "bSno", "bSno2", "a0", "aWea", "aVeg", 
            "Nsite.hat", "Tobs", "Tsim", "p", "psi", "logLik")

# Run the model (increase iteration for publication)           
(jagsOut_BM <- jags(jagsDataL, inits, wanted, "Single_species_single_season_model.jags", DIC=FALSE, 
                    n.chains=3, n.iter=10000, n.adapt=1000, n.burnin=5000, n.thin=2, parallel=TRUE, 
                    codaOnly="logLik"))

####################################
########## PLOTS and MAPS ##########
####################################

# Plot occupancy vs covariate
# E.g., slope
range(slope)
mean.cov <- mean(slope, na.rm=T)
sd.cov <- sd(slope[!is.na(slope)])

nsamp <- jagsOut_BM$mcmc.info$n.samples # nsamp = no. of samples 

# Prepare covariate for prediction
orig.pred.cov <- seq(-640.7336, 475.2769,, 101)
p.cov <- (orig.pred.cov - mean.cov) / sd.cov
p.pred.cov <- plogis(jagsOut_BM$mean$b0 + jagsOut_BM$mean$bSlo * p.cov)

plot(orig.pred.cov, p.pred.cov, main="",
     ylab=expression(paste("Probability of habitat use")), 
     xlab="slope position", 
     ylim=c(0, 1), type="l", lwd=0.01, col="grey95", frame.plot=T, axes=F) 

array.p.pred.cov <- array(NA, dim=c(length(p.cov), nsamp))

for(i in 1:nsamp){
  array.p.pred.cov[,i] <- plogis(jagsOut_BM$sims.list$b0[i] + 
                                   jagsOut_BM$sims.list$bSlo[i] * p.cov)
}

sub.set <- sort(sample(1:nsamp, size = 300)) # size changes the # of lines in the plot

# Plot 95% credible intervals
for (i in sub.set){
  matlines(orig.pred.cov, array.p.pred.cov[,i], type="l",
           lwd=1.2, col=adjustcolor("steelblue", 0.3))
}

# Plot mean
lines(orig.pred.cov, p.pred.cov, type="l", lwd=3, col="navy")

# Add axes
axis(side=1, lwd=1, tcl=-0.5)
axis(side=2, lwd=1, tcl=-0.5)

# Prepare species distribution map (this is the detection-corrected map)
# Calculate occupancy for each pixel
logit_psi <- with(jagsOut_BM$mean, 
                  b0 + bMea * meadRS + bCon * conRS + bCti * ctiRS + bRoa * roadRS + bRou * roughRS + bSlo * sloRS + bSno2 * snowRS)

# Transform into probability scale from logit scale (constrain values between 0 and 1)
psiR <- 1/(1 + exp(-logit_psi))

mapPalette <- colorRampPalette(c("grey", "yellow", "orange", "red"))
plot(psiR, col=mapPalette(100), main="Occupancy probability map of Bhutan takin", axes=FALSE, box=FALSE, legend.shrink=0.35,
     legend.args=list(text="Occupancy probability", side=4, line=-2, font=2))

# Save output for post-processing with GIS software (or you can do it in R!)
writeRaster(psiR,"takin_psi_prediction.tif") 

##############################################################################################
######################################### END ################################################
##############################################################################################
