
#### We will be using a package that is not listed in CRAN anymore. 
#### To install it, uncomment the lines below and run them once.
#### It should be good to go after that

#url = "https://cran.r-project.org/src/contrib/Archive/grofit/grofit_1.1.1-1.tar.gz"
#pkgFile = "grofit_1.1.1-1.tar.gz"

#download.file(url = url, destfile = pkgFile)

#install.packages(pkgs = pkgFile, type = 'source', respo = NULL)

#### Load libraries ####
library(scales); library(nlme); library(lme4); library(bbmle); library(mixsmsn); 
library(lmodel2); library(smatr); library(grofit); library(minpack.lm)
source("./sim/code/Functions for analysing dimorphisms and allometries with intercept.R")
source("./sim/code/Functions for plotting all together.R")
source("./sim/code/scatterhist.R")

##### PHAREICRANAUS MANAUARA #####

## Load the dataset
morpho <- read.csv("./sim/data/manauara-morpho 3.csv",h=T)
males <- morpho[morpho$sex=="male",]

# Sometimes the tusk on one side is broken. So, we made this loop to choose the longest tusk of the individual.

for(i in 1:length(males$left.tusk)){
  ifelse(males$left.tusk[i] >= males$right.tusk[i],
         males$tusk[i] <- males$left.tusk[i],
         males$tusk[i] <- males$right.tusk[i])
}

males2 <- males[ , c(2, 15)]
colnames(males2) <- c("body", "trait")

# This is a custom-built function that will spit out the results of the analyses and two graphs.
# This is the original analyses we made for the species. 
Mixmodels.allometries(bi.trait = males2,
                      label.body = "Dorsal scute length (cm)",
                      label.trait = "Tusk length (cm)")


## This is just a handy-function to calculate AICc. 
myAICc <- function(aic, k, n) {
  aic + ((2*(k^2) + 2*k) / (n - k - 1))
}

## Mixture models and AICc for manauara
SN.SkNorm2 <- smsn.mix(y = na.omit(males2$trait), g=2, nu=6, obs.prob=T, family="Skew.normal") # k = 2 + 2 + 2 + 1 = 7
SN.SkNorm1 <- smsn.mix(y = na.omit(males2$trait), g=1, nu=6, obs.prob=T, family="Skew.normal") # k = 1 + 1 + 1 = 3
AIC2dist <- myAICc(aic = SN.SkNorm2$aic, k = 7, n = SN.SkNorm2$n) # k = 2 + 2 + 2 + 1 = 7
AIC1dist <- myAICc(aic = SN.SkNorm1$aic, k = 3, n = SN.SkNorm1$n) # k = 1 + 1 + 1 = 3


###################################################
### Calculating non-linear static allometries  ####
###################################################

### Our mail goal here is to calculate the best possible allometric fit of the data.
### We will not do any sort of inference of the relationship between the variables.
### Our interest is purely phenomenological because we want to use these non-linear 
### allometries to simulate new trait values based on the allometric relationship.
### Thus, we will choose the best fitted model regardless of how much better a model 
### is when compared to alternative models.

## The Richard's growth equation has four parameters (which we won't go into detail here),
## and we need to provide the initial values for them. Luckily, the package has a built-in 
## function that searches for us.
## We're also using the function nlsLM() for this because we it fits the harvestment better
## than the regular nls()

ric.manauara <- nlsLM(trait ~ richards(body, A, mu, lambda, addpar),
                       start=c(A = 0.502, mu = 2, lambda = 0.75, addpar = 0.1), 
                       data = males2, 
                       na.action=na.exclude) 

## Now we fit the other equation, the self-starting non-linear logistic equation.
## Once again, it has multiple parameters, so we're using the "self-starting" quality
## to search for the parameters

logis.manauara =  nls(trait ~ SSlogis(body, Asym, xmid, scal), data = males2)


## Now we're calling the summary just to be sure no errors creep out of the analyses.
summary(ric.manauara)
summary(logis.manauara)

## AICc calculation time.

AICctab(ric.manauara,logis.manauara,
        base = T, weights = T)
## In this case, the logistic equation was 2.2 AICc units better than the fit with 
## Richard's growth equation. Thus, we're picking the logistic equation

## Now, let's see how this plots.
xv <- seq(min(males2$body), max(males2$body), 0.001) 
yv <- predict(logis.manauara, list(body=xv))
plot(males2$trait ~ males2$body,
     xlab= "Dorsal scute length (cm)",
     ylab= "Tusk length (cm)",
     ylim=c(0, max(males2$trait)),
     xlim=c(min(males2$body), max(males2$body)), frame.plot=F)
lines(xv,yv, col="red", lwd=2)

scatterhist(x = males2$body, y = males2$trait,
            xsize = 2, ylab = "Tusk length (cm)",
            xlab = 'Dorsal scute (cm)', cleanup = F)
lines(xv,yv, col=alpha("red",0.6), lwd=3)
dev.off()
#### Simulating new trait values for manauara ####

## We use the mean, SD and N of the body sizes from the original data
xv2 <- sort(rnorm(n = length(males2$body),
                  mean = mean(males2$body),
                  sd = sd(males2$body)))

## Then, using the allometry curve, we find the trait values for those body sizes 
yv2 <- predict(logis.manauara, list(body=xv2))

## This is how they look
plot(xv2,yv2, pch = 10, col='blue')

  
## Comparing the distributions. On the left column the real data, on the right column 
## the one simulation made here.

png(filename = "./sim/figures/1-manauara-sim-plots.png", res = 600,
    units = 'cm', h = 15, w = 15)
par(mfrow=c(2,2))
hist(males2$body, main = "", xlab = "Real body size data (mm)")
mtext("(A)", side = 3, adj = 0, cex = 1.2)
hist(xv2, main = "", xlab = "Simulated body size")
mtext("(B)", side = 3, adj = 0, cex = 1.2)


hist(males2$trait, main = "", xlab = "Real trait size data (mm)")
mtext("(C)", side = 3, adj = 0, cex = 1.2)
hist(yv2, main = "", xlab = "Simulated trait size")
mtext("(D)", side = 3, adj = 0, cex = 1.2)
dev.off()

#####################################################################
#### Now the mixture models approach with the simulated new dots ####
#####################################################################

hist(yv2, probability=T, xlab="Y (mm)", ylab="Relative frequency density", main=NULL)
SIMYs.den <- density(yv2)
lines(SIMYs.den, lwd=2)

Sim.SN.SkNorm2 <- smsn.mix(y = na.omit(yv2), g=2, nu=6, obs.prob=T, family="Skew.normal") # k = 2 + 2 + 2 + 1 = 7
Sim.SN.SkNorm1 <- smsn.mix(y = na.omit(yv2), g=1, nu=6, obs.prob=T, family="Skew.normal") # k = 1 + 1 + 1 = 3
Sim.AIC2dist <- myAICc(aic = Sim.SN.SkNorm2$aic, k = 7, n = Sim.SN.SkNorm2$n) # k = 2 + 2 + 2 + 1 = 7
Sim.AIC1dist <- myAICc(aic = Sim.SN.SkNorm1$aic, k = 3, n = Sim.SN.SkNorm1$n) # k = 1 + 1 + 1 = 3

Sim.AIC2dist
Sim.AIC1dist

dev.off()

### This is how we got the allometries and simulate data for all species.
### How we actually simulated a 1,000 trait sizes and calculated the AICcs is
### in another code that accompanies this one. It's just a bunch of for loops
### that feel more organized in another code.
### Below, we added the same code you saw in the first step but for the 
### other species. We are not adding a lot of comments below because it is
### similar to what you can see above.

#### FORSTEROPSALIS PUREORA #####

pur <- read.csv("./sim/data/powell-2020.csv",h=T)

# Real data mixture models
pur.SkNorm2 <- smsn.mix(y = na.omit(pur$cheli.length), g=2, nu=6, obs.prob=T, family="Skew.normal") # k = 2 + 2 + 2 + 1 = 7
pur.SkNorm1 <- smsn.mix(y = na.omit(pur$cheli.length), g=1, nu=6, obs.prob=T, family="Skew.normal") # k = 1 + 1 + 1 = 3

# Real data AICc
AIC2.pur = pur.SkNorm2$aic + ((2*(7^2) + 2*7) / (189 - 7 - 1))
AIC1.pur = pur.SkNorm1$aic + ((2*(3^2) + 2*3) / (189 - 3 - 1))

delta.aic.pur = AIC2.pur - AIC1.pur

# Non-linear allometries
ric.mod.pur = nlsLM(cheli.length ~ richards(pro.width, A, mu, lambda, addpar),
                       start=c(A = 14.89, mu = 2, lambda = 0.75, addpar = 0.1), 
                       data = pur, 
                       na.action=na.exclude)

logis.pur =  nls(cheli.length ~ SSlogis(pro.width, Asym, xmid, scal), data = pur)

# Check to see if any errors pop up
summary(ric.mod.pur)
summary(logis.pur)

AICctab(ric.mod.pur, logis.pur,
        weights = T, base = T)

# Checking how the non-linear allometry fits the data
xv <- seq(min(pur$pro.width), max(pur$pro.width), 0.001) 
yv <- predict(logis.pur, list(pro.width=xv))

scatterhist(x = pur$pro.width, y = pur$cheli.length,
            xsize = 2, cleanup =  F)
lines(xv,yv, col="red", lwd=2)
dev.off()

# Simulation of one body and trait data 

xv2 <- sort(rnorm(n = length(pur$pro.width),
                  mean = mean(pur$pro.width),
                  sd = sd(pur$pro.width)))

yv2 <- predict(logis.pur, list(pro.width=xv2))

plot(xv2,yv2, pch = 10, col='blue')

png(filename = "./sim/figures/1-pureora-sim-plots.png", res = 600,
    units = 'cm', h = 15, w = 15)
par(mfrow=c(2,2))
hist(pur$pro.width, main = "", xlab = "Real body size data (mm)")
mtext("(A)", side = 3, adj = 0, cex = 1.2)
hist(xv2, main = "", xlab = "Simulated body size")
mtext("(B)", side = 3, adj = 0, cex = 1.2)

hist(pur$cheli.length, main = "", xlab = "Real trait size data (mm)")
mtext("(C)", side = 3, adj = 0, cex = 1.2)
hist(yv2, main = "", xlab = "Simulated trait size")
mtext("(D)", side = 3, adj = 0, cex = 1.2)

dev.off()


#### COBANIA PICEA ####

cob <- read.csv("./sim/data/Cobania_picea.csv",h=T)

# Real data mixture models
cob.SkNorm2 <- smsn.mix(y = na.omit(cob$C4A), g=2, nu=6, obs.prob=T, family="Skew.normal") # k = 2 + 2 + 2 + 1 = 7
cob.SkNorm1 <- smsn.mix(y = na.omit(cob$C4A), g=1, nu=6, obs.prob=T, family="Skew.normal") # k = 1 + 1 + 1 = 3

# Real data AICc
AIC2.cob = cob.SkNorm2$aic + ((2*(7^2) + 2*7) / (60 - 7 - 1))
AIC1.cob = cob.SkNorm1$aic + ((2*(3^2) + 2*3) / (60 - 3 - 1))

delta.aic.cob = AIC2.cob - AIC1.cob

## We had some issues finding the best starting parameters for the Richard's equation
## That's why we used the trace = T argument. This argument allowed us to see where the 
## algorithm was having issues to find the best solution. We ran the model with a bunch
## of different starting values until we reached the values in the function below.
## The main issue to fit the Richard's growth equation is the asymptote at the large body
## size values, that's why we had to spend extra time trying to fit it well.

ric.mod.cob = nlsLM(C4A ~ richards(Dorsal.scute, A, mu, lambda, addpar),
                  start=c(A = 8.785978 , mu = 3.777412 , lambda = 6.650161 , addpar = 4.536240 ), 
                  data = cob, control = nls.control(maxiter = 500), trace = T,
                  na.action=na.exclude)

logis.cob =  nls(C4A ~ SSlogis(Dorsal.scute, Asym, xmid, scal), data = cob)

# Checking if any errors creep up
summary(ric.mod.cob)
summary(logis.cob)

AICctab(ric.mod.cob, logis.cob,
        weights = T, base = T)

# Viz
xv <- seq(min(cob$Dorsal.scute), max(cob$Dorsal.scute), 0.001) 
yv <- predict(logis.cob, list(Dorsal.scute=xv))

scatterhist(x = cob$Dorsal.scute,
            y = cob$C4A, 
            xsize = 2, cleanup = F)
lines(xv,yv, col="red", lwd=2)

dev.off()

# Simulating data 

xv2 <- sort(rnorm(n = length(cob$Dorsal.scute),
                  mean = mean(cob$Dorsal.scute),
                  sd = sd(cob$Dorsal.scute)))

yv2 <- predict(logis.cob, list(Dorsal.scute=xv2))

plot(xv2,yv2, pch = 10, col='blue')

png(filename = "./sim/figures/1-cobania-sim-plots.png", res = 600,
    units = 'cm', h = 15, w = 15)
par(mfrow=c(2,2))
hist(cob$Dorsal.scute, main = "", xlab = "Real body size data (mm)")
hist(xv2, main = "", xlab = "Simulated body size")

hist(cob$C4A, main = "", xlab = "Real trait size data (mm)")
hist(yv2, main = "", xlab = "Simulated trait size")

dev.off()

#### NEOANCISTROTUS GRACILIS ####

neo <- read.csv("./sim/data/Neoancistrotus_gracilis.csv",h=T)

# Real data mixture models
neo.SkNorm2 <- smsn.mix(y = na.omit(neo$C4A), g=2, nu=6, obs.prob=T, family="Skew.normal") # k = 2 + 2 + 2 + 1 = 7
neo.SkNorm1 <- smsn.mix(y = na.omit(neo$C4A), g=1, nu=6, obs.prob=T, family="Skew.normal") # k = 1 + 1 + 1 = 3

# Real data AICc
AIC2.neo = neo.SkNorm2$aic + ((2*(7^2) + 2*7) / (64 - 7 - 1))
AIC1.neo = neo.SkNorm1$aic + ((2*(3^2) + 2*3) / (64 - 3 - 1))

delta.aic.neo = AIC2.neo - AIC1.neo

# Real data non-linear allometries
ric.mod.neo = nlsLM(C4A ~ richards(Dorsal.scute, A, mu, lambda, addpar),
                     start=c(A = 4.1, mu = 2, lambda = 3, addpar = 1), 
                     data = neo,  control = nls.control(maxiter = 1000))

logis.mod.neo =  nls(C4A ~ SSlogis(Dorsal.scute, Asym, xmid, scal), data = neo)

summary(ric.mod.neo)
summary(logis.mod.neo)

AICctab(logis.mod.neo, ric.mod.neo,
        weights = T, base = T)

## This is a tricky one. Even though the Richard's equation fits better,
## it has a weird asymptote at the end of the distribution. This is isn't 
## a problem for the real data, but it is for the simulated data. 
## So, we decided to use the logistic non-linear model just because it would
## allow us to simulate the data without hiccups. 

# Viz
xv <- seq(min(neo$Dorsal.scute), max(neo$Dorsal.scute), 0.001) 
yv <- predict(logis.mod.neo, list(Dorsal.scute=xv))

scatterhist(x = neo$Dorsal.scute,
            y = neo$C4A,
            xsize = 2, cleanup = F)
lines(xv,yv, col="red", lwd=2)

dev.off()

# Simulating one round

xv2 <- sort(rnorm(n = length(neo$Dorsal.scute),
                  mean = mean(neo$Dorsal.scute),
                  sd = sd(neo$Dorsal.scute)))

yv2 <- predict(logis.mod.neo, list(Dorsal.scute=xv2))

plot(xv2,yv2, pch = 10, col='blue')

png(filename = "./sim/figures/1-gracilis-sim-plots.png", res = 600,
    units = 'cm', h = 15, w = 15)
par(mfrow=c(2,2))
hist(neo$Dorsal.scute, main = "", xlab = "Real body size data (mm)")
mtext("(A)", side = 3, adj = 0, cex = 1.2)
hist(xv2, main = "", xlab = "Simulated body size")
mtext("(B)", side = 3, adj = 0, cex = 1.2)

hist(neo$C4A, main = "", xlab = "Real trait size data (mm)")
mtext("(C)", side = 3, adj = 0, cex = 1.2)
hist(yv2, main = "", xlab = "Simulated trait size")
mtext("(D)", side = 3, adj = 0, cex = 1.2)

dev.off()


#### PROMITOBATES BELLUS ####

promit <- read.csv("./sim/data/Promitobates_bellus.csv",h=T)

# Real data mixture models
promit.SkNorm2 <- smsn.mix(y = na.omit(promit$C4A), g=2, nu=6, obs.prob=T, family="Skew.normal") # k = 2 + 2 + 2 + 1 = 7
promit.SkNorm1 <- smsn.mix(y = na.omit(promit$C4A), g=1, nu=6, obs.prob=T, family="Skew.normal") # k = 1 + 1 + 1 = 3

# Real data AICc
AIC2.promit = promit.SkNorm2$aic + ((2*(7^2) + 2*7) / (54 - 7 - 1))
AIC1.promit = promit.SkNorm1$aic + ((2*(3^2) + 2*3) / (54 - 3 - 1))

delta.aic.promit = AIC2.promit - AIC1.promit

# Non-linear allometries
ric.mod.promit = nlsLM(C4A ~ richards(Dorsal.scute, A, mu, lambda, addpar),
                  start=c(A = 2.155499, mu = 1.850555, lambda = 3.207198, addpar = 17.925252), 
                  data = promit)

logis.mod.promit =  nls(C4A ~ SSlogis(Dorsal.scute, Asym, xmid, scal), data = promit)

summary(ric.mod.promit)
summary(logis.mod.promit)

AICctab(ric.mod.promit, logis.mod.promit, 
         weights = T, base = T)

# Viz
xv <- seq(min(promit$Dorsal.scute), max(promit$Dorsal.scute), 0.001) 
yv <- predict(ric.mod.promit, list(Dorsal.scute=xv))

scatterhist(y = promit$C4A,
            x = promit$Dorsal.scute,
            xsize = 2, cleanup = F)
lines(xv,yv, col="red", lwd=2)

dev.off()

# Simulating

xv2 <- sort(rnorm(n = length(promit$Dorsal.scute),
                  mean = mean(promit$Dorsal.scute),
                  sd = sd(promit$Dorsal.scute)))

yv2 <- predict(ric.mod.promit, list(Dorsal.scute=xv2))

plot(xv2,yv2, pch = 10, col='blue')

png(filename = "./sim/figures/1-promitobates-sim-plots.png", res = 600,
    units = 'cm', h = 15, w = 15)
par(mfrow=c(2,2))
hist(promit$Dorsal.scute, main = "", xlab = "Real body size data (mm)")
mtext("(A)", side = 3, adj = 0, cex = 1.2)
hist(xv2, main = "", xlab = "Simulated body size")
mtext("(B)", side = 3, adj = 0, cex = 1.2)

hist(promit$C4A, main = "", xlab = "Real trait size data (mm)")
mtext("(C)", side = 3, adj = 0, cex = 1.2)
hist(yv2, main = "", xlab = "Simulated trait size")
mtext("(D)", side = 3, adj = 0, cex = 1.2)

dev.off()

### FIGURE ####

png(filename = './sim/figures/2-manauara-allometry-plots.png', res = 600,
    units = 'mm', width = 140, height = 140)
#par(mfrow = c(3,2), las = 1, bty = 'l')

xv <- seq(min(males2$body), max(males2$body), 0.001) 
yv <- predict(logis.manauara, list(body=xv))
scatterhist( y = males2$trait,
             x = males2$body,
             xsize = 2, cleanup = F,
             ylab = "Tusk length (cm)",
             xlab = "Dorsal scute length (cm)")
lines(xv,yv, col=alpha("red",0.5), lwd=3)

dev.off()

png(filename = './sim/figures/2-pureora-allometry-plots.png', res = 600,
    units = 'mm', width = 140, height = 140) 

xv <- seq(min(pur$pro.width), max(pur$pro.width), 0.001) 
yv <- predict(logis.pur, list(pro.width=xv))

scatterhist( y = pur$cheli.length,
             x = pur$pro.width,
             xsize = 2, cleanup = F,
             ylab = "Chelicerae length (cm)",
             xlab = "Prosoma width (cm)")
lines(xv,yv, col=alpha("red",0.5), lwd=5)
dev.off()

png(filename = './sim/figures/2-cobania-allometry-plots.png', res = 600,
    units = 'mm', width = 140, height = 140) 

xv <- seq(min(cob$Dorsal.scute), max(cob$Dorsal.scute), 0.001) 
yv <- predict(logis.cob, list(Dorsal.scute=xv))
scatterhist( y = cob$C4A,
             x = cob$Dorsal.scute,
             xsize = 2, cleanup = F,
             ylab = "Coxa IV Apophysis length (cm)",
             xlab = "Dorsal scute length (cm)")
lines(xv,yv, col=alpha("red",0.5), lwd=5)
dev.off()

png(filename = './sim/figures/2-neo-allometry-plots.png', res = 600,
    units = 'mm', width = 140, height = 140) 
xv <- seq(min(neo$Dorsal.scute), max(neo$Dorsal.scute), 0.001) 
yv <- predict(logis.mod.neo, list(Dorsal.scute=xv))

scatterhist( y = neo$C4A,
             x = neo$Dorsal.scute,
             xsize = 2, cleanup = F,
             ylab = "Coxa IV Apophysis length (cm)",
             xlab = "Dorsal scute length (cm)")
lines(xv,yv, col=alpha("red",0.5), lwd=5)
dev.off()

png(filename = './sim/figures/2-promit-allometry-plots.png', res = 600,
    units = 'mm', width = 140, height = 140) 
xv <- seq(min(promit$Dorsal.scute), max(promit$Dorsal.scute), 0.001) 
yv <- predict(ric.mod.promit, list(Dorsal.scute=xv))

scatterhist( y = promit$C4A,
             x = promit$Dorsal.scute,
             xsize = 2, cleanup = F,
             ylab = "Coxa IV Apophysis length (cm)",
             xlab = "Dorsal scute length (cm)")
lines(xv,yv, col=alpha("red",0.5), lwd=5)
dev.off()
