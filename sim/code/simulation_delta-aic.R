library(parallel); library(doParallel); library(tidyverse)

# This next line will run a previous code entirely, so it might take a few 
# minutes depending on your computer.
source("./sim/code/Simulating sigmoid allometries and detecting dimorphism.R")

## This is the code to simulate a 1,000 delta AICcs.
## You will note that we parallelized the functions to run faster.
## Since we did this in Windows-based systems, the code to parallelize 
## is a bit different from UNIX-based systems. However, we are also 
## providing the data with the simulations already done. That way,
## you don't need to run it yourself if you don't want to. 
## The code is far from optimized, so it will take a while to run each
## species.


#### PHAREICRANAUS MANAUARA ####

aic.sim.manau = function(x) {
  
  require(mixsmsn)
  require(tidyverse)
  
  sim.male.body = data.frame(matrix(NA ,
                                    nrow = length(males2$body),
                                    ncol = 1000))
  sim.male.weap = data.frame(matrix(NA ,
                                    nrow = length(males2$body),
                                    ncol = 1000))
  
  sims2 = vector( mode = 'list', length = 1000)
  sims1 = vector( mode = 'list', length = 1000)
  
  sim.aicc =  data.frame(AIC2 = 0,
                         AIC1 = 0)
  
  delta.sim.aicc = data.frame(matrix(NA ,
                                     nrow = length(males2$body),
                                     ncol = 1000))

  for(i in 1:1000) {
    sim.male.body[i] = sort(rnorm(n = length(males2$body),
                                  mean = mean(males2$body),
                                  sd = sd(males2$body)))
    
    sim.male.weap[i] = predict(logis.manauara, list(body=sim.male.body[,i]))
  }
  
  
  
  for(i in 1:1000) {
    sims2[[i]] = smsn.mix(y = na.omit(sim.male.weap[,i]), g=2, nu=6, obs.prob=T, family="Skew.normal")
    sim.aicc[i,1] = sims2[[i]]$aic + ((2*(7^2) + 2*7) / (125 - 7 - 1))
    
    sims1[[i]] = smsn.mix(y = na.omit(sim.male.weap[,i]), g=1, nu=6, obs.prob=T, family="Skew.normal")
    sim.aicc[i,2] = sims1[[i]]$aic + ((2*(3^2) + 2*3) / (125 - 3 - 1))
  }
  
  
  delta.sim.aicc = sim.aicc[,1] - sim.aicc[,2]
  
  return(sim.aicc)
  return(delta.sim.aicc)
}

# We chose 10 cores based on our computers.
# Choose whatever core value you feel better about.

no_cores = 10 
cl = makeCluster(no_cores)

# THE COMMAND LINE BELOW WILL MAKE THE FUNCTION RUN!!!
# If you want to run, uncomment it.

#system.time( {
#  registerDoParallel(cl)
#  results = foreach(i = 1:10, .combine = 'c') %dopar% aic.sim.manau()
#} )

#

#saveRDS(results,file = './sim/simul_data/manauara_resampling.RDS')

results = readRDS('./sim/simul_data/manauara_resampling.RDS')

# Checking if the list generated is okay.
str(results)

summary(results)

# Getting only the AIC for the mixture model with two distributions
results.AIC2 = results[c(1,3,5,7,9,11,13,15,17,19)]
results.AIC2 = unlist(results.AIC2)

hist(results.AIC2)
length(results.AIC2)

# Now the AIC values for the mixture model with only one distribution
results.AIC1 = results[c(2,4,6,8,10,12,14,16,18,20)]
results.AIC1 = unlist(results.AIC1)

hist(results.AIC1)
length(results.AIC1)

# Generating the AICcs. They're ordered, so we don't need to worry about 
# matching them
delta.sim.aicc = results.AIC2 - results.AIC1

# Getting a sense of median values and whatnot
summary(results.AIC1)
summary(results.AIC2)
summary(delta.sim.aicc)

# Calculating the delta aic from the real dataset
delta.aic.real = AIC2dist - AIC1dist

# First, we will plot the histogram with the generated delta AICcs.
# Then, we will plot where the real AICc lands.
hist(delta.sim.aicc, prob = T, 
     xlab = "Delta AICc", main = "Phareicranaus manauara") 
abline(v = mean(delta.aic.real), lwd = 5)

# Proportion of simulated delta AICcs that are smaller than the 
# real AICc. We used "smaller" instead of "larger" because they're
# mostly negative values. So, we need to control for that.

sum(delta.sim.aicc < delta.aic.real)/length(delta.sim.aicc)
# only 1.4% are smaller. Pretty rare to generate it

# Just getting how many positive values we generated using the simulations
sum(delta.sim.aicc > 0)
sum(delta.sim.aicc > 0)/length(delta.sim.aicc)

### From here on, the idea is identical. We will simulate 10,000 
### AICcs, then calculate the delta AICcs and compare to the observed
### AICc from the dataset. Comments will be kept to a mininum because
### the codes are nearly identical.

#### FORSTEROPSALIS PUREORA ####

aic.sim.pur = function(x) {
  
  require(mixsmsn)
  require(tidyverse)
  
  sim.male.body = data.frame(matrix(NA ,
                                    nrow = length(pur$pro.width),
                                    ncol = 1000))
  sim.male.weap = data.frame(matrix(NA ,
                                    nrow = length(pur$pro.width),
                                    ncol = 1000))
  
  sims2 = vector( mode = 'list', length = 1000)
  sims1 = vector( mode = 'list', length = 1000)
  
  sim.aicc =  data.frame(AIC2 = 0,
                         AIC1 = 0)
  
  delta.sim.aicc = data.frame(matrix(NA ,
                                     nrow = length(pur$pro.width),
                                     ncol = 1000))

  
  for(i in 1:1000) {
    sim.male.body[i] = sort(rnorm(n = length(pur$pro.width),
                                  mean = mean(pur$pro.width),
                                  sd = sd(pur$pro.width)))

    sim.male.weap[i] = predict(logis.pur, list(pro.width=sim.male.body[,i]))
  }
  
  
  
  for(i in 1:1000) {
    sims2[[i]] = smsn.mix(y = na.omit(sim.male.weap[,i]), g=2, nu=6, obs.prob=T, family="Skew.normal")
    sim.aicc[i,1] = sims2[[i]]$aic + ((2*(7^2) + 2*7) / (189 - 7 - 1))

    sims1[[i]] = smsn.mix(y = na.omit(sim.male.weap[,i]), g=1, nu=6, obs.prob=T, family="Skew.normal")
    sim.aicc[i,2] = sims1[[i]]$aic + ((2*(3^2) + 2*3) / (189 - 3 - 1))
  }
  
  return(sim.aicc)
  
  delta.sim.aicc = sim.aicc[,1] - sim.aicc[,2]
  
  return(delta.sim.aicc)
}

no_cores = 10
cl = makeCluster(no_cores)

# THE COMMAND LINE BELOW WILL MAKE THE FUNCTION RUN!!!
# If you want to run, uncomment it.

#system.time( {
#  registerDoParallel(cl)
#  results.pur = foreach(i = 1:10, .combine = 'c') %dopar% aic.sim.pur()
#} )

#saveRDS(results.pur,file = './sim/simul_data/pur_resampling.RDS')

results.pur = readRDS(file = './sim/simul_data/pur_resampling.RDS')

summary(results.pur)

results.AIC2 = results.pur[c(1,3,5,7,9,11,13,15,17,19)]
results.AIC2 = unlist(results.AIC2)

hist(results.AIC2)
length(results.AIC2)

results.AIC1 = results.pur[c(2,4,6,8,10,12,14,16,18,20)]
results.AIC1 = unlist(results.AIC1)

hist(results.AIC1)
length(results.AIC1)

delta.sim.aicc = results.AIC2 - results.AIC1

summary(results.AIC1)
summary(results.AIC2)
summary(delta.sim.aicc)

delta.aic.pur = AIC2.pur - AIC1.pur

hist(delta.sim.aicc, prob = T, 
     xlab = "Delta AICc", main = "Forsteropsalis pureora")
abline(v = mean(delta.aic.pur), lwd = 5)

sum(delta.sim.aicc < delta.aic.pur)/length(delta.sim.aicc)
sum(delta.sim.aicc > 0)
sum(delta.sim.aicc > 0)/length(delta.sim.aicc)


#### COBANIA PICEA ####

aic.sim.cob = function(x) {

  require(grofit)
  require(mixsmsn)
  require(tidyverse)
  
  sim.male.body = data.frame(matrix(NA ,
                                    nrow = length(cob$Dorsal.scute),
                                    ncol = 1000))
  sim.male.weap = data.frame(matrix(NA ,
                                    nrow = length(cob$Dorsal.scute),
                                    ncol = 1000))
  
  sims2 = vector( mode = 'list', length = 1000)
  sims1 = vector( mode = 'list', length = 1000)
  
  sim.aicc =  data.frame(AIC2 = 0,
                         AIC1 = 0)
  
  delta.sim.aicc = data.frame(matrix(NA ,
                                     nrow = length(cob$Dorsal.scute),
                                     ncol = 1000))
  
  
  for(i in 1:1000) {
    sim.male.body[i] = sort(rnorm(n = length(cob$Dorsal.scute),
                                  mean = mean(cob$Dorsal.scute),
                                  sd = sd(cob$Dorsal.scute)))

    sim.male.weap[i] = predict(logis.cob, list(Dorsal.scute=sim.male.body[,i]))
  }
  
  
  
  for(i in 1:1000) {
    sims2[[i]] = smsn.mix(y = na.omit(sim.male.weap[,i]), g=2, nu=6, obs.prob=T, family="Skew.normal")
    sim.aicc[i,1] = sims2[[i]]$aic + ((2*(7^2) + 2*7) / (60 - 7 - 1))
 
    sims1[[i]] = smsn.mix(y = na.omit(sim.male.weap[,i]), g=1, nu=6, obs.prob=T, family="Skew.normal")
    sim.aicc[i,2] = sims1[[i]]$aic + ((2*(3^2) + 2*3) / (60 - 3 - 1))
  }
  
  return(sim.aicc)
  
}
  
no_cores = 10
cl = makeCluster(no_cores)

# THE COMMAND LINE BELOW WILL MAKE THE FUNCTION RUN!!!
# If you want to run, uncomment it.

#system.time( {
#  registerDoParallel(cl)
#  results.cob = foreach(i = 1:10, .combine = 'c') %dopar% aic.sim.cob()
#} )

#saveRDS(results.cob,file = './sim/simul_data/cobania_resampling.RDS')

results.cob = readRDS(file = './sim/simul_data/cobania_resampling.RDS')

summary(results.cob)

results.AIC2 = results.cob[c(1,3,5,7,9,11,13,15,17,19)]
results.AIC2 = unlist(results.AIC2)

hist(results.AIC2)
length(results.AIC2)

results.AIC1 = results.cob[c(2,4,6,8,10,12,14,16,18,20)]
results.AIC1 = unlist(results.AIC1)

hist(results.AIC1)
length(results.AIC1)

delta.sim.aicc = results.AIC2 - results.AIC1

summary(results.AIC1)
summary(results.AIC2)
summary(delta.sim.aicc)

hist(delta.sim.aicc, prob = T, 
     xlab = "Delta AICc", main = "Cobania picea")
abline(v = mean(delta.aic.cob), lwd = 5)

sum(delta.sim.aicc < delta.aic.cob)/length(delta.sim.aicc)
sum(delta.sim.aicc > 0)
sum(delta.sim.aicc > 0)/length(delta.sim.aicc)


#### NEOANCISTROTUS GRACILIS ####

aic.sim.neo = function(x) {
  
  require(mixsmsn)
  require(tidyverse)
  
  sim.male.body = data.frame(matrix(NA ,
                                    nrow = length(neo$Dorsal.scute),
                                    ncol = 1000))
  sim.male.weap = data.frame(matrix(NA ,
                                    nrow = length(neo$Dorsal.scute),
                                    ncol = 1000))
  
  sims2 = vector( mode = 'list', length = 1000)
  sims1 = vector( mode = 'list', length = 1000)
  
  sim.aicc =  data.frame(AIC2 = 0,
                         AIC1 = 0)
  
  delta.sim.aicc = data.frame(matrix(NA ,
                                     nrow = length(neo$Dorsal.scute),
                                     ncol = 1000))
  

  for(i in 1:1000) {
    sim.male.body[i] = sort(rnorm(n = length(neo$Dorsal.scute),
                                  mean = mean(neo$Dorsal.scute),
                                  sd = sd(neo$Dorsal.scute)))

    sim.male.weap[i] = predict(logis.mod.neo, list(Dorsal.scute=sim.male.body[,i]))
  }
  
  
  for(i in 1:1000) {
    sims2[[i]] = smsn.mix(y = na.omit(sim.male.weap[,i]), g=2, nu=6, obs.prob=T, family="Skew.normal")
    sim.aicc[i,1] = sims2[[i]]$aic + ((2*(7^2) + 2*7) / (189 - 7 - 1))

    sims1[[i]] = smsn.mix(y = na.omit(sim.male.weap[,i]), g=1, nu=6, obs.prob=T, family="Skew.normal")
    sim.aicc[i,2] = sims1[[i]]$aic + ((2*(3^2) + 2*3) / (189 - 3 - 1))
  }
  
  return(sim.aicc)
}

no_cores = 10
cl = makeCluster(no_cores)

# THE COMMAND LINE BELOW WILL MAKE THE FUNCTION RUN!!!
# If you want to run, uncomment it.

#system.time( {
#  registerDoParallel(cl)
#  results.neo = foreach(i = 1:10, .combine = 'c') %dopar% aic.sim.neo()
#} )


#saveRDS(results.neo,file = './sim/simul_data/neo_resampling.RDS')

results.neo = readRDS(file = './sim/simul_data/neo_resampling.RDS')

summary(results.neo)

results.AIC2 = results.neo[c(1,3,5,7,9,11,13,15,17,19)]
results.AIC2 = unlist(results.AIC2)

hist(results.AIC2)
length(results.AIC2)

results.AIC1 = results.neo[c(2,4,6,8,10,12,14,16,18,20)]
results.AIC1 = unlist(results.AIC1)

hist(results.AIC1)
length(results.AIC1)

delta.sim.aicc = results.AIC2 - results.AIC1

delta.aic.neo = AIC2.neo - AIC1.neo

summary(results.AIC1)
summary(results.AIC2)
summary(delta.sim.aicc)

hist(delta.sim.aicc, prob = T, 
     xlab = "Delta AICc", main = "Neoancistrotus gracilis")
abline(v = mean(delta.aic.neo), lwd = 5)

sum(delta.sim.aicc < delta.aic.neo)/length(delta.sim.aicc)
sum(delta.sim.aicc > 0)
sum(delta.sim.aicc > 0)/length(delta.sim.aicc)

#### PROMITOBATES BELLUS #####

aic.sim.pro = function(x) {
  
  require(mixsmsn)
  require(tidyverse)
  
  sim.male.body = data.frame(matrix(NA ,
                                    nrow = length(promit$Dorsal.scute),
                                    ncol = 1000))
  sim.male.weap = data.frame(matrix(NA ,
                                    nrow = length(promit$Dorsal.scute),
                                    ncol = 1000))
  
  sims2 = vector( mode = 'list', length = 1000)
  sims1 = vector( mode = 'list', length = 1000)
  
  sim.aicc =  data.frame(AIC2 = 0,
                         AIC1 = 0)
  
  delta.sim.aicc = data.frame(matrix(NA ,
                                     nrow = length(promit$Dorsal.scute),
                                     ncol = 1000))
  
  for(i in 1:1000) {
    sim.male.body[i] = sort(rnorm(n = length(promit$Dorsal.scute),
                                  mean = mean(promit$Dorsal.scute),
                                  sd = sd(promit$Dorsal.scute)))

    sim.male.weap[i] = predict(ric.mod.promit, list(Dorsal.scute=sim.male.body[,i]))
  }
  
  
  
  for(i in 1:1000) {
    sims2[[i]] = smsn.mix(y = na.omit(sim.male.weap[,i]), g=2, nu=6, obs.prob=T, family="Skew.normal")
    sim.aicc[i,1] = sims2[[i]]$aic + ((2*(7^2) + 2*7) / (189 - 7 - 1))

    sims1[[i]] = smsn.mix(y = na.omit(sim.male.weap[,i]), g=1, nu=6, obs.prob=T, family="Skew.normal")
    sim.aicc[i,2] = sims1[[i]]$aic + ((2*(3^2) + 2*3) / (189 - 3 - 1))
  }
  
  return(sim.aicc)
 
}

no_cores = 10
cl = makeCluster(no_cores)

# THE COMMAND LINE BELOW WILL MAKE THE FUNCTION RUN!!!
# If you want to run, uncomment it.

#system.time( {
#  registerDoParallel(cl)
#  results.pro = foreach(i = 1:10, .combine = 'c') %dopar% aic.sim.pro()
#} )

#saveRDS(results.pro,file = './sim/simul_data/promit_resampling.RDS')

results.pro = readRDS(file = './sim/simul_data/promit_resampling.RDS')

results.AIC2 = results.pro[c(1,3,5,7,9,11,13,15,17,19)]
results.AIC2 = unlist(results.AIC2)

hist(results.AIC2)
length(results.AIC2)

results.AIC1 = results.pro[c(2,4,6,8,10,12,14,16,18,20)]
results.AIC1 = unlist(results.AIC1)

hist(results.AIC1)
length(results.AIC1)

delta.sim.aicc = results.AIC2 - results.AIC1

delta.aic.pro = AIC2.promit - AIC1.promit

summary(results.AIC1)
summary(results.AIC2)
summary(delta.sim.aicc)

sum(delta.sim.aicc < delta.aic.pro)/length(delta.sim.aicc)
sum(delta.sim.aicc > 0)
sum(delta.sim.aicc > 0)/length(delta.sim.aicc)

hist(delta.sim.aicc, prob = T, 
     xlab = "Delta AICc", main = "Neoancistrotus gracilis")
abline(v = mean(delta.aic.pro), lwd = 5)

sum(delta.sim.aicc < delta.aic.pro)/length(delta.sim.aicc)
sum(delta.sim.aicc > 0)
sum(delta.sim.aicc > 0)/length(delta.sim.aicc)
