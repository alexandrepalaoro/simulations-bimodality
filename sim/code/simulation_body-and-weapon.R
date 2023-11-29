#### This code was built to generate the simulations graphs in the Supplementary Material
#### It is largely a part of the AICc simulation script, but I cut down several unnecessary parts
#### Thus, it runs much faster.

## This is the main function. It will generate the two data frames - body size simulation
## and weapon size simulation, then merge both. The rest of the code plots each simulation
## as a density function.

# This next line will run a previous code entirely, so it might take a few 
# minutes depending on your computer.
source("Simulating sigmoid allometries and detecting dimorphism.R")

#### PHAREICRANAUS MANAUARA ####

weap.body.sim.manau = function(x) {
  

  sim.male.body = data.frame(matrix(NA ,
                                    nrow = length(males2$body),
                                    ncol = 1000))
  sim.male.weap = data.frame(matrix(NA ,
                                    nrow = length(males2$body),
                                    ncol = 1000))
  
    for(i in 1:1000) {
    sim.male.body[i] = sort(rnorm(n = length(males2$body),
                                  mean = mean(males2$body),
                                  sd = sd(males2$body)))
    sim.male.weap[i] = predict(logis.manauara, list(body=sim.male.body[,i]))
  
    }
  
  male.traits = cbind(body = sim.male.body, weap = sim.male.weap)
  return(male.traits)
}

### This code runs the function

manaura = weap.body.sim.manau()


#tiff(filename = "./figures/manaura-sim.tiff", res = 600, compression = 'lzw',
#     units = 'mm', width = 160, height = 120)

plot(density(manaura[,1]), las = 1, bty = 'l', 
     xlab = "Simulated data (mm)",
     xlim = c(-0.2, 1.3),
     ylim = c(0, 8),
     main = "",
     col = 'grey')

### From column #1 to #1000 we have body size simulations.
### From #1001 to #2001 we have weapon size simulations.
### The for() is simply to plot all columns automatically without 
### the need for a thousand lines. Takes a while to plot though.

## We felt that plotting all 1,000 simulations was making the graph too busy.
## Thus, we changed the code to plot only 100 random samples. 

# This samples 100 random simulations from the entire simulation pool.
sim.to.plot.Body <- sample(c(1:1000), size = 100, replace = F)
sim.to.plot.Weapon <- sample(c(1001:2001), size = 100, replace = F)


# Plotting 100 body size simulated data

for(i in sim.to.plot.Body) {
lines(density(manaura[,i]),
      col = 'grey')
}

lines(density(males2$body), col = 'black', lwd = 2)

# Plotting weapon size simulated data
for(i in sim.to.plot.Weapon) {
  lines(density(manaura[,i]),
        col = 'grey')
}

lines(density(males2$trait), col = 'black', lwd = 2)

#dev.off()


### From here on I won't comment on the code because it is the same thing we did
### for Phareicranaus manaurar. We are only changing the species, but the simulations
### are the same. I'm leaving the species indexed though, in case we need to find them
### faster.


#### FORSTEROPSALIS PUREORA ####

weap.body.sim.pur = function(x) {
  
  
  sim.male.body = data.frame(matrix(NA ,
                                    nrow = length(pur$pro.width),
                                    ncol = 1000))
  sim.male.weap = data.frame(matrix(NA ,
                                    nrow = length(pur$pro.width),
                                    ncol = 1000))
  
  for(i in 1:1000) {
    sim.male.body[i] = sort(rnorm(n = length(pur$pro.width),
                                  mean = mean(pur$pro.width),
                                  sd = sd(pur$pro.width)))
  }
  
  for(i in 1:1000) {
    sim.male.weap[i] = predict(logis.pur, list(pro.width=sim.male.body[,i]))
  }
  
  male.traits = cbind(body = sim.male.body, weap = sim.male.weap)
  return(male.traits)
}

pureora = weap.body.sim.pur()


#tiff(filename = "./figures/pureora-sim.tiff", res = 600, compression = 'lzw',
#     units = 'mm', width = 160, height = 120)

sim.to.plot.Body <- sample(c(1:1000), size = 100, replace = F)
sim.to.plot.Weapon <- sample(c(1001:2001), size = 100, replace = F)

plot(density(pureora[,1]), las = 1, bty = 'l', 
     xlab = "Simulated data (mm)",
     main = "",
     xlim = c(0,16),
     ylim = c(0,1.2),
     col = 'grey')

for(i in sim.to.plot.Body) {
  lines(density(pureora[,i]),
        col = 'grey')
}

lines(density(pur$pro.width),
      col = 'black', lwd = 2)

for(i in sim.to.plot.Weapon) {
  lines(density(pureora[,i]), col = "gray")
}

lines(density(pur$cheli.length),
      lwd = 2)

#dev.off()


#### COBANIA PICEA ####

weap.body.sim.cob = function(x) {
  
  
  sim.male.body = data.frame(matrix(NA ,
                                    nrow = length(cob$Dorsal.scute),
                                    ncol = 1000))
  sim.male.weap = data.frame(matrix(NA ,
                                    nrow = length(cob$Dorsal.scute),
                                    ncol = 1000))
  
  for(i in 1:1000) {
    sim.male.body[i] = sort(rnorm(n = length(cob$Dorsal.scute),
                                  mean = mean(cob$Dorsal.scute),
                                  sd = sd(cob$Dorsal.scute)))
  }
  
  for(i in 1:1000) {
    sim.male.weap[i] = predict(logis.cob, list(Dorsal.scute=sim.male.body[,i]))
  }
  
  
  male.traits = cbind(body = sim.male.body, weap = sim.male.weap)
  return(male.traits)
}


cobania = weap.body.sim.cob()

#tiff(filename = "./figures/cobania-sim.tiff", res = 600, compression = 'lzw',
#     units = 'mm', width = 160, height = 120)
plot(density(cobania[,1]), las = 1, bty = 'l', 
     xlab = "simulated data (mm)",
     main = "",
     xlim = c(-1,14),
     ylim = c(0, 1),
     col = 'grey')

sim.to.plot.Body <- sample(c(1:1000), size = 100, replace = F)
sim.to.plot.Weapon <- sample(c(1001:2001), size = 100, replace = F)

for(i in sim.to.plot.Body) {
  lines(density(cobania[,i]),
        col = 'grey')
}

lines(density(cob$Dorsal.scute), col = 'black', lwd = 2)

for(i in sim.to.plot.Weapon) {
  lines(density(cobania[,i]), col = "gray")
}

lines(density(cob$C4A), lwd = 2)

#dev.off()


#### NEOANCISTROTUS GRACILIS ####

weap.body.sim.neo = function(x) {
  
  
  sim.male.body = data.frame(matrix(NA ,
                                    nrow = length(neo$Dorsal.scute),
                                    ncol = 1000))
  sim.male.weap = data.frame(matrix(NA ,
                                    nrow = length(neo$Dorsal.scute),
                                    ncol = 1000))
  
  for(i in 1:1000) {
    sim.male.body[i] = sort(rnorm(n = length(neo$Dorsal.scute),
                                  mean = mean(neo$Dorsal.scute),
                                  sd = sd(neo$Dorsal.scute)))
  }
  
  for(i in 1:1000) {
    sim.male.weap[i] = predict(logis.mod.neo, list(Dorsal.scute=sim.male.body[,i]))
  }
  
  
  male.traits = cbind(body = sim.male.body, weap = sim.male.weap)
  return(male.traits)
}


neo.sim = weap.body.sim.neo()

#tiff(filename = "./figures/neo-sim.tiff", res = 600, compression = 'lzw',
#     units = 'mm', width = 160, height = 120)
plot(density(neo.sim[,1]), las = 1, bty = 'l', 
     xlab = "simulated data (mm)",
     main = "",
     xlim = c(0, 5),
     ylim = c(0, 3),
     col = 'grey')

sim.to.plot.Body <- sample(c(1:1000), size = 100, replace = F)
sim.to.plot.Weapon <- sample(c(1001:2001), size = 100, replace = F)

for(i in sim.to.plot.Body) {
  lines(density(neo.sim[,i]),
        col = 'grey')
}

lines(density(neo$Dorsal.scute), col = 'black', lwd = 2)

for(i in sim.to.plot.Weapon) {
  lines(density(neo.sim[,i]),
        col = 'grey',)
}

lines(density(neo$C4A), lwd = 2)

#dev.off()


#### PROMITOBATES BELLUS #####

weap.body.sim.promit = function(x) {
  
  
  sim.male.body = data.frame(matrix(NA ,
                                    nrow = length(promit$Dorsal.scute),
                                    ncol = 1000))
  sim.male.weap = data.frame(matrix(NA ,
                                    nrow = length(promit$Dorsal.scute),
                                    ncol = 1000))
  
  for(i in 1:1000) {
    sim.male.body[i] = sort(rnorm(n = length(promit$Dorsal.scute),
                                  mean = mean(promit$Dorsal.scute),
                                  sd = sd(promit$Dorsal.scute)))
  }
  
  for(i in 1:1000) {
    sim.male.weap[i] = predict(ric.mod.promit, list(Dorsal.scute=sim.male.body[,i]))
  }
  
  
  male.traits = cbind(body = sim.male.body, weap = sim.male.weap)
  return(male.traits)
}


promit.sim = weap.body.sim.promit()

#tiff(filename = "./figures/promit-sim.tiff", res = 600, compression = 'lzw',
#     units = 'mm', width = 160, height = 120)

plot(density(promit.sim[,1]), las = 1, bty = 'l', 
     xlab = "simulated data (mm)",
     main = "",
     xlim = c(0, 6),
     ylim = c(0,3.2),
     col = 'grey')

sim.to.plot.Body <- sample(c(1:1000), size = 100, replace = F)
sim.to.plot.Weapon <- sample(c(1001:2001), size = 100, replace = F)

for(i in sim.to.plot.Body) {
  lines(density(promit.sim[,i]),
        col = 'grey')
}

lines(density(promit$Dorsal.scute), lwd = 2, col = 'black')

for(i in sim.to.plot.Weapon) {
  lines(density(promit.sim[,i]), col = 'grey')
}

lines(density(promit$C4A), lwd = 2)

#dev.off()

