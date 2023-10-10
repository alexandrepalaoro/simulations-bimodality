library(parallel)
library(doParallel)
library(tidyverse)
source("./sim/Simulating sigmoid allometries and detecting dimorphism.R")

sim.male.body = data.frame(matrix(NA ,
                                  nrow = length(males2$body),
                                  ncol = 1000))

sim.male.weap = data.frame(matrix(NA ,
                                  nrow = length(males2$body),
                                  ncol = 1000))

sims = vector( mode = 'list', length = 1000)

sim.ef.size = data.frame(med.major = 0,
                         med.minor = 0)

p.value = 0


par.sim = function(x) {

  require(grofit)
  require(mixsmsn)
  require(tidyverse)
  
sim.male.body = data.frame(matrix(NA ,
                                    nrow = length(males2$body),
                                    ncol = 1000))
sim.male.weap = data.frame(matrix(NA ,
                                  nrow = length(males2$body),
                                  ncol = 1000))

sims = vector( mode = 'list', length = 1000)

sim.ef.size = data.frame(med.major = 0,
                         med.minor = 0)

ric = Richards.Model2

for(i in 1:1000) {
  sim.male.body[i] = sort(rnorm(n = length(males2$body),
                                mean = mean(males2$body),
                                sd = sd(males2$body)))
}

for(i in 1:1000) {
  sim.male.weap[i] = predict(Richards.Model2, list(body=sim.male.body[,i]))
}



for(i in 1:1000) {
  sims[[i]] = smsn.mix(y = na.omit(sim.male.weap[,i]), g=2, nu=6, obs.prob=T, family="Skew.normal")
}


for(i in 1:1000) {
  sim.ef.size[i,1] = sort(sims[[i]]$mu, decreasing = T)[1]
  sim.ef.size[i,2] = sort(sims[[i]]$mu, decreasing = T)[2]
}


sim.ef.size = sim.ef.size %>%
  mutate(yi = log(med.major/med.minor))

ef.size = data.frame(med.major = sort(SN.SkNorm2$mu, decreasing = T)[1], 
                     med.minor = sort(SN.SkNorm2$mu, decreasing = T)[2] )

ef.size = ef.size %>%
  mutate(yi = log(med.major/med.minor))

p.value = sum(ef.size$yi < sim.ef.size$yi)/length(sim.ef.size$yi)

return(p.value)

}

no_cores = 15
cl = makeCluster(no_cores)

system.time( {
  registerDoParallel(cl)
  results = foreach(i = 1:100, .combine = 'cbind') %dopar% par.sim()
 } )

results


hist(results, main = "", prob = T, 
     xlab = "'P.values'")
abline(v = mean(results), lwd = 5)


