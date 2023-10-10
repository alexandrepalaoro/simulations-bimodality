#####################
### AICc function ###
#####################
myAICc <- function(aic, k, n) {
  aic + ((2*(k^2) + 2*k) / (n - k - 1))
}

#####################################
### For bivariate with dimorphism ###
#####################################
Mixmodels.allometries <- function(bi.trait, label.body, label.trait) {
  ## Mixture models
  SN.SkNorm2 <- smsn.mix(y = na.omit(bi.trait$trait), g=2, nu=6, obs.prob=T, family="Skew.normal") # k = 2 + 2 + 2 + 1 = 7
  SN.SkNorm1 <- smsn.mix(y = na.omit(bi.trait$trait), g=1, nu=6, obs.prob=T, family="Skew.normal") # k = 1 + 1 + 1 = 3
  AIC2dist <- myAICc(aic = SN.SkNorm2$aic, k = 7, n = SN.SkNorm2$n) # k = 2 + 2 + 2 + 1 = 7
  AIC1dist <- myAICc(aic = SN.SkNorm1$aic, k = 3, n = SN.SkNorm1$n) # k = 1 + 1 + 1 = 3
  
  par(mfrow = c(1, 2))
  hist(bi.trait$trait, prob = T,  xlab = label.trait,
       ylab="Frequency distribution", las=1,breaks=25,
       col=alpha('grey30',0.5),border=alpha('grey30',0.15),main="")
  mix.lines(y=na.omit(bi.trait$trait), model=SN.SkNorm1, lty = 2)
  mix.lines(y=na.omit(bi.trait$trait), model=SN.SkNorm2)
  
  bi.trait$Prob.1 <- SN.SkNorm2$obs.prob[,1]
  bi.trait$Prob.2 <- SN.SkNorm2$obs.prob[,2]
  
  ## WHICH PROB IS MAJOR AND WHICH ONE IS MINOR?
  ifelse(mean(bi.trait$trait[bi.trait$Prob.1 >= 0.8]) >= mean(bi.trait$trait[bi.trait$Prob.2 >= 0.8]),
         colnames(bi.trait) <- c("body", "trait", "Prob.MAJ", "Prob.MIN"),
         colnames(bi.trait) <- c("body", "trait", "Prob.MIN", "Prob.MAJ"))
  
  plot(trait ~ body, data = bi.trait,
       cex = 1.2, ylab = label.trait, xlab = label.body,
       frame.plot = F)
  
  points(trait[bi.trait$Prob.MAJ >= 0.80] ~ body[bi.trait$Prob.MAJ >= 0.80],
         data = bi.trait, cex = 1.2, pch = 19)
  
  points(trait[bi.trait$Prob.MAJ < 0.8 & bi.trait$Prob.MIN < 0.8] ~
           body[bi.trait$Prob.MAJ < 0.8 & bi.trait$Prob.MIN < 0.8],
         data = bi.trait, cex = 1.2, pch = 19, col = "gray")
  
  x.split <- mean(bi.trait$body[bi.trait$Prob.MAJ < 0.8 & bi.trait$Prob.MIN < 0.8])
  abline(v = x.split, lty = 2, col = "gray")
  
  # Allometrical slopes
  minors <-  lmodel2(trait  ~ body, data = bi.trait[bi.trait$Prob.MIN >= 0.8,])
  MIN.slope <- minors$regression.results$Slope[3]
  MIN.int <- minors$regression.results$Intercept[3]
  
  # confidence interval
  MIN.slope2.5 <- minors$confidence.intervals[3,4]
  MIN.slope97.5 <- minors$confidence.intervals[3,5]
  MIN.int2.5 <- minors$confidence.intervals[3,2]
  MIN.int97.5 <- minors$confidence.intervals[3,3]
  
  newx <- data.frame(seq(min(bi.trait$body), max(bi.trait$body[bi.trait$Prob.MIN >= 0.8]), length.out = 100))
  colnames(newx) <- "body"
  newx$newy <- (MIN.slope * newx$body) + minors$regression.results$Intercept[3]
  points(newx$body, newx$newy, type = "l", lty = 2, lwd = 2)
  
  majors <-  lmodel2(trait  ~ body, data = bi.trait[bi.trait$Prob.MAJ >= 0.8,])
  MAJ.slope <- majors$regression.results$Slope[3]
  MAJ.int <- majors$regression.results$Intercept[3]
  
  # confidence interval
  MAJ.slope2.5 <- majors$confidence.intervals[3,4]
  MAJ.slope97.5 <- majors$confidence.intervals[3,5]
  MAJ.int2.5 <- majors$confidence.intervals[3,2]
  MAJ.int97.5 <- majors$confidence.intervals[3,3]
  
  newx2 <- data.frame(seq(min(bi.trait$body[bi.trait$Prob.MAJ >= 0.8]), max(bi.trait$body), length.out = 100))
  colnames(newx2) <- "body"
  newx2$newy2 <- (MAJ.slope * newx2$body) + majors$regression.results$Intercept[3]
  points(newx2$body, newx2$newy2, type = "l", lwd = 2)
  
  ## testing major vs minor allometries
  combined.morphsMAJ <- data.frame(bi.trait$body[bi.trait$Prob.MAJ >= 0.80],
                                bi.trait$trait[bi.trait$Prob.MAJ >= 0.80])
  colnames(combined.morphsMAJ) <- c("body", "trait")
  combined.morphsMAJ$morph <- rep("major", times = length(combined.morphsMAJ[,1]))
  mean.Maj <- mean(combined.morphsMAJ$trait)
  sd.Maj <- sd(combined.morphsMAJ$trait)
  combined.morphsMIN <- data.frame(bi.trait$body[bi.trait$Prob.MIN >= 0.80],
                                   bi.trait$trait[bi.trait$Prob.MIN >= 0.80])
  colnames(combined.morphsMIN) <- c("body", "trait")
  combined.morphsMIN$morph <- rep("minor", times = length(combined.morphsMIN[,1]))
  mean.Min <- mean(combined.morphsMIN$trait)
  sd.Min <- sd(combined.morphsMIN$trait)
  combined.morphs <- rbind(combined.morphsMAJ, combined.morphsMIN)
  test.allo <-  sma(trait  ~ body*morph, data = combined.morphs)
  P.slopetest <- test.allo$commoncoef$p
  
  results <- data.frame(AIC1dist,AIC2dist, x.split,
                   MIN.slope, MIN.slope2.5, MIN.slope97.5,
                   MIN.int, MIN.int2.5, MIN.int97.5,
                   MAJ.slope, MAJ.slope2.5, MAJ.slope97.5,
                   MAJ.int, MAJ.int2.5, MAJ.int97.5, P.slopetest,
                   mean.Maj, sd.Maj, mean.Min, sd.Min)
  print(results)
}

########################################
### For bivariate without dimorphism ###
########################################
Mono.allometries <- function(mono.trait, label.body, label.trait) {
  
  ## Mixture models
  SN.SkNorm2 <- smsn.mix(y = na.omit(mono.trait$trait), g=2, nu=6, obs.prob=T, family="Skew.normal") # k = 2 + 2 + 2 + 1 = 7
  SN.SkNorm1 <- smsn.mix(y = na.omit(mono.trait$trait), g=1, nu=6, obs.prob=T, family="Skew.normal") # k = 1 + 1 + 1 = 3
  AIC2dist <- myAICc(aic = SN.SkNorm2$aic, k = 7, n = SN.SkNorm2$n) # k = 2 + 2 + 2 + 1 = 7
  AIC1dist <- myAICc(aic = SN.SkNorm1$aic, k = 3, n = SN.SkNorm1$n) # k = 1 + 1 + 1 = 3
  
  par(mfrow = c(1, 2))
  hist(mono.trait$trait, prob = T,  xlab = label.trait,
       ylab="Frequency distribution", las=1,breaks=25,
       col=alpha('grey30',0.5),border=alpha('grey30',0.15),main="")
  mix.lines(y=na.omit(mono.trait$trait), model=SN.SkNorm1)
  mix.lines(y=na.omit(mono.trait$trait), model=SN.SkNorm2, lty = 2)
  
  plot(trait ~ body, data = mono.trait,
       cex = 1.2, ylab = label.trait, xlab = label.body,
       frame.plot = F)
  
  # Allometrical slope
  allometry <-  lmodel2(trait  ~ body, data = mono.trait)
  slope <- allometry$regression.results$Slope[3]
  int <- allometry$regression.results$Intercept[3]
  
  # confidence interval
  slope2.5 <- allometry$confidence.intervals[3,4]
  slope97.5 <- allometry$confidence.intervals[3,5]
  int2.5 <- allometry$confidence.intervals[3,2]
  int97.5 <- allometry$confidence.intervals[3,3]
  
  newx <- data.frame(seq(min(mono.trait$body), max(mono.trait$body), length.out = 100))
  colnames(newx) <- "body"
  newx$newy <- (slope * newx$body) + allometry$regression.results$Intercept[3]
  points(newx$body, newx$newy, type = "l", lwd = 2)
  
  results <- cbind(AIC1dist,AIC2dist,
                   slope, slope2.5, slope97.5,
                   int, int2.5, int97.5)
  print(results)
}



######################
### for univariate ###
######################
Mixture.morphs <- function(bi.trait, label.trait) {
  ## Mixture models
  SN.SkNorm2 <- smsn.mix(y = na.omit(bi.trait$trait), g=2, nu=6, obs.prob=T, family="Skew.normal") # k = 2 + 2 + 2 + 1 = 7
  SN.SkNorm1 <- smsn.mix(y = na.omit(bi.trait$trait), g=1, nu=6, obs.prob=T, family="Skew.normal") # k = 1 + 1 + 1 = 3
  AIC2dist <- myAICc(aic = SN.SkNorm2$aic, k = 7, n = SN.SkNorm2$n) # k = 2 + 2 + 2 + 1 = 7
  AIC1dist <- myAICc(aic = SN.SkNorm1$aic, k = 3, n = SN.SkNorm1$n) # k = 1 + 1 + 1 = 3
  
  par(mfrow = c(1, 1))
  hist(bi.trait$trait, prob = T,  xlab = label.trait,
       ylab="Frequency distribution", las=1,breaks=25,
       col=alpha('grey30',0.5),border=alpha('grey30',0.15),main="")
  mix.lines(y=na.omit(bi.trait$trait), model=SN.SkNorm1, lty = 2)
  mix.lines(y=na.omit(bi.trait$trait), model=SN.SkNorm2)
  
  bi.trait$Prob.1 <- SN.SkNorm2$obs.prob[,1]
  bi.trait$Prob.2 <- SN.SkNorm2$obs.prob[,2]
  
  ## WHICH PROB IS MAJOR AND WHICH ONE IS MINOR?
  ifelse(mean(bi.trait$trait[bi.trait$Prob.1 >= 0.8]) >= mean(bi.trait$trait[bi.trait$Prob.2 >= 0.8]),
         colnames(bi.trait) <- c("trait", "Prob.MAJ", "Prob.MIN"),
         colnames(bi.trait) <- c("trait", "Prob.MIN", "Prob.MAJ"))
  mean.Maj <- mean(bi.trait$trait[bi.trait$Prob.MAJ >= 0.8])
  sd.Maj <- sd(bi.trait$trait[bi.trait$Prob.MAJ >= 0.8])
  mean.Min <- mean(bi.trait$trait[bi.trait$Prob.MIN >= 0.8])
  sd.Min <- sd(bi.trait$trait[bi.trait$Prob.MIN >= 0.8])
  
  results <- data.frame(AIC1dist, AIC2dist,
                   mean.Maj, sd.Maj, mean.Min, sd.Min)
  print(results)
}