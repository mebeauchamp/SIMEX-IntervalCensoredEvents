###############################################################################
###                                                                         ###
### Code to run the simulations and preparing results presented in the      ### 
### manuscript:                                                             ###
### "Adapting SIMEX to correct for bias due to interval-censored outcomes   ###
### in survival analysis with time-varying exposure"                        ###
###                                                                         ###
###############################################################################

## Last update: March 10, 2022

rm(list=ls())
gc()

options(warn=1)

## Set the current folder
#setwd("C://path//to//your//directory")


###############################################################################
### Running the simulations 
###############################################################################

### Source functions for simulations ###
source("Functions_for_simulations.R")

### Main simulations with duration-based exposure generation ###

## Scenario 1 
system.time( simulation.fct(n.sims=1000, N=3000, fup=240, incidence.rate=0.333, 
               expoGen="Duration", 
               beta.curUse=log(3), beta.cumDur12w=log(2)/6, 
               beta.sex=0.8, beta.age=0.05, 
               addTVcov=0, beta.TVcov=NULL,
               P=10, move2ndLast.SIMEX=seq(2, 50, by=2), 
               seed=90001) )

## Scenario 2
system.time( simulation.fct(n.sims=1000, N=3000, fup=240, incidence.rate=0.333, 
               expoGen="Duration", 
               beta.curUse=log(2), beta.cumDur12w=log(3)/12, 
               beta.sex=0.8, beta.age=0.05, 
               addTVcov=0, beta.TVcov=NULL,
               P=10, move2ndLast.SIMEX=seq(2, 50, by=2), 
               seed=90002) )

## Scenario 3 
system.time( simulation.fct(n.sims=1000, N=1500, fup=240, incidence.rate=0.333, 
               expoGen="Duration", 
               beta.curUse=log(3), beta.cumDur12w=log(2)/6, 
               beta.sex=0.8, beta.age=0.05, 
               addTVcov=0, beta.TVcov=NULL,
               P=10, move2ndLast.SIMEX=seq(2, 50, by=2), 
               seed=90003) )

## Scenario 4 
system.time( simulation.fct(n.sims=1000, N=1500, fup=240, incidence.rate=0.333, 
               expoGen="Duration", 
               beta.curUse=log(2), beta.cumDur12w=log(3)/12, 
               beta.sex=0.8, beta.age=0.05,  
               addTVcov=0, beta.TVcov=NULL,
               P=10, move2ndLast.SIMEX=seq(2, 50, by=2), 
               seed=90004) )

## Scenario 5 
system.time( simulation.fct(n.sims=1000, N=750, fup=240, incidence.rate=0.333, 
               expoGen="Duration", 
               beta.curUse=log(3), beta.cumDur12w=log(2)/6, 
               beta.sex=0.8, beta.age=0.05, 
               addTVcov=0, beta.TVcov=NULL,
               P=10, move2ndLast.SIMEX=seq(2, 50, by=2), 
               seed=90005) )

## Scenario 6
system.time( simulation.fct(n.sims=1000, N=750, fup=240, incidence.rate=0.333, 
               expoGen="Duration", 
               beta.curUse=log(2), beta.cumDur12w=log(3)/12, 
               beta.sex=0.8, beta.age=0.05, 
               addTVcov=0, beta.TVcov=NULL,
               P=10, move2ndLast.SIMEX=seq(2, 50, by=2), 
               seed=90006) )

### Additional simulations with duration-based exposure generation ###

## Scenario 7 
system.time( simulation.fct(n.sims=1000, N=3000, fup=240, incidence.rate=0.333, 
               expoGen="Duration", 
               beta.curUse=log(3), beta.cumDur12w=log(2)/6, 
               beta.sex=0.8, beta.age=0.05, 
               addTVcov=1, beta.TVcov=log(3),
               P=10, move2ndLast.SIMEX=seq(2, 50, by=2), 
               seed=90007) )

## Scenario 8 
system.time( 
  simulation.fct(n.sims=1000, N=3000, fup=240, incidence.rate=0.333, 
               expoGen="Duration", 
               beta.curUse=log(1), beta.cumDur12w=log(1), 
               beta.sex=0.8, beta.age=0.05, 
               addTVcov=0, beta.TVcov=NULL,
               P=10, move2ndLast.SIMEX=seq(2, 50, by=2), 
               seed=90008) )

### Additional simulations with change-of-status exposure generation ###

## Scenario 9
system.time( simulation.fct(n.sims=1000, N=3000, fup=240, incidence.rate=0.333, 
               expoGen="ChangeStatus", 
               beta.curUse=log(3), beta.cumDur12w=log(2)/6, 
               beta.sex=0.8, beta.age=0.05, 
               addTVcov=0, beta.TVcov=NULL,
               P=10, move2ndLast.SIMEX=seq(2, 50, by=2), 
               seed=90009) )

## Scenario 10 
system.time( simulation.fct(n.sims=1000, N=3000, fup=240, incidence.rate=0.333, 
               expoGen="ChangeStatus", 
               beta.curUse=log(2), beta.cumDur12w=log(3)/12, 
               beta.sex=0.8, beta.age=0.05, 
               addTVcov=0, beta.TVcov=NULL,
               P=10, move2ndLast.SIMEX=seq(2, 50, by=2), 
               seed=90010) )

## Scenario 11
system.time( simulation.fct(n.sims=1000, N=1500, fup=240, incidence.rate=0.333, 
               expoGen="ChangeStatus", 
               beta.curUse=log(3), beta.cumDur12w=log(2)/6, 
               beta.sex=0.8, beta.age=0.05, 
               addTVcov=0, beta.TVcov=NULL,
               P=10, move2ndLast.SIMEX=seq(2, 50, by=2), 
               seed=90011) )

## Scenario 12
system.time( simulation.fct(n.sims=1000, N=1500, fup=240, incidence.rate=0.333, 
               expoGen="ChangeStatus", 
               beta.curUse=log(2), beta.cumDur12w=log(3)/12, 
               beta.sex=0.8, beta.age=0.05, 
               addTVcov=0, beta.TVcov=NULL,
               P=10, move2ndLast.SIMEX=seq(2, 50, by=2), 
               seed=90012) )

## Scenario 13 
system.time( simulation.fct(n.sims=1000, N=750, fup=240, incidence.rate=0.333, 
               expoGen="ChangeStatus", 
               beta.curUse=log(3), beta.cumDur12w=log(2)/6, 
               beta.sex=0.8, beta.age=0.05, 
               addTVcov=0, beta.TVcov=NULL,
               P=10, move2ndLast.SIMEX=seq(2, 50, by=2), 
               seed=90013) )

## Scenario 14 
system.time( simulation.fct(n.sims=1000, N=750, fup=240, incidence.rate=0.333, 
               expoGen="ChangeStatus", 
               beta.curUse=log(2), beta.cumDur12w=log(3)/12, 
               beta.sex=0.8, beta.age=0.05, 
               addTVcov=0, beta.TVcov=NULL,
               P=10, move2ndLast.SIMEX=seq(2, 50, by=2), 
               seed=90014) )


###############################################################################
### Results 
###############################################################################

## Set the current folder
#setwd("C://path//to//your//directory")


#-----------------------------------------------------------------------------#
# FIGURE 2:                                                                   #
#   Diagnostics for the selection of the extrapolating functions, for each of # 
#   the two exposure metrics: a) current drug use U(t) and b) cumulative      #
#   duration of drug use in the past 12 weeks Cum12(t), for scenario 1        #
#   (N=3,000, stronger exposure effects, and duration-based exposure          #
#   generation)                                                               #
#                                                                             #
# and                                                                         #
#                                                                             #
# FIGURE B.1:                                                                 #
#   Diagnostics for the selection of the extrapolating functions, for each of # 
#   the two exposure metrics: a) current drug use U(t) and b) cumulative      #
#   duration of drug use in the past 12 weeks Cum12(t), for scenario 9        #
#   (N=3,000, stronger exposure effects, and change-of-status exposure        #
#   generation)                                                               #
#-----------------------------------------------------------------------------#

################
# Load results #
################

# Load results for 1 scenario at a time, and produce the results for this 
# Table/Figure. Then, delete objects and load results for the next scenario.

rm(list=ls())

# For Figure 2 (scenario 1)
load("Sc1 - N3000 expoDuration BcurUse11 BcumDur12w012 - 1000samples.RData")

# For Figure B.1 (scenario 9)
#load("Sc9 - N3000 expoChangeStatus BcurUse11 BcumDur12w012 - 1000samples.RData")

##########
# Figure #
##########

if (sc.name == "Sc1 - N3000 expoDuration BcurUse11 BcumDur12w012 - 1000samples"){
  fig.no <- '2'
  low.nbr.it.cumDur12w <- 15  # Lower number of iterations in SIMEX procedure
} else if (sc.name == "Sc9 - N3000 expoChangeStatus BcurUse11 BcumDur12w012 - 1000samples"){
  fig.no <- 'B1'
  low.nbr.it.cumDur12w <- 11  # Lower number of iterations in SIMEX procedure
}

png(file=paste0("Figure ", fig.no, ".png"), units="in", height=5, width=9, res=600)
#pdf(file=paste0("Figure ", fig.no, ".pdf"), height=5, width=9)

par(mfrow=c(1,2))

plot(c(0,40), c(0,1), type='n', 
     xlab='Time since imputed event time (weeks)', 
     ylab='Correlation', 
     main='Current drug use')
for (s in 1:10){
  lines(0:40, cor.expo.cases.curUse[s,1:41], col=grey(0.7))
}
lines(0:40, apply(cor.expo.cases.curUse, 2, mean)[1:41], lwd=3, col=1)
it.low <- mean(mean.dist2lastV.curUse[,16]) 
it.max <- mean(mean.dist2lastV.curUse[,26]) 
lines(c(mean(mean.dist2lastV.curUse[,1]),mean(mean.dist2lastV.curUse[,1]))/4, 
      c(0,1), lty=3)
lines(c(it.low,it.low)/4, c(0,1), lty=2)
lines(c(it.max,it.max)/4, c(0,1), lty=2)
legend('topright', c('Individual sample','Mean across 1000 samples'), lty=c(1,1), 
       col=c(grey(0.7),1), lwd=c(1,3), inset=0.01, bg="white")

plot(c(0,40), c(0,1), type='n', 
     xlab='Time since imputed event time (weeks)', 
     ylab='Correlation',  
     main='Cumulative duration of drug use\n in the past 12 weeks')
for (s in 1:10){
  lines(0:40, cor.expo.cases.cumDur12w[s,1:41], col=grey(0.7))
}
lines(0:40, apply(cor.expo.cases.cumDur12w, 2, mean)[1:41], lwd=3, col=1)
it.low <- mean(mean.dist2lastV.cumDur12w[,(low.nbr.it.cumDur12w+1)])
it.max <- mean(mean.dist2lastV.cumDur12w[,26])
lines(c(mean(mean.dist2lastV.cumDur12w[,1]), mean(mean.dist2lastV.cumDur12w[,1]))/4, 
      c(0,7), lty=3)
lines(c(it.low,it.low)/4, c(0,7), lty=2)
lines(c(it.max,it.max)/4, c(0,7), lty=2)
legend('topright', c('Individual sample','Mean across 1000 samples'), 
       lty=c(1,1), col=c(grey(0.7),1), lwd=c(1,3), inset=0.01, bg="white")

dev.off()

# Mean difference between the imputed vs. the true event times 
mean(mean.dist2lastV.curUse[,1])
mean(mean.dist2lastV.curUse[,1])/4

# Mean difference between the last 2 visits for k=15 and k=25 iterations
it.low
it.low/4
it.max
it.max/4


#-----------------------------------------------------------------------------#
# TABLE 1:                                                                    #
#   Comparison of model-specific log HRs for main scenarios with              #
#   duration-based exposure generation                                        #
#                                                                             #
# and                                                                         #
#                                                                             #
# TABLE 2:                                                                    #
#   Comparison of model-specific log HRs for additional scenarios with        #
#   duration-based exposure generation                                        #
#-----------------------------------------------------------------------------#

################
# Load results #
################

# Load results for 1 scenario at a time, and produce the results for this 
# Table/Figure. Then, delete objects and load results for the next scenario.

rm(list=ls())

## For Table 1
load("Sc1 - N3000 expoDuration BcurUse11 BcumDur12w012 - 1000samples.RData")
#load("Sc2 - N3000 expoDuration BcurUse069 BcumDur12w009 - 1000samples.RData")
#load("Sc3 - N1500 expoDuration BcurUse11 BcumDur12w012 - 1000samples.RData")
#load("Sc4 - N1500 expoDuration BcurUse069 BcumDur12w009 - 1000samples.RData")
#load("Sc5 - N750 expoDuration BcurUse11 BcumDur12w012 - 1000samples.RData")
#load("Sc6 - N750 expoDuration BcurUse069 BcumDur12w009 - 1000samples.RData")

## For Table 2
#load("Sc7 - N3000 expoDuration TVcov BcurUse11 BcumDur12w012 - 1000samples.RData")
#load("Sc8 - N3000 expoDuration BcurUse0 BcumDur12w0 - 1000samples.RData")


################################
# Exposure metric: Current use #
################################

### N [mean nbr events] ###

N
#round(mean(n.events))
round(mean(n.events.curUse.mid))

### True log HR ###

exp(beta.curUse)
round(beta.curUse, digits=2)

### Relative bias ###

# Conventional, event middle interval
round((mean(cox.curUse.coef.mid[1,1,]) - beta.curUse) / 
        beta.curUse, digits=3) * 100
# SIMEX, Linear, 15 iterations
round(mean(cox.curUse.expo.SIMEXlin.15it.extrap - beta.curUse) / 
        beta.curUse, digits=3) * 100 
# SIMEX, Linear, 25 iterations
round(mean(cox.curUse.expo.SIMEXlin.25it.extrap - beta.curUse) / 
        beta.curUse, digits=3) * 100

if (sc.name == "Sc8 - N3000 expoDuration BcurUse0 BcumDur12w0 - 1000samples"){
  # Report bias (because true parameter is 0 and then relative bias is Inf)
  # Conventional, event middle interval
  print(round((mean(cox.curUse.coef.mid[1,1,]) - beta.curUse), digits=3))
  # SIMEX, Linear, 15 iterations
  print(round(mean(cox.curUse.expo.SIMEXlin.15it.extrap - beta.curUse), digits=3)) 
  # SIMEX, Linear, 25 iterations
  print(round(mean(cox.curUse.expo.SIMEXlin.25it.extrap - beta.curUse), digits=3))
}

### Standard deviation (SD) ###

# Conventional, event middle interval
round(sd(cox.curUse.coef.mid[1,1,]), digits=3)
# SIMEX, Linear, 15 iterations
round(sd(cox.curUse.expo.SIMEXlin.15it.extrap), digits=3) 
# SIMEX, Linear, 25 iterations
round(sd(cox.curUse.expo.SIMEXlin.25it.extrap), digits=3)

### Root mean squared error (RMSE) ###

rmse.fct <- function(est, true){
  sqrt((mean(est)-true)^2 + var(est))
}
# Conventional, event middle interval
round(rmse.fct(cox.curUse.coef.mid[1,1,], beta.curUse), digits=3)
# SIMEX, Linear, 15 iterations
round(rmse.fct(cox.curUse.expo.SIMEXlin.15it.extrap, beta.curUse), digits=3) 
# SIMEX, Linear, 25 iterations
round(rmse.fct(cox.curUse.expo.SIMEXlin.25it.extrap, beta.curUse), digits=3)

# Reduction in RMSE for conventional vs. SIMEX
if (sc.name == "Sc6 - N750 expoDuration BcurUse069 BcumDur12w009 - 1000samples"){
  print((1 - (rmse.fct(cox.curUse.expo.SIMEXlin.15it.extrap, beta.curUse) / 
          rmse.fct(cox.curUse.coef.mid[1,1,], beta.curUse))) * 100)
}

#########################################################
# Exposure metric: Cumulative duration in past 12 weeks #
#########################################################

### N [mean nbr events] ###

N
#round(mean(n.events))
round(mean(n.events.cumDur12w.mid))

### True log HR ###

round(beta.cumDur12w, digits=2)

### Relative bias ###

# Conventional, event middle interval
round((mean(cox.cumDur12w.coef.mid[1,1,]) - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100
# SIMEX, Linear, 15 iterations
round(mean(cox.cumDur12w.expo.SIMEXlin.15it.extrap - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100 
# SIMEX, Linear, 25 iterations
round(mean(cox.cumDur12w.expo.SIMEXlin.25it.extrap - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100

if (sc.name == "Sc8 - N3000 expoDuration BcurUse0 BcumDur12w0 - 1000samples"){
  # Report bias (because true parameter is 0 and then relative bias is Inf)
  # Conventional, event middle interval
  print(round((mean(cox.cumDur12w.coef.mid[1,1,]) - beta.cumDur12w), digits=3))
  # SIMEX, Linear, 15 iterations
  print(round(mean(cox.cumDur12w.expo.SIMEXlin.15it.extrap - beta.cumDur12w), digits=3)) 
  # SIMEX, Linear, 25 iterations
  print(round(mean(cox.cumDur12w.expo.SIMEXlin.25it.extrap - beta.cumDur12w), digits=3))
}

### Standard deviation (SD) ###

# Conventional, event middle interval
round(sd(cox.cumDur12w.coef.mid[1,1,]), digits=3)
# SIMEX, Linear, 15 iterations
round(sd(cox.cumDur12w.expo.SIMEXlin.15it.extrap), digits=3) 
# SIMEX, Linear, 25 iterations
round(sd(cox.cumDur12w.expo.SIMEXlin.25it.extrap), digits=3)

### Root mean squared error (RMSE) ###

rmse.fct <- function(est, true){
  sqrt((mean(est)-true)^2 + var(est))
}
# Conventional, event middle interval
round(rmse.fct(cox.cumDur12w.coef.mid[1,1,], beta.cumDur12w), digits=3)
# SIMEX, Linear, 15 iterations
round(rmse.fct(cox.cumDur12w.expo.SIMEXlin.15it.extrap, beta.cumDur12w), digits=3) 
# SIMEX, Linear, 25 iterations
round(rmse.fct(cox.cumDur12w.expo.SIMEXlin.25it.extrap, beta.cumDur12w), digits=3)

# Reduction in RMSE for conventional vs. SIMEX
if (sc.name == "Sc6 - N750 expoDuration BcurUse069 BcumDur12w009 - 1000samples"){
  print((1 - (rmse.fct(cox.cumDur12w.expo.SIMEXlin.15it.extrap, beta.cumDur12w) / 
                rmse.fct(cox.cumDur12w.coef.mid[1,1,], beta.cumDur12w))) * 100)
}


#-----------------------------------------------------------------------------#
# TABLE B.1:                                                                  #
#   Comparison of additional model-specific log HRs for main scenarios with   #
#   duration-based exposure generation                                        #
#                                                                             #
# and                                                                         #
#                                                                             #
# TABLE B.3:                                                                  #
#   Comparison of additional model-specific log HRs for additional scenarios  #
#   with duration-based exposure generation                                   #
#                                                                             #
# and                                                                         #
#                                                                             #
# TABLE B.5:                                                                  #
#   Comparison of additional model-specific log HRs for additional scenarios  #
#   with change-of-status exposure generation                                 #
#-----------------------------------------------------------------------------#

################
# Load results #
################

# Load results for 1 scenario at a time, and produce the results for this 
# Table/Figure. Then, delete objects and load results for the next scenario.

rm(list=ls())

## For Table B.1
load("Sc1 - N3000 expoDuration BcurUse11 BcumDur12w012 - 1000samples.RData")
#load("Sc2 - N3000 expoDuration BcurUse069 BcumDur12w009 - 1000samples.RData")
#load("Sc3 - N1500 expoDuration BcurUse11 BcumDur12w012 - 1000samples.RData")
#load("Sc4 - N1500 expoDuration BcurUse069 BcumDur12w009 - 1000samples.RData")
#load("Sc5 - N750 expoDuration BcurUse11 BcumDur12w012 - 1000samples.RData")
#load("Sc6 - N750 expoDuration BcurUse069 BcumDur12w009 - 1000samples.RData")

## For Table B.3
#load("Sc7 - N3000 expoDuration TVcov BcurUse11 BcumDur12w012 - 1000samples.RData")
#load("Sc8 - N3000 expoDuration BcurUse0 BcumDur12w0 - 1000samples.RData")

## For Table B.5
#load("Sc9 - N3000 expoChangeStatus BcurUse11 BcumDur12w012 - 1000samples.RData")
#load("Sc10 - N3000 expoChangeStatus BcurUse069 BcumDur12w009 - 1000samples.RData")
#load("Sc11 - N1500 expoChangeStatus BcurUse11 BcumDur12w012 - 1000samples.RData")
#load("Sc12 - N1500 expoChangeStatus BcurUse069 BcumDur12w009 - 1000samples.RData")
#load("Sc13 - N750 expoChangeStatus BcurUse11 BcumDur12w012 - 1000samples.RData")
#load("Sc14 - N750 expoChangeStatus BcurUse069 BcumDur12w009 - 1000samples.RData")

################################
# Exposure metric: Current use #
################################

### N [mean nbr events] ###

N
#round(mean(n.events))
round(mean(n.events.curUse.mid))

### True log HR ###

exp(beta.curUse)
round(beta.curUse, digits=2)

### Relative bias ###

# Oracle
round((mean(cox.curUse.coef.truth[1,1,]) - beta.curUse) / 
        beta.curUse, digits=3) * 100
# Conventional, event middle interval
round((mean(cox.curUse.coef.mid[1,1,]) - beta.curUse) / 
        beta.curUse, digits=3) * 100
# Conventional, event end interval
round((mean(cox.curUse.coef.end[1,1,]) - beta.curUse) / 
        beta.curUse, digits=3) * 100

if (sc.name == "Sc8 - N3000 expoDuration BcurUse0 BcumDur12w0 - 1000samples"){
  # Report bias (because true parameter is 0 and then relative bias is Inf)
  # Conventional, event middle interval
  print(round((mean(cox.curUse.coef.mid[1,1,]) - beta.curUse), digits=3))
  # Conventional, event middle interval
  print(round((mean(cox.curUse.coef.mid[1,1,]) - beta.curUse), digits=3))
  # Conventional, event end interval
  print(round((mean(cox.curUse.coef.end[1,1,]) - beta.curUse), digits=3))
}

### Standard deviation (SD) ###

# Oracle
round(sd(cox.curUse.coef.truth[1,1,]), digits=3)
# Conventional, event middle interval
round(sd(cox.curUse.coef.mid[1,1,]), digits=3)
# Conventional, event end interval
round(sd(cox.curUse.coef.end[1,1,]), digits=3)

### Root mean squared error (RMSE) ###

rmse.fct <- function(est, true){
  sqrt((mean(est)-true)^2 + var(est))
}
# Oracle
round(rmse.fct(cox.curUse.coef.truth[1,1,], beta.curUse), digits=3)
# Conventional, event middle interval
round(rmse.fct(cox.curUse.coef.mid[1,1,], beta.curUse), digits=3)
# Conventional, event end interval
round(rmse.fct(cox.curUse.coef.end[1,1,], beta.curUse), digits=3)

#########################################################
# Exposure metric: Cumulative duration in past 12 weeks #
#########################################################

### N [mean nbr events] ###

N
#round(mean(n.events))
round(mean(n.events.cumDur12w.mid))

### True log HR ###

round(beta.cumDur12w, digits=2)

### Relative bias ###

# Oracle
round((mean(cox.cumDur12w.coef.truth[1,1,]) - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100
# Conventional, event middle interval
round((mean(cox.cumDur12w.coef.mid[1,1,]) - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100
# Conventional, event end interval
round((mean(cox.cumDur12w.coef.end[1,1,]) - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100

if (sc.name == "Sc8 - N3000 expoDuration BcurUse0 BcumDur12w0 - 1000samples"){
  # Report bias (because true parameter is 0 and then relative bias is Inf)
  # Conventional, event middle interval
  print(round((mean(cox.cumDur12w.coef.mid[1,1,]) - beta.cumDur12w), digits=3))
  # Conventional, event middle interval
  print(round((mean(cox.cumDur12w.coef.mid[1,1,]) - beta.cumDur12w), digits=3))
  # Conventional, event end interval
  print(round((mean(cox.cumDur12w.coef.end[1,1,]) - beta.cumDur12w), digits=3))
}
### Standard deviation (SD) ###

# Oracle
round(sd(cox.cumDur12w.coef.truth[1,1,]), digits=3)
# Conventional, event middle interval
round(sd(cox.cumDur12w.coef.mid[1,1,]), digits=3)
# Conventional, event end interval
round(sd(cox.cumDur12w.coef.end[1,1,]), digits=3)

### Root mean squared error (RMSE) ###

rmse.fct <- function(est, true){
  sqrt((mean(est)-true)^2 + var(est))
}
# Oracle
round(rmse.fct(cox.cumDur12w.coef.truth[1,1,], beta.cumDur12w), digits=3)
# Conventional, event middle interval
round(rmse.fct(cox.cumDur12w.coef.mid[1,1,], beta.cumDur12w), digits=3)
# Conventional, event end interval
round(rmse.fct(cox.cumDur12w.coef.end[1,1,], beta.cumDur12w), digits=3)


#-----------------------------------------------------------------------------#
# TABLE B.2:                                                                  #
#  Comparison of SIMEX-based model-specific log HRs for additional            #
#  extrapolating functions for scenario 1, with duration-based exposure       #
#  generation                                                                 #
#-----------------------------------------------------------------------------#

rm(list=ls())

## Load results for scenario 1

load("Sc1 - N3000 expoDuration BcurUse11 BcumDur12w012 - 1000samples.RData")

################################
# Exposure metric: Current use #
################################

### N [mean nbr events] ###

N
#round(mean(n.events))
round(mean(n.events.curUse.mid))

### True log HR ###

exp(beta.curUse)
round(beta.curUse, digits=2)

### Relative bias ###

# Conventional, event middle interval
round((mean(cox.curUse.coef.mid[1,1,]) - beta.curUse) / 
        beta.curUse, digits=3) * 100
# SIMEX, Linear, 15 iterations
round(mean(cox.curUse.expo.SIMEXlin.15it.extrap - beta.curUse) / 
        beta.curUse, digits=3) * 100 
# SIMEX, Linear, 25 iterations
round(mean(cox.curUse.expo.SIMEXlin.25it.extrap - beta.curUse) / 
        beta.curUse, digits=3) * 100
# SIMEX, Log, 15 iterations
round(mean(cox.curUse.expo.SIMEXlog2.15it.extrap - beta.curUse) / 
        beta.curUse, digits=3) * 100 
# SIMEX, Log, 25 iterations
round(mean(cox.curUse.expo.SIMEXlog2.25it.extrap - beta.curUse) / 
        beta.curUse, digits=3) * 100
# SIMEX, 1-df quadratic, 15 iterations
round(mean(cox.curUse.expo.SIMEXquad1df.15it.extrap - beta.curUse) / 
        beta.curUse, digits=3) * 100 
# SIMEX, 1-df quadratic, 25 iterations
round(mean(cox.curUse.expo.SIMEXquad1df.25it.extrap - beta.curUse) / 
        beta.curUse, digits=3) * 100

### Standard deviation (SD) ###

# Conventional, event middle interval
round(sd(cox.curUse.coef.mid[1,1,]), digits=3)
# SIMEX, Linear, 15 iterations
round(sd(cox.curUse.expo.SIMEXlin.15it.extrap), digits=3) 
# SIMEX, Linear, 25 iterations
round(sd(cox.curUse.expo.SIMEXlin.25it.extrap), digits=3)
# SIMEX, Log, 15 iterations
round(sd(cox.curUse.expo.SIMEXlog2.15it.extrap), digits=3) 
# SIMEX, Log, 25 iterations
round(sd(cox.curUse.expo.SIMEXlog2.25it.extrap), digits=3)
# SIMEX, 1-df quadratic, 15 iterations
round(sd(cox.curUse.expo.SIMEXquad1df.15it.extrap), digits=3) 
# SIMEX, 1-df quadratic, 25 iterations
round(sd(cox.curUse.expo.SIMEXquad1df.25it.extrap), digits=3)

### Root mean squared error (RMSE) ###

rmse.fct <- function(est, true){
  sqrt((mean(est)-true)^2 + var(est))
}
# Conventional, event middle interval
round(rmse.fct(cox.curUse.coef.mid[1,1,], beta.curUse), digits=3)
# SIMEX, Linear, 15 iterations
round(rmse.fct(cox.curUse.expo.SIMEXlin.15it.extrap, beta.curUse), digits=3) 
# SIMEX, Linear, 25 iterations
round(rmse.fct(cox.curUse.expo.SIMEXlin.25it.extrap, beta.curUse), digits=3)
# SIMEX, Log, 15 iterations
round(rmse.fct(cox.curUse.expo.SIMEXlog2.15it.extrap, beta.curUse), digits=3) 
# SIMEX, Log, 25 iterations
round(rmse.fct(cox.curUse.expo.SIMEXlog2.25it.extrap, beta.curUse), digits=3)
# SIMEX, 1-df quadratic, 15 iterations
round(rmse.fct(cox.curUse.expo.SIMEXquad1df.15it.extrap, beta.curUse), digits=3) 
# SIMEX, 1-df quadratic, 25 iterations
round(rmse.fct(cox.curUse.expo.SIMEXquad1df.25it.extrap, beta.curUse), digits=3)

#########################################################
# Exposure metric: Cumulative duration in past 12 weeks #
#########################################################

### N [mean nbr events] ###

N
#round(mean(n.events))
round(mean(n.events.cumDur12w.mid))

### True log HR ###

round(beta.cumDur12w, digits=2)

### Relative bias ###

# Conventional, event middle interval
round((mean(cox.cumDur12w.coef.mid[1,1,]) - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100
# SIMEX, Linear, 15 iterations
round(mean(cox.cumDur12w.expo.SIMEXlin.15it.extrap - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100 
# SIMEX, Linear, 25 iterations
round(mean(cox.cumDur12w.expo.SIMEXlin.25it.extrap - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100
# SIMEX, Log, 15 iterations
round(mean(cox.cumDur12w.expo.SIMEXlog2.15it.extrap - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100 
# SIMEX, Log, 25 iterations
round(mean(cox.cumDur12w.expo.SIMEXlog2.25it.extrap - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100
# SIMEX, 1-df quadratic, 15 iterations
round(mean(cox.cumDur12w.expo.SIMEXquad1df.15it.extrap - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100 
# SIMEX, 1-df quadratic, 25 iterations
round(mean(cox.cumDur12w.expo.SIMEXquad1df.25it.extrap - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100

### Standard deviation (SD) ###

# Conventional, event middle interval
round(sd(cox.cumDur12w.coef.mid[1,1,]), digits=3)
# SIMEX, Linear, 15 iterations
round(sd(cox.cumDur12w.expo.SIMEXlin.15it.extrap), digits=3) 
# SIMEX, Linear, 25 iterations
round(sd(cox.cumDur12w.expo.SIMEXlin.25it.extrap), digits=3)
# SIMEX, Log, 15 iterations
round(sd(cox.cumDur12w.expo.SIMEXlog2.15it.extrap), digits=3) 
# SIMEX, Log, 25 iterations
round(sd(cox.cumDur12w.expo.SIMEXlog2.25it.extrap), digits=3)
# SIMEX, 1-df quadratic, 15 iterations
round(sd(cox.cumDur12w.expo.SIMEXquad1df.15it.extrap), digits=3) 
# SIMEX, 1-df quadratic, 25 iterations
round(sd(cox.cumDur12w.expo.SIMEXquad1df.25it.extrap), digits=3)

### Root mean squared error (RMSE) ###

rmse.fct <- function(est, true){
  sqrt((mean(est)-true)^2 + var(est))
}
# Conventional, event middle interval
round(rmse.fct(cox.cumDur12w.coef.mid[1,1,], beta.cumDur12w), digits=3)
# SIMEX, Linear, 15 iterations
round(rmse.fct(cox.cumDur12w.expo.SIMEXlin.15it.extrap, beta.cumDur12w), digits=3) 
# SIMEX, Linear, 25 iterations
round(rmse.fct(cox.cumDur12w.expo.SIMEXlin.25it.extrap, beta.cumDur12w), digits=3)
# SIMEX, Log, 15 iterations
round(rmse.fct(cox.cumDur12w.expo.SIMEXlog2.15it.extrap, beta.cumDur12w), digits=3) 
# SIMEX, Log, 25 iterations
round(rmse.fct(cox.cumDur12w.expo.SIMEXlog2.25it.extrap, beta.cumDur12w), digits=3)
# SIMEX, 1-df quadratic, 15 iterations
round(rmse.fct(cox.cumDur12w.expo.SIMEXquad1df.15it.extrap, beta.cumDur12w), digits=3) 
# SIMEX, 1-df quadratic, 25 iterations
round(rmse.fct(cox.cumDur12w.expo.SIMEXquad1df.25it.extrap, beta.cumDur12w), digits=3)


#-----------------------------------------------------------------------------#
# TABLE B.4:                                                                  #
#  Comparison of model-specific log HRs for main scenarios with               #
#  change-of-status exposure generation                                       #
#-----------------------------------------------------------------------------#

################
# Load results #
################

# Load results for 1 scenario at a time, and produce the results for this 
# Table/Figure. Then, delete objects and load results for the next scenario.

rm(list=ls())

load("Sc9 - N3000 expoChangeStatus BcurUse11 BcumDur12w012 - 1000samples.RData")
#load("Sc10 - N3000 expoChangeStatus BcurUse069 BcumDur12w009 - 1000samples.RData")
#load("Sc11 - N1500 expoChangeStatus BcurUse11 BcumDur12w012 - 1000samples.RData")
#load("Sc12 - N1500 expoChangeStatus BcurUse069 BcumDur12w009 - 1000samples.RData")
#load("Sc13 - N750 expoChangeStatus BcurUse11 BcumDur12w012 - 1000samples.RData")
#load("Sc14 - N750 expoChangeStatus BcurUse069 BcumDur12w009 - 1000samples.RData")

################################
# Exposure metric: Current use #
################################

### N [mean nbr events] ###

N
#round(mean(n.events))
round(mean(n.events.curUse.mid))

### True log HR ###

exp(beta.curUse)
round(beta.curUse, digits=2)

### Relative bias ###

# Conventional, event middle interval
round((mean(cox.curUse.coef.mid[1,1,]) - beta.curUse) / 
        beta.curUse, digits=3) * 100
# SIMEX, Linear, 15 iterations
round(mean(cox.curUse.expo.SIMEXlin.15it.extrap - beta.curUse) / 
        beta.curUse, digits=3) * 100 
# SIMEX, Linear, 25 iterations
round(mean(cox.curUse.expo.SIMEXlin.25it.extrap - beta.curUse) / 
        beta.curUse, digits=3) * 100
# SIMEX, Log, 15 iterations
round(mean(cox.curUse.expo.SIMEXlog2.15it.extrap - beta.curUse) / 
        beta.curUse, digits=3) * 100 
# SIMEX, Log, 25 iterations
round(mean(cox.curUse.expo.SIMEXlog2.25it.extrap - beta.curUse) / 
        beta.curUse, digits=3) * 100

### Standard deviation (SD) ###

# Conventional, event middle interval
round(sd(cox.curUse.coef.mid[1,1,]), digits=3)
# SIMEX, Linear, 15 iterations
round(sd(cox.curUse.expo.SIMEXlin.15it.extrap), digits=3) 
# SIMEX, Linear, 25 iterations
round(sd(cox.curUse.expo.SIMEXlin.25it.extrap), digits=3)
# SIMEX, Log, 15 iterations
round(sd(cox.curUse.expo.SIMEXlog2.15it.extrap), digits=3) 
# SIMEX, Log, 25 iterations
round(sd(cox.curUse.expo.SIMEXlog2.25it.extrap), digits=3)

### Root mean squared error (RMSE) ###

rmse.fct <- function(est, true){
  sqrt((mean(est)-true)^2 + var(est))
}
# Conventional, event middle interval
round(rmse.fct(cox.curUse.coef.mid[1,1,], beta.curUse), digits=3)
# SIMEX, Linear, 15 iterations
round(rmse.fct(cox.curUse.expo.SIMEXlin.15it.extrap, beta.curUse), digits=3) 
# SIMEX, Linear, 25 iterations
round(rmse.fct(cox.curUse.expo.SIMEXlin.25it.extrap, beta.curUse), digits=3)
# SIMEX, Log, 15 iterations
round(rmse.fct(cox.curUse.expo.SIMEXlog2.15it.extrap, beta.curUse), digits=3) 
# SIMEX, Log, 25 iterations
round(rmse.fct(cox.curUse.expo.SIMEXlog2.25it.extrap, beta.curUse), digits=3)


#########################################################
# Exposure metric: Cumulative duration in past 12 weeks #
#########################################################

### N [mean nbr events] ###

N
#round(mean(n.events))
round(mean(n.events.cumDur12w.mid))

### True log HR ###

round(beta.cumDur12w, digits=2)

### Relative bias ###

# Conventional, event middle interval
round((mean(cox.cumDur12w.coef.mid[1,1,]) - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100
# SIMEX, Linear, 11 iterations
round(mean(cox.cumDur12w.expo.SIMEXlin.11it.extrap - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100 
# SIMEX, Linear, 25 iterations
round(mean(cox.cumDur12w.expo.SIMEXlin.25it.extrap - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100
# SIMEX, Log, 15 iterations
round(mean(cox.cumDur12w.expo.SIMEXlog2.15it.extrap - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100 
# SIMEX, Log, 25 iterations
round(mean(cox.cumDur12w.expo.SIMEXlog2.25it.extrap - beta.cumDur12w) / 
        beta.cumDur12w, digits=3) * 100

### Standard deviation (SD) ###

# Conventional, event middle interval
round(sd(cox.cumDur12w.coef.mid[1,1,]), digits=3)
# SIMEX, Linear, 11 iterations
round(sd(cox.cumDur12w.expo.SIMEXlin.11it.extrap), digits=3) 
# SIMEX, Linear, 25 iterations
round(sd(cox.cumDur12w.expo.SIMEXlin.25it.extrap), digits=3)
# SIMEX, Log, 15 iterations
round(sd(cox.cumDur12w.expo.SIMEXlog2.15it.extrap), digits=3) 
# SIMEX, Log, 25 iterations
round(sd(cox.cumDur12w.expo.SIMEXlog2.25it.extrap), digits=3)

### Root mean squared error (RMSE) ###

rmse.fct <- function(est, true){
  sqrt((mean(est)-true)^2 + var(est))
}
# Conventional, event middle interval
round(rmse.fct(cox.cumDur12w.coef.mid[1,1,], beta.cumDur12w), digits=3)
# SIMEX, Linear, 11 iterations
round(rmse.fct(cox.cumDur12w.expo.SIMEXlin.11it.extrap, beta.cumDur12w), digits=3) 
# SIMEX, Linear, 25 iterations
round(rmse.fct(cox.cumDur12w.expo.SIMEXlin.25it.extrap, beta.cumDur12w), digits=3)
# SIMEX, Linear, 15 iterations
round(rmse.fct(cox.cumDur12w.expo.SIMEXlog2.15it.extrap, beta.cumDur12w), digits=3) 
# SIMEX, Linear, 25 iterations
round(rmse.fct(cox.cumDur12w.expo.SIMEXlog2.25it.extrap, beta.cumDur12w), digits=3)

