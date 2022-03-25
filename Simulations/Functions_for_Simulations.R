###############################################################################
###                                                                         ###
### Functions used in the simulations presented in the manuscript:          ###
### "Adapting SIMEX to correct for bias due to interval-censored outcomes   ###
### in survival analysis with time-varying exposure"                        ###
###                                                                         ###
###############################################################################

## Last update: March 10, 2022

#-----------------------------------------------------------------------------#
# genVisits.fct:                                                              #
#   Generation of visit times for each subject                                #
#-----------------------------------------------------------------------------#

## Arguments:
  # fup: maximum follow-up, in weeks.  
  # P: average number of visits generated during maximum follow-up time

genVisits.fct <- function(fup, P){
  
  # Expected distance between visits for the current subject
  d <- max(rmultinom(n=1, size=1, prob=c(0.2,0.2,0.2,0.2,0.2)) * 
             c(0.5,0.75,1,1.25,1.5)) * (fup/P)
  
  # Observed visit times until study end
  Vt.tmp <- round(cumsum(runif(n = (fup/(d * 0.5)), 0.5, 1.5) * d)) 
    # n = (fup/(d*0.5)) is the maximum number of visits possible
  Vt <- c(0, Vt.tmp[Vt.tmp <= fup])  # Always a visit at t0
  
  return(Vt)
}

#-----------------------------------------------------------------------------#
# dataPrepSims.fct:                                                           #       
#   Data preparation for each subject, in the simulations, with event times   #
#   imputed at either the middle or end of interval between last physician    #
#   visit before and first visit after the true event time                    #
#-----------------------------------------------------------------------------#

## Arguments:
  # dat: data for one subject (generated data with true event times). 
  # Vt: visit times of the subject.
  # full.TVvariables: matrix of full history of time-varying variables for the   
  #                   subject, beyond true event time, with exposure in 1st 
  #                   column and additional time-varying covariate, if 
  #                   applicable, in the 2nd column.
  # impute: event time imputation at 'mid' or 'end' of interval.
  # addTVcov: simulation setting for additional time-varying covariate

dataPrepSims.fct <- function(dat, Vt, full.TVvariables, impute, addTVcov){
    
  # Remove observations beyond the last physician visit
  dat.tmp <- dat[dat$stop <= max(Vt),]
  
  # Add columns indicating if suject is case (0 or 1), 2nd last for cases (NA  
  # for controls), last visit, and follow-up
  dat.tmp$case <- rep(dat.tmp$event[nrow(dat.tmp)], times=nrow(dat.tmp))  
  dat.tmp$vis.2ndLast <- ifelse(dat.tmp$event[nrow(dat.tmp)] == 1,
                                max(Vt[Vt < dat.tmp$stop[nrow(dat.tmp)]]), NA)
  dat.tmp$vis.Last <- min(Vt[Vt >= dat.tmp$stop[nrow(dat.tmp)]])
  dat.tmp$fup <- rep(dat.tmp$stop[nrow(dat.tmp)], times=nrow(dat.tmp))  

  # If subjects is a control or a case with event after last visit: follow-up
  # is censored at last visit
  if (dat.tmp$event[nrow(dat.tmp)] == 0){
    return(dat.tmp)  
    
  # If subjects is a case
  } else {
    
    # Determine the imputed event time
    if (impute == 'mid'){
      imp.et <- ceiling((max(Vt[Vt < dat.tmp$stop[nrow(dat.tmp)]]) + 
                         min(Vt[Vt >= dat.tmp$stop[nrow(dat.tmp)]]))/2)
    } else if (impute == 'end'){
      imp.et <- min(Vt[Vt >= dat.tmp$stop[nrow(dat.tmp)]])        
    }

    # If imputed time < true time: delete observations after imputed event time  
    if (imp.et <= tail(dat.tmp$stop, n=1)){
      dat.imput <- dat.tmp[dat.tmp$stop <= imp.et,]
      dat.imput$event[nrow(dat.imput)] <- 1
          
    # If imputed time > true time: add observations after true event time     
    } else {
      nlines <- imp.et - tail(dat.tmp$stop, n=1) + 1
      data.toadd <- cbind(rep(dat.tmp$id[1], times=nlines),          # id
                          c(rep(0, times=nlines-1), 1),              # event
                          (nrow(dat.tmp)-1):(imp.et-1),              # start
                          nrow(dat.tmp):imp.et,                      # stop
                          full.TVvariables[nrow(dat.tmp):imp.et, 1], # exposure
                          rep(dat.tmp$sex[1], times=nlines),         # sex
                          rep(dat.tmp$age.m40[1], times=nlines))     # age

      # Add additional time-varying covariate, if applicable
      if (addTVcov == 1){
        data.toadd <- cbind(data.toadd, full.TVvariables[nrow(dat.tmp):imp.et, 2])          
      }
      
      data.toadd <- cbind(data.toadd, 
                          rep(dat.tmp$case[1], times=nlines),        # case
                          rep(dat.tmp$vis.2ndLast[1], times=nlines), # vis.2ndLast
                          rep(dat.tmp$vis.Last[1], times=nlines),    # vis.Last
                          rep(dat.tmp$fup[1], times=nlines))         # fup

      colnames(data.toadd) <- colnames(dat.tmp)
      if ((nrow(dat.tmp)-1) > 0){
        dat.imput <- rbind(dat.tmp[1:(nrow(dat.tmp)-1),], data.toadd)
      } else {
        dat.imput <- as.data.frame(data.toadd)
      }
    }
    
    # Update fup after imputation of event time
    dat.imput$fup <- rep(dat.imput$stop[nrow(dat.imput)], times=nrow(dat.imput))
    return(dat.imput)
  }
} 

#-----------------------------------------------------------------------------#
# seqExpoCases.fct:                                                           #
#   Returns the sequence of exposure values at the time of event and the 40   #
#   preceding values for a subject.                                           #
#   Used for calculating the correlation of exposures helping to chose the    #
#   extrapolation function.                                                   #
#-----------------------------------------------------------------------------#

## Arguments:
  # dat: data for one subject. 

seqExpoCases.fct <- function(dat){   
  
  # Controls are ignored
  if (dat$event[nrow(dat)]==0){
    return(NULL)
  
  } else if (nrow(dat) > 40){  
    return(dat[nrow(dat):(nrow(dat)-40), 5])  # 5th column is exposure
  } else {
    return(rev(dat[, 5]))  # 5th column is exposure
  }
} 

#-----------------------------------------------------------------------------#
# SIMEX_dataPrepSims.fct:                                                     #
#   Data preparation for SIMEX method in the simulations: adjust follow-up    #    
#   according to changes in visits                                            #
#-----------------------------------------------------------------------------#

## Arguments:
  # dat: a dataset, with event times imputed at mid-point between last 2 visits
  #      for cases, or censored at last visit for controls.
  # back.units: number of units by which 2nd last visit is moved backward.  

SIMEX_dataPrepSims.fct <- function(dat, back.units){

  # Amount of time to cut at the end of follow-up of the subject:
    # case: event time imputed at mid-point between last (unchanged) and 2nd last 
    #       (moved ahead by 'back.units') visits; so, cut back.units/2.
    # control: censored at last visit (moved ahead by 'back.units'); so, cut
    #          back.units.
  fup.cut <- ifelse(dat$case==1, floor(back.units/2), back.units)
  
  # For cases, 2nd last visit cannot be moved earlier than t0. Apply the 
  # restriction to t0 if needed.
  fup.cut[dat$case==1 & back.units > dat$vis.2ndLast] <- 
    floor(dat$vis.2ndLast/2)[dat$case==1 & back.units > dat$vis.2ndLast]  

  # For controls, we cannot cut follow-up by more than dat$fup
  fup.cut[dat$case==0 & fup.cut > dat$fup] <- dat$fup[dat$case==0 & fup.cut > dat$fup]  
  
  # Cut follow-up
  dat$fupSIMEX <- dat$fup - fup.cut
  dat2 <- dat[dat$stop <= dat$fupSIMEX,]       
  dat2$event <- ifelse(dat2$case==1 & dat2$stop==dat2$fupSIMEX, 1, 0)
  
  dat2$vis.2ndLastSIMEX <- dat2$vis.2ndLast - back.units
  dat2$vis.2ndLastSIMEX[dat2$vis.2ndLastSIMEX < 0] <- 0  
  
  return(dat2)
}

#-----------------------------------------------------------------------------#
# simulation.fct:                                                             #
#   Simulation function to generate data, estimate models, and save results   #
#-----------------------------------------------------------------------------#

## Arguments:
  # n.sims: number of simulated samples per scenario.
  # N: sample size.
  # fup: maximum follow-up in weeks.
  # incidence.rate: proportion of events during follow-up.
  # expoGen: generation of time-varying exposure, with possible options:
  #          "Duration" = Duration-based exposure generation
  #          "ChangeStatus" = Change-of-status exposure generation
  # beta.curUse: true parameter value for "current use" exposure in true model 
  #              for data generation.
  # beta.cumDur12w: true parameter value for "cumulative duration in past 12   
  #                 weeks" exposure in true model for data generation.
  # beta.sex: true parameter value for sex in true model for data generation.
  # beta.age: true parameter value for age in true model for data generation.
  # addTVcov: binary indicator (1 or 0) of whether an additional time-varying 
  #           covariate.
  # beta.TVcov: true parameter value for additional time-varying covariate in  
  #             true model for data generation.
  # P: average number of visits generated during maximum follow-up time for 
  #    each subject.
  # move2ndLast.SIMEX: vector of values of the number of units by which the 2nd
  #                    last visits (cases) or last visits (controls) are moved
  #                    backward in step 1 of SIMEX-like proposed method.
  # seed: seed for the random number generator (to allow reproducibility of 
  #       results).
  
simulation.fct <- function(n.sims, N, fup, incidence.rate, expoGen,  
                           beta.curUse, beta.cumDur12w, beta.sex, beta.age, 
                           addTVcov, beta.TVcov = NULL, 
                           P, move2ndLast.SIMEX, seed){
  
  ### Additional parameters for exposure generation ###
  
  if (expoGen=="Duration"){
    # Initial probability of being exposed
    Pt0expo <- 0.5 
    # Probability of being exposed for subsequent periods
    Ptexpo <- 0.5 
    # Lower/upper bound of meta-parameter for subject's mean duration of use episodes 
    DurON.lb <- 4      
    DurON.ub <- 12    
    # Lower/upper bound of meta-parameter for subject's mean duration of non-use episodes
    DurOFF.lb <- 4     
    DurOFF.ub <- 16  
    
  } else if (expoGen=="ChangeStatus"){ 
    # Initial probability of being exposed
    Pt0expo <- 0.5     
    # Lower/upper bound of meta-parameter for subject's probability continuing drug use
    Pcont.lb <- 0.4     
    Pcont.ub <- 0.9 
    # Lower/Upper bound of meta-parameter for subject's probability starting new drug use period 
    Prestart.lb <- 0   
    Prestart.ub <- 0.3  
  }

  #######################
  # Objects to be saved #
  #######################
  
  ## Number of events
  n.events <- rep(NA, n.sims)
  n.events.curUse.mid <- rep(NA, n.sims)
  n.events.cumDur12w.mid <- rep(NA, n.sims)
  
  ## Correlation of exposures at different times for choosing extrapolation function  
  cor.expo.cases.curUse <- matrix(NA, nrow=n.sims, ncol=41)
  cor.expo.cases.cumDur12w <- matrix(NA, nrow=n.sims, ncol=41) 
  
  ## Coefficients and SEs for conventional Cox PH models 
  cox.curUse.coef.truth <- array(NA, dim=c((3+addTVcov),2,n.sims))
  cox.cumDur12w.coef.truth <- array(NA, dim=c((3+addTVcov),2,n.sims))
  cox.curUse.coef.mid <- array(NA, dim=c((3+addTVcov),2,n.sims))
  cox.cumDur12w.coef.mid <- array(NA, dim=c((3+addTVcov),2,n.sims))
  cox.curUse.coef.end <- array(NA, dim=c((3+addTVcov),2,n.sims))
  cox.cumDur12w.coef.end <- array(NA, dim=c((3+addTVcov),2,n.sims))
  
  if (addTVcov == 0){
    covariates.names <- c('sex','age.m40')
  } else if (addTVcov == 1){
    covariates.names <- c('sex','age.m40','curTVcov')
  }
  rownames(cox.curUse.coef.truth) <- c('curUse',covariates.names)
  rownames(cox.curUse.coef.mid) <- c('curUse',covariates.names)
  rownames(cox.curUse.coef.end) <- c('curUse',covariates.names)
  rownames(cox.cumDur12w.coef.truth) <- c('cumDur12w',covariates.names)
  rownames(cox.cumDur12w.coef.mid) <- c('cumDur12w',covariates.names)
  rownames(cox.cumDur12w.coef.end) <- c('cumDur12w',covariates.names) 
  colnames(cox.curUse.coef.truth) <- colnames(cox.cumDur12w.coef.truth) <- 
    colnames(cox.curUse.coef.mid) <- colnames(cox.cumDur12w.coef.mid) <- 
    colnames(cox.curUse.coef.end) <- colnames(cox.cumDur12w.coef.end) <- 
    c('coef','se(coef)')
  
  ## Mean distance between last 2 vistis for cases; 1 col for each SIMEX iteration
  mean.dist2lastV.curUse <- matrix(rep(NA, times=n.sims*(1+length(move2ndLast.SIMEX))), 
                            nrow=n.sims, ncol=1+length(move2ndLast.SIMEX), byrow=F)
  mean.dist2lastV.cumDur12w <- matrix(rep(NA, times=n.sims*(1+length(move2ndLast.SIMEX))), 
                               nrow=n.sims, ncol=1+length(move2ndLast.SIMEX), byrow=F)
  
  ## Coefficients for Cox PH models at each SIMEX iteration
  cox.curUse.move2ndLast.SIMEX <- array(NA, dim=c((3+addTVcov), 
                                                  length(move2ndLast.SIMEX), n.sims))
  cox.cumDur12w.move2ndLast.SIMEX <- array(NA, dim=c((3+addTVcov), 
                                                     length(move2ndLast.SIMEX), n.sims))
  rownames(cox.curUse.move2ndLast.SIMEX) <- c('curUse',covariates.names)
  rownames(cox.cumDur12w.move2ndLast.SIMEX) <- c('cumDur12w',covariates.names)
  colnames(cox.curUse.move2ndLast.SIMEX) <- 
    colnames(cox.cumDur12w.move2ndLast.SIMEX) <- 
    paste0('Ext', move2ndLast.SIMEX)
  
  ## Coefficients from SIMEX extrapolation step
  cox.curUse.expo.SIMEXlin.15it.extrap <- rep(NA, times=n.sims)
  cox.cumDur12w.expo.SIMEXlin.15it.extrap <- rep(NA, times=n.sims)
  cox.curUse.expo.SIMEXlin.25it.extrap <- rep(NA, times=n.sims)
  cox.cumDur12w.expo.SIMEXlin.25it.extrap <- rep(NA, times=n.sims)
  
  cox.curUse.expo.SIMEXlog2.15it.extrap <- rep(NA, times=n.sims)
  cox.cumDur12w.expo.SIMEXlog2.15it.extrap <- rep(NA, times=n.sims)
  cox.curUse.expo.SIMEXlog2.25it.extrap <- rep(NA, times=n.sims)
  cox.cumDur12w.expo.SIMEXlog2.25it.extrap <- rep(NA, times=n.sims)
  
  cox.curUse.expo.SIMEXquad1df.15it.extrap <- rep(NA, times=n.sims)
  cox.cumDur12w.expo.SIMEXquad1df.15it.extrap <- rep(NA, times=n.sims)
  cox.curUse.expo.SIMEXquad1df.25it.extrap <- rep(NA, times=n.sims)
  cox.cumDur12w.expo.SIMEXquad1df.25it.extrap <- rep(NA, times=n.sims)
  
  cox.cumDur12w.expo.SIMEXlin.11it.extrap <- rep(NA, times=n.sims)

  ###################
  # Simulation loop #
  ###################
  
  library(survival)
  library(PermAlgo)
  library(rlist)
  
  options(warn=1)
  
  set.seed(seed)
  
  cat("Progression of simulations:", "\n")
  for (s in 1:n.sims){
    cat(" s=", s, "\n")
      
    ###################
    # Data generation #
    ###################
    
    ### Generate time-varying exposure ###
  
    if (expoGen == "Duration"){
      curUse <- NULL
      for (i in 1:N){
        # Meta-parameters for subject's mean duration of use/non-use episodes 
        MeanDurON <- runif(1, DurON.lb, DurON.ub)
        MeanDurOFF <- runif(1, DurOFF.lb, DurOFF.ub)
        curUsei <- NULL  
        step <- 1
        while (step < fup){
          curUseit <- rbinom(1, 1, Pt0expo)
          curUsei <- c(curUsei, 
                       rep(curUseit, 
                           times=ifelse(curUseit==1, 
                                        round(runif(1, 0.5*MeanDurON, 2*MeanDurON)), 
                                        round(runif(1, 0.5*MeanDurOFF, 2*MeanDurOFF)))))
          step <- length(curUsei)
        }
        curUse <- rbind(curUse, curUsei[1:fup])
      }
      cumUse <- t(apply(curUse, 1, cumsum))
      cumDur12w <- cbind(cumUse[,1:12], t(apply(cumUse, 1, diff, 12))) 
  
    } else if (expoGen == "ChangeStatus"){
      curUse <- matrix(c(rbinom(N,1,Pt0expo), rep(NA,(N*fup-N))), nrow=N, ncol=fup)
      for (i in 1:N){
        # Meta-parameters for subject's probability to continue, or restart, drug use
        Pcont <- runif(1, Pcont.lb, Pcont.ub)
        Prestart <- runif(1, Prestart.lb, Prestart.ub)
        for (t in 2:fup){
          curUse[i,t] <- ifelse(curUse[i,t-1]==1, 
                                rbinom(1,1,Pcont), 
                                rbinom(1,1,Prestart))
        }
      }
      cumUse <- t(apply(curUse, 1, cumsum))
      cumDur12w <- cbind(cumUse[,1:12], t(apply(cumUse, 1, diff, 12))) 
    }

    ### Generate covariates ###
      
    sex <- rbinom(N, 1, 0.7)   # 70% women
    age <- round(runif(N, 40, 80)) - 40
  
    ###  Generate time-varying covariate (same generation process as for exposure) ###
      
    if (addTVcov == 1 & (expoGen == "Duration")){
      curTVcov <- NULL
      for (i in 1:N){
        # Meta-parameters for subject's mean duration of use/non-use episodes 
        MeanDurON <- runif(1, DurON.lb, DurON.ub)
        MeanDurOFF <- runif(1, DurOFF.lb, DurOFF.ub)
        curTVcovi <- NULL  
        step <- 1
        while (step < fup){
          curTVcovit <-  rbinom(1, 1, Pt0expo)
          curTVcovi <-  c(curTVcovi, 
                          rep(curTVcovit, 
                              times=ifelse(curTVcovit==1, 
                                           round(runif(1, 0.5*MeanDurON, 2*MeanDurON)), 
                                           round(runif(1, 0.5*MeanDurOFF, 2*MeanDurOFF)))))
          step <- length(curTVcovi)
        }
        curTVcov <- rbind(curTVcov, curTVcovi[1:fup])
      }  
    }        
  
    ### Generate event times with permutational algorithm ###
    
    # All censoring times at end of study
    censortimes <- rep(fup, N)  
    
    # Event times followed exponential distribution 
    eventtimes <- ceiling(rexp(N, -log(1-incidence.rate)/fup))  
    
    if (addTVcov == 0){
      covariates.betas <- c(beta.sex, beta.age)
      covariates.Xmat <- apply(t(cbind(sex, age)), 1, rep, each=fup)
    
    } else if (addTVcov == 1){
      covariates.betas <- c(beta.sex, beta.age, beta.TVcov) 
      covariates.Xmat <- cbind(apply(t(cbind(sex, age)), 1, rep, each=fup), c(t(curTVcov)))
    }
  
    # Use of permutational algorithm  
    dat.curUse <- permalgorithm(numSubjects=N, maxTime=fup, 
                                Xmat=cbind(c(t(curUse)), covariates.Xmat),
                                eventRandom=eventtimes, 
                                censorRandom=censortimes,
                                betas=c(beta.curUse, covariates.betas))
    dat.curUse <- dat.curUse[,c(1,2,4:(8+addTVcov))]
    colnames(dat.curUse) <- c('id','event','start','stop','curUse',covariates.names)
    
    dat.cumDur12w <- permalgorithm(numSubjects=N, maxTime=fup, 
                                   Xmat=cbind(c(t(cumDur12w)), covariates.Xmat),
                                   eventRandom=eventtimes, censorRandom=censortimes, 
                                   betas=c(beta.cumDur12w, covariates.betas)) 
    dat.cumDur12w <- dat.cumDur12w[,c(1,2,4:(8+addTVcov))]
    colnames(dat.cumDur12w) <- c('id','event','start','stop','cumDur12w',covariates.names)
    
    ### Generate physician visits ###
      
    visits <- sapply(1:N, function(i) genVisits.fct(fup=fup, P=P))
  
    ### Create datasets with event times imputed at mid or end of interval ###
    ### between last 2 visits                                              ###
  
    # Split datasets by subjects    
    dat.curUse.i <- split(dat.curUse, dat.curUse$id)
    dat.cumDur12w.i <- split(dat.cumDur12w, dat.cumDur12w$id)
  
    # Preparation of full history time-varying variables 
    if (addTVcov == 0){
      full.TVvariables.curUse <- lapply(1:N, function(i) 
                                        matrix(curUse[i,], ncol=1))
      full.TVvariables.cumDur12w <- lapply(1:N, function(i) 
                                            matrix(cumDur12w[i,], ncol=1))
    } else if (addTVcov == 1){
      full.TVvariables.curUse <- lapply(1:N, function(i) 
                                        cbind(curUse[i,], curTVcov[i,]))
      full.TVvariables.cumDur12w <- lapply(1:N, function(i) 
                                            cbind(cumDur12w[i,], curTVcov[i,]))
    }

    dat.curUse.mid <- do.call("rbind", lapply(1:N, function(i) 
                        dataPrepSims.fct(dat = dat.curUse.i[[i]], 
                                          Vt = visits[[i]], 
                                          full.TVvariables = full.TVvariables.curUse[[i]],
                                          impute = 'mid', 
                                          addTVcov = addTVcov))) 
  
    dat.curUse.end <- do.call("rbind", lapply(1:N, function(i) 
                        dataPrepSims.fct(dat = dat.curUse.i[[i]], 
                                          Vt = visits[[i]],
                                          full.TVvariables = full.TVvariables.curUse[[i]],
                                          impute = 'end', 
                                          addTVcov = addTVcov))) 
     
    dat.cumDur12w.mid <- do.call("rbind", lapply(1:N, function(i) 
                           dataPrepSims.fct(dat = dat.cumDur12w.i[[i]], 
                                             Vt = visits[[i]], 
                                             full.TVvariables = full.TVvariables.cumDur12w[[i]],
                                             impute = 'mid',
                                             addTVcov = addTVcov)))
    
    dat.cumDur12w.end <- do.call("rbind", lapply(1:N, function(i) 
                           dataPrepSims.fct(dat = dat.cumDur12w.i[[i]], 
                                             Vt = visits[[i]], 
                                             full.TVvariables = full.TVvariables.cumDur12w[[i]],
                                             impute = 'end', 
                                             addTVcov = addTVcov)))  
     
    n.events[s] <- sum(dat.curUse$event)
    n.events.curUse.mid[s] <- sum(dat.curUse.mid$event)
    n.events.cumDur12w.mid[s] <- sum(dat.cumDur12w.mid$event)
  
    ### Distance between last 2 visits for case (before and after true event time) ###
  
    mean.dist2lastV.curUse[s,1] <- 
      mean(dat.curUse.mid[dat.curUse.mid$event==1,'vis.Last'] - 
           dat.curUse.mid[dat.curUse.mid$event==1,'vis.2ndLast'])
    mean.dist2lastV.cumDur12w[s,1] <- 
      mean(dat.cumDur12w.mid[dat.cumDur12w.mid$event==1,'vis.Last'] - 
      dat.cumDur12w.mid[dat.cumDur12w.mid$event==1,'vis.2ndLast'])
  
    ##############################################################################
    # Correlation of exposures at different times for choosing the extrapolation #
    # function                                                                   #
    ##############################################################################
  
    # Split datasets by subjects    
    dat.curUse.mid.i <- split(dat.curUse.mid, dat.curUse.mid$id)
    dat.cumDur12w.mid.i <- split(dat.cumDur12w.mid, dat.cumDur12w.mid$id)
    
    tmp.curUse <- lapply(1:N, function(i) 
                            seqExpoCases.fct(dat = dat.curUse.mid.i[[i]])) 
    seq.expo.cases.curUse <- tmp.curUse[-which(sapply(tmp.curUse, is.null))]
    tmp.cumDur12w <- lapply(1:N, function(i) 
                            seqExpoCases.fct(dat = dat.cumDur12w.mid.i[[i]])) 
    seq.expo.cases.cumDur12w <- tmp.cumDur12w[-which(sapply(tmp.cumDur12w, is.null))]   
    
    for (j in 1:41){
          
      x.expo.cases.curUse <- NULL
      y.expo.cases.curUse <- NULL
      for (i in 1:n.events.curUse.mid[s]){
        if (length(seq.expo.cases.curUse[[i]]) >= j ){
          x.expo.cases.curUse <- 
            c(x.expo.cases.curUse, seq.expo.cases.curUse[[i]][1])
          y.expo.cases.curUse <- 
            c(y.expo.cases.curUse, seq.expo.cases.curUse[[i]][j])
        }
      }
      cor.expo.cases.curUse[s,j] <- 
        cor(x.expo.cases.curUse, y.expo.cases.curUse)
      
      x.expo.cases.cumDur12w <- NULL
      y.expo.cases.cumDur12w <- NULL
      for (i in 1:n.events.cumDur12w.mid[s]){
        if (length(seq.expo.cases.cumDur12w[[i]]) >= j ){
          x.expo.cases.cumDur12w <- 
            c(x.expo.cases.cumDur12w, seq.expo.cases.cumDur12w[[i]][1])
          y.expo.cases.cumDur12w <- 
            c(y.expo.cases.cumDur12w, seq.expo.cases.cumDur12w[[i]][j])
        }
      }
      cor.expo.cases.cumDur12w[s,j] <- 
        cor(x.expo.cases.cumDur12w, y.expo.cases.cumDur12w)
    }
  
    #########################################
    # Estimation of conventional Cox models #
    #########################################
      
    if (addTVcov == 0){ 
      covariates <- 'sex + age.m40'
    } else if (addTVcov == 1){
      covariates <- 'sex + age.m40 + curTVcov'
    }
    
    ### Oracle ###
      
    cox.curUse.truth <- coxph(
      as.formula(paste0('Surv(start, stop, event) ~ curUse +', covariates)), 
      data=dat.curUse)
    cox.curUse.coef.truth[,,s] <- summary(cox.curUse.truth)$coef[,c(1,3)]
    
    cox.cumDur12w.truth <- coxph(
      as.formula(paste0('Surv(start, stop, event) ~ cumDur12w +', covariates)), 
      data=dat.cumDur12w)  
    cox.cumDur12w.coef.truth[,,s] <- summary(cox.cumDur12w.truth)$coef[,c(1,3)]  
    
    ### Event times imputed at middle interval ###
    
    cox.curUse.coef.mid[,,s] <- summary(coxph(
      as.formula(paste0('Surv(start, stop, event) ~ curUse +', covariates)), 
      data=dat.curUse.mid))$coef[,c(1,3)]
    cox.cumDur12w.coef.mid[,,s] <- summary(coxph(
      as.formula(paste0('Surv(start, stop, event) ~ cumDur12w +', covariates)), 
      data=dat.cumDur12w.mid))$coef[,c(1,3)] 
    
    ### Event times imputed at end of interval ###
    
    cox.curUse.coef.end[,,s] <- summary(coxph(
      as.formula(paste0('Surv(start, stop, event) ~ curUse +', covariates)), 
      data=dat.curUse.end))$coef[,c(1,3)]
    cox.cumDur12w.coef.end[,,s] <- summary(coxph(
      as.formula(paste0('Surv(start, stop, event) ~ cumDur12w +', covariates)), 
      data=dat.cumDur12w.end))$coef[,c(1,3)]
   
    #############################################
    # SIMEX-like proposed method: Steps 1 and 2 #
    #############################################
     
    for (j in 1:length(move2ndLast.SIMEX)){
      #cat("  j=", j, "\n")
      
      # Step 1: Gradually increase mean time between visits for cases by moving 
      #         backward 2nd last visits. Move backward last visits for controls.
      dat.curUse.mid.modif <- SIMEX_dataPrepSims.fct(dat = dat.curUse.mid, 
                                                 back.units = move2ndLast.SIMEX[j])
      dat.cumDur12w.mid.modif <- SIMEX_dataPrepSims.fct(dat = dat.cumDur12w.mid, 
                                                    back.units = move2ndLast.SIMEX[j])
        
      # Step 2: Refit Cox Ph model on the modified dataset
      cox.curUse.move2ndLast.SIMEX[,j,s] <- coxph(
        as.formula(paste0('Surv(start, stop, event) ~ curUse +', covariates)), 
        data=dat.curUse.mid.modif)$coef    
      cox.cumDur12w.move2ndLast.SIMEX[,j,s] <- coxph(
        as.formula(paste0('Surv(start, stop, event) ~ cumDur12w +', covariates)), 
        data=dat.cumDur12w.mid.modif)$coef   
  
      # Modified mean time between visits for cases 
      mean.dist2lastV.curUse[s,1+j] <- mean(
        dat.curUse.mid.modif[dat.curUse.mid.modif$event==1,'vis.Last'] - 
        dat.curUse.mid.modif[dat.curUse.mid.modif$event==1,'vis.2ndLastSIMEX'])
      mean.dist2lastV.cumDur12w[s,1+j] <- mean(
        dat.cumDur12w.mid.modif[dat.cumDur12w.mid.modif$event==1,'vis.Last'] - 
        dat.cumDur12w.mid.modif[dat.cumDur12w.mid.modif$event==1,'vis.2ndLastSIMEX']) 
    }
    rm(dat.curUse, dat.cumDur12w, dat.curUse.mid, dat.cumDur12w.mid, 
       dat.curUse.end, dat.cumDur12w.end, visits)
      
    #############################################
    # SIMEX-like proposed method: Steps 3 and 4 #
    #############################################
  
    ### Linear extrapolation function ###
    
    # Mean difference between the (modified) times of the two relevant visits for cases
    x.lin.curUse <- mean.dist2lastV.curUse[s,]
    x.lin.cumDur12w <- mean.dist2lastV.cumDur12w[s,]
    
    # Coefficients from Cox PH models for modified events and censoring times 
    cox.curUse.expo.SIMEX <- c(cox.curUse.coef.mid[1,1,s], 
                               cox.curUse.move2ndLast.SIMEX[1,,s])
    cox.cumDur12w.expo.SIMEX <- c(cox.cumDur12w.coef.mid[1,1,s], 
                                  cox.cumDur12w.move2ndLast.SIMEX[1,,s])
    
    # Step 3: Linear regression
      # 15 iterations
    curUse.expo.SIMEXlin.15it <- lm(cox.curUse.expo.SIMEX[1:16] ~ x.lin.curUse[1:16])
    cumDur12w.expo.SIMEXlin.15it <- lm(cox.cumDur12w.expo.SIMEX[1:16] ~ x.lin.cumDur12w[1:16])
      # 25 iterations  
    curUse.expo.SIMEXlin.25it <- lm(cox.curUse.expo.SIMEX ~ x.lin.curUse)
    cumDur12w.expo.SIMEXlin.25it <- lm(cox.cumDur12w.expo.SIMEX ~ x.lin.cumDur12w)
    
    # Step 4: Extrapolation
      # 15 iterations
    cox.curUse.expo.SIMEXlin.15it.extrap[s] <- curUse.expo.SIMEXlin.15it$coef[1] + 
                                                1*curUse.expo.SIMEXlin.15it$coef[2]  
    cox.cumDur12w.expo.SIMEXlin.15it.extrap[s] <- cumDur12w.expo.SIMEXlin.15it$coef[1] + 
                                                  1*cumDur12w.expo.SIMEXlin.15it$coef[2] 
      # 25 iterations
    cox.curUse.expo.SIMEXlin.25it.extrap[s] <- curUse.expo.SIMEXlin.25it$coef[1] + 
                                                1*curUse.expo.SIMEXlin.25it$coef[2]  
    cox.cumDur12w.expo.SIMEXlin.25it.extrap[s] <- cumDur12w.expo.SIMEXlin.25it$coef[1] + 
                                                1*cumDur12w.expo.SIMEXlin.25it$coef[2] 
    
    ### Alternative extrapolation function (relevant for some simulation ###
    ### scenarios): Log in base 2                                        ###
    
    x.log2.curUse <- log2(x.lin.curUse)
    x.log2.cumDur12w <- log2(x.lin.cumDur12w)  
  
    curUse.expo.SIMEXlog2.15it <- lm(cox.curUse.expo.SIMEX[1:16] ~ x.log2.curUse[1:16])
    cumDur12w.expo.SIMEXlog2.15it <- lm(cox.cumDur12w.expo.SIMEX[1:16] ~ x.log2.cumDur12w[1:16])
    curUse.expo.SIMEXlog2.25it <- lm(cox.curUse.expo.SIMEX ~ x.log2.curUse)
    cumDur12w.expo.SIMEXlog2.25it <- lm(cox.cumDur12w.expo.SIMEX ~ x.log2.cumDur12w)
  
    cox.curUse.expo.SIMEXlog2.15it.extrap[s] <- curUse.expo.SIMEXlog2.15it$coef[1] + 
                                                log2(1)*curUse.expo.SIMEXlog2.15it$coef[2] 
    cox.cumDur12w.expo.SIMEXlog2.15it.extrap[s] <- cumDur12w.expo.SIMEXlog2.15it$coef[1] + 
                                                   log2(1)*cumDur12w.expo.SIMEXlog2.15it$coef[2] 
    cox.curUse.expo.SIMEXlog2.25it.extrap[s] <- curUse.expo.SIMEXlog2.25it$coef[1] + 
                                                log2(1)*curUse.expo.SIMEXlog2.25it$coef[2] 
    cox.cumDur12w.expo.SIMEXlog2.25it.extrap[s] <- cumDur12w.expo.SIMEXlog2.25it$coef[1] + 
                                                   log2(1)*cumDur12w.expo.SIMEXlog2.25it$coef[2] 
  
    ### Alternative extrapolation function (relevant for some simulation ###
    ### scenarios): 1-df quadratic                                       ###
    
    x.quad1df.curUse <- x.lin.curUse*x.lin.curUse
    x.quad1df.cumDur12w <- x.lin.cumDur12w*x.lin.cumDur12w  
  
    curUse.expo.SIMEXquad1df.15it <- lm(cox.curUse.expo.SIMEX[1:16] ~ x.quad1df.curUse[1:16])
    cumDur12w.expo.SIMEXquad1df.15it <- lm(cox.cumDur12w.expo.SIMEX[1:16] ~ x.quad1df.cumDur12w[1:16])  
    curUse.expo.SIMEXquad1df.25it <- lm(cox.curUse.expo.SIMEX ~ x.quad1df.curUse)
    cumDur12w.expo.SIMEXquad1df.25it <- lm(cox.cumDur12w.expo.SIMEX ~ x.quad1df.cumDur12w)
  
    cox.curUse.expo.SIMEXquad1df.15it.extrap[s] <- curUse.expo.SIMEXquad1df.15it$coef[1] + 
                                                   1*curUse.expo.SIMEXquad1df.15it$coef[2]  
    cox.cumDur12w.expo.SIMEXquad1df.15it.extrap[s] <- cumDur12w.expo.SIMEXquad1df.15it$coef[1] + 
                                                      1*cumDur12w.expo.SIMEXquad1df.15it$coef[2]  
    cox.curUse.expo.SIMEXquad1df.25it.extrap[s] <- curUse.expo.SIMEXquad1df.25it$coef[1] + 
                                                   1*curUse.expo.SIMEXquad1df.25it$coef[2] 
    cox.cumDur12w.expo.SIMEXquad1df.25it.extrap[s] <- cumDur12w.expo.SIMEXquad1df.25it$coef[1] + 
                                                      1*cumDur12w.expo.SIMEXquad1df.25it$coef[2] 
    
    ### Alternative extrapolation function (relevant for some simulation ###
    ### scenarios): Linear with 11 iterations                            ###
    
    cumDur12w.expo.SIMEXlin.11it <- lm(cox.cumDur12w.expo.SIMEX[1:12] ~ x.lin.cumDur12w[1:12])
    
    cox.cumDur12w.expo.SIMEXlin.11it.extrap[s] <- cumDur12w.expo.SIMEXlin.11it$coef[1] + 
                                                  1*cumDur12w.expo.SIMEXlin.11it$coef[2] 
  }

  ################
  # Save results #
  ################
  
  ### Scenario number and description ###

  if (N==3000 & expoGen=="Duration" & 
      beta.curUse==log(3) & beta.cumDur12w==log(2)/6 & 
      addTVcov==0){ sc.no <- 1
  } else if (N==3000 & expoGen=="Duration" & 
      beta.curUse==log(2) & beta.cumDur12w==log(3)/12 & 
      addTVcov==0){ sc.no <- 2
  } else if (N==1500 & expoGen=="Duration" & 
      beta.curUse==log(3) & beta.cumDur12w==log(2)/6 & 
      addTVcov==0){ sc.no <- 3
  } else if (N==1500 & expoGen=="Duration" & 
      beta.curUse==log(2) & beta.cumDur12w==log(3)/12 & 
      addTVcov==0){ sc.no <- 4
  } else if (N==750 & expoGen=="Duration" & 
      beta.curUse==log(3) & beta.cumDur12w==log(2)/6 & 
      addTVcov==0){ sc.no <- 5
  } else if (N==750 & expoGen=="Duration" & 
      beta.curUse==log(2) & beta.cumDur12w==log(3)/12 & 
      addTVcov==0){ sc.no <- 6

  } else if (N==3000 & expoGen=="Duration" & 
      beta.curUse==log(3) & beta.cumDur12w==log(2)/6 & 
      addTVcov==1){ sc.no <- 7
  } else if (N==3000 & expoGen=="Duration" & 
      beta.curUse==log(1) & beta.cumDur12w==log(1) & 
      addTVcov==0){ sc.no <- 8

  } else if (N==3000 & expoGen=="ChangeStatus" & 
      beta.curUse==log(3) & beta.cumDur12w==log(2)/6 & 
      addTVcov==0){ sc.no <- 9
  } else if (N==3000 & expoGen=="ChangeStatus" & 
      beta.curUse==log(2) & beta.cumDur12w==log(3)/12 & 
      addTVcov==0){ sc.no <- 10
  } else if (N==1500 & expoGen=="ChangeStatus" & 
      beta.curUse==log(3) & beta.cumDur12w==log(2)/6 & 
      addTVcov==0){ sc.no <- 11
  } else if (N==1500 & expoGen=="ChangeStatus" & 
      beta.curUse==log(2) & beta.cumDur12w==log(3)/12 & 
      addTVcov==0){ sc.no <- 12
  } else if (N==750 & expoGen=="ChangeStatus" & 
      beta.curUse==log(3) & beta.cumDur12w==log(2)/6 & 
      addTVcov==0){ sc.no <- 13
  } else if (N==750 & expoGen=="ChangeStatus" & 
      beta.curUse==log(2) & beta.cumDur12w==log(3)/12 & 
      addTVcov==0){ sc.no <- 14
  }

  if (addTVcov==0){
    sc.name <- paste0("Sc", sc.no,
                      " - N", N,
                      " expo", expoGen,
                      " BcurUse", gsub("\\.", "", round(beta.curUse, 2)),
                      " BcumDur12w", gsub("\\.", "", round(beta.cumDur12w, 2)),
                      " - ", n.sims, "samples")
    
  } else if (addTVcov==1){
    sc.name <- paste0("Sc", sc.no,
                      " - N", N,
                      " expo", expoGen,
                      " TVcov",
                      " BcurUse", gsub("\\.", "", round(beta.curUse, 2)),
                      " BcumDur12w", gsub("\\.", "", round(beta.cumDur12w, 2)),
                      " - ", n.sims, "samples")
  } 

  ## Parameters to save depending in parameter for exposure generation
  expopara <- NULL
  if (expoGen == "ChangeStatus"){
    expopara <- c('Pt0expo', 'Pcont.lb', 'Pcont.ub', 'Prestart.lb', 'Prestart.ub')
  } else if (expoGen == "Duration"){
    expopara <- c('Pt0expo','Ptexpo','DurON.lb','DurON.ub','DurOFF.lb','DurOFF.ub') 
  } 
  
  ### Save results ###

  save(n.sims, fup, N, incidence.rate,
       expoGen, expopara, 
       list=c(expopara),
       beta.curUse, beta.cumDur12w, beta.sex, beta.age, beta.TVcov,
       P, move2ndLast.SIMEX, 
       sc.name, 
  
       n.events, n.events.curUse.mid, n.events.cumDur12w.mid,
       mean.dist2lastV.curUse, mean.dist2lastV.cumDur12w,
       cor.expo.cases.curUse, cor.expo.cases.cumDur12w,
       
       cox.curUse.coef.truth, cox.cumDur12w.coef.truth,
       cox.curUse.coef.mid, cox.cumDur12w.coef.mid,
       cox.curUse.coef.end, cox.cumDur12w.coef.end,
       
       cox.curUse.move2ndLast.SIMEX, 
       cox.cumDur12w.move2ndLast.SIMEX,
       cox.curUse.expo.SIMEXlin.15it.extrap, 
       cox.cumDur12w.expo.SIMEXlin.15it.extrap,
       cox.curUse.expo.SIMEXlog2.15it.extrap, 
       cox.cumDur12w.expo.SIMEXlog2.15it.extrap,
       cox.curUse.expo.SIMEXquad1df.15it.extrap, 
       cox.cumDur12w.expo.SIMEXquad1df.15it.extrap,
       cox.curUse.expo.SIMEXlin.25it.extrap,
       cox.cumDur12w.expo.SIMEXlin.25it.extrap,
       cox.curUse.expo.SIMEXlog2.25it.extrap, 
       cox.cumDur12w.expo.SIMEXlog2.25it.extrap,
       cox.curUse.expo.SIMEXquad1df.25it.extrap, 
       cox.cumDur12w.expo.SIMEXquad1df.25it.extrap,
       cox.cumDur12w.expo.SIMEXlin.11it.extrap,
       
       file=paste0(sc.name, ".RData"))
}

