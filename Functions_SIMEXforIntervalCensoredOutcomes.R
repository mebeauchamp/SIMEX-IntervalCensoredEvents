###############################################################################
###                                                                         ###
### Functions to implement the proposed methods for a given dataset. For    ###
### details, see the manuscript:                                            ###
### "Adapting SIMEX to correct for bias due to interval-censored outcomes   ###
### in survival analysis with time-varying exposure"                        ###
###                                                                         ###
### Code by Marie-Eve Beauchamp                                             ###
### Last update: March 25, 2022                                             ###
###                                                                         ###
###############################################################################

library(survival)

options(warn=1)

#-----------------------------------------------------------------------------#
# imputeEventMid:                                                             #
#   Imputes event times at the middle of intervals between the last 2         #
#   physician visits for cases, and leaves controls censored at their last    #
#   physician visit.                                                          #
#-----------------------------------------------------------------------------#

## Arguments:
  # data.end: a data frame, with event times imputed at end of interval between 
  #     last 2 physician visits for cases (i.e. subjects with an event), and
  #     controls (i.e. subjects without an event during follow-up) censored at  
  #     their last physician visit. Must have the counting process format with  
  #     each line corresponding to one and only one time unit for a given 
  #     individual.
  # id: name of column in data.end for subject unique identifiers, e.g. "id".
  # time: name of column in data.end indicating the beginning of time  
  #     intervals, e.g. "start".
  # time2: name of column in data.end indicating the end of time intervals,  
  #     e.g. "stop".
  # event: name of column in data.end for the binary status indicator, e.g. 
  #     "event".
  # visit.2ndLast: name of column in data.end for the time of 2nd last 
  #     physician visit for cases (NA for controls). Constant value over all  
  #     rows of a subject.
  # visit.Last: name of column in data.end for the time of last physician visit
  #     for cases and controls. Constant value over all rows of a subject.  

## Details:
  # The function requires that the dataset passed to it has the format: 
  #   i) cases have their event times imputed at the end of interval between   
  #      last 2 physician visits, and                                              
  #   ii) controls are censored at their last physician visit.                
  # If your data settings are different, adapt the functions imputeEventMid and 
  # imputeEventMid.1id accordingly.

imputeEventMid <- function(data.end, id, time, time2, event, visit.2ndLast, 
                           visit.Last){
  
  # Check for errors in arguments for functions imputeEventMid and imputeEventMid.1id 
  if (max(data.end[,time2] - data.end[,time]) > 1){
    stop("Argument data.end must have each line corresponding to one time unit.")
  }
  if (sum(is.na(data.end[data.end[,event] == 1, visit.2ndLast])) > 0 | 
       sum(is.na(data.end[data.end[,event] == 1, visit.Last])) > 0){
    stop("NA for the time of 2nd last or last visit for some cases.  
    They must be known for all cases. If there is no visit before a 
    true event time, you can assign t0 as the 2nd last visit time.")
  }
  
  # Impute event time for each subject
  N <- length(unique(data.end[,id]))
  data.i <- split(data.end, data.end$id)
  data.middle <- do.call("rbind", lapply(1:N, function(i) 
                  imputeEventMid.1id(data.end.1id = data.i[[i]], time = time, 
                                         time2 = time2, event = event, 
                                         visit.2ndLast = visit.2ndLast, 
                                         visit.Last = visit.Last)))   
  return(data.middle)
}

#-----------------------------------------------------------------------------#
# imputeEventMid.1id:                                                         #
#   Internal function performing event imputation for 1 individual.           #
#-----------------------------------------------------------------------------#

## Arguments:
  # data.end.1id: a data frame including data for 1 subject, with event times 
  #     imputed at the end of interval between last 2 physician visits for cases
  #     (i.e. subjects with an event), and controls (i.e. subjects without an
  #     event during follow-up) censored at their last physician visit. Must 
  #     have the counting process format with each line corresponding to one   
  #     and only one time unit for a given individual.
  # time: name of column in data.end.1id indicating the beginning of time  
  #     intervals, e.g. "start".
  # time2: name of column in data.end.1id indicating the end of time intervals,  
  #     e.g. "stop".
  # event: name of column in data.end.1id for the binary status indicator, e.g. 
  #     "event".
  # visit.2ndLast: name of column in data.end.1id for the time of 2nd last 
  #     physician visit for cases (NA for controls). Constant value over all 
  #     rows of a subject.
  # visit.Last: name of column in data.end.1id for the time of last physician  
  #     visit for cases and controls. Constant value over all rows of a subject.

imputeEventMid.1id <- function(data.end.1id, time, time2, event,  
                                     visit.2ndLast, visit.Last){
  
  # If the subject is a control: no changes to data
  if (data.end.1id[nrow(data.end.1id), event] == 0){
    return(data.end.1id)  
    
  # If the subject is a case
  } else {
    
    # Determine the new imputed event time (middle of interval between 2nd last 
      # and last visits)
    imp.et <- ceiling((data.end.1id[nrow(data.end.1id), visit.2ndLast] + 
                         data.end.1id[nrow(data.end.1id), visit.Last]) / 2)
    
    # Delete observations after the new imputed event time 
    if (imp.et > data.end.1id[nrow(data.end.1id), time2]){ 
      stop("The middle of interval is later than the follow-up data available
      for at least one of the cases.")
    }
    
    data.imput <- data.end.1id[data.end.1id[,time2] <= imp.et,]
    data.imput[nrow(data.imput), event] <- 1
    return(data.imput)
  }
} 

#-----------------------------------------------------------------------------#
# SIMEX_steps1_2:                                                             #
#   Applies steps 1 and 2 of the SIMEX-like proposed method (see section 2.2  #
#   of the manuscript for details), i.e.:                                     #
#   - (step 1) to gradually increase the mean time between last 2 visits for  #
#     cases in j=1,... iterations, and                                        # 
#   - (step 2) to refit the Cox PH model of interest using the modified times #
#     of events for each iteration j                                          #
#-----------------------------------------------------------------------------#

## Arguments:
  # data.middle: a data frame, with event times imputed at the middle of  
  #     interval between the last 2 physician visits for cases (i.e. subjects  
  #     with an event), and controls (i.e. subjects without an event during    
  #     follow-up) censored at their last physician visit. Must have the    
  #     counting process format with each line corresponding to one and only  
  #     one time unit for a given individual.
  # id: name of column in data.middle for subject unique identifiers, e.g. "id".
  # time: name of column in data.middle indicating the beginning of time  
  #     intervals, e.g. "start".
  # time2: name of column in data.middle indicating the end of time intervals,  
  #     e.g. "stop".
  # event: name of column in data.middle for the binary status indicator, e.g. 
  #     "event".
  # exposure: name of column in data.middle for the time-varying exposure  
  #     variable to which the SIMEX-like method is applied, e.g. "cumDur12w".
  # covariates: vector of names of covariates, other than exposure, to include
  #     in the models, e.g. c("sex", "age").
  # visit.2ndLast: name of column in data.middle for the time of the 2nd last 
  #     physician visit for cases (NA for controls). Constant value over all 
  #      rows of a subject.
  # visit.Last: name of column in data.middle for the time of the last  
  #     physician visit for cases and controls. Constant value over all rows of
  #     a subject.
  # move2ndLast: vector of the number of time units by which the 2nd last visit
  #     for cases is moved backward at subsequent SIMEX iterations j. Values  
  #     must be greater than 0. E.g., seq(2, 50, by=2).

## Details: 
  # In this implementation of the SIMEX-like method, the frequency of visits at
  # step 1 is decreased by moving backward the time of the 2nd last visit for 
  # cases (algorithm (i) in section A.1 of the Supporting Information).  
  # Algorithm (ii), which implies to delete the 2nd last visit, was not 
  # implemented in this program.

## Value:
  # cox.coef.SIMEXit: coefficients (log HR) of exposure and covariates at each
  #     SIMEX iteration j=0,1,2,...
  # mean.dist2lastV.SIMEXit: mean distance between the last 2 visits for cases
  #     in the (modified) data at each SIMEX iteration j=0,1,2,...
  # arguments: list of arguments passed to the function SIMEX_steps1_2.
  # data.middle: list including the number of rows (n.rows), of individuals  
  #     (n.ids), and of events (n.events) in the argument data.middle. 

SIMEX_steps1_2 <- function(data.middle, id, time, time2, event, exposure, 
                               covariates, visit.2ndLast, visit.Last, 
                               move2ndLast){

  # Arguments passed to the function SIMEX_steps1_2 (for output)
  arguments <- list(data.middle = deparse(substitute(data.middle)), 
                    id = id, time = time, time2 = time2, 
                    event = event, exposure = exposure, 
                    covariates = covariates, 
                    visit.2ndLast = visit.2ndLast, 
                    visit.Last = visit.Last, 
                    move2ndLast = move2ndLast)
  
  # Add a column indicating if subjects is a case
  fup.eachId <- by(data.middle[,event], data.middle[,id], length)
  data.middle$case <- rep(by(data.middle[,event], data.middle[,id], max), 
                          times = fup.eachId)
  
  # Check for errors in arguments for functions SIMEX_steps1_2 and SIMEX_dataModif
  if (max(data.middle[,time2] - data.middle[,time2]) > 1){
    stop("Argument data.middle must have each line corresponding to one time unit.")
  }
  if (sum(is.na(data.middle[,visit.Last])) > 0) {
    stop("Cannot have NA in visit.Last.")
  }
  if (max(by(data.middle[,visit.Last], data.middle[,id], max) -
      by(data.middle[,visit.Last], data.middle[,id], min)) > 0){
    stop("visit.Last must take a constant value across all rows of a subject.")
  }
  data.tmp <- data.middle[data.middle$case == 1,]
  if (sum(is.na(data.tmp[,visit.2ndLast])) > 0) {
    stop("Cannot have NA in visit.2ndLast for cases.")
  }
  if (max(by(data.tmp[,visit.2ndLast], data.tmp[,id], max) -
      by(data.tmp[,visit.2ndLast], data.tmp[,id], min)) > 0){
    stop("visit.2ndLast must take a constant value across all rows of a case.")
  }
  if (min(move2ndLast) <= 0){
    stop("All values in move2ndLast must be greater than 0.")
  }
  
  # Number of rows, subjects, and events in data.middle (for output) 
  n.rows <- nrow(data.middle)
  n.ids <- length(unique(data.middle[,id]))
  n.events <- sum(data.middle[,event])
  
  # Add a column representing total follow-up of each subject in original 
  # (unmodified) data
  data.middle$fupORI <- rep(by(data.middle[,time2], data.middle[,id], max), 
                        times = fup.eachId)

  # Formula for the Cox PH model to be estimated
  formula <- paste0("Surv(", time, ", ", time2, ", ", event, ") ~ ", exposure)
  for (i in 1:length(covariates)){
    formula <- paste0(formula, " + ", covariates[i])  
  }

  # Mean distance between the last 2 visits for cases at SIMEX iterations  
    # j=0,1,...; to be filled later 
  mean.dist2lastV <- rep(NA, times = (length(move2ndLast)+1))
  names(mean.dist2lastV) <- paste0("it.", 0:length(move2ndLast))

  # Mean distance between the last 2 visits in original (unmodified) data (for j=0)
  mean.dist2lastV[1] <- mean(data.middle[data.middle[,event] == 1, visit.Last] -
                                 data.middle[data.middle[,event] == 1,
                                            visit.2ndLast])
  
  # Coefficient from the Cox PH model at iterations j=0,1,...; to be filled later    
  cox.coef.SIMEXit <- matrix(NA, nrow = 3, ncol = (length(move2ndLast)+1))
  rownames(cox.coef.SIMEXit) <- c(exposure, covariates)
  colnames(cox.coef.SIMEXit) <- paste0("it.", 0:length(move2ndLast))

  # Coefficient from the Cox model for the original analysis on unmodified data
    # (for j=0)
  cox.coef.SIMEXit[,1] <- coxph(eval(parse(text = formula)), 
                                data = data.middle)$coef

  for (j in 1:length(move2ndLast)){
    #cat("  j=", j, "\n")

    # SIMEX step 1: Gradually increase mean time between last 2 visits for
      # cases by moving backward the 2nd last visits. Move backward last visits
      # for controls.
    data.modif <- SIMEX_dataModif(data.middle = data.middle, id = id,
                                  time = time, time2 = time2, event = event,
                                  exposure = exposure, covariates = covariates,
                                  visit.2ndLast = visit.2ndLast,
                                  visit.Last = visit.Last,
                                  back.units = move2ndLast[j])

    # SIMEX step 2: Refit the Cox Ph model on the modified data
    cox.coef.SIMEXit[,1+j] <- coxph(eval(parse(text = formula)),
                                    data = data.modif)$coef

    # Modified mean time between the last 2 visits for cases
    mean.dist2lastV[1+j] <- mean(data.modif[data.modif$event == 1, visit.Last] -
                                  data.modif[data.modif$event == 1,
                                             "visit.2ndLastSIMEX"])
  }
   
  return(list(cox.coef.SIMEXit = cox.coef.SIMEXit,
              mean.dist2lastV.SIMEXit = mean.dist2lastV,
              arguments = arguments,
              data.middle = list(n.rows = n.rows, n.ids = n.ids, 
                                 n.events = n.events)))
}

#-----------------------------------------------------------------------------#
# SIMEX_dataModif:                                                            #
#   Internal function for data modification at SIMEX step 1 to adjust         #
#   follow-up of each subject according to changes in visit times.            #
#-----------------------------------------------------------------------------#

## Arguments:
  # data.middle: a data frame, with event times imputed at the middle of  
  #     interval between the last 2 physician visits for cases (i.e. subjects  
  #     with an event), and controls (i.e. subjects without an event during 
  #     follow-up) censored at their last physician visit. Must have the
  #     counting process format with each line corresponding to one and only 
  #     one time unit for a given individual.
  # id: name of column in data.middle for subject unique identifiers, e.g. "id".
  # time: name of column in data.middle indicating the beginning of time  
  #     intervals, e.g. "start".
  # time2: name of column in data.middle indicating the end of time intervals,  
  #     e.g. "stop".
  # event: name of column in data.middle for the binary status indicator, e.g. 
  #     "event".
  # exposure: name of column in data.middle for the time-varying exposure  
  #     variable to which the SIMEX-like method is applied, e.g. "cumDur12w".
  # covariates: vector of names of covariates, other than exposure, to include
  #     in the models, e.g. c("sex", "age").
  # visit.2ndLast: name of column in data.middle for the time of the 2nd last 
  #     physician visit for cases (NA for controls). Constant value over all 
  #      rows of a subject.
  # visit.Last: name of column in data.middle for the time of last physician 
  #     visit for cases and controls. Constant value over all rows of a subject.
  # back.units: the number of time units by which the 2nd last visit for cases
  #     is moved backward (one value).

## Details: 
  # In this implementation of the SIMEX-like method, the frequency of visits at
  # step 1 is decreased by moving backward the time of the 2nd last visit
  # (algorithm (i) in section A.1 of the Supporting Information). Algorithm 
  # (ii), which implies to delete the 2nd last visit, was not implemented in this 
  # program.

SIMEX_dataModif <- function(data.middle, id, time, time2, event, exposure, 
                               covariates, visit.2ndLast, visit.Last, 
                               back.units){

  # Check for errors in arguments
  if (length(back.units) > 1){ 
    stop("Argument back.units must take only 1 value.")
  }
  
  # Amount of time to cut at the end of follow-up of each subject:
    # - case: event time imputed at mid-point between last (unchanged) and 2nd  
    #         last (moved ahead by back.units) visits; therefore, cut 
    #         back.units/2.
    # - control: censored at last visit (moved ahead by back.units); therefore,
    #         cut back.units.
  fup.cut <- ifelse(data.middle$case == 1, floor(back.units/2), back.units)
  
  # For cases, the 2nd last visit cannot be moved earlier than t0. Apply the 
  # restriction to t0 if needed.
  fup.cut[data.middle$case == 1 & back.units > data.middle[,visit.2ndLast]] <- 
    floor(data.middle[,visit.2ndLast]/2)[data.middle$case==1 & 
                                        back.units > data.middle[,visit.2ndLast]]  
  
  # For controls, follow-up cannot be cut by more than data.middle$fupORI
  fup.cut[data.middle$case == 0 & fup.cut > data.middle$fupORI] <- 
    data.middle$fupORI[data.middle$case == 0 & fup.cut > data.middle$fupORI]  
  
  # New total follow-up for each subject in the modified data
  data.middle$fupSIMEX <- data.middle$fupORI - fup.cut  
  
  # Cut follow-up and add event to the new last line of cases
  data2 <- data.middle[data.middle[,time2] <= data.middle$fupSIMEX,]       
  data2[,event] <- ifelse(data2$case == 1 & data2[,time2] == data2$fupSIMEX, 1, 0)
  
  # Adjust time of the 2nd last visit
  data2$visit.2ndLastSIMEX <- data2[,visit.2ndLast] - back.units
  data2$visit.2ndLastSIMEX[data2$visit.2ndLastSIMEX < 0] <- 0  
  
  return(data2)
}

#-----------------------------------------------------------------------------#
# SIMEX_steps3_4:                                                             #
#   Applies steps 3 and 4 of the SIMEX-like proposed method (see section 2.2  #
#   of the manuscript for details), i.e.:                                     #
#   - (step 3) to regress the k+1 point estimates from step 2, for SIMEX      #
#     iterations j=0,...,k, and                                               # 
#   - (step 4) to extrapolate to a hypothetical ideal scenario where the mean #
#     time between last 2 visits for cases is 1                               #
#-----------------------------------------------------------------------------#

## Arguments:
  # SIMEX_steps1_2.out: output from the function SIMEX_steps1_2.
  # k: the selected number of SIMEX iterations.
  # extrapFUN: extrapolating function f(.) used in SIMEX steps 3 and 4. By
  #     default "linear", alternative options are "log", "log2" (log in base 2)
  #     and "1df-quadratic" (1 degree of freedom quadratic).

## Details:
  # The diagnostic plot implemented in the function plot.diag can help to   
  # select the extrapolating function f(.) and the number of SIMEX iterations k.                                                          

## Value:
  # expo.coef.corrected: SIMEX-corrected coefficient (log HR) estimate for
  #     exposure.
  # linear.reg.step3: output from the linear regression at step 3.
  # SIMEX_steps1_2.out: name of output from the function SIMEX_steps1_2 used.
  # n.iterations.k: number of SIMEX iterations used.
  # extrapolation.FUN: extrapolating function used.

SIMEX_steps3_4 <- function(SIMEX_steps1_2.out, k, extrapFUN = "linear"){
  
  cox.coef.SIMEXit <- SIMEX_steps1_2.out$cox.coef.SIMEXit
  mean.dist2lastV.SIMEXit <- SIMEX_steps1_2.out$mean.dist2lastV.SIMEXit
  
  if (extrapFUN == "linear"){  
    
    ### Step 3: Linear regression on k+1 point estimates ###
    f.mean.Delta.t <- mean.dist2lastV.SIMEXit[1:(k+1)]
    linReg.coef <- summary(lm(cox.coef.SIMEXit[1,1:(k+1)] ~ f.mean.Delta.t))$coefficients
    
    ### Step 4: Extrapolation at mean dela.t of 1 ###
    extrap <- linReg.coef[1,1] + 1*linReg.coef[2,1] 

  } else if (extrapFUN == "log"){  

    ### Step 3: Linear regression on k+1 point estimates ###
    f.mean.Delta.t <- log(mean.dist2lastV.SIMEXit[1:(k+1)])
    linReg.coef <- summary(lm(cox.coef.SIMEXit[1,1:(k+1)] ~ f.mean.Delta.t))$coefficients
    
    ### Step 4: Extrapolation at mean dela.t of 1 ###
    extrap <- linReg.coef[1,1] + log(1)*linReg.coef[2,1] 
    
  } else if (extrapFUN == "log2"){  
    
    ### Step 3: Linear regression on k+1 point estimates ###
    f.mean.Delta.t <- log2(mean.dist2lastV.SIMEXit[1:(k+1)])
    linReg.coef <- summary(lm(cox.coef.SIMEXit[1,1:(k+1)] ~ f.mean.Delta.t))$coefficients
    
    ### Step 4: Extrapolation at mean dela.t of 1 ###
    extrap <- linReg.coef[1,1] + log2(1)*linReg.coef[2,1] 
    
  } else if (extrapFUN == "1df-quadratic"){  
    
    ### Step 3: Linear regression on k+1 point estimates ###
    f.mean.Delta.t <- mean.dist2lastV.SIMEXit[1:(k+1)] * mean.dist2lastV.SIMEXit[1:(k+1)]
    linReg.coef <- summary(lm(cox.coef.SIMEXit[1,1:(k+1)] ~ f.mean.Delta.t))$coefficient
    
    ### Step 4: Extrapolation at mean dela.t of 1 ###
    extrap <- linReg.coef[1,1] + 1*linReg.coef[2,1] 
  }
  
  return(list(expo.coef.corrected = extrap,
              linear.reg.step3 = linReg.coef,
              SIMEX_steps1_2.out = deparse(substitute(SIMEX_steps1_2.out)),
              n.iterations.k = k, 
              extrapolation.FUN = extrapFUN))
}    

#-----------------------------------------------------------------------------#
# plot.diag:                                                                  #
#   Produces a diagnostic plot to select the functional form for the          #
#   extrapolating function and the number of iterations k for the proposed    #
#   SIMEX-like method (see section 2.3 of the manuscript for details)         #
#-----------------------------------------------------------------------------#

## Arguments:
  # data.middle: a data frame, with event times imputed at the middle of  
  #     interval between the last 2 physician visits for cases (i.e. subjects  
  #     with an event), and controls (i.e. subjects without an event during 
  #     follow-up) censored at their last physician visit. Must have the 
  #     counting process format with each line corresponding to one and only 
  #     one time unit for a given individual. Must be the same data as used by
  #     the function SIMEX_steps1_2.  
  # SIMEX_steps1_2.out: output from the function SIMEX_steps1_2 call for which
  #     the diagnostic plot is produced.
  # nPrevTimes: number of previous time units for which the correlation is 
  #     calculated between exposure value at event time with values at previous
  #     times. Must be greater than 1/4 of the mean distance between the last 2 
  #     visits for cases in observed data, because the upper bound p'·delta.t  
  #     to be selected must be beyond this point.

## Details:
  # A dotted vertical line is added on the plot at the x axis value 
  # corresponding to the mean distance between the last 2 physician visits of 
  # cases in observed data, i.e. at 1/4 of this mean distance (Delta.t bar).   
  # The upper bound p'·delta.t to be selected, to determine the number k of   
  # SIMEX iterations, must be beyond this point.  

  # A dotdashed vertical line is added on the plot at the x axis value 
  # corresponding to the mean distance between the last 2 physician visits of 
  # cases in the modified data for the last SIMEX iteration included in the  
  # object SIMEX_steps1_2.out, i.e. at 1/4 of this mean distance (Delta.t bar).    
  # If the desired upper bound p'·delta.t to be selected is beyond this line,   
  # increase the number of SIMEX iterations in SIMEX_steps1_2.out.

plot.diag <- function(data.middle, SIMEX_steps1_2.out, nPrevTimes){

  id <- SIMEX_steps1_2.out$argument$id
  event <- SIMEX_steps1_2.out$argument$event 
  exposure <- SIMEX_steps1_2.out$argument$exposure
  mean.dist2lastV <- SIMEX_steps1_2.out$mean.dist2lastV
  
  # Check for errors in arguments
  if (deparse(substitute(data.middle)) != 
      SIMEX_steps1_2.out$argument$data.middle){
    stop("The argument data.middle differs than the one used for the function 
    SIMEX_steps1_2. The same data must be used by functions plot.diag and 
    SIMEX_steps1_2.")
  }
  if (nrow(data.middle) != SIMEX_steps1_2.out$data.middle$n.rows){
    stop("The number of rows of argument data.middle differs than the one used
    for the function SIMEX_steps1_2. The same data must be used by functions 
    plot.diag and SIMEX_steps1_2.")
  }
  if (length(unique(data.middle[,id])) != 
      SIMEX_steps1_2.out$data.middle$n.ids){
    stop("The number of subjects in argument data.middle differs than the one 
    used for the function SIMEX_steps1_2. The same data must be used by  
    functions plot.diag and SIMEX_steps1_2.")
  }  
  if (sum(data.middle[,event]) != SIMEX_steps1_2.out$data.middle$n.events){
    stop("The number of events in argument data.middle differs than the one 
    used for the function SIMEX_steps1_2. The same data must be used by  
    functions plot.diag and SIMEX_steps1_2.")
  } 
  if (nPrevTimes < SIMEX_steps1_2.out$mean.dist2lastV[1]/4){
    stop("Argument nPrevTimes must be greater than 1/4 of the mean distance  
    between the last 2 visits for cases in observed data, because the upper 
    bound p'·delta.t to be selected must be beyond this point.")  
  }
  if (nPrevTimes < tail(SIMEX_steps1_2.out$mean.dist2lastV, 1)/4){
    warning("Argument nPrevTimes is smaller than 1/4 of the mean distance   
    between the last 2 visits for cases in the modified data at last SIMEX
    iteration.")  
  }

  ### Calculate correlation between exposure at event time with previous times ###

  # Split dataset by subjects    
  data.middle.i <- split(data.middle, data.middle[,id])
  
  # Create the list of sequence of exposure values from event time and the  
  # preceeding values
  N <- length(unique(data.middle[,id]))
  seq.expo.cases.tmp <- lapply(1:N, function(i) seq.expoCases(
                                                  data.1id = data.middle.i[[i]], 
                                                  event = event, 
                                                  exposure = exposure,
                                                  nPrevTimes = nPrevTimes)) 
  # Remove NULL entries of controls in the list seq.expo.cases.tmp
  seq.expo.cases <- seq.expo.cases.tmp[-which(sapply(seq.expo.cases.tmp, is.null))]   
  
  # Correlation between the exposure at event time and with each of preceeding values
  cor.expoCases <- rep(NA, nPrevTimes + 1)
  for (j in 1:(nPrevTimes + 1)){
    x.expo.cases <- NULL
    y.expo.cases <- NULL
    n.events <- sum(data.middle[,event])
    for (i in 1:n.events){
      if (length(seq.expo.cases[[i]]) >= j){
        x.expo.cases <- c(x.expo.cases, seq.expo.cases[[i]][1])
        y.expo.cases <- c(y.expo.cases, seq.expo.cases[[i]][j])
      }
    }
    cor.expoCases[j] <- cor(x.expo.cases, y.expo.cases)
  }
  
  ### Diagnostic plot ###
  
  plot(0:(length(cor.expoCases) - 1), cor.expoCases, type="l", lwd=2,
       xlim=c(0, max(nPrevTimes, tail(mean.dist2lastV, 1)/4)), ylim=c(0,1),
       xlab="Time since imputed event time", 
       ylab="Correlation",  
       main="Diagnostic plot to select the extrapolating function\nand the number k of SIMEX iterations")
  
  # Add a dotted vertical line at the x axis value corresponding to the mean 
  # distance Delta.t between relevant visits in observed data, i.e. 1/4 of 
  # Delta.t  
  lines(c(mean.dist2lastV[1]/4, mean.dist2lastV[1]/4), c(0,1), lty=3)
  
  # Add a dotdashed vertical line at the x axis value corresponding to 1/4 of 
  # maximum Delta.t for the in the last iteration of SIMEX steps 1-2.
  lines(c(tail(mean.dist2lastV, 1)/4, tail(mean.dist2lastV, 1)/4), c(0,1), lty=4)
  
  legend("topright", 
        c('Correlation', 
          expression(paste('¼·', italic(bar(paste(Delta, 't'))), ' in observed data')),
          expression(paste('¼·', italic(bar(paste(Delta, 't'))), ' at last SIMEX iteration'))), 
         lty = c(1,3,4), lwd = c(2,1,1), cex = 0.9)
}

#-----------------------------------------------------------------------------#
# seq.expoCases:                                                              #
#   Internal function that returns the sequence of exposure values at the     #
#   time of event and at the preceding values for a case (NULL is returned    #
#   for a control). This function is used when calculating the correlation of #
#   exposure values for the diagnostic plot helping to choose the             #
#   extrapolating function.                                                   #
#-----------------------------------------------------------------------------#

## Arguments:
  # data.1id: data for one subject. Must have the counting process format with
  #     each line corresponding to one and only one time unit for a given 
  #     individual. 
  # nPrevTimes: number of previous time units for which the correlation  
  #     between exposure at event time with previous times is calculated.

seq.expoCases <- function(data.1id, event, exposure, nPrevTimes){  
  
  # Controls are ignored
  if (data.1id[,event][nrow(data.1id)] == 0){
    return(NULL)
    
  # For a case, return exposure values from event time and the preceding values  
  } else {
    return(data.1id[nrow(data.1id):max(1, nrow(data.1id) - nPrevTimes), exposure])  
  }
} 

#-----------------------------------------------------------------------------#
# SIMEX_bootCI:                                                               #
#   Bootstrap confidence interval (CI) for the SIMEX-corrected log HR         #
#-----------------------------------------------------------------------------#

## Arguments:
  # data.middle: a data frame, with event times imputed at the middle of  
  #     interval between the last 2 physician visits for cases (i.e. subjects  
  #     with an event), and controls (i.e. subjects without an event during 
  #     follow-up) censored at their last physician visit. Must have the
  #     counting process format with each line corresponding to one and only 
  #     one time unit for a given individual.
  # id: name of column in data.middle for subject unique identifiers, e.g. "id".
  # time: name of column in data.middle indicating the beginning of time  
  #     intervals, e.g. "start".
  # time2: name of column in data.middle indicating the end of time intervals,  
  #     e.g. "stop".
  # event: name of column in data.middle for the binary status indicator, e.g. 
  #     "event".
  # exposure: name of column in data.middle for the time-varying exposure  
  #     variable to which the SIMEX-like method is applied, e.g. "cumDur12w".
  # covariates: vector of names of covariates, other than exposure, to include
  #     in the models, e.g. c("sex", "age").
  # visit.2ndLast: name of column in data.middle for the time of the 2nd last 
  #     physician visit for cases (NA for controls). Constant value over all 
  #      rows of a subject.
  # visit.Last: name of column in data.middle for the time of the last  
  #     physician visit for cases and controls. Constant value over all rows 
  #     of a subject.
  # move2ndLast: vector of the number of time units by which the 2nd last visit
  #     for cases is moved backward at subsequent SIMEX iterations j. Values  
  #     must be greater than 0. E.g., seq(2, 50, by=2).
  # k: the selected number of SIMEX iterations.
  # extrapFUN: extrapolating function f(.) used in SIMEX steps 3 and 4. By
  #     default "linear", alternative options are "log", "log2" (log in base 2)
  #     and "1df-quadratic" (1 degree of freedom quadratic).
  # B: number of bootstrap resamples. 
  # conf.level: value of the confidence level of the confidence interval. The  
  #     default value is 95 (for 95% CI).

## Details:
  # The number of resamples should be at least B = 100, but a higher value would 
  # be preferable, e.g. B = 300 or even B = 1000.

## Value:
  # expo.coef.corrected: SIMEX-corrected coefficient (log HR) estimate for
  #     exposure.
  # expo.coef.corrected.bootCI: bootstrap CI for the SIMEX-corrected 
  #     coefficient (log HR) estimate for exposure.
  # arguments: list of arguments passed to the function SIMEX_bootCI.

SIMEX_bootCI <- function(data.middle, id, time, time2, event, exposure,  
                          covariates, visit.2ndLast, visit.Last, move2ndLast,
                          k, extrapFUN = "linear", B, conf.level = 95){

  # Arguments passed to the function SIMEX_bootCI (for output)
  arguments <- list(data.middle = deparse(substitute(data.middle)), 
                    id = id, time = time, time2 = time2, 
                    event = event, exposure = exposure, 
                    covariates = covariates, 
                    visit.2ndLast = visit.2ndLast, 
                    visit.Last = visit.Last, 
                    move2ndLast = move2ndLast,
                    n.iterations.k = k, 
                    extrapolation.FUN = extrapFUN,
                    B = B,
                    conf.level = conf.level)
  
  ### Apply the SIMEX method to original data ###
  
  SIMEX_12.out <- SIMEX_steps1_2(data.middle = data.middle, id = id,
                                  time = time, time2 = time2, 
                                  event = event, exposure = exposure,
                                  covariates = covariates, 
                                  visit.2ndLast = visit.2ndLast,
                                  visit.Last = visit.Last,
                                  move2ndLast = move2ndLast[1:k])
  SIMEX_34.out <- SIMEX_steps3_4(SIMEX_steps1_2.out = SIMEX_12.out, k = k,
                                  extrapFUN = extrapFUN)
  
  # Declare objects that will store the extrapolation for each resample
  b.extrap <- rep(NA, times=B)
  
  ### Bootstrap loop ###
  
  for (b in 1:B){
    cat(" b=", b, "\n")
    
    ids <- (unique(data.middle[,id]))
    max.id <- max(ids)
    
    ## Create the resample  
    
    # Resample with replacement subjects 
    ids.resamp <- sort(sample(ids, replace = TRUE))
    
    # Assemble the bootstrap resample (duplicated id ignored at this stage)
    bdata.middle <- data.middle[data.middle[,id] %in% ids.resamp, ] 
    # id variable in the resample (duplicated ids will need to be changed)
    bdata.middle$bid <- bdata.middle[,id]  
    
    # Deal with duplicated ids
    step <- 1
    repeat {
      # Select duplicated ids only
      ids.resamp <- ids.resamp[duplicated(ids.resamp) == TRUE]
      
      # Break loop when no more duplicated ids
      if (length(ids.resamp) == 0) break  
      dup.mid <- data.middle[data.middle[,id] %in% ids.resamp, ]
      # Change id of duplicates
      dup.mid$bid <- max.id * step + dup.mid$id
      bdata.middle <- rbind(bdata.middle, dup.mid)
      step <- step+1
    }
    
    ## Apply the SIMEX method to the resample    
    b.SIMEX_12.out <- SIMEX_steps1_2(data.middle = bdata.middle, id = "bid",
                                     time = time, time2 = time2, event = event,
                                     exposure = exposure, 
                                     covariates = covariates, 
                                     visit.2ndLast = visit.2ndLast, 
                                     visit.Last = visit.Last,
                                     move2ndLast = move2ndLast[1:k])
    b.SIMEX_34.out <- SIMEX_steps3_4(SIMEX_steps1_2.out = b.SIMEX_12.out, 
                                     k = k, extrapFUN = extrapFUN)   
    b.extrap[b] <- b.SIMEX_34.out$expo.coef.corrected
  }
  
  ### Bootstrap 95% CI ###
  
  alpha <- (1 - (conf.level/100))
  cox.SIMEX.bootCI <- quantile(b.extrap, prob=c(alpha/2, (1 - (alpha/2))))         
  
  return(list(expo.coef.corrected = SIMEX_34.out$expo.coef.corrected,
              expo.coef.corrected.bootCI = cox.SIMEX.bootCI,
              arguments = arguments))
}

