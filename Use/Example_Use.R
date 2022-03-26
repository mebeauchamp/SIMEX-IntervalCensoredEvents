###############################################################################
###                                                                         ###
### Example of use on a dataset of the methods presented in the manuscript: ### 
### "Adapting SIMEX to correct for bias due to interval-censored outcomes   ###
### in survival analysis with time-varying exposure"                        ###
###                                                                         ###
### Code by Marie-Eve Beauchamp                                             ###
### Last update: March 25, 2022                                             ###
###                                                                         ###
###############################################################################

### Description: This program illustrates with a dataset how to:
  #  - Implement the proposed SIMEX-like method to correct the effect for a 
  #    time-varying exposure (TVE)
  #  - Produce the proposed diagnostic plot to select the extrapolating 
  #    function and the number k of SIMEX iterations  
  #  - Calculate the bootstrap confidence interval for the corrected TVE effect  
  #    from the SIMEX-like method

# Note: The program Functions_SIMEXforIntervalCensoredOutcomes.R contains the
  # functions to implement the proposed methods for a given dataset. See 
  # Functions_SIMEXforIntervalCensoredOutcomes.R for the description and 
  # details about the functions.   

rm(list=ls())
gc()

options(warn=1)

## Set the current folder
setwd("C://path//to//your//directory")

#-----------------------------------------------------------------------------#
# Load the data                                                               #
#-----------------------------------------------------------------------------#

load("Data_for_Example.RData")
ls()

########################
# The dataset: dat.end #
########################

## Details: 
  # A data frame including a time-varying exposure (TVE). 

  # The functions in the program Functions_SIMEXforIntervalCensoredOutcomes.R
  # require that the data have the counting process format with each line 
  # corresponding to one and only one time unit, delimited by (start, stop], 
  # for a given individual.

  # In the dataset dat.end the event times are imputed at the end of interval 
  # between the 2nd last and last physician visits. However, your data could 
  # have instead event times imputed at the middle of these intervals. Then,  
  # the step of data transformation to impute event times at the interval 
  # middle could be skipped.
  
head(dat.end)

## Description of columns in dat.end:
  # id: unique identifier of the subjects.
  # event: binary indicator of an event for this observation (1=event, 0=no event).
  # start: beginning time of each time interval.
  # stop: end time of each time interval.
  # cumDur12w: cumulative duration of drug exposure in the last 12 weeks. 
  #     This is the time-varying exposure of interest.
  # sex: sex of the subject (1=male, 0=female). This is a covariate for 
  #     model adjustment.
  # agem40: age of the subject minus 40. This is a covariate for model 
  #     adjustment.
  # vis2ndLast: time of 2nd last physician visit if the current subject is a 
  #     case (i.e. subject with an event). NA for controls (i.e. subject 
  #     censored without event). Constant value over all rows of a subject.
  #     If some cases do not have a visit before their true event time, assign 
  #     t0 for their 2nd last visit.
  # visLast: time of last physician visit for the current subject. Constant   
  #     value over all rows of a subject.

# Shows beginning and end of follow-up for a case (i.e. subject with an event)
cases.id <- unique(dat.end$id[dat.end$event == 1])
head(dat.end[dat.end$id == cases.id[1],])
tail(dat.end[dat.end$id == cases.id[1],])

# Shows beginning and end of follow-up for a control (i.e. subject censored)
controls.id <- setdiff(unique(dat.end$id), cases.id)
head(dat.end[dat.end$id == controls.id[1],])
tail(dat.end[dat.end$id == controls.id[1],])

###############
# Visit times #
###############

is.list(visits)
length(visits)
visits[[cases.id[1]]]    # Visit times for a case
visits[[controls.id[1]]] # Visit times for a control

## Details: 
  # Each item in the list visits corresponds to all visit times for 1 subject.

  # Note: in the current implementation of the SIMEX-like method, the frequency
  # of visits at step 1 of the SIMEX procedure are decreased by moving backward
  # the time of the 2nd last visit for cases (algorithm (i) in section A.1 of
  # the Supporting Information). The information on 2nd last and last visits is
  # already included in the object dat.end. Therefore, the object visits will 
  # not be used here in the application of the SIMEX-like, but is provided for
  # information.

#-----------------------------------------------------------------------------#
# Source the functions to implement the proposed methods                      #
#-----------------------------------------------------------------------------#

source("Functions_SIMEXforIntervalCensoredOutcomes.R")

#-----------------------------------------------------------------------------#
# SIMEX-like proposed method: STEPS 1 and 2                                   #
# (see section 2.2 in the manuscript for details)                             # 
#-----------------------------------------------------------------------------#

## Data transformation to impute event times at the middle of intervals        
  # between the 2nd last and last physician visits for cases                  

dat.mid <- imputeEventMid(data = dat.end, id = "id", time = "start", 
                          time2 = "stop", event = "event", 
                          visit.2ndLast = "vis2ndLast", visit.Last = "visLast")

# Confirm that the number of rows decreases from dat.end t dat.mid, but that 
  # the number of events and subjects remain identical
dim(dat.end)
dim(dat.mid)
sum(dat.end$event)
sum(dat.mid$event)
length(unique(dat.end$id))
length(unique(dat.mid$id))

## Apply steps 1 and 2 of the proposed SIMEX-like method for the time-varying
  # exposure "cumDur12w", by moving backward the 2nd last physician visit of  
  # cases by 2,4,6,...,50 time units at each subsequent SIMEX iteration.
  # This results in j=0,1,...,25 SIMEX iterations, with j=0 for the Cox PH 
  # model estimation on the original (unmodified) data.

# Note: in this implementation of the SIMEX-like method, the frequency of    
  # visits at step 1 are decreased by moving backward the time of the 2nd   
  # last visit (algorithm (i) in section A.1 of the Supporting Information).
  # Algorithm (ii), which implies to delete the previous visit, was not 
  # implemented.

SIMEX12.out <- SIMEX_steps1_2(data.middle = dat.mid, id = "id", time = "start",
                              time2 = "stop", event = "event", 
                              exposure = "cumDur12w", 
                              covariates = c("sex","agem40"), 
                              visit.2ndLast = "vis2ndLast", visit.Last="visLast",
                              move2ndLast = seq(2, 50, by=2))
# Full output 
SIMEX12.out
  
# Mean distance between last 2 visits for cases at each SIMEX iteration 
  # j=0,1,...,25, and difference with mean distance at the previous iteration
round(SIMEX12.out$mean.dist2lastV.SIMEXit, 1)

# Difference with mean distance at the previous iteration
round(diff(SIMEX12.out$mean.dist2lastV.SIMEXit), 1)

  # Note: the increase in mean distance between last 2 visits for cases is 
  # smaller than 2, despite that move2ndLast = seq(2, 50, by=2), because some 
  # events occurred between t0 and the 1st visit after t0 for some cases.
  # Therefore, the 2nd last visit (t0 for those cases) could not be moved 
  # earlier at each SIMEX step 1 iteration.

  # The 2nd last visit can never be moved earlier than t0, which explains why
  # the increase in mean distance becomes smaller with higher iteration  
  # numbers, as t0 is reached for more cases when moving backward their 2nd 
  # last visit.

#-----------------------------------------------------------------------------#
# DIAGNOSTIC PLOT: to select the extrapolating function and the number k of   #
# iterations for the proposed SIMEX-like method                               #
# (see section 2.3 in the manuscript for details)                             # 
#-----------------------------------------------------------------------------#

## Step a): Visual assessment of the pattern of changes in correlation between
  # exposure value at event time with values at the previous times.

plot.diag(data.middle = dat.mid, SIMEX_steps1_2.out = SIMEX12.out, 
          nPrevTimes = 40)

## Step b): Check for major departure from linearity for some interval between
  # 0 and an upper bound p'·delta.t beyond the dotted vertical line.
  # => Linear function is selected (default option) 

## Step c): If the linear function is rejected in step b), select the 
  # appropriate function among log, log in base 2, or 1-df quadratic function.
  # => Step skipped

## Step d): Identify an approximate upper bound p'·delta.t, beyond the vertical
  # dotted line, over which the function chosen in step b) or c) mimics the 
  # corresponding correlations.

# Add a dashed bold vertical line at the selected upper bound (p'·delta.t) = 14
lines(c(14, 14), c(0,1), lty=2, lwd=2)

# Select the number k of SIMEX iteration corresponding approximately to 14*4 = 56
SIMEX12.out$mean.dist2lastV.SIMEXit
  # At position 19 is approximately 56
SIMEX12.out$mean.dist2lastV.SIMEXit[19] 
  # => Select k=18 SIMEX iterations 

#-----------------------------------------------------------------------------#
# SIMEX-like proposed method: STEPS 3 and 4                                   #
# (after the extrapolating function and number of k iterations are selected;  #
# see section 2.2 in the manuscript for details)                              # 
#-----------------------------------------------------------------------------#

SIMEX34.k18.out <- SIMEX_steps3_4(SIMEX_steps1_2.out = SIMEX12.out, 
                                  k = 18, extrapFUN = "linear")

# Full output
SIMEX34.k18.out

# SIMEX-corrected log HR estimate for the exposure
SIMEX34.k18.out$expo.coef.corrected

# SIMEX-corrected HR estimate for the exposure
exp(SIMEX34.k18.out$expo.coef.corrected)

#-----------------------------------------------------------------------------#
# Bootstrap 95% confidence interval (CI) for the SIMEX-corrected log hazard   # 
# ratio (HR) for exposure                                                     #
#-----------------------------------------------------------------------------#

SIMEXbootCI.k18.out <- SIMEX_bootCI(data.middle = dat.mid, id = "id", 
                                time = "start", time2 = "stop", 
                                event = "event", exposure = "cumDur12w", 
                                covariates = c("sex","agem40"), 
                                visit.2ndLast = "vis2ndLast", 
                                visit.Last = "visLast",
                                move2ndLast=seq(2, 50, by=2), k = 18,
                                extrapFUN = "linear", B = 10)

# Note: B = 10 is used for illustration only. A value of at least B = 100 or
  # should be used, e.g. B = 300 or even B = 1000.

# Full output
SIMEXbootCI.k18.out

# SIMEX-corrected log HR and 95% boostrap CI for exposure
round(SIMEXbootCI.k18.out$expo.coef.corrected, 3)
round(SIMEXbootCI.k18.out$expo.coef.corrected.bootCI, 3)

# SIMEX-corrected HR and 95% boostrap CI for 1 week increase in duration of 
  # drug use
round(exp(SIMEXbootCI.k18.out$expo.coef.corrected), 3)
round(exp(SIMEXbootCI.k18.out$expo.coef.corrected.bootCI), 3)

#-----------------------------------------------------------------------------#
# Comparison with conventional Cox models                                     #
#-----------------------------------------------------------------------------#

## Event times imputed at the middle of intervals
cox.mid <- coxph(Surv(start, stop, event) ~ cumDur12w + sex + agem40, 
                 data = dat.mid)
#summary(cox.mid)

## Event times imputed at the last visit
cox.end <- coxph(Surv(start, stop, event) ~ cumDur12w + sex + agem40, 
                 data = dat.end)
#summary(cox.end)

# HR from function SIMEX_steps3_4 (without CI)
round(exp(SIMEX34.k18.out$expo.coef.corrected), 3)

# HR from function SIMEX_bootCI (with CI)
round(exp(SIMEXbootCI.k18.out$expo.coef.corrected), 3)
round(exp(SIMEXbootCI.k18.out$expo.coef.corrected.bootCI), 3)

# Conventional Cox with event imputed at mid-intervals
round(summary(cox.mid)$conf.int[1, c(1,3,4)], 3)

# Conventional Cox with event imputed at end of intervals
round(summary(cox.end)$conf.int[1, c(1,3,4)], 3)

