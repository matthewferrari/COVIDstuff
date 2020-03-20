##############################################################################
# Consider two control scenarios for control:
# - there are 3 levels of control: mild, medium, severe, which reduce the mean-field transmission rate
# 
# Scenario 1: initial trigger results in move first to "mild" control, subsequent increases in outbreak severity (accel_triggers) result in progressive acceleration to "medium" and "severe". If outbreak subsides (return_trigger_accel) then control can decelerate.  
# Scenario 2: initial trigger results in move first to "severe" control, subsequent increases in outbreak severity (decel_triggers) result in progressive acceleration to "medium" and "mild". If outbreak increases (return_trigger_decel) then control can accelerate.  
##############################################################################
##############################################################################
# Intial observations -- 
# Scenario 1 tends to result in lowest cost outcome; mildest control for the longest period (area of grey). 
# Scenario 2 tends to result in smallest overall size of outbreak and (in general) lowest overall magnitude (peak height)
#  -- when R0 is high, overall outbreak size is similar, but peak hight much higher in scenario 1
##############################################################################

mild <- .9  # 1 - reduction in transmission in mild control
medium <- .5  # 1 - reduction in transmission in medium control
severe <- .01  # 1 - reduction in transmission in severe control

par(mfrow=c(1,2))
duration <- 500 # simulation duration
R0 <- 2.5 # as you make this higher, severity of initial outbreak changes dramatically 
population_size <- 1000
initial_cases <- 5  # number of initial cases
transmission <- R0/population_size  
silent_spread <- 7 # number of days of silent spread before any control can start

observation_lag <- 1  # values 1-6; number of days lag in information for triggers
accel_triggers <- c( 50, 100, 200, 400) # for levels of accelerating -- move up a level if you hit this case level
decel_triggers <- c(50, 25, 10, 5, 1) # note there are 5 levels; one for start, 4 for decelerating - move down a level if you hit this case level
return_trigger_accel <- 25  # what is the trigger for decelerating control in accelerating scenario
return_trigger_decel <- 25  # what is the trigger for accelerating control in decelerating scenario
change_delay <- 7 # must wait this many days to change a control level

#############################################################################################
# Here is a short script for a stochastic SIR model 
# this genrates cases per day and then aggregates into a series of observation intervals; here 13 biweeks, or half a year
# 

gilstep<-function(S.t,I.t,I.p,beta,births){
  # S.t = current susceptibles
  # I = current infecteds
  # I.p = infecteds one infectious period previous (note, this is analogous to your fixed infectious period)
  # beta = transmission rate
  # births = births per time step
  #browser()
  ifelse(I.t>0,nSI<-rpois(1,max(rnorm(1,beta,.2*beta),1e-5)*S.t*I.t),nSI<-0) # number of new cases in time step t
#  ifelse(I.t>0,nSI<-rpois(1,rgamma(1,beta,1)*S.t*I.t),nSI<-0) # number of new cases in time step t
  nIR<-I.p #rpois(1,I.p)                               # number of cases recovering in time step t
  I.new<-max(0,I.t+nSI-nIR)                       # update I
  S.new<-max(0,S.t-nSI+births)                    # update S
  return(list(S=S.new,I=I.new,new.case=nSI))      # store S, I, and new cases
}


################################################################
################################################################
# new stuff
################################################################
################################################################
# accelerating control
gilrun<-function(T,S.init,I.init,beta_levels,beta_triggers_high,beta_triggers_low,obs_lag,inf.per,births=0){
  # T = the total time of the simulation
  # S.init = the initial susceptible population
  # I.init = the initial infectious population
  # beta = the transmission rate
  # inf.per = the infectious period
  # births = the births per time step
  # 
  
  S<-S.init
  I<-I.init
  trig <- 1
  trig.stor <- 1
  delay_counter <- 0
  new.case<-0

    for(i in 2:T){
    #browser()
    ifelse(i>(inf.per+1),I.p<-new.case[i-inf.per],I.p<-0) # If time is earlier than the first infectious period, nobody can recover

        # count cases in beds (I) set control and beta
    #
    # Don't trigger unless cases get at least to X
    if(i > silent_spread){
      if(I[i-obs_lag] > beta_triggers_high[trig] & delay_counter > change_delay){
      trig <- trig + 1
      delay_counter <- 0
      }
      if(I[i-obs_lag] < beta_triggers_low & delay_counter > change_delay){
        trig <- max(trig - 1,1)
        delay_counter <- 0
      }
    }
    delay_counter <- delay_counter + 1  # increment delay counter
    beta <- beta_levels[trig] 
    #
if(length(I[i-1])==0){browser()}
    out<-gilstep(S[i-1],I[i-1],I.p,beta,births)           # generate new infectious and recoveries 
    I<-c(I,out$I)                                         # Append new I
    new.case<-c(new.case,out$new.case)                    # Append new cases
    S<-c(S,out$S)                                         # Append new S
    trig.stor <- c(trig.stor,trig)
  }
  return(list(I=I,S=S,new.case=new.case,trig.stor =trig.stor))
}

out<-gilrun(duration,population_size,initial_cases,
            beta_levels = c( 1, mild, medium, severe) * transmission/14, #beta_levels -- accelerating
            beta_triggers_high = accel_triggers, # beta_triggers -- accelerating
            beta_triggers_low = return_trigger_accel, # beta_triggers -- accelerating
            obs_lag = observation_lag,
            14,0)
plot(NA,xlim=c(1,200),ylim=c(0,1000),ylab="cases",xlab="time") # plot panels
title("accelerating control")

#barplot(out$trig.stor*100,col=grey(.25),add=T)
width <- .45  # width of bar for plotting control levels
severity <- (out$trig.stor-1) * 100  # scale so it shows up
severity[1:silent_spread] <- 0 # don't plot control before it starts
rect((1:length(out$I)) - width, rep(0,length(out$I)), (1:length(out$I)) + width,severity, col = gray(.75), border=grey(.75)) # plot bars
points(out$I,col="red",ylim=c(0,1e3)) # plot cases 
points(out$S,col="blue")  # plot susceptibles
#points(out$new.case,col=grey(.5))



################################################################
# decelerating control
gilrun<-function(T,S.init,I.init,beta_levels,beta_triggers_high, beta_triggers_low, obs_lag,inf.per,births=0){
  # T = the total time of the simulation
  # S.init = the initial susceptible population
  # I.init = the initial infectious population
  # beta = the transmission rate
  # inf.per = the infectious period
  # births = the births per time step
  # 
  
  S<-S.init
  I<-I.init
  trig <- 1
  trig.stor <- 1
  delay_counter <- 0
  new.case<-0
  for(i in 2:T){
    #browser()
    ifelse(i>(inf.per+1),I.p<-new.case[i-inf.per],I.p<-0) # If time is earlier than the first infectious period, nobody can recover
    
    # count cases in beds (I) set control and beta
    #
    # Don't trigger unless cases get at least to X
    if(trig ==1){
     if(i >  silent_spread) {
       #browser()
       if(I[i-obs_lag] > beta_triggers_low[trig] & delay_counter > change_delay){ # force change_delay days silent spread
      trig <- min(trig + 1,length(beta_triggers_low))
      delay_counter <- 0
       }
     } 
    }
    if(trig >1){
      if(i >  silent_spread){
        if(I[i-obs_lag] < beta_triggers_low[trig] & delay_counter > change_delay){ # force change_delay days silent spread
        trig <- min(trig + 1,length(beta_triggers_low))
        delay_counter <- 0
        }
        if(trig >=4 & I[i-obs_lag] > beta_triggers_high & delay_counter > change_delay){ # force change_delay days silent spread
          trig <- max(trig - 1,2)
          delay_counter <- 0
        }
      } 
    }
    delay_counter <- delay_counter + 1  # increment delay counter
    beta <- beta_levels[trig] 
    #cat(beta,"-",trig,".\n")
    #
    #browser()
    out<-gilstep(S[i-1],I[i-1],I.p,beta,births)           # generate new infectious and recoveries 
    I<-c(I,out$I)                                         # Append new I
    new.case<-c(new.case,out$new.case)                    # Append new cases
    S<-c(S,out$S)                                         # Append new S
    trig.stor <- c(trig.stor,trig)
  }
  return(list(I=I,S=S,new.case=new.case,trig.stor =trig.stor))
}

out<-gilrun(duration,population_size,initial_cases,
            beta_levels = c(1,severe,medium,mild ,1) * transmission/14, #beta_levels -- accelerating
            beta_triggers_high = c(return_trigger_decel), # beta_triggers -- accelerating
            beta_triggers_low = decel_triggers, # beta_triggers -- accelerating
            obs_lag = observation_lag,
            14,0)


plot(NA,xlim=c(1,200),ylim=c(0,1000),ylab="cases",xlab="time")
title("decelerating control")

#barplot(out$trig.stor*100,col=grey(.25),add=T)
width <- .45 # width of bar for plotting control levels
severity <- c(0,3,2,1,0)[out$trig.stor] * 100 # scale so it shows up
severity[which(out$trig.stor==1)] <- 0 # don't plot control before it starts
rect((1:length(out$I)) - width, rep(0,length(out$I)), (1:length(out$I)) + width,severity, col = gray(.75), border=grey(.75)) # plot bars
points(out$I,col="red",ylim=c(0,1e3)) # plot cases 
points(out$S,col="blue")  # plot susceptibles
#points(out$new.case,col=grey(.5))


legend(0,600,legend=c("susceptible", "infected","control level"),pch=c(1,1,15),col=c(4,2,grey(.75)),bty="n")


################################################################
################################################################

