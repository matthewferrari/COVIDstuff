# SIR model to track isolation and quarantine room usage for COVID-like parameters
# triggers NPIs that reduce transmission rate if a threshold capacity is reached
# 

library(mc2d)
sir_step <- function (sims, S, E, I1, I2, R, N, newSympt1, newSympt2, beta, theta,gamma_I1I2,gamma_I2R, delta.t=1, tests,ppn_sympt=1/20,contacts=7, CT=TRUE, ...) {
  # adapted from Aaron King's code
  #sims - number of stochastic simulations
  #S vector of susceptibles, length=sims
  #E vector of exposed, length=sims
  #I1 vector of pre symptomatic infecteds, length=sims
  #I2 vector of possibly symptomatic infecteds, length=sims
  #R vector of recoverdeds, length=sims 
  #N vector of population size, length=sims
  #newSympt1  counter for new possibly symptomatic individuals, to allow 2 day health seeking delay
  #newSympt2  counter for new possibly symptomatic individuals, to allow 2 day health seeking delay
  #beta  transmission rate
  #theta rate from exposed to infected 
  #gamma_I1I2 rate from pre-symptomatic infecteds to possibly symptomatic
  #gamma_I2R rate from possibly symptomatic to recovered
  #delta.t time step length default is 1 day timestep
  #tests number of asymptomatic tests per day
  #ppn_sympt proportion of possibly symptomatic individuals who are symptomatic
  #contacts  average number of contaccts per individual
  #CT TRUE == very efficient contact tracing, implies contact tracing disproportionately finds infected indivduals
  
  # transitions between classes
  dN_SE <- rbinom(n=sims,size=S,prob=1-exp(-beta*(I1+I2)/N*delta.t)) + rbinom(sims,1,.1) # add random introductions
  dN_EI1 <- rbinom(n=sims,size=E,prob=1-exp(-theta*delta.t))
  dN_I1I2 <- rbinom(n=sims,size=I1,prob=1-exp(-gamma_I1I2*delta.t))
  dN_I2R <- rbinom(n=sims,size=I2,prob=1-exp(-gamma_I2R*delta.t))

  # update classes
  S <- S - dN_SE 
  E <- E + dN_SE - dN_EI1 
  I1 <- I1 + dN_EI1 - dN_I1I2
  I2 <- I2 + dN_I1I2 - dN_I2R
  newSympt2 <- newSympt1
  newSympt1 <- dN_I1I2
  newSymptReported <- rbinom(sims,newSympt2,ppn_sympt) # randomly draw symtomatic individuals
  R <- R + dN_I2R
  
  out <- cbind(S ,  E,  I1,  I2, R, dN_I1I2 ) # assume that I1->I2 is when cases become detectable
  atests <- rmultinomial(sims,rep(tests,sims),out[,1:5]) # radomly draw the indivdiuals tested via asymptomatic tests
  #browser()
  # randomly draw the contacts from the different classes
  ifelse(CT==FALSE,contacts <- rmultinomial(sims,rep(rpois(sims,contacts)*(newSymptReported + apply(atests[,3:4], 1, sum)),sims),matrix(c(1,1,1,1,1),nr=sims,nc=5,byrow=T)*out[,1:5]), # tracing preferentially identifies infected individuals
  contacts <- rmultinomial(sims,rep(contacts*(newSymptReported + apply(atests[,3:4], 1, sum)),sims),log(out[,1:5]+1))) # tracing preferentially identifies infected individuals
  #browser()
  
  atests.isolate <- atests # holder for which tests will be positive that need to be isolated 
  atests.isolate[,c(1,2,5)] <- 0 # set non-infected classes to 0
  sympt.isolate <- matrix(0,nr=sims,nc=5) # stoarge for symoptomatic cases to isolate
  sympt.isolate[,4] <- newSymptReported # fill this in with the randomly drawn symptomatics
  out[,1:5] <- pmax(out[,1:5] - sympt.isolate - atests.isolate - contacts,0) # remove asymptomatic tests and contacts from population
  out <- cbind(out,atests, newSympt1, newSympt2, newSymptReported, contacts) # store all states -- SIR states plus tested, reported, contacts
}


############################################################################################################
############################################################################################################
############################################################################################################

plot_metrics <- function(introductions, ppn_sympt){
# wrapper function to plot figures
  
library(RColorBrewer)
pal <- brewer.pal(4,"Dark2")
par(mfcol = c(4,12))
par(mar=c(1,2,2,2)+0.1)

ct <- c(FALSE, TRUE)  # effectiveness of contact tracing
thresh <-c(10000,150) # threshold to trigger imposition of NPIs (very high number is eqivalent to no thershold)
                      # this is in terms of isolation and quarantine rooms occupied
tst <- c(0,100,300) # daily asymptomatic tests

for(cc in 1:2){
for(tt in 1:2){
for(hh in 1:3){

T <- 100        # time to simulate over, we only care about start
sims <- 1000    # number of simulations
inf <- matrix(NA,T,sims)  # storage for infectious class
case <- matrix(NA,T,sims) # storage for daily cases
distancing_reduction <- 0.5 # if trigger is crossed, NPIs are imposed and transmission is reduced by this fraction

R0 <- 2.5   # R0
RE <- R0  # RE assuming some fraction of population is already immune
beta.normal <- RE * (1/9)  # calculate Beta
beta.outbreak <- RE * (1/9) * (1-distancing_reduction)   # calculate Beta
beta_vec <- rep(beta.normal,sims)
theta <- 1/5  # 5 days from infection to infectious
gamma_I1I2 <- 1/2 # 2 days asymptomatic infectious
gamma_I2R <- 1/7 # 7 days infectious (this is probably too short)

tests <- tst[hh]           # assume this many asymptomatic tests per day


# storage 
S <- matrix(round(50000*0.85),1,sims)  # start with 50K * 0.85 susceptible
E <- matrix(introductions,1,sims)
I1 <- matrix(0,1,sims)                 # start with 5 asymptomatic infectious
I2 <- matrix(0,1,sims)
R <- matrix(50000 - S,1,sims)
sympt1 <- matrix(0,1,sims)
sympt2 <- matrix(0,1,sims)
symptrep <- matrix(0,1,sims)
new_contacts <- matrix(0,1,sims)
isolation <- matrix(0,1,sims)
quarantine <- matrix(0,1,sims)
new_cases <- matrix(0,1,sims)
true_test_positives <- matrix(0,1,sims)
N <- S+E+I1+I2+R        # population size

qi_trigger <- numeric(sims)
for(ts in 2:T){
  
  out <- sir_step(sims, S[ts-1,], E[ts-1,], I1[ts-1,], I2[ts-1,], R[ts-1,], N[ts-1,], sympt1[ts-1,], sympt2[ts-1,], beta_vec, theta,gamma_I1I2,gamma_I2R, delta.t=1, tests, CT=ct[tt], ppn_sympt = ppn_sympt) # call to SIR step function above
  #browser()
  S <- rbind(S,out[,1])  # update state
  E <- rbind(E,out[,2])  # update state
  I1 <- rbind(I1,out[,3])  # update state
  I2 <- rbind(I2,out[,4])  # update state
  R <- rbind(R,out[,5])  # update state
  N <- rbind(N,N[ts-1])  # update state
  sympt1 <- rbind(sympt1,out[,12])  # update state
  sympt2 <- rbind(sympt2,out[,13])  # update state
  symptrep <- rbind(symptrep,out[,14])  # update state
  new_contacts <- rbind(new_contacts,apply(out[,15:19],1,sum))  # update state
  true_test_positives <- rbind(true_test_positives,apply(out[,9:10],1,sum))  # update state
  new_cases <- rbind(new_cases,out[,6])  # update state
  ####################################################################################
  # Total in Isolation/Qurantine
  isolation <- rbind(isolation,apply(true_test_positives[(max(1,ts-10)):ts,],2,sum) + apply(symptrep[(max(1,ts-10)):ts,],2,sum)) # isolate for 10 days
  quarantine <- rbind(quarantine, apply(new_contacts[(max(1,ts-14)):ts,],2,sum) ) # quarantine for 14 days
  #browser()
  ####################################################################################
  # is the number in isolation/quarantine greater than your qi_trigger value? (or has it ever been)
  qi_trigger <- qi_trigger + ((isolation[ts,] + quarantine[ts,]) > thresh[cc] | qi_trigger)
  # note that qi_trigger also marks the time step at which the trigger is reached
  
  #test_trigger <- true_test_positives[ts,] > 1 # can set alternate triggers based on # positive tests, prevalence, etc.
  #trigger <- pmax(trigger,(apply(out[,9:10],1,sum)>0),)
  #if ... 
  beta_vec <- rep(beta.normal,sims)
  beta_vec[which(qi_trigger>0)] <- beta.outbreak
  #
  # add counter to identify when trigger happens
  #if(ts == 20){browser()}
  
}

inf <- I1+I2   # total infectious
case <- new_cases # daily cases

ifelse(hh==1 & cc==1 & tt == 1, lb<-TRUE, lb<-FALSE)
inf.bd <- apply(inf,1,quantile, probs=c(.05,.95))
plot(NA,xlim=c(0,100),ylim=c(0,4000),xlab="days from start",ylab="active infections",axes=F)
polygon(c(1:100,100:1),c(inf.bd[1,],rev(inf.bd[2,])),col=pal[1],border=pal[1])
axis(1,labels=FALSE)
axis(2,labels=lb)

sympt.bd <- apply(symptrep,1,quantile, probs=c(.05,.95))
plot(NA,xlim=c(0,100),ylim=c(0,20),xlab="days from start",ylab="symptomatic cases",axes=F)
polygon(c(1:100,100:1),c(sympt.bd[1,],rev(sympt.bd[2,])),col=pal[2],border=pal[2])
axis(1,labels=FALSE)
axis(2,labels=lb)

test.bd <- apply(true_test_positives,1,quantile, probs=c(.05,.95))
plot(NA,xlim=c(0,100),ylim=c(0,20),xlab="days from start",ylab="asymptomatic test pos",axes=F)
polygon(c(1:100,100:1),c(test.bd[1,],rev(test.bd[2,])),col=pal[3],border=pal[3])
axis(1,labels=FALSE)
axis(2,labels=lb)

quar.bd <- apply(quarantine+isolation,1,quantile, probs=c(.05,.95))
plot(NA,xlim=c(0,100),ylim=c(0,1500),xlab="days from start",ylab="students in isol/quar",axes=F)
polygon(c(1:100,100:1),c(quar.bd[1,],rev(quar.bd[2,])),col=pal[4],border=pal[4])
axis(1,labels=TRUE)
axis(2,labels=lb)

abline(h=300) # add threshold at 300 isolation/quarantine rooms
# add hashmark on top boder for each simulation that crossed the 150 room trigger
rug((100-qi_trigger),side=3,ticksize = .3,col="grey")
}
}
}  
par(mar=c(5, 4, 4, 2) + 0.1)
}

# run and plot setting with 0 cases on day 1 and 1/20 infected individuals are symptomatic
plot_metrics(0, 1/20)

# run and plot setting with 50 cases on day 1 and 1/20 infected individuals are symptomatic
plot_metrics(50, 1/20)

# run and plot setting with 5 cases on day 1 and 1/100 infected individuals are symptomatic
plot_metrics(5, 1/100)

############################################################################################
############################################################################################
# alternate version with negative binomial introductions

sir_step <- function (sims, S, E, I1, I2, R, N, newSympt1, newSympt2, beta, theta,gamma_I1I2,gamma_I2R, delta.t=1, tests,ppn_sympt=1/20,contacts=7, CT=TRUE, ...) {
  # adapted from Aaron King's code
  #sims - number of stochastic simulations
  #S vector of susceptibles, length=sims
  #E vector of exposed, length=sims
  #I1 vector of pre symptomatic infecteds, length=sims
  #I2 vector of possibly symptomatic infecteds, length=sims
  #R vector of recoverdeds, length=sims 
  #N vector of population size, length=sims
  #newSympt1  counter for new possibly symptomatic individuals, to allow 2 day health seeking delay
  #newSympt2  counter for new possibly symptomatic individuals, to allow 2 day health seeking delay
  #beta  transmission rate
  #theta rate from exposed to infected 
  #gamma_I1I2 rate from pre-symptomatic infecteds to possibly symptomatic
  #gamma_I2R rate from possibly symptomatic to recovered
  #delta.t time step length default is 1 day timestep
  #tests number of asymptomatic tests per day
  #ppn_sympt proportion of possibly symptomatic individuals who are symptomatic
  #contacts  average number of contaccts per individual
  #CT TRUE == very efficient contact tracing, implies contact tracing disproportionately finds infected indivduals
  
  # transitions between classes
  dN_SE <- rbinom(n=sims,size=S,prob=1-exp(-beta*(I1+I2)/N*delta.t)) + rnbinom(sims,mu=1,size=.1) # add random introductions
  dN_EI1 <- rbinom(n=sims,size=E,prob=1-exp(-theta*delta.t))
  dN_I1I2 <- rbinom(n=sims,size=I1,prob=1-exp(-gamma_I1I2*delta.t))
  dN_I2R <- rbinom(n=sims,size=I2,prob=1-exp(-gamma_I2R*delta.t))
  
  # update classes
  S <- S - dN_SE 
  E <- E + dN_SE - dN_EI1 
  I1 <- I1 + dN_EI1 - dN_I1I2
  I2 <- I2 + dN_I1I2 - dN_I2R
  newSympt2 <- newSympt1
  newSympt1 <- dN_I1I2
  newSymptReported <- rbinom(sims,newSympt2,ppn_sympt) # randomly draw symtomatic individuals
  R <- R + dN_I2R
  
  out <- cbind(S ,  E,  I1,  I2, R, dN_I1I2 ) # assume that I1->I2 is when cases become detectable
  atests <- rmultinomial(sims,rep(tests,sims),out[,1:5]) # radomly draw the indivdiuals tested via asymptomatic tests
  #browser()
  # randomly draw the contacts from the different classes
  ifelse(CT==FALSE,contacts <- rmultinomial(sims,rep(rpois(sims,contacts)*(newSymptReported + apply(atests[,3:4], 1, sum)),sims),matrix(c(1,1,1,1,1),nr=sims,nc=5,byrow=T)*out[,1:5]), # tracing preferentially identifies infected individuals
         contacts <- rmultinomial(sims,rep(contacts*(newSymptReported + apply(atests[,3:4], 1, sum)),sims),log(out[,1:5]+1))) # tracing preferentially identifies infected individuals
  #browser()
  
  atests.isolate <- atests # holder for which tests will be positive that need to be isolated 
  atests.isolate[,c(1,2,5)] <- 0 # set non-infected classes to 0
  sympt.isolate <- matrix(0,nr=sims,nc=5) # stoarge for symoptomatic cases to isolate
  sympt.isolate[,4] <- newSymptReported # fill this in with the randomly drawn symptomatics
  out[,1:5] <- pmax(out[,1:5] - sympt.isolate - atests.isolate - contacts,0) # remove asymptomatic tests and contacts from population
  out <- cbind(out,atests, newSympt1, newSympt2, newSymptReported, contacts) # store all states -- SIR states plus tested, reported, contacts
}

