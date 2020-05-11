sir_step <- function (S, E, I1, I2, R, N, beta, theta,gamma_I1I2,gamma_I2R, delta.t, ...) {
# adapted from Aaron King's code
  
  dN_SE <- rbinom(n=1,size=S,prob=1-exp(-beta*(I1+I2)/N*delta.t))
  dN_EI1 <- rbinom(n=1,size=E,prob=1-exp(-theta*delta.t))
  dN_I1I2 <- rbinom(n=1,size=I1,prob=1-exp(-gamma_I1I2*delta.t))
  dN_I2R <- rbinom(n=1,size=I2,prob=1-exp(-gamma_I2R*delta.t))
  #browser()
  S <- S - dN_SE
  E <- E + dN_SE - dN_EI1
  I1 <- I1 + dN_EI1 - dN_I1I2
  I2 <- I2 + dN_I1I2 - dN_I2R
  R <- R + dN_I2R
  c(S ,  E,  I1,  I2, R, dN_I1I2 ) # assume that I1->I2 is when cases become detectable
}


############################################################################################################
############################################################################################################
############################################################################################################
T <- 100        # time to simulate over, we only care about start
sims <- 1000    # number of simulations
inf <- matrix(NA,T,sims)  # storage for infectious class
case <- matrix(NA,T,sims) # storage for daily cases
date.asypmt.detected.p <- numeric(sims)  # date at which asymptomatic cases are detected with probabilty P
#case.aseympt.detected.p<- numeric(sims) # cass on date at which asymptomatic cases are detected with probabilty P
Pdetect <- .9 # we want Pdetect certainty that we will have detected >=1 asymptomatic case
date.asympt.detected <- numeric(sims)  # date at which asymptomatic cases are detected (directly comparable with sympt)
date.sympt.detected <- numeric(sims) # date at which symptomatic cases
false_neg <- matrix(NA,T,sims)
false_pos <- matrix(NA,T,sims)

R0 <- 2   # R0
RE <- R0*.85  # RE assuming some fraction of population is already immune
beta <- R0 * (1/9)  # calculate Beta
theta <- 1/5  # 5 days from infection to infectious
gamma_I1I2 <- 1/2 # 2 days asymptomatic infectious
gamma_I2R <- 1/7 # 7 days infectious (this is probably too short)

ppn_sympt <- 1/20 # proportion of infectious cases that become symptomatic -- 1/20 is low for population, but might be high for student population
tests <- 1000           # assume this many asymptomatic tests per day
sensitivity <- 0.95  # test sensitivity
specificity <- 0.95 # test specificity


for(k in 1:sims){
  S <- round(50000*0.85)  # start with 50K * 0.85 susceptible
  E <- 0
  I1 <- 5                 # start with 5 asymptomatic infectious
  I2 <- 0
  R <- 50000 - S
  new_cases <- 0
  N <- S+E+I1+I2+R        # population size

  for(ts in 2:T){
  
  out <- sir_step(S[ts-1], E[ts-1], I1[ts-1], I2[ts-1], R[ts-1], N, beta, theta,gamma_I1I2,gamma_I2R, delta.t=1) # call to SIR step function above
  
  S <- c(S,out[1])  # update state
  E <- c(E,out[2])  # update state
  I1 <- c(I1,out[3])  # update state
  I2 <- c(I2,out[4])  # update state
  R <- c(R,out[5])  # update state
  new_cases <- c(new_cases,out[6])  # update state
}

cprod <- 1 - cumprod(dbinom(0,tests,(I1+I2) / N)) # cumulative probability of detecting at least 1 case from conducting asymptomatic tests at current prevalence 
sympt_cases <- c(0,0,rbinom(T,new_cases,ppn_sympt)) # assume all symptomatic cases are seen, but don't present for 2 days after symptom onset

inf[,k] <- I1+I2    # total infectious
case[,k] <- new_cases # daily cases
date.asypmt.detected.p[k] <- min(which(cprod>Pdetect))  # date at by which >=1 asymptomatic case is detected with probabiliyt Pdetect
#case.aseympt.detected.p[k] <- cumsum(case)[min(which(cprod>Pdetect))] # cases on date above

true_test_positives <- rbinom(T,tests,(I1+I2) / N)
date.asympt.detected[k] <- min(which(true_test_positives>0))  # first date on which an asymptomatic case is detected
date.sympt.detected[k] <- min(which(sympt_cases>0)) # first date on which a symptomatic case is detected

############################################################################
#False positives and False negatives
neg <- tests - true_test_positives
true_pos <- rbinom(T, true_test_positives, sensitivity) 
false_neg[,k] <- true_test_positives - true_pos

true_neg <- rbinom(T, neg, specificity)
false_pos[,k] <- neg - true_neg
############################################################################
}

layout(matrix(c(1,1,2,2),2,2,byrow=T))
par(mar=c(0,4,4,2)+0.1)
hh <- hist(date.asympt.detected,
           breaks=seq(0,100,by=5),col=rgb(0,0,0,.25),
           main="days until first case detected",axes=F,xlab=NA)
axis(2)
hist(date.sympt.detected,add=T,col=rgb(0,0,1,.25),breaks=seq(0,100,by=5))
legend(60,.5*max(hh$counts),pch=c(15,15),legend =c("asymptomatic testing","symptomatic testing"),col=c(rgb(0,0,0,.25),rgb(0,0,1,.25)),bty="n")

par(mar=c(5,4,0,2)+0.1)
matplot(case,type="l",col=rgb(1,0,0,.1),xlab="days from start")  # reported cases
matplot(inf,type="l",col=rgb(0,1,0,.05),add=T) # infectious class
rug(date.asympt.detected,col=rgb(0,0,0,.1),side = 3,lwd=2)
rug(date.sympt.detected,col=rgb(0,0,1,.1),side = 1,lwd=2)

# abline(v=date.asympt.detected,col=rgb(0,0,0,.25)) # dates of detection of 
# abline(v=date.sympt.detected,col=rgb(0,0,1,.25)) # dates of detection of 

legend(0,.75*max(case),lty=c(1,1),col=c(3,2),legend=c("total infected","new cases"),box.col="white")

par(mar=c(5, 4, 4, 2) + 0.1)
