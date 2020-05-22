##########################################################################################
# Check that the order of the rows in the EVPI table are correct
##########################################################################################

require(shiny)
#require(fields)
#require(viridis)

#load stuff
source("sir_step_vec.R")

shinyServer(function(input,output) {
pop_size<-reactive({input$pop_size})
R0<-reactive({input$R0})
ppn_immune<-reactive({input$ppn_immune})
tests<-reactive({input$tests})
care_seeking_delay<-reactive({input$care_seeking_delay})
ppn_sympt<-reactive({input$ppn_sympt})
sensitivity<-reactive({input$sensitivity})
specificity<-reactive({input$specificity})


output$Time_to_detection_plot <- renderPlot({
  ############################################################################################################
  T <- 100        # time to simulate over, we only care about start (longer will do full outbreak, and take longer)
  sims <- 1000    # number of simulations, can set this as a slider so that things run faster while testing scenarios
  inf <- matrix(NA,T,sims)  # storage for infectious class
  case <- matrix(NA,T,sims) # storage for daily cases
  date.asypmt.detected.p <- numeric(sims)  # date at which asymptomatic cases are detected with probabilty P
  #case.aseympt.detected.p<- numeric(sims) # cass on date at which asymptomatic cases are detected with probabilty P
  #Pdetect <- .9 # we want Pdetect certainty that we will have detected >=1 asymptomatic case
  date.asympt.detected <- numeric(sims)  # date at which asymptomatic cases are detected (directly comparable with sympt)
  date.sympt.detected <- numeric(sims) # date at which symptomatic cases
  false_neg <- matrix(NA,T,sims)
  false_pos <- matrix(NA,T,sims)
  
  #R0 <- 2   # R0
  RE <- R0()*(1-ppn_immune())  # RE assuming some fraction of population is already immune
  beta <- RE * (1/9)  # calculate Beta
  theta <- 1/5  # 5 days from infection to infectious
  gamma_I1I2 <- 1/2 # 2 days asymptomatic infectious
  gamma_I2R <- 1/7 # 7 days infectious (this is probably too short)
  
  # ppn_sympt <- 1/20 # proportion of infectious cases that become symptomatic -- 1/20 is low for population, but might be high for student population
  # care_seeking_delay <- 2 # number of days before seeking care after onset of symptoms
  # tests <- 300           # assume this many asymptomatic tests per day
  # sensitivity <- 0.99  # test sensitivity
  # specificity <- 0.99 # test specificity
  
  
  #
  S <- matrix(round(pop_size()*(1-ppn_immune())),1,sims)  # start with 50K * 0.85 susceptible
  E <- matrix(0,1,sims)
  I1 <- matrix(5,1,sims)                 # start with 5 asymptomatic infectious
  I2 <- matrix(0,1,sims)
  R <- matrix(pop_size() - S,1,sims)
  new_cases <- matrix(0,1,sims)
  N <- S+E+I1+I2+R        # population size
  
  for(ts in 2:T){
    
    out <- sir_step(sims, S[ts-1,], E[ts-1,], I1[ts-1,], I2[ts-1,], R[ts-1,], N[ts-1,], beta, theta,gamma_I1I2,gamma_I2R, delta.t) # call to SIR step function above
    #browser()
    S <- rbind(S,out[,1])  # update state
    E <- rbind(E,out[,2])  # update state
    I1 <- rbind(I1,out[,3])  # update state
    I2 <- rbind(I2,out[,4])  # update state
    R <- rbind(R,out[,5])  # update state
    N <- rbind(N,N[ts-1,])  # update state
    new_cases <- rbind(new_cases,out[,6])  # update state
  }
  #browser()
  inf <- I1+I2   # total infectious
  case <- new_cases # daily cases
  
  cprod <- 1 - apply(dbinom(0,tests(),inf / N),2,cumprod) # 100x1000 cumulative probability of detecting at least 1 case from conducting asymptomatic tests at current prevalence 
  sympt_cases <- rbind(matrix(0,care_seeking_delay(),sims,byrow=T),matrix(rbinom(T*sims,new_cases,ppn_sympt()),T,sims)) # assume all symptomatic cases are seen, but don't present for 2 days after symptom onset
  
  #date.asypmt.detected.p <- apply(cprod<=Pdetect,2,sum)  # date at by which >=1 asymptomatic case is detected with probabiliyt Pdetect
  #case.aseympt.detected.p[k] <- cumsum(case)[min(which(cprod>Pdetect))] # cases on date above
  
  true_test_positives <- matrix(rbinom(T*sims,tests(),inf / N),T,sims)
  date.asympt.detected <- apply(apply(true_test_positives,2,cumsum)==0,2,sum)  # first date on which an asymptomatic case is detected
  date.sympt.detected <- apply(apply(sympt_cases,2,cumsum)==0,2,sum) # first date on which a symptomatic case is detected
  
  ############################################################################
  #False positives and False negatives
  neg <- tests() - true_test_positives
  true_pos <- matrix(rbinom(T*sims, true_test_positives, sensitivity()),T,sims) 
  false_neg <- true_test_positives - true_pos
  
  true_neg <- matrix(rbinom(T*sims, neg, specificity()),T,sims)
  false_pos <- neg - true_neg
  ############################################################################
  layout(matrix(c(1,1,2,2),2,2,byrow=T))
  par(mar=c(0,4,4,2)+0.1)
  hh <- hist(pmin(date.asympt.detected,100),
             breaks=seq(0,100,by=5),col=rgb(0,0,0,.25),
             main="days until first case detected",axes=F,xlab=NA,)
  axis(2)
  hist(pmin(date.sympt.detected,100),add=T,col=rgb(0,.65,.5,.25),breaks=seq(0,100,by=5))
  legend(60,.5*max(hh$counts),pch=c(15,15),legend =c("asymptomatic testing","symptomatic testing"),col=c(rgb(0,0,0,.25),rgb(0,.65,.5,.25)),bty="n")
  
  par(mar=c(5,4,0,2)+0.1)
  matplot(case,type="l",col=rgb(.9,.6,0,.1),xlab="days from start")  # reported cases
  matplot(inf,type="l",col=rgb(.35,.7,.9,.1),add=T) # infectious class
  rug(date.asympt.detected,col=rgb(0,0,0,.1),side = 3,lwd=2)
  rug(date.sympt.detected,col=rgb(0,.65,.5,.1),side = 1,lwd=2)
  
  # abline(v=date.asympt.detected,col=rgb(0,0,0,.25)) # dates of detection of 
  # abline(v=date.sympt.detected,col=rgb(0,0,1,.25)) # dates of detection of 
  
  legend(0,.75*max(case),lty=c(1,1),col=c(rgb(.35,.7,.9,1),rgb(.9,.6,0,1)),legend=c("total infected","new cases"),box.col="white")
  
  par(mar=c(5, 4, 4, 2) + 0.1)
})

output$False_positives_plot <- renderPlot({
  par(mfrow=c(2,1))
  plot(apply(false_neg,1,mean),pch=19,col=1,ylim=c(0,50),ylab="false negatives per day", xlab="day")
  plot(apply(false_pos,1,mean),pch=19,col=1,ylim=c(0,50),ylab="false positives per day", xlab="day")
    })



})
