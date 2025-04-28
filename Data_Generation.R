#simulate from the SIR model using Gillespie’s algorithm,
#originally formulated in Gillespie (1976)

#We're using an improvement of the original algorithm where
#we have to sample only once by treating the min of 2 exp var
#as an exp variable

#Inputs: Population size, transmission rate beta, recovery
#rate gamma
#Outputs: Filled vectors S,I and R (no: of these at each time step)

SIR<-function(N,beta, gamma ){

  tmax=1000;
  T=rep(NA,tmax);S=rep(NA,tmax);I=rep(NA,tmax);R=rep(NA,tmax);
  S[1]=N-1;I[1]=1;R[1]=0;T[1]=0;
  k=1;

  while(I[k]!=0 && k<tmax){

    r_inf<-beta*S[k]*I[k]/N;
    r_rec<-gamma*I[k];
    tot<-r_inf+r_rec;

    dt<-rexp(1,tot);
    T[k+1]=T[k]+dt;
    u<-runif(1);

    if(u<r_inf/tot){
      S[k+1]=S[k]-1;
      I[k+1]=I[k]+1;
      R[k+1]=R[k];
    }

    else{
      S[k+1]=S[k];
      I[k+1]=I[k]-1;
      R[k+1]=R[k]+1;
    }
    k=k+1
  }
  return(list(time=T[1:k], S=S[1:k], I=I[1:k], R=R[1:k]))
}


#--------------------------------------------------------
#--------------------------------------------------------
#for MCMC sampling and the likelihood function
#we need the times of recovery and infection as well
#the following function implements the same algorithm
#as above but also returns the {i}, {r} vectors that give infection
#and recovery times

SIR2 <- function(N, beta, gamma) {

  tmax <- 1000

  T<-rep(NA, tmax); S<-rep(NA, tmax)
  I<-rep(NA, tmax);R<-rep(NA, tmax)

  T[1]<-0; S[1]<- N-1
  I[1]<-1; R[1] <-0

  k <-1

  infected <-rep(Inf,N)
  recovered <-rep(Inf,N)
  infected[1] <-0
  recovered[1] <-rexp(1,gamma)

  next_id <- 2
  I_set <- c(1)

  events <- data.frame(
    time = 0,
    type = "infection",
    person = 1,
    stringsAsFactors = FALSE
  )

  while (I[k]>0 && k<tmax) {

    #next possible infection time

    rate_inf <- beta*S[k]*I[k]/N

    if (rate_inf > 0) {
      next_inf_time <- T[k]+rexp(1,rate_inf)
    } else{
      next_inf_time <- Inf
    }

    #next recovery: find the earliest recovery among infected
    next_rec_times <- recovered[I_set]
    next_rec_time <- min(next_rec_times)
    next_event_time <- min(next_inf_time,next_rec_time)

    if (next_event_time==Inf) break  #end of epidemic

    #updating time and state
    T[k+1] <- next_event_time

    if (next_inf_time<next_rec_time) {
      #infection event
      if (next_id>N) break  # no one left to infect

      S[k+1] <- S[k]-1
      I[k+1] <- I[k]+1
      R[k+1] <- R[k]

      infected[next_id] <- T[k+1]
      recovered[next_id] <- T[k+1]+rexp(1,gamma)
      I_set <- c(I_set,next_id)

      events <- rbind(events, data.frame(
        time = T[k+1], type="infection", person=next_id))

      next_id <- next_id+1

    } else {

      #recovery event

      recovering <- I_set[which.min(next_rec_times)]

      S[k+1] <- S[k]
      I[k+1] <- I[k]-1
      R[k+1] <- R[k]+1

      I_set <- setdiff(I_set,recovering)

      events <- rbind(events, data.frame(
        time = T[k+1], type = "recovery", person=recovering))
    }

    k <- k+1
  }

  return(list(time = T[1:k],S= S[1:k],I= I[1:k],R= R[1:k],
              i=infected,
              r=recovered,
              events=events))}

#------------------------------------------------------------
#------------------------------------------------------------

#let i be the set of infection event times, r be the set of recovery event times, and t be the set of
#both infection and recovery event times- order the recovery times such that r1 ≤ r2 ≤ · · · ≤ rN ,
#and let ij be the infection time for the individual who recovers at time rj- start the model at the
#moment of the first infection, so i1 = 0.
#If n < N, then not all individuals become infective, so let
#in+1, . . . , iN , rn+1, . . . , rN = ∞



sim<-SIR2(100,2,0.6)


#get ordered recovery set, corresponding infection set
#and time of all events
getsetsSIR<-function(N,beta,gamma){

  sim<-SIR2(N,beta,gamma)
  #get i and r vectors
  i<- sim$i
  r<- sim$r

  # i1 = 0 and r is sorted
  ord<-order(r)
  rs<-r[ord]
  is<-i[ord]

  #all event times
  ts<-sort(unique(c(i, r)))
  return(list("is"=is,"rs"=rs,"ts"=ts,"S"=sim$S,"I"=sim$I,"R"=sim$R))}



mysim <- SIR2(N = 100, beta=1, gamma = 0.2)
plot(mysim$time, mysim$I, type='s', col='red', ylab="count", xlab="Time")
lines(mysim$time,mysim$S, type='s', col='blue')
lines(mysim$time, mysim$R, type='s', col='green')
legend("right", legend=c("S","I","R"), col = c("blue", "red", "green"), lty=1)
