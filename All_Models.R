

#------------------------------------------------------------
# MCMC FOR BETA, GAMMA | i, r (by Gibbs)
#------------------------------------------------------------

#for MCMC, we need to find the conditional posterior
#distributions for beta and gamma
#we calculate the joint posterior distributions in paper
#based on this, the conditional post. dists are


#---------------------------------------------------
#function sigma(i,r) to use for sampling beta
#summation j (1->n), k(1->n) {min(rj,ik)-min(ij,ik)}
#--------------------------------------------------

#derived from individual j can
#only being able to infect individual k for the length of time
#where both j is infective (ij < t < rj ) and k
#has not yet been infected (t < ik), summed over all j and k.

sigma<-function(i,r){
  ij_mat<- matrix(i, nrow=length(i), ncol=length(i), byrow=TRUE)
  ik_mat<- matrix(i, nrow=length(i), ncol=length(i), byrow=FALSE)
  rj_mat<- matrix(r, nrow=length(r), ncol=length(r), byrow=TRUE)

  term1<- pmin(rj_mat,ik_mat)
  term2<- pmin(ij_mat,ik_mat)

  return(sum(term1-term2))}

#dgamma x shape rate

#conditional distriburion of beta | i, r, gamma


#function to remove infinity values from i and r
clean<-function(x){
  return(x[is.finite(x)])
}

betac<-function(a,b,i,r,gamma,N){
  i=clean(i);r=clean(r)
  n<-length(i)
  sample<-rgamma(1,shape=n+a-1,rate=sigma(i,r)/N +b)
  return(sample)
}

#conditional distribution of gamma | i, r, beta

gammac<-function(a,b,i,r,beta,N){
  i=clean(i);r=clean(r)
  n<-length(i)
  sample<-rgamma(1,shape=n+a,rate=(sum(r-i)+b))
}

#Gibbs sampler for beta, gamma
#for SIR Epidemic Model

GibSIR<-function(N,M,i,r,beta_init,gamma_init){
  betav <- numeric(M)
  gammav <- numeric(M)
  betav[1] <- beta_init
  gammav[1] <- gamma_init
  Ba=1;Ga=1;
  Bb=10^(-4);Gb=10^(-4)

  for (k in 2:M) {
    #2:Generate βk directly from π(β|ik−1, r, γk−1).
    sb=betac(Ba,Bb,i,r,gammav[k-1],N)
    sg=gammac(Ga,Gb,i,r,sb,N)
    betav[k]=sb
    gammav[k]=sg
  }
  return(list("Beta"=betav,"Gamma"=gammav))
}


m=GibSIR(100,1000,sets$is,sets$rs,0.1,0.01)
#function to test the accuracy of our GibSIR

testg<-function(beta, gamma, N, M){
  finalsim<-getsetsSIR(N,1,0.2)
  finalsim<-sets
  N=100
  M=2000
  values<-GibSIR(N,M,finalsim$is, finalsim$rs,0,0)
  #lines(values$Beta)
  betaf<-tail(values$Beta,n=1)
  gammaf<-tail(values$Gamma,n=1)
  return(list("Your Beta:"=betaf, "Actual Beta:"=beta,
              "Your Gamma:"=gammaf, "Actual Gamma:"=gamma))
}


#------------------------------------------------------------
# MCMC FOR BETA, GAMMA, i | r (by Gibbs+M-H)
#------------------------------------------------------------

#we need a metropolis step for infection times unobserved as
#we cant directly sample from the conditional distribution


#π(i|r, β, γ) ∝  prod j(2->n){(I(ij ) − 1)}
#exp(-beta/N σ(i, r) − γ*sum(r-i))


infcon<-function(i,r,t,I,beta,gamma,N){
  i=clean(i);r=clean(r)
  n=length(i)

  matching_indices <- which(t %in% i)
  p=prod(I[matching_indices[-1]]-1)
  term=-(beta/N)*sigma(i,r)-gamma*sum(r-i)
  fin=p*exp(term)

  return(fin)
}

log_infcon <- function(i, r, t, I, beta, gamma, N) {
  i <- clean(i); r <- clean(r)
  n <- length(i)

  matching_indices <- which(t %in% i)

  log_p <- sum(log(I[matching_indices[-1]] - 1))
  term <- -(beta/N) * sigma(i,r) - gamma * sum(r - i)

  return(log_p+term)
}



#function f_X(x) where x is the infection duration
f_X <- function(x, gamma) {
  #calculate the exponential density
  return(gamma*exp(-gamma * x))}



MCMC_MH_GIB <- function(N, M,I, init_i, r,t, beta_init, gamma_init) {
  #vectors for beta, gamma, and infection times (i)
  #N=100
  #M=100
  #I=sets$I
  #init_i=sets$rs/2
  #r=sets$rs
  #beta_init=0
  #gamma_init=0
  betav <- numeric(M)
  gammav <- numeric(M)

  r <- clean(r)
  n <- length(r)

  main_i <- matrix(nrow = M, ncol = n)  #for infection times

  betav[1] <- beta_init
  gammav[1] <- gamma_init
  Ba <- 1; Ga <- 1
  Bb <- 10^(-4); Gb <- 10^(-4)

  i<-clean(init_i)  #initial infection times
  main_i[1,]<-i


  #Main MCMC loop
  for (k in 2:M) {
    #generate β_k from the posterior distribution π(β|i_k-1, r, γ_k-1)
    #k=2
    sb <- betac(Ba, Bb, i, r, gammav[k-1], N)
    sg <- gammac(Ga, Gb, i, r, sb, N)

    betav[k] <- sb
    gammav[k] <- sg

    #loop over each individual for the Metropolis-Hastings step
    for (j in 1:n) {
      #generate a new infection duration X from Exp(γ_k)
      #j=2
      x <- rexp(1, gammav[k])  #exponentially distributed with rate γ_k

      #set the proposed infection time i_star
      i_star <- main_i[k-1,]

      repeat {
        x <- rexp(1, gammav[k])
        i_proposed <- r[j] - x
        if (i_proposed >= 0) break
      }

      i_star[j] <- r[j] - x  #propose new infection time for individual j

      #compute acceptance ratio α
      likelihood_current <- infcon(i, r,t, I, betav[k], gammav[k], N)  # π(i|r, β, γ)
      likelihood_proposed <- infcon(i_star, r, t,I, betav[k], gammav[k], N)  # π(i_star|r, β, γ)
      likeratio=likelihood_proposed/likelihood_current

      #calculate f_X(r_j - i_j) and f_X(r_j - i_star_j)
      f_current <- f_X(r[j] - i[j], gammav[k])
      f_proposed <- f_X(r[j] - i_star[j], gammav[k])
      f_ratio=f_current/f_proposed

      #acceptance probability
      alpha <-min(1, likeratio*f_ratio)

      main_i[k,j ]<-main_i[k-1,j ]
      #accept or reject the proposal
      if (runif(1) < alpha) {
        #i[j] <- i_star[j]
        main_i[k, j] <- i_star[j]
      } else {
        main_i[k, j] <- i[j]
      }


    }
    # Update the working vector for infection times to the new iteration's state
    i <- main_i[k, ]



  }

  #get the results
  return(list("Beta" = betav, "Gamma" = gammav, "infection"=main_i))
}


#log-density for infection times i (conditional density)
loginf <- function(i, r, beta, gamma, n, N) {
  #xxclude the earliest infection (which is fixed)
  a <- which.min(i)
  sum <- 0
  for(j in 1:n) {
    if(j != a) {
      #the infectives count: positive counts before taking log.
      inf_val <- (sum(i<i[j]) - sum(r<i[j]))
      if(inf_val <= 0) {
        #return a very low likelihood if the count is zero.
        sum <- sum-1e6
      } else {
        sum <- sum+log(inf_val)
      }
    }
  }
  value <- sum - (beta*sigma(i,r)/N) - gamma*sum(r-i)
  return(value)
}


logMCMC_MH_GIB <- function(N, M, init_i, r, beta_init, gamma_init) {
  betav <- numeric(M)
  gammav <- numeric(M)
  r <- clean(r)
  n <- length(r)
  main_i <- matrix(nrow = M, ncol = n)

  #hyperparameters
  Ba <- 1; Bb <- 10
  Ga <- 1; Gb <- 10

  i <- clean(init_i)
  main_i[1, ] <- i
  betav[1] <- beta_init
  gammav[1] <- gamma_init

  for (k in 2:M) {
    i_current <- main_i[k - 1, ]
    #sample beta and gamma from their conditional posteriors using the latest i_current.
    sb <- betac(Ba, Bb, i_current, r, gammav[k-1], N)
    sg <- gammac(Ga, Gb, i_current, r, sb, N)
    betav[k] <- sb
    gammav[k] <- sg

    #Metropolis–Hastings updates for each infection time.
    for (j in 1:n) {
      #propose a new infection time for individual j, ensuring it's nonnegative.
      repeat {
        x <- rexp(1, gammav[k])
        i_proposed_j <- r[j] - x
        if (i_proposed_j >= 0) break
      }

      i_star <- i_current
      i_star[j] <- i_proposed_j

      #compute the log-density using the log_i function.
      log_lik_current <- loginf(i_current, r, betav[k], gammav[k], n, N)
      log_lik_proposed <- loginf(i_star, r, betav[k], gammav[k], n, N)

      #proposal density corrections (for an Exp(gamma) distribution)
      log_f_current <- log(f_X(r[j] - i_current[j], gammav[k]))
      log_f_proposed <- log(f_X(r[j] - i_star[j], gammav[k]))

      log_alpha <- min(0, log_lik_proposed + log_f_current - log_lik_current - log_f_proposed)

      if (!is.nan(log_alpha) && log(runif(1)) < log_alpha) {
        i_current[j] <- i_star[j]  # Accept the proposal.
      }
      #the updated infection time for individual j.
      main_i[k, j] <- i_current[j]
    }
    #update the working infection vector for the next iteration.
    i <- i_current
    #print("iter")
    #print(k)
  }

  return(list("Beta" = betav, "Gamma" = gammav, "infection" = main_i))
}


sets<-getsetsSIR(100,1,0.2)
sets2<-getsetsSIR(100,2,1)

#MCMC
results <- MCMC_MH_GIB(N = 100, M = 1000, I = sets$I, init_i = sets$rs/2, r = sets$rs,sets$ts, beta_init = 0, gamma_init = 0)
results4 <- MCMC_MH_GIB(N = 100, M = 1000, I = sets2$I, init_i = sets2$rs/2, r = sets2$rs,sets2$ts, beta_init = 0, gamma_init = 0)

results2 <- logMCMC_MH_GIB(N = 100, M = 1000, I = sets2$I, init_i = sets2$rs/2, r = sets2$rs,sets2$ts, beta_init = 0, gamma_init = 0)
results3 <- logMCMC_MH_GIB(N = 100, M = 2000, init_i = sets$rs/2, r = sets$rs, beta_init = 0, gamma_init = 0)

#function to remove burn-in

remove_burnin <- function(mcmc_matrix, burnin) {
  if (burnin >= nrow(mcmc_matrix)) {
    stop("Burn-in is larger than or equal to the number of MCMC samples.")
  }
  return(mcmc_matrix[(burnin + 1):nrow(mcmc_matrix), , drop = FALSE])
}


#PLOTTING ALL THE RESULTS TO VISUALIZE MIXING

#MCMC traces
individual_55_infection_times <- results3$infection[, 55]
beta_chain <- results3$Beta
gamma_chain <- results3$Gamma
burn_inf<-remove_burnin(results3$infection,0)



#plot infection times for individuals 30 and 60

par(mfrow = c(2, 1))

plot(burn_inf[, 30], type = "l",
     main = "Individual 30",
     xlab = "Iteration", ylab = "Infection Time")
grid(col = "grey", lty = 2,lwd=1) # Draw light grid
lines(burn_inf[, 30], col = "black", lwd = 1)  # Now draw the actual trace
abline(h = sets$is[30], col = "red", lwd = 3, lty=2)


plot(burn_inf[, 60], type = "l",
     main = "Individual 60",
     xlab = "Iteration", ylab = "Infection Time")
grid(col = "grey", lty = 2,lwd=1) # Draw light grid
lines(burn_inf[, 60], col = "black", lwd = 1)  # Now draw the actual trace
abline(h = sets$is[60], col = "red", lwd = 3, lty=2)


plot(burn_inf[, 59], type = "l",
     main = "Infection Time (Ind. 50)",
     xlab = "Iteration", ylab = "Infection Time")
abline(h = sets$is[59], col = "red", lwd = 2, lty=2)


beta_chain<-remove_burnin(beta_chain,50)
gamma_chain<-remove_burnin(gamma_chain,50)


par(mfrow = c(2, 1))

plot(beta_chain[1:length(beta_chain)], type = "l",
     main = expression("Trace plot of " * beta),
     xlab = "Iteration", ylab = expression(beta))
grid(col = "grey", lty = 2,lwd=1) # Draw light grid
lines(beta_chain, col = "black", lwd = 1)  # Now draw the actual trace

abline(h = 1, col = "red", lwd = 2,lty=2)

plot(gamma_chain[1:length(gamma_chain)], type = "l",
     main = expression("Trace plot of " * gamma),
     xlab = "Iteration", ylab = expression(gamma))
grid(col = "grey", lty = 2,lwd=1) # Draw light grid
lines(gamma_chain, col = "black", lwd = 1)  # Now draw the actual trace

abline(h = 0.2, col = "red", lwd = 2,lty=2)

#legend("topright", legend = "True value", col = "darkred", lty = 2, lwd = 1.5, bty = "n", cex = 0.8)








