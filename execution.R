
#helper function to remove burnin
remove_burnin <- function(mcmc_matrix, burnin) {
  if (burnin >= nrow(mcmc_matrix)) {
    stop("Burn-in is larger than or equal to the number of MCMC samples.")
  }
  return(mcmc_matrix[(burnin + 1):nrow(mcmc_matrix), , drop = FALSE])
}

#get data
sets<-getsetsSIR(100,1,0.2)


#########################################################
#########################################################
#########################################################
############GIBBS-SAMPLER################################


#get results for gibbs sampler
results1=GibSIR(100,1000,sets$is,sets$rs,0.1,0.01)


#betagammaplotting for gibbs sampler
beta_chain1 <- results1$Beta
gamma_chain1 <- results1$Gamma

par(mfrow = c(2, 1))

plot(beta_chain1[1:length(beta_chain1)], type = "l",
     main = expression("Trace plot of " * beta),
     xlab = "Iteration", ylab = expression(beta))
grid(col = "grey", lty = 2,lwd=1) # Draw light grid
lines(beta_chain1, col = "black", lwd = 1)  # Now draw the actual trace

abline(h = 1, col = "red", lwd = 2,lty=2)

plot(gamma_chain1[1:length(gamma_chain1)], type = "l",
     main = expression("Trace plot of " * gamma),
     xlab = "Iteration", ylab = expression(gamma))
grid(col = "grey", lty = 2,lwd=1) # Draw light grid
lines(gamma_chain1, col = "black", lwd = 1)  # Now draw the actual trace

abline(h = 0.2, col = "red", lwd = 2,lty=2)


#########################################################
#########################################################
#########################################################
############METROPOLIS-GIBBS#############################



#metropolis-gibbs results
results3 <- logMCMC_MH_GIB(N = 100, M = 2000, init_i = sets$rs/2, r = sets$rs, beta_init = 0, gamma_init = 0)

#betagammaplotting for m-h
beta_chain <- results3$Beta
gamma_chain <- results3$Gamma

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


#plot infection times for individuals 30 and 60 with 0 burnin
burn_inf<-remove_burnin(results3$infection,0)

par(mfrow = c(2, 1))

plot(burn_inf[, 30], type = "l",
     main = "Individual 30",
     xlab = "Iteration", ylab = "Infection Time")
grid(col = "grey", lty = 2,lwd=1)
lines(burn_inf[, 30], col = "black", lwd = 1)
abline(h = sets$is[30], col = "red", lwd = 3, lty=2)


plot(burn_inf[, 60], type = "l",
     main = "Individual 60",
     xlab = "Iteration", ylab = "Infection Time")
grid(col = "grey", lty = 2,lwd=1) # Draw light grid
lines(burn_inf[, 60], col = "black", lwd = 1)  # Now draw the actual trace
abline(h = sets$is[60], col = "red", lwd = 3, lty=2)


#########################################################
#########################################################
#########################################################
############RJMCMC#######################################

#function to get the first half of a vector
gethalf<-function(vec){
  return(clean((vec)[1:(length(vec)/2)]))
}

beta2=0.9
gamma2=0.2
sets2<-getsetsSIR(100,beta2,gamma2)

#partial recovery vector
prec<-gethalf(sets2$rs)

#get rjmcmc sampler results
results4<-rjmcmc_diagnostics_adapt_old(prec, N=100, T_max=max(prec),M=5000,z=55,L_move = 15, p_dim=0.4, prob = 0.5)



par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))  # 1 row, 2 plots; tidy margins

# trace plot for m (number of infections) ---
m_vals <- as.numeric(results4$m_mcmc)
plot(m_vals, type = "l", col = "black", lwd = 1,
     main = "Trace plot of inf size",
     xlab = "Iteration", ylab = expression(m))
grid(col = "grey", lty = 2,lwd=1) # Draw light grid
lines(m_vals, col = "black", lwd = 1)  # Now draw the actual trace

abline(h = sum(sort(sets2$is) <= max(prec)),
       col = "red", lwd = 2, lty = 2)

#beta-gamma plots for rjmcmc
beta_chain3 <- as.numeric(results4$beta_mcmc)
gamma_chain3 <- as.numeric(results4$gamma_mcmc)

par(mfrow = c(2,1), mar = c(4, 4, 2, 1))  # Adjust margins for better spacing
plot(beta_chain3[1:5000], type = "l", col = "black", lwd = 0.6,
     main = expression("Trace plot of " * beta),
     xlab = "Iteration", ylab = expression(beta))
grid(col = "grey", lty = 2) # Draw light grid
lines(beta_chain3, col = "black", lwd = 0.6)
abline(h = 0.9, col = "red", lwd = 3, lty=2)

# --- Plot for gamma ---
plot(gamma_chain3[1:5000], type = "l", col = "black", lwd = 1,
     main = expression("Trace plot of " * gamma),
     xlab = "Iteration", ylab = expression(gamma))
grid(col = "grey", lty = 2)
lines(gamma_chain3, col = "black", lwd = 0.3)
abline(h = 0.2, col = "red", lwd = 3)
