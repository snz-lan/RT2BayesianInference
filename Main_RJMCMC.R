valid_path2 <- function(inf, rec, N, T_max) {
  inf <- sort(inf)
  rec <- sort(rec)

  #trim to [0,T_max]
  inf <- inf[inf >= 0 & inf <= T_max]
  rec <- rec[rec >= 0 & rec <= T_max]

  #index‐wise: jth infection < jth recovery
  k <- min(length(inf), length(rec))
  if (k>0 && any(inf[1:k] >= rec[1:k])) return(FALSE)

  #aggregate: at each event time, #inf >= #rec
  events <- sort(unique(c(inf, rec)))
  for (t in events) {
    if (sum(inf<t) < sum(rec<t))
      return(FALSE)
  }

  TRUE
}


log_density2 <- function(infections, recoveries, beta, gamma, N, T_obs) {
  infections <- clean(infections)
  recoveries <- clean(recoveries)

  #build the X,Y path
  paths <- get_XY_paths(infections,recoveries,N)
  X_t   <- paths$X
  Y_t   <- paths$Y
  times <- paths$times

  #recovery term
  rec_mask <- times %in% recoveries
  #only sum when Y>0
  Ys_rec   <- Y_t[rec_mask]
  Ys_rec   <- Ys_rec[Ys_rec > 0]
  log_rec  <- length(Ys_rec) * log(gamma) + sum(log(Ys_rec))

  #infection term (skip the first “fixed” infection)
  inf_term <- 0
  if (length(infections) > 1) {
    inf_mask <- times %in% infections
    X_inf    <- X_t[inf_mask][-1]
    Y_inf    <- Y_t[inf_mask][-1]
    ok       <- (X_inf > 0 & Y_inf > 0)
    m        <- sum(ok)
    # each adds log(beta/N)+log X + log Y
    inf_term <- m * (log(beta) - log(N)) +
      sum(log(X_inf[ok]) + log(Y_inf[ok]))
  }

  #integral term via trapezoid rule, with β/N
  dt  <- diff(times)
  bXY <- (beta / N) * (X_t[-length(X_t)] * Y_t[-length(Y_t)] +
                         X_t[-1]          * Y_t[-1]) / 2
  gY  <- gamma * (Y_t[-length(Y_t)] + Y_t[-1]) / 2
  int_term <- sum(dt * (bXY + gY))

  #account for [last event, T_obs]
  t_last <- times[length(times)]
  if (t_last < T_obs) {
    delta_t <- T_obs - t_last
    bXY_last <- (beta / N) * X_t[length(X_t)] * Y_t[length(Y_t)]
    gY_last  <- gamma * Y_t[length(Y_t)]
    int_term <- int_term + delta_t * (bXY_last + gY_last)
  }
  #print("log_rec")
  #print(log_rec)
  #print("inf_term")
  #print(inf_term)
  #print("int_term")
  #print(int_term)

  return(log_rec + inf_term - int_term)
}


#this is to get the number of infectives and susceptibles right before
#the time of event
#KEEP IN MIND HERE THAT X AND Y WILL
#HAVE THE SAME LENGTH..IE..NO: OF EVENTS THAT HAPPENED UNTILL
#LAST OBSERVED RECOVERY EVENT
#DURING THE CURRENT BAYESIAN UPDATE

get_XY_paths <- function(infections, recoveries, N) {
  #recoveries<-rec
  #infections<-rec/2
  infections<-clean(infections)
  recoveries<-clean(recoveries)
  #infections<-recoveries/2

  #recoveries<-c(recoveries,6.9,7.3,9.9,9.0)
  #recoveries<-sort(recoveries)
  T_max<-max(recoveries)
  infections  <- infections[ infections <= T_max ]
  recoveries  <- recoveries[ recoveries <= T_max ]
  #events <- sort(unique(c(infs, recs)))
  events <- sort(unique(c(infections, recoveries)))
  X_vals <- integer(length(events))
  Y_vals <- integer(length(events))

  for (i in seq_along(events)) {
    #cat(i)
    t <- events[i]
    #infections and recoveries BEFORE time t (i.e., t^-)
    num_infected_before <- sum(infections < t)
    num_recovered_before <- sum(recoveries < t)

    X_vals[i] <- N - num_infected_before
    Y_vals[i] <- num_infected_before - num_recovered_before
  }

  return(list(times = events, X = X_vals, Y = Y_vals))
}



log_density <- function(infections, recoveries, beta, gamma, N) {

  #beta=1
  #gamma=0.2
  #infections<-prop_inf#mysim$is
  #recoveries<-rec#mysim$rs[-length(mysim$rs)]

  infections<-clean(infections)
  recoveries<-clean(recoveries)
  N<-100
  paths <- get_XY_paths(infections, recoveries, N)
  X_t <- paths$X
  Y_t <- paths$Y
  times <- paths$times
  T_max<-tail(times, n=1)

  recovery_mask <- times %in% recoveries
  valid_recovery_Y <- Y_t[recovery_mask] > 0
  log_recovery <- sum(log(gamma) + log(Y_t[recovery_mask][valid_recovery_Y]))

  # Terms
  #log_recovery <- sum(log(gamma) + log(Y_t[times %in% recoveries])

  log_infection <- 0
  if (length(infections) > 1) {
    #I FIXED IT
    #I FIXEDDD ITTTTT!!!
    #infection_mask <- times %in% infections
    infection_mask <- times %in% infections
    infect_X <- X_t[infection_mask][-1]
    infect_Y <- Y_t[infection_mask][-1]
    valid_infect <- infect_X > 0 & infect_Y > 0
    log_infection <- sum(log(beta) + log(infect_X[valid_infect]) + log(infect_Y[valid_infect]))
  }

  # Integral term using rectangular approximation
  integral_term <- 0
  # for (i in 1:(length(times) - 1)) {
  #   delta_t <- times[i + 1] - times[i]
  #   integral_term <- integral_term +
  #     delta_t * (beta * X_t[i] * Y_t[i] + gamma * Y_t[i])
  # }
  for (i in 1:(length(times) - 1)) {
    delta_t <- times[i + 1] - times[i]
    beta_term <- beta * (X_t[i] * Y_t[i] + X_t[i+1] * Y_t[i+1]) / 2
    gamma_term <- gamma * (Y_t[i] + Y_t[i+1]) / 2
    integral_term <- integral_term + delta_t * (beta_term + gamma_term)
  }


  # Add final bit if last event < T_max
  if (tail(times, 1) < T_max) {
    delta_t <- T_max - tail(times, 1)
    integral_term <- integral_term +
      delta_t * (beta * tail(X_t, 1) * tail(Y_t, 1) + gamma * tail(Y_t, 1))
  }
  # cat("log_recovery =", log_recovery, "\n")
  # cat("log_infection =", log_infection, "\n")
  # cat("integral_term =", integral_term, "\n")
  # cat("log_density =", log_recovery + log_infection - integral_term, "\n")
  #
  return(log_recovery + log_infection - integral_term)
}

update_beta <- function(paths, lambda_beta, nu_beta) {
  t   <- paths$times
  X   <- paths$X
  Y   <- paths$Y
  T0  <- max(t)                      # length of observation window

  # non-dimensional ∫0^1 X(u)Y(u) du
  dt_star <- diff(t) / T0
  int_XY  <- sum(dt_star * (X[-length(X)] * Y[-length(Y)]))

  m       <- length(t) - 1           # total infections
  #m <- sum(t %in% infections) - 1

  shape   <- lambda_beta + int_XY
  rate    <- (m - 1) + nu_beta
  beta_star <- rgamma(1, shape = shape, rate = rate)

  # convert back into real‐time units:
  return(beta_star/T0)
}

update_gamma <- function(paths, lambda_gamma, nu_gamma) {
  t   <- paths$times
  Y   <- paths$Y
  T0  <- max(t)

  # non-dimensional ∫0^1 Y(u) du
  dt_star <- diff(t) / T0
  int_Y   <- sum(dt_star * Y[-length(Y)])

  n       <- length(t) - 1           # total removals
  shape   <- lambda_gamma + int_Y
  rate    <- n + nu_gamma
  gamma_star <- rgamma(1, shape = shape, rate = rate)

  return(gamma_star/T0)
}


rjmcmc_diagnostics_adapt_old <- function(
    rec, N, T_max, z,M=2000,
    L_move    = 20,
    p_dim     = 0.3,
    prob,
    lambda_beta = 1, nu_beta = 1e-4,
    lambda_gamma= 1, nu_gamma= 1e-4,
    target_acc = 0.3,   # target local‐move acceptance
    adapt_rate = 0.05   # how strongly to adapt w
){
  rec <- sort(clean(rec))
  inf <- rec/2
  inf[1]<-0
  beta <- gamma <- 1
  nob<-length(rec)

  #initial RW window:
  w <- T_max/10

  # storage
  beta_chain  <- numeric(M)
  gamma_chain <- numeric(M)
  m_chain     <- integer(M)
  I_chain     <- vector("list",M)
  i20_chain   <- numeric(M)  # infection #20
  add_prop <- add_acc <- 0L
  rem_prop <- rem_acc <- 0L

  for(k in seq_len(M)){
    ##Gibbs‐update beta gamma
    paths <- get_XY_paths(inf, rec, N)
    beta  <- update_beta (paths, lambda_beta, nu_beta)
    gamma <- update_gamma(paths, lambda_gamma,nu_gamma)

    ##local 1‐D updates (adapt w)
    local_prop <- local_acc <- 0L
    #if(runif(1) > p_dim){
    for(l in seq_len(L_move)){


      j <- sample(2:length(inf),1)

      repeat {
        y <- inf[j] + rnorm(1, 0, w)
        #print("w")
        #print(w)
        if (y > 0 && y < T_max) break
      }
      prop <- inf
      prop[j] <- y
      prop<-sort(prop)            #pmax(0, pmin(T_max, prop)))
      if(!valid_path2(prop, rec, N, T_max)) next
      local_prop <- local_prop + 1L
      log0 <- log_density2(inf,  rec, beta, gamma, N,T_max)
      log1 <- log_density2(prop, rec, beta, gamma, N,T_max)
      if(log(runif(1)) < (log1-log0)){
        inf <- prop
        local_acc <- local_acc + 1L
      }

    }
    if(local_prop > 0){
      acc_rate <- local_acc / local_prop
      # adapt window up/down
      w <- w * exp(adapt_rate * (acc_rate - target_acc))
      #w <- w * exp(adapt_rate * (acc_rate - target_acc))
      w <- max(w, 1e-3)          #don’t let it collapse to 0
      w <- min(w, T_max / 2)     #don’t let it blow up larger than time horizon
      #w <- min(w, T_max/2)
      if(acc_rate < 0.3) w <- w*0.5
    }
    #}
    ## Add/Remove

    if(runif(1)<p_dim){
      #else{
      m <- length(inf)
      if(runif(1) < prob && m > 1){
        ## REMOVE move
        rem_prop <- rem_prop + 1L
        l=0
        repeat{
          l=l+1
          j   <- sample(2:m,1)
          prop<- inf[-j]
          if(valid_path2(prop[-j], rec, N, T_max) || l>100 ) break}

        #}
        if(valid_path2(prop, rec, N, T_max) && length(prop)>=length(rec)){

          if (length(prop) < 1){next}
          log0 <- log_density2(inf,  rec, beta, gamma, N,T_max)
          log1 <- log_density2(prop, rec, beta, gamma, N,T_max)
          log_alpha <- (log1 - log0)+log(length(inf)) - log(T_max)#+ (log_q_bwd - log_q_fwd)
          #print("log alpha remove")
          #print(log_alpha)
          #print(log_q_bwd)
          if(log(runif(1)) < log_alpha){
            inf<-prop
            rem_acc <- rem_acc+1L
          }
        }

      }
      else{

        # add_prop <- add_prop + 1L
        # j_rec    <- sample(seq_along(rec), 1)
        # x_new    <- rexp(1, rate = gamma)
        # i_new    <- rec[j_rec] - x_new
        # prop_inf <- sort(c(inf, i_new))
        #
        # if (i_new >= 0 && valid_path2(prop_inf, rec, N, T_max)) {
        #   # forward proposal density
        #   log_q_fwd <- -log(length(rec)) + (log(gamma) - gamma * x_new)
        #   # backward (would remove this one) density
        #   log_q_bwd <- -log(length(inf) + 1L)
        #
        #   log0 <- log_density(inf,  rec, beta, gamma, N)
        #   log1 <- log_density(prop_inf, rec, beta, gamma, N)
        #   log_alpha <- (log1 - log0) + (log_q_bwd - log_q_fwd)
        #
        #   if (log(runif(1)) < log_alpha) {
        #     inf     <- prop_inf
        #     add_acc <- add_acc + 1L
        #   }
        # }

        ## ADD move
        add_prop <- add_prop + 1L
        i_new <- runif(1,inf[2],T_max)
        prop  <- sort(c(inf, i_new))

        #if(length(prop)>=N){next}

        if(i_new > 0 && valid_path2(prop, rec, N, T_max)&& length(prop)<=N){
          #log_q_fwd <- -log(length(rec)) + (log(gamma) - gamma * x)
          #log_q_bwd <- -log(m + 1)
          log0 <- log_density2(inf,  rec, beta, gamma, N,T_max)
          log1 <- log_density2(prop, rec, beta, gamma, N,T_max)
          delta_ll <- log1 - log0
          log_alpha <- (log1-log0) +log(T_max)-log(length(inf)+1)
          if(log(runif(1)) < log_alpha){
            inf <- prop
            add_acc <- add_acc + 1L
          }
        }

      }
    }

    #record
    beta_chain[k]  <- beta
    gamma_chain[k] <- gamma
    m_chain[k]     <- length(inf)
    I_chain[[k]]   <- inf
    i20_chain[k]   <- if(length(inf)>=z) inf[z] else NA


    #every 100 its, print summary
    if(k %% 100 == 0){

      cat(sprintf(
        "[%4d] w=%.4f  local acc=%.2f  add=%.2f  rem=%.2f\n",
        k, w,
        if(local_prop>0) local_acc/local_prop else NA,
        if(add_prop>0) add_acc/add_prop else NA,
        if(rem_prop>0) rem_acc/rem_prop else NA
      ))
    }
  }

  library(coda)
  list(
    beta_mcmc  = mcmc(beta_chain),
    gamma_mcmc = mcmc(gamma_chain),
    m_mcmc     = mcmc(m_chain),
    infections = I_chain,
    i20_mcmc   = mcmc(na.omit(i20_chain))
  )
}


every<-function(rec, M, z, sets,b,g){

  out <- rjmcmc_diagnostics_adapt_old(rec, N=100, T_max=max(rec),M=M,z=z,L_move = 15, p_dim=0.4, prob = 0.5)
  m_vals <- as.numeric(out$m_mcmc)
  plot(m_vals, type = "l", main = "Trace plot w/ true vector size")
  abline(h = sum(sort(sets$is) <= max(rec)), col = "red", lwd = 2)

  #plot(out$m_mcmc);   #effectiveSize(res$m_mcmc)
  ploty(0,samples=out,b,g)

  burnin=0
  mcmc_vals <- as.numeric(out$i20_mcmc)
  mcmc_valsb <- mcmc_vals[(burnin+1):length(mcmc_vals)]
  plot(mcmc_valsb, type = "l", main = "Trace plot w/ zth true inf value")
  abline(h = sort(sets$is)[z], col = "red", lwd = 2)
  return(out)
}
