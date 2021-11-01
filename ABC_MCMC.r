
prior_assort = function(x, sigma = 1, a = 0.05, b = 20){ # Uniform prior for assortativity
  return(dunif(x, min = a, max = b))
}

prior_b = function(x, mu = 1, sigma = 3){
  return(dnorm(x, mu, sigma))
}

prior = function(theta, calibrate_assort = F){
  if(calibrate_assort){
    assort = theta[2]
    b = theta[1]
    return(prior_assort(x = assort) * prior_b(x = b))
  }
  else{
    return(prior_b(x = theta))
  }
}

perform_sim = function(){
  assort = runif(1, min = 0.05, max = 20)
  b = rnorm(1, mean = 1.2, sd = 3)
  while(b < 0){
    b = rnorm(1, mean = 1.2, sd = 3)
  }
  error = compare_cases(c(assort, b), N_observed_infected_A, N_observed_infected_B)
  return(c(assort, b, error))
  
}

compare_cases = function(theta_p, 
                         N_obs_inf_A,
                         N_obs_inf_B){
  b = theta_p[1]
  assort = theta_p[2]
  
  
  ### Simulate cases with input parameters ###
  res <- run_4x4(matrix(c(assort,1,1, assort), nrow=2),
                 matrix(c(1,1, 1, 1), nrow=2),
                 1.3,
                 c(1, 1.2, 1 * b, 1.2 * b),
                 rep(1, 4), n=1)
  Total_infected = res$R  
  
  N_sim_inf_A = Total_infected[1] + Total_infected[2]
  N_sim_inf_B = Total_infected[3] + Total_infected[4]
  
  
  Dist_A = abs(N_obs_inf_A - N_sim_inf_A)/90000
  Dist_B = abs(N_obs_inf_B - N_sim_inf_B)/10000
  
  return(c(Dist_A, Dist_B))
}

ABC_MCMC = function(N_A,
                    N_B, # observed counts in A and B
                    start_val, # start point of algorithm 
                    var_proposal_dist, # variance of proposal distribution, can be tuned for acceptance rate 
                    epsilon, # threshold, larger value gives faster convergence, smaller value gives more accurate result
                    N, # Length of chain 
                    calibrate_assort = F, # Boolean for whether or not the assortativity is unknown 
                    assort = NULL) # True assortativity if known
{                 
  thetas = c()
  theta_0 = start_val
  for(i in 1:N){
    if(calibrate_assort){
      theta_p = c(rnorm(1, mean = theta_0[1], sd = sqrt(var_proposal_dist[1])), 
                  rnorm(1, mean = theta_0[2], sd = sqrt(var_proposal_dist[2])))
      while(any(theta_p < 0)){
        theta_p = c(rnorm(1, mean = theta_0[1], sd = sqrt(var_proposal_dist[1])), 
                    rnorm(1, mean = theta_0[2], sd = sqrt(var_proposal_dist[2])))
      }
      params = theta_p
    }
    else{
      theta_p = rnorm(1, mean = theta_0, sd = sqrt(var_proposal_dist))
      params = c(theta_p, assort)
      while(theta_p < 0){
        theta_p = rnorm(1, mean = theta_0, sd = sqrt(var_proposal_dist))
        params = c(theta_p, assort)        
      }
    }
      
    
    dists = compare_cases(params, N_obs_inf_A = N_A, N_obs_inf_B = N_B)
    if(sum(dists) < epsilon){
      accept_prob = min(1, prior(theta_p, calibrate_assort) * dnorm(theta_0, mean = theta_p, sd = sqrt(var_proposal_dist)) / (prior(theta_0, calibrate_assort) * dnorm(theta_p, mean = theta_0, sd = sqrt(var_proposal_dist))))
      rand = runif(1)
      if(rand < accept_prob){
        print("accepted a particle!")
        thetas = rbind(thetas, theta_p)
        theta_0 = theta_p
      }
      else{
        thetas = rbind(thetas, theta_0)
      }
    }
    else{
      thetas = rbind(thetas, theta_0)
    }
  }
  return(thetas)
}  
