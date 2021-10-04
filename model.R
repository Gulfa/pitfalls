odin_model <- odin.dust::odin_dust("odinmodel.R")




run_model <- function(mixing_matrix, beta_day, N, t, I_ini,
                      n_particles, n_threads=1, beta_norm=NULL,
                      susceptibility=NULL,
                      transmisibility=NULL, dt=0.2, print_index=FALSE, gamma=1/3){
  if(is.null(beta_norm)) beta_norm <- N
  if(is.null(susceptibility)) susceptibility <- rep(1, dim(mixing_matrix)[1])
  if(is.null(transmisibility)) transmisibility <- rep(1, dim(mixing_matrix)[1])
  params <- list(
    n=dim(mixing_matrix)[1],
    S_ini=N - I_ini,
    I_ini=I_ini,
    mixing_matrix=mixing_matrix,
    beta_day=rep(beta_day, each=1/dt),
    beta_norm=beta_norm,
    susceptibility=susceptibility,
    transmisibility=transmisibility,
    N_steps=t/dt,
    dt=dt,
    gamma=gamma
  )
  dust_model <- odin_model$new(pars = params,
                               step = 1,
                               n_particles = n_particles,
                               n_threads = n_threads
                               )
  if(print_index){
    print(dust_model$info()$index)
  }
  res <- dust_model$simulate(1:t/dt)
  n <- dim(mixing_matrix)[1]
  R <- res[(2+2*n):(1+3*n),,t]
  fractions <- (R-I_ini) / (N + 1e-9)

  return(list(fractions=fractions,
              R=R,
              full_results=res))
}

plot_2x2_RR <- function(fractions){

  RR <- fractions[2,]  / fractions[1,]

  ggplot(data.frame(RR=RR)) + geom_histogram(aes(x=RR))
  
  
}

plot_history <- function(hist, var="I"){
  n <- (dim(hist)[1] - 1)/3
  sims <- dim(hist)[2]
  all_dfs <- list()
  for(i in 1:n){
    group <- hist[c((1 + i), (1+n+i), (1+2*n+i)),,]
    dfs <- list()
    for(sim in 1:sims){
      dfs[[sim]] <- data.frame(S=group[1, sim,],
                       I=group[2, sim,],
                       R=group[3, sim,],
                       t=1:dim(hist)[3],
                       sim=sim,
                       group=i)
    }
    
    all_dfs[[i]] <- rbindlist(dfs)
  }
  df <- rbindlist(all_dfs)
  df[, factor_group:=paste(sim, group)]
  
  ggplot(df) + geom_line(aes(x=t, y=get(var), group=factor_group, color=factor(group))) + ylab(var)

}

cij_NGM <- function(c_ij, N, susceptibility, transmisibility, gamma=1/3, norm_contacts=NULL){
  if(is.null(norm_contacts)){
    norm <- c_ij %*% N/sum(N)
  }else{
    N_conts <- as.numeric(norm_contacts %*% N)
    norm <- c_ij %*% N/sum(N)/N_conts*sum(N_conts)
  }

  c_ij <- c_ij/as.numeric(norm)
  NGM <- c_ij %*% diag(transmisibility)*N/sum(N)*susceptibility #NGM = suscept_i * S_i/N c_ij *trans_j
#  print(NGM)
  beta_R <- Re(eigen(NGM, only.values=T)$values[1]/gamma)
  return(list(c_ij=c_ij,
              NGM=NGM,
              beta_R=beta_R))
}


run_4x4 <- function(cimat, crmat, R0, susceptibility, transmisibility, gamma=1/3, n=100){

  N_N <- 90000
  N_I <- 10000
  N_IH <- 5000
  N_IL <- N_I - N_IH
  N_NH <- 0.1*N_N
  N_NL <- N_N - N_NH
  N <- c(N_NL, N_NH, N_IL, N_IH)  
  #cimat <- 
  cimat <- cimat/rowSums(cimat)
  #crmat <- matrix(c(1, 1, 1, 1), nrow=2)
  #crmat <- crmat/rowSums(crmat)
  
  c_ij <- matrix(1, nrow=4, ncol=4)
  c_ij[1:2, 1:2] <- cimat[1,1]*crmat
  c_ij[1:2, 3:4] <- cimat[1,2]*crmat
  c_ij[3:4, 1:2] <- cimat[2,1]*crmat
  c_ij[3:4, 3:4] <- cimat[2,2]*crmat

  norm_contacts <- cbind(rbind(crmat, crmat),rbind(crmat, crmat))

  input_mats <- cij_NGM(c_ij, N, susceptibility, transmisibility, norm_contacts=norm_contacts)
  beta <- R0/input_mats$beta_R
#  print("final c_ij")
#  print(input_mats$c_ij)
  res <- run_model(input_mats$c_ij, rep(beta, 200), N, 200, 0.001*N, n, 1,
                   susceptibility=susceptibility, transmisibility = transmisibility, gamma=gamma)
  res$N <- N
  return(res)

}


run_regs <- function(cimat, crmat, R0, susceptibility, transmisibility, on_mean=FALSE,
                     n=100, gamma=1/3){
  
  res <- run_4x4(cimat,
                 crmat,
                 R0,
                 susceptibility,
                 transmisibility,
                 n=n,
                 gamma=gamma)
  N <- res$N

  if(on_mean){
    data <- data.frame(N=N, I=rowMeans(res$full_results[10:13, , 200]), ethnicity=factor(c("N", "N", "I", "I"), levels=c("N", "I")), risk=c("L", "H", "L", "H"))
    return(glm(I~offset(log(N)) + ethnicity + risk, data=data, family=poisson))
  }else{
    sums <- list()
    for(i in 1:n){
      data <- data.frame(N=N, I=res$full_results[10:13,i, 200], ethnicity=c("N", "N", "I", "I"), risk=c("L", "H", "L", "H"))
      m <- glm(I~offset(log(N)) + ethnicity + risk, data=data, family=poisson)
      sums[[i]] <- m
    }
    return(sums)
  }
}



add_theme <- function(q){
  q + theme_bw() + theme(text = element_text(size=8))+ scale_size_identity()
}
