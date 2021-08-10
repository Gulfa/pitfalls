library(odin.dust)
library(ggplot2)
library(data.table)
odin_model <- odin.dust::odin_dust("odinmodel.R")




run_model <- function(mixing_matrix, beta_day, N, t, I_ini,
                      n_particles, n_threads=1){
  print(dim(mixing_matrix)[1])
  params <- list(
    n=dim(mixing_matrix)[1],
    S_ini=N - I_ini,
    I_ini=I_ini,
    mixing_matrix=mixing_matrix,
    beta_day=beta_day,
    N_steps=t
  )
  dust_model <- odin_model$new(pars = params,
                               step = 1,
                               n_particles = n_particles,
                               n_threads = n_threads
                               )
  
  res <- dust_model$simulate(1:t)
  n <- dim(mixing_matrix)[1]
  R <- res[(2+2*n):(1+3*n),,t]
  fractions <- R / (N + 1e-9)

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

res <- run_model(diag(c(1/90000, 1.2/10000)), rep(0.2, 200), c(90000, 10000), 200, c(90,10), 200, 1)

plot_history(res$full_results, var="R")

plot_2x2_RR(res$fractions)


res <- run_model(diag(c(1/10000, 1.1/10000, 1.2/10000, 1.3/10000)), rep(0.2, 200), c(10000, 10000, 10000, 10000), 200, c(10,10,10,10), 200, 1)

plot_history(res$full_results, var="I")
