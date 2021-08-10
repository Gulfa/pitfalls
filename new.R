library(odin.dust)

odin_model <- odin.dust::odin_dust("odinmodel.R")

params <- list(
  n=4,
  S_ini=c(90000, 0, 10000, 0),
  I_ini=c(1,0,1,0),
  mixing_matrix=diag(c(1,0,1.2,0)),
  beta_day=rep(0.3, 300),
  N_steps=300
)

  

dust_model <- odin_model$new(pars = params,
                             step = 1,
                             n_particles = 100,
                             n_threads = 1
                             )

res <- dust_model$simulate(1:300)
dust_index <- names(dust_model$info()$index)



