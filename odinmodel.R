dt <- user(1)
initial(time) <- 0
update(time) <- (step + 1) * dt
N_steps <- user()

beta <- beta_day[step]
dim(beta_day) <- N_steps
beta_day[] <- user()

update(S[]) <- S[i] - n_SI[i]
update(I[]) <- I[i] + n_SI[i] - n_IR[i]
update(R[]) <- R[i] + n_IR[i]

## Individual probabilities of transition:
p_SI[] <- 1 - exp(-sum(lambda_ij[i,])* dt) # S to I
p_IR <- 1 - exp(-gamma * dt) # I to R

## Draws from binomial distributions for numbers changing between
## compartments:
n_SI[] <- rbinom(S[i], p_SI[i])
n_IR[] <- rbinom(I[i], p_IR)


lambda_ij[,] <- beta * mixing_matrix[i,j]*I[j]/sum(N)*susceptibility[i]*transmisibility[j]



## Total population size
N[] <- S[i] + I[i] + R[i]

## Initial states:
initial(S[]) <- S_ini[i]
initial(I[]) <- I_ini[i]
initial(R[]) <- 0

dim(S) <- n
dim(I) <- n
dim(R) <- n
dim(N) <- n
dim(p_SI) <- n
dim(n_SI) <- n
dim(n_IR) <- n
dim(lambda_ij) <- c(n,n)
dim(beta_norm) <- n
dim(susceptibility) <- n
dim(transmisibility) <- n
## User defined parameters - default in parentheses:
gamma <- user(0.1)

n <- user(4)
S_ini[] <- user()
I_ini[] <- user()
mixing_matrix[,] <- user()
dim(mixing_matrix) <- c(n,n)
dim(S_ini) <- n
dim(I_ini) <- n
beta_norm[] <- user()
susceptibility[] <- user()
transmisibility[] <- user()
              
