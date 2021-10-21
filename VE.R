library(odin.dust)
library(ggplot2)
library(data.table)
library(dplyr)
source("model.R")


# Checking final size
fs <- c()
N <- c(500000, 500000)
input_mats <- cij_NGM(matrix(1,nrow=2, ncol=2), N, c(1,0.5), c(1,1), gamma=1/3)
input_mats$c_ij
beta <- 1.3/input_mats$beta_R
res <- run_model(input_mats$c_ij, rep(beta, t), N, t,
                 c(500,500), 1000, 3, susceptibility = c(1,0.5), print_index=T)


RR1 = colMeans(res$full_results[9,,]/res$full_results[8,,])


input_mats <- cij_NGM(diag(c(1,1)), N, c(1,0.5), c(1,1), gamma=1/3)
input_mats$c_ij
beta <- 1.3/input_mats$beta_R
res <- run_model(input_mats$c_ij, rep(beta, t), N, t,
                 c(500,500), 1000, 3, susceptibility = c(1,0.5), print_index=T)


RR2 = colMeans(res$full_results[9,,]/res$full_results[8,,])




input_mats <- cij_NGM(matrix(c(2,1,1,2),nrow=2, ncol=2), N, c(1,0.5), c(1,1), gamma=1/3)
input_mats$c_ij
beta <- 1.3/input_mats$beta_R
res <- run_model(input_mats$c_ij, rep(beta, t), N, t,
                 c(500,500), 1000, 3, susceptibility = c(1,0.5), print_index=T)

RR3 = colMeans(res$full_results[9,,]/res$full_results[8,,])


df <- rbind(data.frame(t=1:100,
                       RR=RR1,
                       Model="Random mixing"),
            data.frame(t=1:100,
                       RR=RR2,
                       Model="No contact"),
            data.frame(t=1:100,
                       RR=RR3,
                       Model="Assortative mixing"))

add_theme(ggplot(df) + geom_line(aes(x=t, y=RR, color=Model), size=2) + ylab("Incidence ratio") + xlab("Days"))
ggsave("RR_over_time.png")
