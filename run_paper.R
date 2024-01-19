library(odin.dust)
library(ggplot2)
library(data.table)
library(dplyr)
library(gridExtra)
library(ggpubr)
source("model.R")

# Script for producing the figures for the paper

# Molde validation to check final sizes
fs <- c()
Rs <- seq(0.8, 2, by=0.05)
for(R in Rs){
  N <- c(90000, 10000)
  input_mats <- cij_NGM(matrix(1,nrow=2, ncol=2), N, c(1,1), c(1,1), gamma=1/3)
  input_mats$c_ij
  beta <- R/input_mats$beta_R
  res <- run_model(input_mats$c_ij, rep(beta, t), N, t,
                   c(90,10), 1000, 3, susceptibility = c(1,1))
  fs <- c(fs, mean(res$fractions))
  
}


# Producing fig 1 - How relative risk changes with R0 and time

out <- data.frame()
for (R in seq(1.1, 3, by=0.01)){
  for(t in c(20,50, 100)){
    input_mats <- cij_NGM(diag(c(1,9)), c(90000,10000), c(1,1.2), c(1,1))
    beta <- R/input_mats$beta_R
    res <- run_model(diag(c(1, 9)), rep(beta, t), c(90000, 10000), t, c(90,10), 2000, n_threads=3,
                     susceptibility = c(1,1.2))
    RR <- quantile(res$fractions[2,]  / res$fractions[1,], probs=c(0.025, 0.5, 0.975))
    out <- rbind(out, data.frame(min=RR[1],
                                 med=RR[2],
                                 max=RR[3],
                                 R=R,
                                 time=t))
  }
}
add_theme(ggplot(out ,aes(x=R, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(time)), alpha = 0.5) + geom_line(aes(color=factor(time)))+ scale_y_continuous(trans='log10') + ylab("Relative risk in B compared to A") + xlab(expression(paste("Reproduction number ", R[0], sep = " "))) + scale_fill_brewer("Time (days)", palette = "Dark2")+ scale_color_brewer("Time (days)", palette = "Dark2"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))))
ggsave("fig1_RR.png")
ggsave("fig1_RR.eps", device = cairo_ps, width = 4, height = 2.5)
ggsave("fig1_RR.pdf", device = "pdf", width = 4, height = 2.5)
ggsave("fig1_RR.tiff", device = "tiff", width = 4, height = 2.5, dpi = 700)



# Producing fig 2a -  How well does a linear model predict the total number of cases
out <- data.frame()
N <- c(90000, 10000)
input_mats <- cij_NGM(matrix(1,nrow=2, ncol=2), N, c(1,1), c(1,1))
beta <- 1.2/input_mats$beta_R
actual_out <- data.frame()

for(sus in seq(1, 4, by=1)){
  t <- 400
  b <- 0.05
  c_ij <- matrix(c(1,1,1,1), nrow=2)
  
  N <- c(90000, 10000)
  input_mats <- cij_NGM(c_ij, N, c(1,sus), c(1,1))
  res <- run_model(input_mats$c_ij, rep(beta, t), N, t,
                   c(90,10), 100, 3, susceptibility = c(1,sus))


  actual <- run_model(input_mats$c_ij, rep(beta, t), c(sum(N), 0), t,
                      c(100,0), 100, 3, susceptibility = c(1,sus))
  
  tot_inf <- quantile(actual$full_results[6,,t], probs=c(0.0275, 0.5, 0.975))
  actual_out <- rbind(actual_out,
                      data.frame(
                        min=tot_inf[1],
                        med=tot_inf[2],
                        max=tot_inf[3],
                        sus=sus))
  n <- 100
  N_sim <- 10000
  for(i in 1:n){
    data <- data.frame(N=N, I=res$full_results[6:7,i, 400],  risk=factor(c("L", "H"), levels=c("L", "H")))
    mod <-  glm(I~offset(log(N)) + risk, family=poisson, data=data)
    cs <- summary(mod)$coefficients[1,1:2]
    coefs <- rnorm(N_sim, mean=cs[1], sd=cs[2])
    vals <- rpois(N_sim,lambda=sum(N)*exp(coefs))
    preds <- quantile(vals, probs=c(0.025, 0.5, 0.975))
    out <- rbind(out, data.frame(min_p=preds[1],
                                 med_p=preds[2],
                                 max_p=preds[3],
                                 sim_med=actual$full_results[6,i,t],
                                 sim=i,
                                 sus=sus)
                 )
  }
}

f2a = add_theme(ggplot(out, aes(x=sus, y=med_p, ymin=min_p, ymax=max_p, color="Prediction", size = 0.001)) + geom_pointrange( position="jitter", size = 0.0001) + ylab("Total cases") + xlab("a") +geom_point(aes(x=sus, y=sim_med, color="Simulation"), size=0.5, position="jitter") + scale_color_brewer("Model", palette="Dark2")) + ggtitle("a)") + theme(plot.title = element_text(size = 8, face = "plain", margin=margin(0,0,0,0))) +  theme(legend.position="top",
                                                                                                                                                                                                                                                                                                                                                                                                                      legend.box.margin=margin(-0.38,-0.38,-0.38,-0.38, unit = "cm"))
ggsave(f2a, file = "fig_2a.png")
ggsave(f2a, file = "fig2_a.eps", device = "eps", width = 3, height = 2.5)
ggsave(f2a, file = "fig2_a.pdf", device = "pdf", width = 3, height = 2.5)



# Producing fig 2b -  How infections in the low risk group changes when changing risk in high risk group
out <- data.frame()
N <- c(90000, 10000)
input_mats <- cij_NGM(matrix(1,nrow=2, ncol=2), N, c(1,1), c(1,1))
beta <- 0.9/input_mats$beta_R
for(sus in seq(1, 10, by=0.2)){
  t <- 400
  b <- 0.05
  c_ij <- matrix(c(1,1,1,1), nrow=2)
  
  N <- c(90000, 10000)
  input_mats <- cij_NGM(c_ij, N, c(1,sus), c(1,1))
  res <- run_model(input_mats$c_ij, rep(beta, t), N, t,
                   c(90,10), 1000, 3, susceptibility = c(1,sus))
  tot_inf <- quantile((res$full_results[6,,t]-90)/90000, probs=c(0.025, 0.5, 0.975))
  mean(tot_inf)
  out <- rbind(out, data.frame(min=tot_inf[1],
                               med=tot_inf[2],
                               max=tot_inf[3],
                               sus=sus))
}
f2b = add_theme(ggplot(out, aes(x=sus, y=med, ymin=min, ymax=max)) + geom_line() + geom_ribbon(alpha=0.3) + ylab("Fraction infected in low-risk group") + xlab(expression(a)) + scale_y_continuous(labels = scales::percent_format(accuracy = 1))) + ggtitle("b)") + theme(plot.title = element_text(size = 8, face = "plain"))
ggsave(f2b, file = "fig2_b.png")
ggsave(f2b, file = "fig2_b.eps", device = "eps", width = 3, height = 2.5)
ggsave(f2b, file = "fig2_b.pdf", device = "pdf", width = 3, height = 2.5)

# Comnbing fig 2a and 2b
f2 = grid.arrange(f2a, f2b, nrow = 1)
ggsave(f2, file = "fig2.eps", device = cairo_ps, width = 6, height = 2.5)
ggsave(f2, file = "fig2.pdf", device = "pdf", width = 6, height = 2.5)




# Figure 3a- Unexplained relative risk as assocaited mixing increases
out <- data.frame()
for( sus in seq(1, 3, by=0.1)){
  for(assort in c(1, 2, 5, 10)){
    res <- run_4x4(matrix(c(assort,1,1, assort), nrow=2),
                    matrix(c(1,1, 1, 1), nrow=2),
                    1.3,
                    c(1,sus,1,sus),
                   rep(1,4), n=500)
    RR_HL <- ((res$full_results[11,, 200] + res$full_results[13,, 200])/15000)/((res$full_results[10,, 200] + res$full_results[12,, 200])/85000)
    RR_eth <- ((res$full_results[12,, 200] + res$full_results[13,, 200])/10000)/((res$full_results[10,, 200] + res$full_results[11,, 200])/90000)
    explained <- (0.5 + 0.5*RR_HL) / (0.9 + 0.1*RR_HL)
    unexplained <-quantile(RR_eth - explained, probs=c(0.025, 0.5, 0.975))
    out <- rbind(out,data.frame(min=unexplained[1],
                                med=unexplained[2],
                                max=unexplained[3],
                                sus=sus,
                                assort=assort))
                       
                 
  }
}

f3a = add_theme(ggplot(out ,aes(x=sus, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(assort)), alpha=0.5) + geom_line(aes(color=factor(assort)))+   ylab("Unexplained relative risk") + xlab(expression("a")) + scale_fill_brewer("Assortative mixing", palette = "Dark2")+ scale_color_brewer("Assortative mixing", palette = "Dark2")+ scale_y_continuous(labels = scales::percent_format(accuracy = 1))) + ggtitle("a)") + theme(plot.title = element_text(size = 8, face = "plain"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))
ggsave(f3a, "fig3_1_uo.png")
ggsave(f3a, file = "fig3_a.eps", device = cairo_ps, width = 5, height = 2.5)
ggsave(f3a, file = "fig3_a.pdf", device = "pdf", width = 5, height = 2.5)


# Producing figure 3b - Regression coefficients chaning with associate mixing where high risk group is more susceptible
out <- data.frame()
for( sus in c(2)){
  for(assort in c(1, 1.5, 2, 5, 10)){
    res <- run_regs(matrix(c(assort,1,1, assort), nrow=2),
                    matrix(c(1,1, 1, 1), nrow=2),
                    1.3,
                    c(1,sus,1,sus),
                    rep(1,4), n=100)
    
    for(r in res){
      mid <- exp(coefficients(r)[2])
      ci <- exp(confint(r)[2,])
      out <- rbind(out, data.frame(
                          mid=mid,
                          lower=ci[1],
                          upper=ci[2],
                          assort=assort,
                          sus=sus))
    }
  }
}
f3b = add_theme(ggplot(out, aes(x=factor(assort), y=mid, ymin=lower, ymax=upper)) +  geom_pointrange( position="jitter", size = 0.05, fatten = 0.1, fill = "pink", shape = 21) + ylab(expression(paste("Regression coefficent ethnicity ", beta[e], " "))) + xlab("Assortative mixing")) + ggtitle("b)") + theme(plot.title = element_text(size = 8, face = "plain")) 
ggsave(f3b, "fig3_1_beta.png")
ggsave(f3b, file = "fig3_b.eps", device = "eps", width = 5, height = 2.5)
ggsave(f3b, file = "fig3_b.pdf", device = "pdf", width = 5, height = 2.5)


# Producing figure 3c - Regressision coefficients changing with associate mixing where high risk group have more contacts
out <- data.frame()
for( hc in c(3)){
  for(assort in c(1, 1.5, 2, 5, 10)){
    res <- run_regs(matrix(c(assort,1,1, assort), nrow=2),
                    matrix(c(1,1, 1, hc), nrow=2),
                    1.3,
                    c(1,1,1,1),
                    rep(1,4), n=100)
    
    for(r in res){
      mid <- exp(coefficients(r)[2])
      ci <- exp(confint(r)[2,])
      out <- rbind(out, data.frame(
                          mid=mid,
                          lower=ci[1],
                          upper=ci[2],
                          assort=assort,
                          sus=sus))
    }
  }
}


f3c = add_theme(ggplot(out, aes(x=factor(assort), y=mid, ymin=lower, ymax=upper))  +  geom_pointrange( position="jitter", size = 0.05, fatten = 0.1, shape = 21) + ylab(expression(paste("Regression coefficent ethnicity ", beta[e], " "))) + xlab("Assortative mixing")) + ggtitle("c)") + theme(plot.title = element_text(size = 8, face = "plain"))
ggsave(f3c,"fig3_1_beta.png")
ggsave(f3c, file = "fig3_c.eps", device = "eps", width = 5, height = 2.5)
ggsave(f3c, file = "fig3_c.pdf", device = "pdf", width = 5, height = 2.5)

# Combing fig 3a, 3b and 3c
f3 = grid.arrange(f3a, f3b, f3c, nrow = 3)
ggsave(f3, file = "fig3.eps", device = cairo_ps, width = 5, height = 8)
ggsave(f3, file = "fig3.pdf", device = "pdf", width = 5, height = 8)

### Figure 4 - Incidence ratio over time for different types of mixing
# Random mixing
t <- 100
fs <- c()
N <- c(500000, 500000)
input_mats <- cij_NGM(matrix(1,nrow=2, ncol=2), N, c(1,0.5), c(1,1), gamma=1/3)
input_mats$c_ij
beta <- 1.3/input_mats$beta_R
res <- run_model(input_mats$c_ij, rep(beta, t), N, t,
                 c(500,500), 1000, 3, susceptibility = c(1,0.5), print_index=T)
RR1 = colMeans(res$full_results[9,,]/res$full_results[8,,])

# No contact
input_mats <- cij_NGM(diag(c(1,1)), N, c(1,0.5), c(1,1), gamma=1/3)
input_mats$c_ij
beta <- 1.3/input_mats$beta_R
res <- run_model(input_mats$c_ij, rep(beta, t), N, t,
                 c(500,500), 1000, 3, susceptibility = c(1,0.5), print_index=T)
RR2 = colMeans(res$full_results[9,,]/res$full_results[8,,])


#Assortative mixing
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
ggsave("RR_over_time.pdf")
ggsave("RR_over_time.eps")

