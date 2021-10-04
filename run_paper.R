library(odin.dust)
library(ggplot2)
library(data.table)
library(dplyr)
library(gridExtra)
source("model.R")

# Checking final size
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


#add_theme(ggplot(out ,aes(x=R, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(col=factor(time)), fill = "NA", linetype="dashed") + geom_line(aes(color=factor(time)))+ scale_y_continuous(trans='log10') + ylab("Relative risk") + xlab(expression(paste("Reproduction number ", R[0], sep = " "))) + scale_fill_brewer("Time (days)", palette = "Dark2")+ scale_color_brewer("Time (days)", palette = "Dark2"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))))
add_theme(ggplot(out ,aes(x=R, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(time)), alpha = 0.5) + geom_line(aes(color=factor(time)))+ scale_y_continuous(trans='log10') + ylab("Relative risk") + xlab(expression(paste("Reproduction number ", R[0], sep = " "))) + scale_fill_brewer("Time (days)", palette = "Dark2")+ scale_color_brewer("Time (days)", palette = "Dark2"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))))

ggsave("fig1_RR.png")

ggsave("fig1_RR.eps", device = cairo_ps, width = 4, height = 2.5)

ggsave("fig1_RR.pdf", device = "pdf", width = 4, height = 2.5)

ggsave("fig1_RR.tiff", device = "tiff", width = 4, height = 2.5, dpi = 700)



out <- data.frame()

N <- c(90000, 10000)
input_mats <- cij_NGM(matrix(1,nrow=2, ncol=2), N, c(1,1), c(1,1))
beta <- 0.9/input_mats$beta_R



#for(bc in c(0.5, 1, 1.5)){
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

#f2b = add_theme(ggplot(out, aes(x=sus, y=med, ymin=min, ymax=max)) + geom_ribbon(col="black", linetype = "dashed", fill = "light grey") + geom_line() + ylab("Fraction infected in low-risk group") + xlab(expression(a)) + scale_y_continuous(labels = scales::percent_format(accuracy = 1))) + ggtitle("b)") + theme(plot.title = element_text(size = 8, face = "plain"))

f2b = add_theme(ggplot(out, aes(x=sus, y=med, ymin=min, ymax=max)) + geom_line() + geom_ribbon(alpha=0.3) + ylab("Fraction infected in low-risk group") + xlab(expression(a)) + scale_y_continuous(labels = scales::percent_format(accuracy = 1))) + ggtitle("b)") + theme(plot.title = element_text(size = 8, face = "plain"))

ggsave(f2b, file = "fig2_b.png")
ggsave(f2b, file = "fig2_b.eps", device = "eps", width = 3, height = 2.5)
ggsave(f2b, file = "fig2_b.pdf", device = "pdf", width = 3, height = 2.5)



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
                                                                                                                                                                                                                                                                                                                                                                                                                       legend.margin=margin(0,0,0,0),
                                                                                                                                                                                                                                                                                                                                                                                                                       legend.spacing.x = unit(0, 'cm'),
                                                                                                                                                                                                                                                                                                                                                                                                                       legend.box.margin=margin(-0.38,-0.38,-0.38,-0.38, unit = "cm"))
ggsave(f2a, file = "fig_2a.png")
ggsave(f2a, file = "fig2_a.eps", device = "eps", width = 3, height = 2.5)
ggsave(f2a, file = "fig2_a.pdf", device = "pdf", width = 3, height = 2.5)


f2 = grid.arrange(f2a, f2b, nrow = 1)
ggsave(f2, file = "fig2.eps", device = cairo_ps, width = 6, height = 2.5)
ggsave(f2, file = "fig2.pdf", device = "pdf", width = 6, height = 2.5)


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

f3a = add_theme(ggplot(out ,aes(x=sus, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(assort)), alpha=0.5) + geom_line(aes(color=factor(assort)))+   ylab("Unexplained overrepresentation") + xlab(expression("a")) + scale_fill_brewer("Assortative mixing", palette = "Dark2")+ scale_color_brewer("Assortative mixing", palette = "Dark2")+ scale_y_continuous(labels = scales::percent_format(accuracy = 1))) + ggtitle("a)") + theme(plot.title = element_text(size = 8, face = "plain"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))
# f3a = add_theme(ggplot(out ,aes(x=sus, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(col=factor(assort)), fill = "NA", linetype = "dashed") + geom_line(aes(color=factor(assort)), size = 1.5)+   ylab("Unexplained overrepresentation") + xlab(expression("a")) + scale_fill_brewer("Assortative mixing", palette = "Dark2")+ scale_color_brewer("Assortative mixing", palette = "Dark2")+ scale_y_continuous(labels = scales::percent_format(accuracy = 1))) + ggtitle("a)") + theme(plot.title = element_text(size = 8, face = "plain"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))

ggsave(f3a, "fig3_1_uo.png")
ggsave(f3a, file = "fig3_a.eps", device = cairo_ps, width = 5, height = 2.5)

ggsave(f3a, file = "fig3_a.pdf", device = "pdf", width = 5, height = 2.5)






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
ggsave(f3c,"fig4_1_beta.png")
ggsave(f3c, file = "fig3_c.eps", device = "eps", width = 5, height = 2.5)
ggsave(f3c, file = "fig3_c.pdf", device = "pdf", width = 5, height = 2.5)

#hlay = rbind(c(NA, 1, 1, NA), c(2, 2, 3, 3))
f3 = grid.arrange(f3a, f3b, f3c, nrow = 3)
ggsave(f3, file = "fig3.eps", device = cairo_ps, width = 5, height = 8)
ggsave(f3, file = "fig3.pdf", device = "pdf", width = 5, height = 8)

out <- data.frame()
for( hc in seq(1, 4, by=0.2)){
  for(assort in c(1, 2, 5, 10)){
    res <- run_4x4(matrix(c(assort,1,1, assort), nrow=2),
                    matrix(c(1,1, 1, hc), nrow=2),
                    1.3,
                    c(1,1,1,1),
                   rep(1,4), n=500)
    RR_HL <- ((res$full_results[11,, 200] + res$full_results[13,, 200])/15000)/((res$full_results[10,, 200] + res$full_results[12,, 200])/85000)
    RR_eth <- ((res$full_results[12,, 200] + res$full_results[13,, 200])/10000)/((res$full_results[10,, 200] + res$full_results[11,, 200])/90000)
    explained <- (0.5 + 0.5*RR_HL) / (0.9 + 0.1*RR_HL)
    unexplained <-quantile(RR_eth - explained, probs=c(0.025, 0.5, 0.975))
    out <- rbind(out,data.frame(min=unexplained[1],
                                med=unexplained[2],
                                max=unexplained[3],
                                hc=hc,
                                assort=assort))
                       
                 
  }
}

add_theme(ggplot(out ,aes(x=hc, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(assort)), alpha=0.5) + geom_line(aes(color=factor(assort)))+   ylab("Unexplained overrepresentation") + xlab("High-high contact")+ scale_fill_brewer("Assortative mixing", palette = "Dark2")+ scale_color_brewer("Assortative mixing", palette = "Dark2")+ scale_y_continuous(labels = scales::percent_format(accuracy = 1)))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))
ggsave("fig4_1_uo.png")


## ABC-analysis ##

### Here we compute the "truth" from the SIR-model which we would like to recover ###
### We need observed count infected in A, observed count infected in B ###
source("ABC_MCMC.r")
out <- data.frame()
sus = 1.2
assort = 10
b = 1.05
res <- run_4x4(matrix(c(assort,1,1, assort), nrow=2),
                   matrix(c(1,1, 1, 1), nrow=2),
                   1.3,
                   c(1,sus,1 * b,sus * b),
                   rep(1,4), n=1000)
Total_infected = rowMeans(res$R)  
  
N_observed_infected_A = Total_infected[1] + Total_infected[2]
N_observed_infected_B = Total_infected[3] + Total_infected[4]

params1 = ABC_MCMC(N_observed_infected_A, N_observed_infected_B, start_val = c(1.2, 5), var_proposal_dist = c(0.2, 4), epsilon = 100/100000, N = 400000, calibrate_assort = T)
params2 = ABC_MCMC(N_observed_infected_A, N_observed_infected_B, start_val = 1.2, var_proposal_dist = 0.2, epsilon = 100/100000, N = 400000, calibrate_assort = F, assort = 10)

save(params1, params2, file = "params2.Rdata")



N_cores = 120
sus = 1.2
assort = 10
b = 1.05
res_be_neq1 <- run_4x4(matrix(c(assort,1,1, assort), nrow=2),
               matrix(c(1,1, 1, 1), nrow=2),
               1.3,
               c(1, sus ,1 * b, sus * b),
               rep(1,4), n=100)

b = 1.0
res_be_1 <- run_4x4(matrix(c(assort,1,1, assort), nrow=2),
               matrix(c(1,1, 1, 1), nrow=2),
               1.3,
               c(1, sus, 1 * b, sus * b),
               rep(1,4), n=100)
N_chain = 100000
burn_in = 50000

perform_one_MCMC = function(res1, res2){
  N_observed_infected_A = res1[1] + res1[2]
  N_observed_infected_B = res1[3] + res1[4] 
  params1 = ABC_MCMC(N_observed_infected_A, N_observed_infected_B, start_val = c(1.2, 5), var_proposal_dist = c(0.2, 4), epsilon = 100/100000, N = N_chain, calibrate_assort = T)
  params2 = ABC_MCMC(N_observed_infected_A, N_observed_infected_B, start_val = 1.2, var_proposal_dist = 0.2, epsilon = 100/100000, N = N_chain, calibrate_assort = F, assort = 10)
  
  
  N_observed_infected_A = res2[1] + res2[2]
  N_observed_infected_B = res2[3] + res2[4] 
  params3 = ABC_MCMC(N_observed_infected_A, N_observed_infected_B, start_val = c(1.2, 5), var_proposal_dist = c(0.2, 4), epsilon = 100/100000, N = N_chain, calibrate_assort = T)
  params4 = ABC_MCMC(N_observed_infected_A, N_observed_infected_B, start_val = 1.2, var_proposal_dist = 0.2, epsilon = 100/100000, N = N_chain, calibrate_assort = F, assort = 10)
  
  return(list(params1 = params1, params2 = params2, params3 = params3, params4 = params4))
}


combres = parallel::mclapply(1:N, function(i){perform_one_MCMC(res_be_neq1$R[, i], res_be_1$R[, i])}, mc.cores=N_cores)
save(combres, file = "combres.Rdata")

all_be1a = c()
all_be1b = c()
all_be2a = c()
all_be2b = c()
all_assort1 = c()
all_assort2 = c()
mid1a = c()
lower1a = c()
upper1a = c()
mid1b = c()
lower1b = c()
upper1b = c()
mid2a = c()
lower2a = c()
upper2a = c()
mid2b = c()
lower2b = c()
upper2b = c()
for (res in combres){
  be1a = res$params1[burn_in:N_chain, 1]
  be1b = res$params2[burn_in:N_chain, 1]
  be2a = res$params3[burn_in:N_chain, 1]
  be2b = res$params4[burn_in:N_chain, 1]
  
  all_be1a = c(all_be1a, be1a)
  all_be1b = c(all_be1b, be1b)
  all_be2a = c(all_be2a, be2a)
  all_be2b = c(all_be2b, be2b)
  all_assort1 = c(all_assort1, params1[burn_in:N_chain, 2])
  all_assort2 = c(all_assort2, params3[burn_in:N_chain, 2])
  
  mid1a <- c(mid1a, mean(be1a))
  lower1a <- c(lower1a, quantile(be1a, c(0.025, 0.975)))
  upper1a <- c(upper1a, quantile(be1a, c(0.025, 0.975)))
  mid1b <- c(mid1b, mean(be1b))
  lower1b <- c(lower1b, quantile(be1b, c(0.025, 0.975)))
  upper1b <- c(upper1b, quantile(be1b, c(0.025, 0.975)))
  mid2a <- c(mid2a, mean(be2a))
  lower2a <- c(lower2a, quantile(be2a, c(0.025, 0.975)))
  upper2a <- c(upper2a, quantile(be2a, c(0.025, 0.975)))
  mid2b <- c(mid2b, mean(be2b))
  lower2b <- c(lower2b, quantile(be2b, c(0.025, 0.975)))
  upper2b <- c(upper2b, quantile(be2b, c(0.025, 0.975)))

}

par(mfrow = c(3, 2))
hist(all_be1a, xlab = expression(beta[e]), main = paste("True ", expression(beta[e] = 1.05), " assortativity unknown"))
hist(all_be1b, xlab = expression(beta[e]), main = paste("True ", expression(beta[e] = 1.05), " assortativity known"))
hist(all_be2a, xlab = expression(beta[e]), main = paste("True ", expression(beta[e] = 1), " assortativity unknown"))
hist(all_be2b, xlab = expression(beta[e]), main = paste("True ", expression(beta[e] = 1), " assortativity known"))

hist(all_assort1, xlab = "Assortativity", main = paste("True ", expression(beta[e] = 1.05), " assortativity unknown"))
hist(all_assort2, xlab = "Assortativity", main = paste("True ", expression(beta[e] = 1), " assortativity unknown"))

dev.off()

df1a = data.frame("mid" = mid1a, "lower" = lower1a, "upper" = upper1a)
df1b = data.frame("mid" = mid1b, "lower" = lower1b, "upper" = upper1b)
df2a = data.frame("mid" = mid2a, "lower" = lower2a, "upper" = upper2a)
df2b = data.frame("mid" = mid2b, "lower" = lower2b, "upper" = upper2b)

supp1a = add_theme(ggplot(df1a, aes(x = 1:100, y=mid, ymin=lower, ymax=upper))  +  geom_pointrange(size = 0.05, fatten = 0.1, shape = 21) + ylab(expression(paste("Regression coefficent ethnicity ", beta[e], " "))) + ggtitle(paste("True ", expression(beta[e] = 1.05), " assortativity unknown")) + theme(plot.title = element_text(size = 8, face = "plain")) + scale_x_discrete(labels = NULL, breaks = NULL) + xlab(""))
supp1b = add_theme(ggplot(df1b, aes(x = 1:100, y=mid, ymin=lower, ymax=upper))  +  geom_pointrange(size = 0.05, fatten = 0.1, shape = 21) + ylab(expression(paste("Regression coefficent ethnicity ", beta[e], " "))) + ggtitle(paste("True ", expression(beta[e] = 1.05), " assortativity known")) + theme(plot.title = element_text(size = 8, face = "plain")) + scale_x_discrete(labels = NULL, breaks = NULL) + xlab(""))
supp2a = add_theme(ggplot(df2a, aes(x = 1:100, y=mid, ymin=lower, ymax=upper))  +  geom_pointrange(size = 0.05, fatten = 0.1, shape = 21) + ylab(expression(paste("Regression coefficent ethnicity ", beta[e], " "))) + ggtitle(paste("True ", expression(beta[e] = 1), " assortativity unknown")) + theme(plot.title = element_text(size = 8, face = "plain")) + scale_x_discrete(labels = NULL, breaks = NULL) + xlab(""))
supp2b = add_theme(ggplot(df2b, aes(x = 1:100, y=mid, ymin=lower, ymax=upper))  +  geom_pointrange(size = 0.05, fatten = 0.1, shape = 21) + ylab(expression(paste("Regression coefficent ethnicity ", beta[e], " "))) + ggtitle(paste("True ", expression(beta[e] = 1), " assortativity known")) + theme(plot.title = element_text(size = 8, face = "plain")) + scale_x_discrete(labels = NULL, breaks = NULL) + xlab(""))


fsup = grid.arrange(supp1a, supp1b, supp2a, supp2b, nrow = 2)
ggsave(fsup, file = "fig_supp.eps", device = cairo_ps, width = 5, height = 8)

ggsave(fsup, file = "fig_supp.pdf", device = pdf, width = 5, height = 8)
