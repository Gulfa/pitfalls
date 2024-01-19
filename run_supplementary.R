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


#  Running the ABC
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
N_chain = 1000000
burn_in = 100000

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

# Extracting the results for plotting
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
  all_assort1 = c(all_assort1, res$params1[burn_in:N_chain, 2])
  all_assort2 = c(all_assort2, res$params3[burn_in:N_chain, 2])
  
  mid1a <- c(mid1a, mean(be1a))
  lower1a <- c(lower1a, quantile(be1a, c(0.025)))
  upper1a <- c(upper1a, quantile(be1a, c(0.975)))
  mid1b <- c(mid1b, mean(be1b))
  lower1b <- c(lower1b, quantile(be1b, c(0.025)))
  upper1b <- c(upper1b, quantile(be1b, c(0.975)))
  mid2a <- c(mid2a, mean(be2a))
  lower2a <- c(lower2a, quantile(be2a, c(0.025)))
  upper2a <- c(upper2a, quantile(be2a, c(0.975)))
  mid2b <- c(mid2b, mean(be2b))
  lower2b <- c(lower2b, quantile(be2b, c(0.025)))
  upper2b <- c(upper2b, quantile(be2b, c(0.975)))

}

# Printing key results
print(c(mean(all_be1a), quantile(all_be1a, c(0.025, 0.975))))
print(c(mean(all_be1b), quantile(all_be1b, c(0.025, 0.975))))
print(c(mean(all_be2a), quantile(all_be2a, c(0.025, 0.975))))
print(c(mean(all_be2b), quantile(all_be2b, c(0.025, 0.975))))

print(c(mean(all_assort1), quantile(all_assort1, c(0.025, 0.975))))
print(c(mean(all_assort2), quantile(all_assort2, c(0.025, 0.975))))

print(cor(all_be1a, all_assort1))
print(cor(all_be2a, all_assort2))

# Parameter histograms

pdf("histograms.pdf", width = 5, height = 8)
par(mfrow = c(3, 2))
hist(all_be1a, xlab = expression(beta[e]), main = expression(paste(beta[e], "=1.05 assortativity unknown")), font.main = 1, cex.lab = 1, cex.main = 1)
hist(all_be1b, xlab = expression(beta[e]), main = expression(paste(beta[e], "=1.05 assortativity known")), font.main = 1, cex.lab = 1, cex.main = 1)
hist(all_be2a, xlab = expression(beta[e]), main = expression(paste(beta[e], "=1 assortativity unknown")), font.main = 1, cex.main = 1, cex.lab = 1)
hist(all_be2b, xlab = expression(beta[e]), main = expression(paste(beta[e], "=1 assortativity known")), font.main = 1, cex.main = 1, cex.lab = 1)

hist(all_assort1, xlab = "Assortativity", main = expression(paste(beta[e], "=1.05 assortativity unknown")), font.main = 1, cex.lab = 1, cex.main=1)
hist(all_assort2, xlab = "Assortativity", main = expression(paste(beta[e], "=1 assortativity unknown")), font.main = 1, cex.lab=1, cex.main=1)

dev.off()

# Plot of regression coefficients
df1a = data.frame("mid" = mid1a, "lower" = lower1a, "upper" = upper1a)
df1b = data.frame("mid" = mid1b, "lower" = lower1b, "upper" = upper1b)
df2a = data.frame("mid" = mid2a, "lower" = lower2a, "upper" = upper2a)
df2b = data.frame("mid" = mid2b, "lower" = lower2b, "upper" = upper2b)


supp1a = add_theme(ggplot(df1a, aes(x = 1:dim(df1a)[1], y=mid, ymin=lower, ymax=upper))  +  geom_pointrange(size = 0.4, fatten = 0.9) + ylab(expression(paste("Regression coefficent ethnicity ", beta[e], " "))) + ggtitle(expression(paste(beta[e], "=1.05 assortativity unknown"))) + theme(plot.title = element_text(size = 8, face = "plain")) + scale_x_discrete(labels = NULL, breaks = NULL) + xlab("")) + ylim(0.9, 1.25)
supp1b = add_theme(ggplot(df1b, aes(x = 1:dim(df1b)[1], y=mid, ymin=lower, ymax=upper))  +  geom_pointrange(size = 0.4, fatten = 0.9) + ylab(expression(paste("Regression coefficent ethnicity ", beta[e], " "))) + ggtitle(expression(paste(beta[e], "=1.05 assortativity known"))) + theme(plot.title = element_text(size = 8, face = "plain")) + scale_x_discrete(labels = NULL, breaks = NULL) + xlab("")) + ylim(0.9, 1.25)
supp2a = add_theme(ggplot(df2a, aes(x = 1:dim(df2a)[1], y=mid, ymin=lower, ymax=upper))  +  geom_pointrange(size = 0.4, fatten = 0.9) + ylab(expression(paste("Regression coefficent ethnicity ", beta[e], " "))) + ggtitle(expression(paste(beta[e], "=1 assortativity unknown"))) + theme(plot.title = element_text(size = 8, face = "plain")) + scale_x_discrete(labels = NULL, breaks = NULL) + xlab("")) + ylim(0.9, 1.25)
supp2b = add_theme(ggplot(df2b, aes(x = 1:dim(df2b)[1], y=mid, ymin=lower, ymax=upper))  +  geom_pointrange(size = 0.4, fatten = 0.9) + ylab(expression(paste("Regression coefficent ethnicity ", beta[e], " "))) + ggtitle(expression(paste(beta[e], "=1 assortativity known"))) + theme(plot.title = element_text(size = 8, face = "plain")) + scale_x_discrete(labels = NULL, breaks = NULL) + xlab("")) + ylim(0.9, 1.25)


fsup = grid.arrange(supp1a, supp1b, supp2a, supp2b, nrow = 2)
ggsave(fsup, file = "fig_supp.eps", device = cairo_ps, width = 5, height = 6)

ggsave(fsup, file = "fig_supp.pdf", device = pdf, width = 5, height = 6)




## Figure S1 - Time serier of I(t) for different R-values

R = 1.1
input_mats <- cij_NGM(diag(c(1,9)), c(90000,10000), c(1,1.2), c(1,1))

beta <- R/input_mats$beta_R

res <- run_model(diag(c(1, 9)), rep(beta, t), c(90000, 10000), t, c(90,10), 2000, n_threads=3,
                 susceptibility = c(1,1.2))
tmp11A = t(apply(res$full_results[4, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
tmp11B = t(apply(res$full_results[5, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
colnames(tmp11A) = c("min", "med", "max")
colnames(tmp11B) = c("min", "med", "max")
out11 = data.frame(cbind(rbind(tmp11A, tmp11B), group = c(rep("A", t), rep("B", t))), Day = c(1:t, 1:t))
out11$min = as.numeric(out11$min)
out11$med = as.numeric(out11$med)
out11$max = as.numeric(out11$max)
S1a = add_theme(ggplot(out11 ,aes(x=Day, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(group)), alpha = 0.5) + geom_line(aes(color=factor(group))) + ylab("I(t), R = 1.1") + xlab("Day") + scale_fill_brewer("Group", palette = "Dark2")+ scale_color_brewer("Group", palette = "Dark2"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))))

R = 1.5 
input_mats <- cij_NGM(diag(c(1,9)), c(90000,10000), c(1,1.2), c(1,1))

beta <- R/input_mats$beta_R

res <- run_model(diag(c(1, 9)), rep(beta, t), c(90000, 10000), t, c(90,10), 2000, n_threads=3,
                 susceptibility = c(1,1.2))
tmp15A = t(apply(res$full_results[4, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
tmp15B = t(apply(res$full_results[5, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
colnames(tmp15A) = c("min", "med", "max")
colnames(tmp15B) = c("min", "med", "max")
out15 = data.frame(cbind(rbind(tmp15A, tmp15B), group = c(rep("A", t), rep("B", t))), Day = c(1:t, 1:t))
out15$min = as.numeric(out15$min)
out15$med = as.numeric(out15$med)
out15$max = as.numeric(out15$max)
S1b = add_theme(ggplot(out15 ,aes(x=Day, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(group)), alpha = 0.5) + geom_line(aes(color=factor(group))) + ylab("I(t), R = 1.5") + xlab("Day") + scale_fill_brewer("Group", palette = "Dark2")+ scale_color_brewer("Group", palette = "Dark2"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))))

R = 2.0 
input_mats <- cij_NGM(diag(c(1,9)), c(90000,10000), c(1,1.2), c(1,1))

beta <- R/input_mats$beta_R

res <- run_model(diag(c(1, 9)), rep(beta, t), c(90000, 10000), t, c(90,10), 2000, n_threads=3,
                 susceptibility = c(1,1.2))
tmp2A = t(apply(res$full_results[4, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
tmp2B = t(apply(res$full_results[5, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
colnames(tmp2A) = c("min", "med", "max")
colnames(tmp2B) = c("min", "med", "max")
out2 = data.frame(cbind(rbind(tmp2A, tmp2B), group = c(rep("A", t), rep("B", t))), Day = c(1:t, 1:t))
out2$min = as.numeric(out2$min)
out2$med = as.numeric(out2$med)
out2$max = as.numeric(out2$max)
S1c = add_theme(ggplot(out2 ,aes(x=Day, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(group)), alpha = 0.5) + geom_line(aes(color=factor(group))) + ylab("I(t), R = 2.0") + xlab("Day") + scale_fill_brewer("Group", palette = "Dark2")+ scale_color_brewer("Group", palette = "Dark2"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))))


R = 3.0 
input_mats <- cij_NGM(diag(c(1,9)), c(90000,10000), c(1,1.2), c(1,1))

beta <- R/input_mats$beta_R

res <- run_model(diag(c(1, 9)), rep(beta, t), c(90000, 10000), t, c(90,10), 2000, n_threads=3,
                 susceptibility = c(1,1.2))
tmp3A = t(apply(res$full_results[4, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
tmp3B = t(apply(res$full_results[5, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
colnames(tmp3A) = c("min", "med", "max")
colnames(tmp3B) = c("min", "med", "max")
out3 = data.frame(cbind(rbind(tmp3A, tmp3B), group = c(rep("A", t), rep("B", t))), Day = c(1:t, 1:t))
out3$min = as.numeric(out3$min)
out3$med = as.numeric(out3$med)
out3$max = as.numeric(out3$max)
S1d = add_theme(ggplot(out3 ,aes(x=Day, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(group)), alpha = 0.5) + geom_line(aes(color=factor(group))) + ylab("I(t), R = 3.0") + xlab("Day") + scale_fill_brewer("Group", palette = "Dark2")+ scale_color_brewer("Group", palette = "Dark2"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))))

S1 = ggarrange(S1a, S1b, S1c, S1d, nrow = 2, ncol = 2, common.legend = T, legend = "bottom")
ggsave(S1, file = "S1.eps", device = cairo_ps, width = 5, height = 5)
ggsave(S1, file = "S1.pdf", device = "pdf", width = 5, height = 5)



### FIGURE S2 - Time series of I(t) for different alpha-values
sus = 1 
t <- 400
b <- 0.05
c_ij <- matrix(c(1,1,1,1), nrow=2)
  
N <- c(90000, 10000)
input_mats <- cij_NGM(c_ij, N, c(1,sus), c(1,1))
res <- run_model(input_mats$c_ij, rep(beta, t), N, t,
                  c(90,10), 100, 3, susceptibility = c(1,sus))
tmp1A = t(apply(res$full_results[4, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
tmp1B = t(apply(res$full_results[5, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
colnames(tmp1A) = c("min", "med", "max")
colnames(tmp1B) = c("min", "med", "max")
out1 = data.frame(cbind(rbind(tmp1A, tmp1B), group = c(rep("A", t), rep("B", t))), Day = c(1:t, 1:t))
out1$min = as.numeric(out1$min)
out1$med = as.numeric(out1$med)
out1$max = as.numeric(out1$max)
S2a = add_theme(ggplot(out1 ,aes(x=Day, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(group)), alpha = 0.5) + geom_line(aes(color=factor(group))) + ylab("I(t), a = 1.0") + xlab("Day") + scale_fill_brewer("Group", palette = "Dark2")+ scale_color_brewer("Group", palette = "Dark2"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))))


sus = 2 
t <- 400
b <- 0.05
c_ij <- matrix(c(1,1,1,1), nrow=2)

N <- c(90000, 10000)
input_mats <- cij_NGM(c_ij, N, c(1,sus), c(1,1))
res <- run_model(input_mats$c_ij, rep(beta, t), N, t,
                 c(90,10), 100, 3, susceptibility = c(1,sus))
tmp2A = t(apply(res$full_results[4, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
tmp2B = t(apply(res$full_results[5, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
colnames(tmp2A) = c("min", "med", "max")
colnames(tmp2B) = c("min", "med", "max")
out2 = data.frame(cbind(rbind(tmp2A, tmp2B), group = c(rep("A", t), rep("B", t))), Day = c(1:t, 1:t))
out2$min = as.numeric(out2$min)
out2$med = as.numeric(out2$med)
out2$max = as.numeric(out2$max)
S2b = add_theme(ggplot(out2 ,aes(x=Day, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(group)), alpha = 0.5) + geom_line(aes(color=factor(group))) + ylab("I(t), a = 2.0") + xlab("Day") + scale_fill_brewer("Group", palette = "Dark2")+ scale_color_brewer("Group", palette = "Dark2"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))))



sus = 3 
t <- 400
b <- 0.05
c_ij <- matrix(c(1,1,1,1), nrow=2)

N <- c(90000, 10000)
input_mats <- cij_NGM(c_ij, N, c(1,sus), c(1,1))
res <- run_model(input_mats$c_ij, rep(beta, t), N, t,
                 c(90,10), 100, 3, susceptibility = c(1,sus))
tmp3A = t(apply(res$full_results[4, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
tmp3B = t(apply(res$full_results[5, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
colnames(tmp3A) = c("min", "med", "max")
colnames(tmp3B) = c("min", "med", "max")
out3 = data.frame(cbind(rbind(tmp3A, tmp3B), group = c(rep("A", t), rep("B", t))), Day = c(1:t, 1:t))
out3$min = as.numeric(out3$min)
out3$med = as.numeric(out3$med)
out3$max = as.numeric(out3$max)
S2c = add_theme(ggplot(out3 ,aes(x=Day, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(group)), alpha = 0.5) + geom_line(aes(color=factor(group))) + ylab("I(t), a = 3.0") + xlab("Day") + scale_fill_brewer("Group", palette = "Dark2")+ scale_color_brewer("Group", palette = "Dark2"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))))


sus = 4 
t <- 400
b <- 0.05
c_ij <- matrix(c(1,1,1,1), nrow=2)

N <- c(90000, 10000)
input_mats <- cij_NGM(c_ij, N, c(1,sus), c(1,1))
res <- run_model(input_mats$c_ij, rep(beta, t), N, t,
                 c(90,10), 100, 3, susceptibility = c(1,sus))
tmp4A = t(apply(res$full_results[4, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
tmp4B = t(apply(res$full_results[5, , ], 2, quantile, probs = c(0.025, 0.5, 0.975)))
colnames(tmp4A) = c("min", "med", "max")
colnames(tmp4B) = c("min", "med", "max")
out4 = data.frame(cbind(rbind(tmp4A, tmp4B), group = c(rep("A", t), rep("B", t))), Day = c(1:t, 1:t))
out4$min = as.numeric(out4$min)
out4$med = as.numeric(out4$med)
out4$max = as.numeric(out4$max)
S2d = add_theme(ggplot(out4 ,aes(x=Day, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(group)), alpha = 0.5) + geom_line(aes(color=factor(group))) + ylab("I(t), a = 4.0") + xlab("Day") + scale_fill_brewer("Group", palette = "Dark2")+ scale_color_brewer("Group", palette = "Dark2"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))))


S2 = ggarrange(S2a, S2b, S2c, S2d, nrow = 2, ncol = 2, common.legend = T, legend = "bottom")
ggsave(S2, file = "S2.eps", device = cairo_ps, width = 5, height = 5)
ggsave(S2, file = "S2.pdf", device = "pdf", width = 5, height = 5)

