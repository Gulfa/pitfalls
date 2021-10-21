library(odin.dust)
library(ggplot2)
library(data.table)
library(dplyr)
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
    RR <- quantile(res$fractions[2,]  / res$fractions[1,], probs=c(0.2, 0.5, 0.8))
    
    out <- rbind(out, data.frame(min=RR[1],
                                 med=RR[2],
                                 max=RR[3],
                                 R=R,
                                 time=t))
  }
}


add_theme(ggplot(out ,aes(x=R, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(time)), alpha=0.5) + geom_line(aes(color=factor(time)))+ scale_y_continuous(trans='log10') + ylab("Relative Risk") + xlab("Reproduction Number (R0)")+ scale_fill_brewer("Time (days)", palette = "Dark2")+ scale_color_brewer("Time (days)", palette = "Dark2"))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))))

ggsave("fig1_RR.eps", device="cairo_ps")






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
  tot_inf <- quantile((res$full_results[6,,t]-90)/90000, probs=c(0.2, 0.5, 0.8))
  mean(tot_inf)
  out <- rbind(out, data.frame(min=tot_inf[1],
                               med=tot_inf[2],
                               max=tot_inf[3],
                               bc=bc,
                               sus=sus))
}

add_theme(ggplot(out, aes(x=sus, y=med, ymin=min, ymax=max)) + geom_line() + geom_ribbon(alpha=0.3) + ylab("Fraction infected in low-risk group") + xlab("Susceptibility in High-Risk group") + scale_y_continuous(labels = scales::percent_format(accuracy = 1)))
ggsave("fig2_b.png")



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
  
  tot_inf <- quantile(actual$full_results[6,,t], probs=c(0.1, 0.5, 0.8))
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
    preds <- quantile(vals, probs=c(0.2, 0.5, 0.8))
    out <- rbind(out, data.frame(min_p=preds[1],
                                 med_p=preds[2],
                                 max_p=preds[3],
                                 sim_med=actual$full_results[6,i,t],
                                 sim=i,
                                 sus=sus)
                 )
  }
}

add_theme(ggplot(out, aes(x=sus, y=med_p, ymin=min_p, ymax=max_p, color="Prediction")) + geom_pointrange( position="jitter") + ylab("Total cases") + xlab("Susceptibility") +geom_point(aes(x=sus, y=sim_med, color="Simulation"), size=2, position="jitter") + scale_color_brewer("Model", palette="Dark2"))
ggsave("fig_2a.png")




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






add_theme(ggplot(out, aes(x=factor(assort), y=mid, ymin=lower, ymax=upper)) + geom_point(position="jitter") + geom_pointrange( position="jitter") + ylab("Regression coefficent ethnicity (beta_e)") + xlab("Assortativ miksing"))
ggsave("fig3_1_beta.png")



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
    unexplained <-quantile(RR_eth - explained, probs=c(0.2, 0.5, 0.8))
    out <- rbind(out,data.frame(min=unexplained[1],
                                med=unexplained[2],
                                max=unexplained[3],
                                sus=sus,
                                assort=assort))
                       
                 
  }
}

add_theme(ggplot(out ,aes(x=sus, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(assort)), alpha=0.5) + geom_line(aes(color=factor(assort)))+   ylab("Unexplained Overrepresentation") + xlab("Suceptibility in High Risk Group")+ scale_fill_brewer("Assortative Mixing", palette = "Dark2")+ scale_color_brewer("Assortative Mixing", palette = "Dark2")+ scale_y_continuous(labels = scales::percent_format(accuracy = 1)))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))
ggsave("fig3_1_uo.png")






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






add_theme(ggplot(out, aes(x=factor(assort), y=mid, ymin=lower, ymax=upper)) + geom_point(position="jitter") + geom_pointrange( position="jitter") + ylab("Regression coefficent ethnicity (beta_e)") + xlab("Assortativ miksing"))

ggsave("fig4_1_beta.png")



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
    unexplained <-quantile(RR_eth - explained, probs=c(0.2, 0.5, 0.8))
    out <- rbind(out,data.frame(min=unexplained[1],
                                med=unexplained[2],
                                max=unexplained[3],
                                hc=hc,
                                assort=assort))
                       
                 
  }
}

add_theme(ggplot(out ,aes(x=hc, y=med, ymin=min, ymax=max)) + geom_ribbon(aes(fill=factor(assort)), alpha=0.5) + geom_line(aes(color=factor(assort)))+   ylab("Unexplained Overrepresentation") + xlab("High-High contact")+ scale_fill_brewer("Assortative Mixing", palette = "Dark2")+ scale_color_brewer("Assortative Mixing", palette = "Dark2")+ scale_y_continuous(labels = scales::percent_format(accuracy = 1)))# + scale_x_continuous(expand = expand_scale(mult = c(0, 0.1))
ggsave("fig4_1_uo.png")
