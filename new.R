library(odin.dust)
library(ggplot2)
library(data.table)
library(dplyr)
source("model.R")

out <- data.frame()
for (b in seq(0.05, 0.6, by=0.025)){
  for(t in c(50, 100, 200)){
    res <- run_model(diag(c(1, 9)), rep(b, t), c(90000, 10000), t, c(90,10), 1000, 1,
                     susceptibility = c(1,1.2))
    RR <- quantile(res$fractions[2,]  / res$fractions[1,], probs=c(0.2, 0.5, 0.8))
    
    out <- rbind(out, data.frame(min=RR[1],
                                 med=RR[2],
                                 max=RR[3],
                                 beta=b,
                                 time=t))
  }
}
  
ggplot(out ,aes(x=beta, y=med, ymin=min, ymax=max, fill=factor(time), color=factor(time))) + geom_ribbon(alpha=0.3) + geom_line()+ scale_y_continuous(trans='log10') + ylab("RR (B/A)")#+ geom_line(data=data.frame(x=c(0.08, 0.3), y=c(1.2, 1.2)), aes(x=x, y=y))
ggsave("fig1_RR.png")


out <- data.frame()

N <- c(90000, 10000)
input_mats <- cij_NGM(matrix(1,nrow=2, ncol=2), N, c(1,1), c(1,1))
beta <- 1.3/input_mats$beta_R

for(bc in c(0.1, 0.5, 1)){
  for(sus in seq(1, 10)){
    t <- 400
    b <- 0.05
    c_ij <- matrix(c(1,bc,bc,1), nrow=2)
    N <- c(90000, 10000)
    input_mats <- cij_NGM(c_ij, N, c(1,sus), c(1,1))
    res <- run_model(input_mats$c_ij, rep(beta, t), N, t,
                     c(90,10), 400, 1, susceptibility = c(1,sus))
    tot_inf <- quantile(res$full_results[6,,t]/90000, probs=c(0.2, 0.5, 0.8))
    mean(tot_inf)
    out <- rbind(out, data.frame(min=tot_inf[1],
                                 med=tot_inf[2],
                                 max=tot_inf[3],
                                 bc=bc,
                                 sus=sus))
  }
}

ggplot(out, aes(x=sus, y=med, ymin=min, ymax=max, fill=factor(bc), color=factor(bc))) + geom_line() + geom_ribbon(alpha=0.3) + ylab("Fraction infected in low-risk group") + xlab("Increased susceptibility in high-risk group")
ggsave("fig2_RR.png")





out <- data.frame()
for( sus in c(2)){
  for(assort in c(0.1,0.5, 1)){
    res <- run_regs(matrix(c(1,assort, assort, 1), nrow=2),
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

ggplot(out, aes(x=assort, y=mid, ymin=lower, ymax=upper, color=factor(sus))) + geom_point(position="jitter") + geom_pointrange( position="jitter") + ylab("RR(N)") + xlab("Off-diagonal element in matrix")
ggsave("fig3_betas.png")
