library(data.table)
library(ggplot2)
N_B = 10000 # Population size ethnicity group B.
N_A = 90000 # Population size ethnicity group A.
N_tot = N_B + N_A

cmat1 = matrix(c(1, 1, 1, 1), nrow = 2) # Contact matrix for random mixing.
cmat1 = cmat1/rowSums(cmat1)

cmat125 = matrix(c(1.25, 1, 1, 1.25), nrow = 2) # Contact matrix for 25% assortative mixing.
cmat125 = cmat125/rowSums(cmat125)

cmat15 = matrix(c(1.5, 1, 1, 1.5), nrow = 2) # Contact matrix for 50% assortative mixing.
cmat15 = cmat15/rowSums(cmat15)

cmat175 = matrix(c(1.75, 1, 1, 1.75), nrow = 2) # Contact matrix for 75% assortative mixing. 
cmat175 = cmat175/rowSums(cmat175)

cmats = list(cmat1, cmat125, cmat15, cmat175)




run_model <- function(cimat, # Contact structure ethnicity, B first row, A second row.
                      crmat, # Contact stucture risk groups, high first row, low second row. 
                      R=NULL, # Either provide reproduction number R, or beta for the low-risk group (beta_l). 
                      beta=NULL, # Transmissibility in low-risk group (beta_l). 
                      P_Bh=0.5, # Proportion of ethnicity group B in high-risk group.
                      P_Ah = 0.1, # Proportion of ethnicity group A in low-risk group. 
                      rel_beta_h=2.0,# Relative transmissibility of high-risk group.
                      delta_t=0.1, # time step. 
                      gamma=1/3.0, # 1/duration of infectious period.
                      run_lm=T # boolean for whether a regression model should be ran, TRUE/FALSE.
){
  
  
  
  
  # Initialisation. Start with 1 infectious high-risk individual in each ethnnicity group. 
  Sprev_Bh = N_B * P_Bh
  N_Bh = Sprev_Bh
  Sprev_Bl = (N_B - Sprev_Bh) 
  N_Bl = Sprev_Bl
  Iprev_Bh = 1 
  Sprev_Bh = (Sprev_Bh - 1) 
  Iprev_Bl = 0


  Sprev_Ah = N_A * P_Ah 
  N_Ah = Sprev_Ah
  Sprev_Al = (N_A - Sprev_Ah)
  N_Al = Sprev_Al
  Iprev_Ah = 1 
  Sprev_Ah = (Sprev_Ah - 1) 
  Iprev_Al = 0

  # Compute beta from R if R is provided instead of beta.
  if(!is.null(R) & is.null(beta)){
    NGM <- matrix(1, nrow=4, ncol=4)
    NGM[1:2, 1:2] <- cimat[1,1] * crmat
    NGM[1:2, 3:4] <- cimat[1,2] * crmat
    NGM[3:4, 1:2] <- cimat[2,1] * crmat
    NGM[3:4, 3:4] <- cimat[2,2] * crmat
    NGM[1, ] <- NGM[1, ] * rel_beta_h
    NGM[3, ] <- NGM[3, ] * rel_beta_h
    NGM[,1] <- NGM[,1] * P_Ah
    NGM[,2] <- NGM[,2] * (1-P_Ah)
    NGM[,3] <- NGM[,3] * P_Bh
    NGM[,4] <- NGM[,4] * (1-P_Bh)
    
    R_beta1 <- Re(eigen(NGM, only.values=T)$values[1]/gamma)
    beta = R/R_beta1
  }
  print(beta)
  betal = beta  # beta in low-risk group.
  betah = rel_beta_h * betal # beta in high-risk group. 
  
  St_Ah = Sprev_Ah
  St_Al = Sprev_Al
  St_Bh = Sprev_Bh
  St_Bl = Sprev_Bl
  It_Ah = Iprev_Ah
  It_Al = Iprev_Al
  It_Bh = Iprev_Bh
  It_Bl = Iprev_Bl

  # Simulate the SIR equations using the forward Euler-method.
  while(Iprev_Bh + Iprev_Bl + Iprev_Bh + Iprev_Bl > 1){
    Iprev_I = Iprev_Bh + Iprev_Bl
    Iprev_N = Iprev_Ah + Iprev_Al
    
    S_Ah = Sprev_Ah - Sprev_Ah * delta_t * betah * (cimat[2, 1] * (crmat[1,1]*Iprev_Bh + crmat[1,2]*Iprev_Bl) / N_B + cimat[2, 2] * (crmat[1,1] * Iprev_Ah + crmat[1,2] * Iprev_Al) / N_A)
    S_Al = Sprev_Al - Sprev_Al * delta_t * betal * (cimat[2, 1] * (crmat[2,1]*Iprev_Bh + crmat[2,2]*Iprev_Bl) / N_B + cimat[2, 2] * (crmat[2,1] * Iprev_Ah + crmat[2,2] * Iprev_Al) / N_A)
    S_Bh = Sprev_Bh - Sprev_Bh * delta_t * betah * (cimat[1, 1] * (crmat[1,1]*Iprev_Bh + crmat[1,2]*Iprev_Bl) / N_B + cimat[1, 2] * (crmat[1,1] * Iprev_Ah + crmat[1,2] * Iprev_Al) / N_A)
    S_Bl = Sprev_Bl - Sprev_Bl * delta_t * betal * (cimat[1, 1] * (crmat[2,1]*Iprev_Bh + crmat[2,2]*Iprev_Bl) / N_B + cimat[1, 2] * (crmat[2,1] * Iprev_Ah + crmat[2,2] * Iprev_Al) / N_A)
    I_Ah = Iprev_Ah + Sprev_Ah * delta_t * betah * (cimat[2, 1] * (crmat[1,1]*Iprev_Bh + crmat[1,2]*Iprev_Bl) / N_B + cimat[2, 2] * (crmat[1,1] * Iprev_Ah + crmat[1,2] * Iprev_Al) / N_A) - delta_t * gamma * Iprev_Ah
    I_Al = Iprev_Al + Sprev_Al * delta_t * betal * (cimat[2, 1] * (crmat[2,1]*Iprev_Bh + crmat[2,2]*Iprev_Bl) / N_B + cimat[2, 2] * (crmat[2,1] * Iprev_Ah + crmat[2,2] * Iprev_Al) / N_A) - delta_t * gamma * Iprev_Al
    I_Bh = Iprev_Bh + Sprev_Bh * delta_t * betah * (cimat[1, 1] * (crmat[1,1]*Iprev_Bh + crmat[1,2]*Iprev_Bl) / N_B + cimat[1, 2] * (crmat[1,1] * Iprev_Ah + crmat[1,2] * Iprev_Al) / N_A) - delta_t * gamma * Iprev_Bh
    I_Bl = Iprev_Bl + Sprev_Bl * delta_t * betal * (cimat[1, 1] * (crmat[2,1]*Iprev_Bh + crmat[2,2]*Iprev_Bl) / N_B + cimat[1, 2] * (crmat[2,1] * Iprev_Ah + crmat[2,2] * Iprev_Al) / N_A) - delta_t * gamma * Iprev_Bl
    
    St_Ah = c(St_Ah, S_Ah)
    St_Al = c(St_Al, S_Al)
    St_Bh = c(St_Bh, S_Bh)
    St_Bl = c(St_Bl, S_Bl)
    It_Ah = c(It_Ah, I_Ah)
    It_Al = c(It_Al, I_Al)
    It_Bh = c(It_Bh, I_Bh)
    It_Bl = c(It_Bl, I_Bl)
  
    Sprev_Ah = S_Ah
    Sprev_Al = S_Al
    Sprev_Bh = S_Bh
    Sprev_Bl = S_Bl
    Iprev_Ah = I_Ah
    Iprev_Al = I_Al
    Iprev_Bh = I_Bh
    Iprev_Bl = I_Bl
  }
  
  R_Bh = N_Bh - Sprev_Bh # Total number of infected high-risk individuals of ethnicity group B. 
  R_Bl = N_Bl - Sprev_Bl # Total number of infected low-risk individuals of ethnicity group B. 

  R_Ah = N_Ah - Sprev_Ah # Total number of infected high-risk individuals of ethnicity group A. 
  R_Al = N_Al - Sprev_Al # Total number of infected low-risk individuals of ethnicity group A. 
  
  fit <- NULL
  if(run_lm){
  # Make data set for analysis.  
    infection = c(rep(1, round(R_Bh) + round(R_Bl)), rep(0, round(Sprev_Bh) + round(Sprev_Bl)), rep(1, round(R_Ah) + round(R_Al)), rep(0, round(Sprev_Ah) + round(Sprev_Al)))
    ethnicity = c(rep("B", N_B), rep("A", N_A))
    risk = c(rep("High", round(R_Bh)), rep("Low", round(R_Bl)), rep("High", round(Sprev_Bh)), rep("Low", round(Sprev_Bl)), rep("High", round(R_Ah)), rep("Low", round(R_Al)), rep("High", round(Sprev_Ah)), rep("Low", round(Sprev_Al)))
    
    data = data.frame("infection" = infection, "ethnicity" = ethnicity, "risk" = risk)
    # Fit linear regression model
    fit = lm(infection~factor(risk, levels = c("Low", "High")) + factor(ethnicity, levels = c("A", "B")), data = data)
  }

  Rt_Bl <-N_Bl - St_Bl
  Rt_Bh <-N_Bh - St_Bh
  Rt_Al <-N_Al - St_Al
  Rt_Ah <-N_Ah - St_Ah
  return(
    list(
      fractions=
        data.frame("Ah"=R_Ah/N_Ah,
                   "Al"=R_Al/N_Al,
                   "Bh"=R_Bh/N_Bh,
                   "Bl"=R_Bl/N_Bl,
                   "A"=(R_Ah+R_Al)/(N_Ah+N_Al),
                   "B"=(R_Bh+R_Bl)/(N_Bh+N_Bl),
                   "H"=(R_Ah+R_Bh)/(N_Ah+N_Bh),
                   "L"=(R_Al+R_Bl)/(N_Al+N_Bl)),
      fit=fit,
      full_results=data.frame(
        St_Ah = St_Ah,
        St_Al = St_Al,
        St_Bh = St_Bh,
        St_Bl = St_Bl,
        It_Ah = It_Ah/N_Ah,
        It_Al = It_Al/N_Al,
        It_Bh = It_Bh/N_Bh,
        It_Bl = It_Bl/N_Bl,
        B = (Rt_Bl + Rt_Bh) / (N_Bh+N_Bl),
        A = (Rt_Al + Rt_Ah) / (N_Ah+N_Al),
        H = (Rt_Bh + R_Ah) / (N_Bh+N_Ah),
        L = (Rt_Bl + R_Al) / (N_Bl+N_Al)
      )
    )
  )
}



add_theme <- function(q){
  q + fhiplot::theme_fhi_lines()+ theme(text = element_text(size=10))+
    scale_size_identity()
  
}


### Example runs of the model ###

#r <- run_model(matrix(1, nrow=2, ncol=2), cmats[[1]], beta=0.35, rel_beta_h = 2, P_Bh=0.5, P_Ah=0.2, run_lm=T)

#r2 <- run_model(matrix(1, nrow=2, ncol=2), cmats[[3]], beta=0.35, rel_beta_h = 1, P_Bh=0.5, P_Ah=0.2, run_lm=T)



# Figure 1a
results <- list()
for(cmat in cmats[1:4]){
     for(rel_beta_h in c(1, 1.5, 2, 2.5, 3)){
      crmat <- copy(cmat)
      r <- run_model(matrix(1, nrow=2, ncol=2), crmat, beta=0.35, rel_beta_h = rel_beta_h, P_Bh=0.5, P_Ah=0.1, run_lm=T)
      r2 <- run_model(matrix(1, nrow=2, ncol=2), crmat, beta=0.35, rel_beta_h = 1, P_Bh=0.5, P_Ah=0.1, run_lm=F) # rel_beta_h = 1, dvs. alle lik oppforsel, ingen i hoyrisiko

      
      results[[length(results) + 1]] <- data.frame(
        rel_beta_h=rel_beta_h,
        assortative=glue::glue("{(cmat[1,1]/cmat[1,2]-1)*100} %"),
        actual=r2$fractions[["H"]],
        expected= coef(r$fit)[1]
      )
    }
}
res <- rbindlist(results)
res$diff <- (res$actual - res$expected)/res$expected
options(OutDec = ".")
q1 <- add_theme(ggplot(res) + geom_line(aes(x=rel_beta_h, y=diff, color=assortative), size=1.2) + xlab(expression(beta[h]/beta[l])) + scale_y_continuous("Relative difference", labels = scales::percent)) 

#Figure 1b
results <- list()
for(cmat in cmats[1:4]){
    for(c in c(1,2,3,4,5)){
      crmat <- copy(cmat)
      crmat[1,1] <- crmat[1,1] * c
      r <- run_model(matrix(1, nrow=2, ncol=2), crmat, R=NULL, beta=0.2, rel_beta_h = 2, P_Bh=0.5, P_Ah=0.1, run_lm=F)
      results[[length(results) + 1]] <- data.frame(
        c=c,
        assortative=glue::glue("{(cmat[1,1]/cmat[1,2]-1)*100} %"),
        RR=r$fractions[["H"]]/r$fractions[["L"]],
        frac_L = r$fractions[["L"]]
      )
  }
}
res <- rbindlist(results)
q2 <- add_theme(ggplot(res) + geom_line(aes(x=c, y=frac_L, color=assortative, group=assortative), size=1.2) + xlab(expression(paste("High-high contacts ", c[hh]^{r}, sep = " "))) + scale_y_continuous("Prop. infected in low-risk group ", labels = scales::percent_format(decimal.mark=".")))

# Figure 1c
results <- list()
for(cmat in cmats[1:4]){
  for(rel_beta_h in c(1, 2, 3, 4, 5)){
    r <- run_model(cmat, matrix(c(1,1,1,1), nrow=2, ncol=2), 1.05,
                   rel_beta_h = rel_beta_h, P_Bh=0.5, P_Ah=0.1, run_lm=F)
    RR <- r$fractions[["H"]]/r$fractions[["L"]]
    explained <- (0.5 + 0.5*RR) / (0.9 + 0.1*RR)
    results[[length(results) + 1]] <- data.frame(
      unexplained=r$fractions[["B"]]/r$fractions[["A"]] - explained,
      rel_beta_h=rel_beta_h,
      assortative=glue::glue("{(cmat[1,1]/cmat[1,2]-1)*100} %")
    )
  }
}
res <- rbindlist(results)
options(OutDec = ".")
q3 <- add_theme(ggplot(res) + geom_line(aes(x=rel_beta_h, y=unexplained, color=assortative, group=assortative), size=1.2) + scale_y_continuous("Unexplained overrep. (UO)", labels = scales::percent_format(decimal.mark=".")) + xlab(expression(beta[h]/beta[l])))

#Figure 1d
results <- list()
for(cmat in cmats[1:4]){
  for(c in c(1, 2, 3, 4, 5)){
    r <- run_model(cmat, matrix(c(c,1,1,1), nrow=2, ncol=2), R=1.05, rel_beta_h = 1, P_Bh=0.5, P_Ah=0.1, run_lm=F)
    RR <- r$fractions[["H"]]/r$fractions[["L"]]
    explained <- (0.5 + 0.5*RR) / (0.9 + 0.1*RR)
    results[[length(results) + 1]] <- data.frame(
      unexplained=r$fractions[["B"]]/r$fractions[["A"]] - explained,
      high_contact=c,
      assortative=glue::glue("{(cmat[1,1]/cmat[1,2]-1)*100} %")
    )
  }
}

res <- rbindlist(results)
options(OutDec = ".")
q4 <- add_theme(ggplot(res) + geom_line(aes(x=high_contact, y=unexplained, color=assortative, group=assortative), size=1.2) + scale_y_continuous("Unexplained overrep. (UO)", labels = scales::percent_format(decimal.mark=".")) + xlab(expression(paste("High-high contacts ", c[hh]^{r}, sep = " "))))


#Combine plots
legend <- cowplot::get_legend(
  # Create some space to the left of the legend
  q1 + labs(color="Assortative mixing") + 
    theme(legend.position = "bottom")
)

plots <- cowplot::plot_grid(q1 + theme(legend.position="none"),
                   q2 + theme(legend.position="none"),
                   q3 + theme(legend.position="none"),
                   q4 + theme(legend.position="none"),
                   labels=c("a", "b", "c", "d"), label_y = 1.01, label_size = 12)


cowplot::plot_grid(plots, legend, nrow=2, rel_heights = c(1, 0.05))

ggsave("fig1.png", width = 7, height = 6)


#Create Table 1
results <- list()
# Model 1
for(cmat in cmats[1:4]){
  r <- run_model(cmat, matrix(c(1,1,1,1), nrow=2, ncol=2), beta = 0.35,
                 rel_beta_h = 2, P_Bh=0.5, P_Ah=0.1, run_lm=T)
    results[[length(results) + 1]] <- data.frame(
      pvalue_e = summary(r$fit)$coefficients[3,4],
      beta_e = summary(r$fit)$coefficients[3,1],
      pvalue_r = summary(r$fit)$coefficients[2,4],
      beta_r = summary(r$fit)$coefficients[2,1],
      model="model1",
      assortative=glue::glue("{(cmat[1,1]/cmat[1,2]-1)*100} %")
    )
}
# Model 2
for(cmat in cmats[1:4]){
    r <- run_model(cmat, matrix(c(2,1,1,1), nrow=2, ncol=2), beta = 0.35, rel_beta_h = 1, P_Bh=0.5, P_Ah=0.1, run_lm=T)
    results[[length(results) + 1]] <- data.frame(
      pvalue_e = summary(r$fit)$coefficients[3,4],
      beta_e = summary(r$fit)$coefficients[3,1],
      pvalue_r = summary(r$fit)$coefficients[2,4],
      beta_r = summary(r$fit)$coefficients[2,1],
      model="model2",
      assortative=glue::glue("{(cmat[1,1]/cmat[1,2]-1)*100} %")
    )
}
tabel1 <- rbindlist(results)
print(tabel1)

