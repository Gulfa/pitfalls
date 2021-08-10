N_I = 10000 # Antall innvandrere
N_N = 90000 # Antall etniske nordmenn
N_tot = N_I + N_N

cmat05 = matrix(c(0.5, 1, 1, 0.5), nrow = 2)
cmat05 = cmat05/rowSums(cmat05)

cmat1 = matrix(c(1, 1, 1, 1), nrow = 2) # Kontaktmatrise for assortativ miksing
cmat1 = cmat1/rowSums(cmat1)

cmat125 = matrix(c(1.25, 1, 1, 1.25), nrow = 2) # Kontaktmatrise for assortativ miksing
cmat125 = cmat125/rowSums(cmat125)

cmat15 = matrix(c(1.5, 1, 1, 1.5), nrow = 2) # Kontaktmatrise for assortativ miksing
cmat15 = cmat15/rowSums(cmat15)

cmat175 = matrix(c(1.75, 1, 1, 1.75), nrow = 2)
cmat175 = cmat175/rowSums(cmat175)

cmats = list(cmat05, cmat1, cmat125, cmat15, cmat175)




run_model <- function(cimat,
                      crmat,
                      R=NULL,
                      beta=NULL,
                      P_Ih=0.5,
                      P_Nh = 0.2,
                      rel_beta_h=2.0,
                      delta_t=0.1,
                      gamma=1/3.0,
                      run_lm=T
                      ){
  
  
  
  
                                        # Initialisering. Starter med 1 smittsom etnisk nordmann i høyrisikogruppe, og 1 smittsom innvandrer i høyrisikogruppe
  Sprev_Ih = N_I * P_Ih
  N_Ih = Sprev_Ih
  Sprev_Il = (N_I - Sprev_Ih) 
  N_Il = Sprev_Il
  Iprev_Ih = 1 
  Sprev_Ih = (Sprev_Ih - 1) 
  Iprev_Il = 0


  Sprev_Nh = N_N * P_Nh 
  N_Nh = Sprev_Nh
  Sprev_Nl = (N_N - Sprev_Nh)
  N_Nl = Sprev_Nl
  Iprev_Nh = 1 
  Sprev_Nh = (Sprev_Nh - 1) 
  Iprev_Nl = 0

  if(!is.null(R) & is.null(beta)){
    NGM <- matrix(1, nrow=4, ncol=4)
    NGM[1:2, 1:2] <- cimat[1,1]*crmat
    NGM[1:2, 3:4] <- cimat[1,2]*crmat
    NGM[3:4, 1:2] <- cimat[2,1]*crmat
    NGM[3:4, 3:4] <- cimat[2,2]*crmat
    NGM[1, ] <- NGM[1, ]*rel_beta_h
    NGM[3, ] <- NGM[3, ]*rel_beta_h
    NGM[,1] <- NGM[,1]*P_Nh
    NGM[,2] <- NGM[,2]*(1-P_Nh)
    NGM[,3] <- NGM[,3]*P_Ih
    NGM[,4] <- NGM[,4]*(1-P_Ih)
    
    R_beta1 <- Re(eigen(NGM, only.values=T)$values[1]/gamma)
    beta = R/R_beta1
  }
  print(beta)
  betal = beta
  betah = rel_beta_h * betal # beta i høyrisikogruppe. 
  
  St_Nh = Sprev_Nh
  St_Nl = Sprev_Nl
  St_Ih = Sprev_Ih
  St_Il = Sprev_Il
  It_Nh = Iprev_Nh
  It_Nl = Iprev_Nl
  It_Ih = Iprev_Ih
  It_Il = Iprev_Il
  
  I_Ih_l = list()
  I_Il_l = list()
  I_Nh_l = list()
  I_Nl_l = list()
  
  
  Sprev_Ih = N_I * P_Ih
  N_Ih = Sprev_Ih
  Sprev_Il = (N_I - Sprev_Ih) 
  N_Il = Sprev_Il
  Iprev_Ih = 1 
  Sprev_Ih = (Sprev_Ih - 1) 
  Iprev_Il = 0
  
  
  Sprev_Nh = N_N * P_Nh 
  N_Nh = Sprev_Nh
  Sprev_Nl = (N_N - Sprev_Nh)
  N_Nl = Sprev_Nl
  Iprev_Nh = 1 
  Sprev_Nh = (Sprev_Nh - 1) 
  Iprev_Nl = 0
  
  St_Nh = Sprev_Nh
  St_Nl = Sprev_Nl
  St_Ih = Sprev_Ih
  St_Il = Sprev_Il
  It_Nh = Iprev_Nh
  It_Nl = Iprev_Nl
  It_Ih = Iprev_Ih
  It_Il = Iprev_Il
  
  while(Iprev_Ih + Iprev_Il + Iprev_Ih + Iprev_Il > 1){
    Iprev_I = Iprev_Ih + Iprev_Il
    Iprev_N = Iprev_Nh + Iprev_Nl
    
    S_Nh = Sprev_Nh - Sprev_Nh * delta_t * betah * (cimat[2, 1] * (crmat[1,1]*Iprev_Ih + crmat[1,2]*Iprev_Il) / N_I + cimat[2, 2] * (crmat[1,1]*Iprev_Nh + crmat[1,2]*Iprev_Nl) / N_N)
    S_Nl = Sprev_Nl - Sprev_Nl * delta_t * betal * (cimat[2, 1] * (crmat[2,1]*Iprev_Ih + crmat[2,2]*Iprev_Il) / N_I + cimat[2, 2] * (crmat[2,1]*Iprev_Nh + crmat[2,2]*Iprev_Nl) / N_N)
    S_Ih = Sprev_Ih - Sprev_Ih * delta_t * betah * (cimat[1, 1] * (crmat[1,1]*Iprev_Ih + crmat[1,2]*Iprev_Il) / N_I + cimat[1, 2] *(crmat[1,1]*Iprev_Nh + crmat[1,2]*Iprev_Nl)  / N_N)
    S_Il = Sprev_Il - Sprev_Il * delta_t * betal * (cimat[1, 1] * (crmat[2,1]*Iprev_Ih + crmat[2,2]*Iprev_Il)  / N_I + cimat[1, 2] *  (crmat[2,1]*Iprev_Nh + crmat[2,2]*Iprev_Nl)/ N_N)
    I_Nh = Iprev_Nh + Sprev_Nh * delta_t * betah * (cimat[2, 1] * (crmat[1,1]*Iprev_Ih + crmat[1,2]*Iprev_Il) / N_I + cimat[2, 2] * (crmat[1,1]*Iprev_Nh + crmat[1,2]*Iprev_Nl) / N_N) - delta_t * gamma * Iprev_Nh
    I_Nl = Iprev_Nl + Sprev_Nl * delta_t * betal * (cimat[2, 1] * (crmat[2,1]*Iprev_Ih + crmat[2,2]*Iprev_Il) / N_I + cimat[2, 2] * (crmat[2,1]*Iprev_Nh + crmat[2,2]*Iprev_Nl) / N_N) - delta_t * gamma * Iprev_Nl
    I_Ih = Iprev_Ih + Sprev_Ih * delta_t * betah * (cimat[1, 1] * (crmat[1,1]*Iprev_Ih + crmat[1,2]*Iprev_Il) / N_I + cimat[1, 2] * (crmat[1,1]*Iprev_Nh + crmat[1,2]*Iprev_Nl) / N_N) - delta_t * gamma * Iprev_Ih
    I_Il = Iprev_Il + Sprev_Il * delta_t * betal * (cimat[1, 1] * (crmat[2,1]*Iprev_Ih + crmat[2,2]*Iprev_Il) / N_I + cimat[1, 2] * (crmat[2,1]*Iprev_Nh + crmat[2,2]*Iprev_Nl) / N_N)  - delta_t * gamma * Iprev_Il
    
    St_Nh = c(St_Nh, S_Nh)
    St_Nl = c(St_Nl, S_Nl)
    St_Ih = c(St_Ih, S_Ih)
    St_Il = c(St_Il, S_Il)
    It_Nh = c(It_Nh, I_Nh)
    It_Nl = c(It_Nl, I_Nl)
    It_Ih = c(It_Ih, I_Ih)
    It_Il = c(It_Il, I_Il)
  
    Sprev_Nh = S_Nh
    Sprev_Nl = S_Nl
    Sprev_Ih = S_Ih
    Sprev_Il = S_Il
    Iprev_Nh = I_Nh
    Iprev_Nl = I_Nl
    Iprev_Ih = I_Ih
    Iprev_Il = I_Il
  }
  
#  I_Ih_l[[i]] = It_Ih/N_Ih
#  I_Il_l[[i]] = It_Il/N_Il
#  I_Nh_l[[i]] = It_Nh/N_Nh
#  I_Nl_l[[i]] = It_Nl/N_Nl


  R_Ih = N_Ih - Sprev_Ih # Totalt antall smittede høyrisikoinnvandrere
  R_Il = N_Il - Sprev_Il # Totalt antall smittede lavrisikoinnvandrere

  R_Nh = N_Nh - Sprev_Nh # Totalt antall smittede høyrisiko-etniske nordmenn
  R_Nl = N_Nl - Sprev_Nl # Totalt antall smittede lavrisiko-etniske nordmenn
  
  fit <- NULL
  if(run_lm){
  # Lag datasett for analyse 
    infection = c(rep(1, round(R_Ih) + round(R_Il)), rep(0, round(Sprev_Ih) + round(Sprev_Il)), rep(1, round(R_Nh) + round(R_Nl)), rep(0, round(Sprev_Nh) + round(Sprev_Nl)))
    ethnicity = c(rep("Immigrant", N_I), rep("Norwegian", N_N))
    risk = c(rep("High", round(R_Ih)), rep("Low", round(R_Il)), rep("High", round(Sprev_Ih)), rep("Low", round(Sprev_Il)), rep("High", round(R_Nh)), rep("Low", round(R_Nl)), rep("High", round(Sprev_Nh)), rep("Low", round(Sprev_Nl)))
    

    data = data.frame("infection" = infection, "ethnicity" = ethnicity, "risk" = risk)
    
                                        # Tilpasning lineær regresjonsmodell
    fit = lm(infection~factor(risk, levels = c("Low", "High")) + factor(ethnicity, levels = c("Norwegian", "Immigrant")), data = data)
    print(summary(fit))
    
    fit_log = glm(infection~factor(risk, levels = c("Low", "High")) + factor(ethnicity, levels = c("Norwegian", "Immigrant")), family = "binomial")
  }
  
 # fits[[i]] = fit
                                        #fits_log[[i]] = fit_log

  Rt_Il <-N_Il - St_Il
  Rt_Ih <-N_Ih - St_Ih
  Rt_Nl <-N_Nl - St_Nl
  Rt_Nh <-N_Nh - St_Nh
  return(
    list(
      fractions=
        data.frame("Nh"=R_Nh/N_Nh,
                   "Nl"=R_Nl/N_Nl,
                   "Ih"=R_Ih/N_Ih,
                   "Il"=R_Il/N_Il,
                   "N"=(R_Nh+R_Nl)/(N_Nh+N_Nl),
                   "I"=(R_Ih+R_Il)/(N_Ih+N_Il),
                   "H"=(R_Nh+R_Ih)/(N_Nh+N_Ih),
                   "L"=(R_Nl+R_Il)/(N_Nl+N_Il)),
      fit=fit,
      full_results=data.frame(
        St_Nh = St_Nh,
        St_Nl = St_Nl,
        St_Ih = St_Ih,
        St_Il = St_Il,
        It_Nh = It_Nh/N_Nh,
        It_Nl = It_Nl/N_Nl,
        It_Ih = It_Ih/N_Ih,
        It_Il = It_Il/N_Il,
        I = (Rt_Il + Rt_Ih) / (N_Ih+N_Il),
        N = (Rt_Nl + Rt_Nh) / (N_Nh+N_Nl),
        H = (Rt_Ih + R_Nh) / (N_Ih+N_Nh),
        L = (Rt_Il + R_Nl) / (N_Il+N_Nl)
      )
    )
  )
}



add_theme <- function(q){
  q + fhiplot::theme_fhi_lines()+ theme(text = element_text(size=22))+
    scale_size_identity()
  
}



results <- list()
for(cmat in cmats[2:5]){
  for(R in c(1.05)){
    for(rel_beta_h in c(1,1.5,2, 2.5, 3)){
                                        #    cmat[1,1] <- cmat[1,1]*1.5
      crmat <- copy(cmat)
      r <- run_model(matrix(1, nrow=2, ncol=2), crmat,beta=0.34, rel_beta_h = rel_beta_h, P_Ih=0.5, P_Nh=0.2, run_lm=F)
      results[[length(results) + 1]] <- data.frame(
        R=R,
        rel_beta_h=rel_beta_h,
        assortative=glue::glue("{(cmat[1,1]/cmat[1,2]-1)*100} %"),
        RR=r$fractions[["H"]]/r$fractions[["L"]],
        frac_L = r$fractions[["L"]],
        frac_H = r$fractions[["H"]]
      )
    }
  }
}
res <- rbindlist(results)
res$assortative <- factor(res$assortative)

q1 <- add_theme(ggplot(res) + geom_line(aes(x=rel_beta_h, y=RR, color=assortative, group=assortative), size=2) + ylab("Risk ratio") + xlab("Relative beta") + geom_line(data=data.frame(x=1:3, y=1:3), aes(x=x, y=y)))
q1
ggsave("1a.png")
results <- list()
for(cmat in cmats[2:5]){
  for(R in c(1.05)){
    for(c in c(1,2,3,4,5)){
                                        #    cmat[1,1] <- cmat[1,1]*1.5
      crmat <- copy(cmat)
      crmat[1,1] <- crmat[1,1] * c
      r <- run_model(matrix(1, nrow=2, ncol=2), crmat,R=NULL, beta=0.2, rel_beta_h = 2, P_Ih=0.5, P_Nh=0.2, run_lm=F)
      results[[length(results) + 1]] <- data.frame(
        R=R,
        c=c,
        assortative=glue::glue("{(cmat[1,1]/cmat[1,2]-1)*100} %"),
        RR=r$fractions[["H"]]/r$fractions[["L"]],
        frac_L = r$fractions[["L"]]
      )
    }
  }
}
res <- rbindlist(results)
res$assortative <- factor(res$assortative)

q2 <- add_theme(ggplot(res) + geom_line(aes(x=c, y=frac_L, color=assortative, group=assortative), size=2) + xlab("Antall høy-høy kontakter") + scale_y_continuous("Andel smittet lav-risiko", labels = scales::percent))


results <- list()
for(cmat in cmats[2:5]){
  for(rel_beta_h in c(1, 2, 3, 4, 5)){
    r <- run_model(cmat, matrix(c(1,1,1,1), nrow=2, ncol=2), 1.05,
                   rel_beta_h = rel_beta_h, P_Ih=0.5, P_Nh=0.2, run_lm=F)
    RR <- r$fractions[["H"]]/r$fractions[["L"]]
    explained <- (0.5 + 0.5*RR) / (0.8 + 0.2*RR)
    results[[length(results) + 1]] <- data.frame(
      unexplained=r$fractions[["I"]]/r$fractions[["N"]] - explained,
      rel_beta_h=rel_beta_h,
      assortative=glue::glue("{(cmat[1,1]/cmat[1,2]-1)*100} %")
    )
  }
}
res <- rbindlist(results)
res$assortative <- factor(res$assortative)
q3 <- add_theme(ggplot(res) + geom_line(aes(x=rel_beta_h, y=unexplained, color=assortative, group=assortative), size=2) + scale_y_continuous("Uforklart overrepresentasjon", labels = scales::percent) + xlab("Relative beta høy risiko"))

results <- list()
for(cmat in cmats[2:5]){
  for(c in c(1, 2, 3, 4, 5)){
    r <- run_model(cmat, matrix(c(c,1,1,1), nrow=2, ncol=2), R=1.05, rel_beta_h = 1, P_Ih=0.5, P_Nh=0.2, run_lm=F)
    RR <- r$fractions[["H"]]/r$fractions[["L"]]
    explained <- (0.5 + 0.5*RR) / (0.8 + 0.2*RR)
    results[[length(results) + 1]] <- data.frame(
      unexplained=r$fractions[["I"]]/r$fractions[["N"]] - explained,
      high_contact=c,
      assortative=glue::glue("{(cmat[1,1]/cmat[1,2]-1)*100} %")
    )
  }
}

res <- rbindlist(results)
res$assortative <- factor(res$assortative)
q4 <- add_theme(ggplot(res) + geom_line(aes(x=high_contact, y=unexplained, color=assortative, group=assortative), size=2) + scale_y_continuous("Uforklart overrepresentasjon", labels = scales::percent) + xlab("Antall høy-høy kontakter"))

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  q1 + labs(color="Assortativ miksing") + 
    theme(legend.position = "bottom")
)

plots <- cowplot::plot_grid(q1 + theme(legend.position="none"),
                   q2 + theme(legend.position="none"),
                   q3 + theme(legend.position="none"),
                   q4 + theme(legend.position="none"),
                   labels=c("a", "b", "c", "d"))


cowplot::plot_grid(plots, legend, nrow=2, rel_heights = c(1, 0.05))

ggsave("fig1.png")


results <- list()
for(cmat in cmats[2:5]){
  r <- run_model(cmat, matrix(c(1,1,1,1), nrow=2, ncol=2), 1.05,
                 rel_beta_h = rel_beta_h, P_Ih=0.5, P_Nh=0.2, run_lm=T)
    results[[length(results) + 1]] <- data.frame(
      pvalue_e = summary(r$fit)$coefficients[3,4],
      beta_e = summary(r$fit)$coefficients[3,1],
      pvalue_r = summary(r$fit)$coefficients[2,4],
      beta_r = summary(r$fit)$coefficients[2,1],
      model="model1",
      assortative=glue::glue("{(cmat[1,1]/cmat[1,2]-1)*100} %")
    )
}
for(cmat in cmats[2:5]){
    r <- run_model(cmat, matrix(c(3,1,1,1), nrow=2, ncol=2), R=1.05, rel_beta_h = 1, P_Ih=0.5, P_Nh=0.2, run_lm=T)
    results[[length(results) + 1]] <- data.frame(
      pvalue_e = summary(r$fit)$coefficients[3,4],
      beta_e = summary(r$fit)$coefficients[3,1],
      pvalue_r = summary(r$fit)$coefficients[2,4],
      beta_r = summary(r$fit)$coefficients[2,1],
      model="model2",
      assortative=glue::glue("{(cmat[1,1]/cmat[1,2]-1)*100} %")
    )
}

rbindlist(results)

## pvals = c()
## coefs = c()
## for(r in results){
##   pvals = c(pvals, summary(r$fit)$coefficients[3, 4])
##   coefs = c(coefs, summary(r$fit)$coefficients[3, 1])
## }


## tab <- list()
## for(r in results){
##   hl = r$fractions[["H"]]/r$fractions[["L"]]
##   row <- data.frame(rel_beta_h=r$rel_beta_h,
##                     high_risk_contancts=r$high_risk_contacts,
##                     assortative=glue::glue("{(r$assortative - 1)*100} %"),
##                     norwegian=glue::glue('{round(r$fractions[["N"]], 2)*100} %'),
##                     Imigrant=glue::glue('{round(r$fractions[["I"]], 2)*100} %'),
##                     high_risk=glue::glue('{round(r$fractions[["H"]], 2)*100} %'),
##                     low_risk=glue::glue('{round(r$fractions[["L"]], 2)*100} %'),
##                     high_low_RR=round(hl,2),
##                     over_representation=glue::glue('{round(r$fractions[["I"]] / r$fractions[["N"]] - 1,2)*100} %'),
##                     explained_overrepresentation=glue::glue('{round((0.5+0.5*hl)/(0.8+0.2*hl) - 1,2)*100} %'),
##                     unexplained_overrepresentation=glue::glue('{round((r$fractions[["I"]] / r$fractions[["N"]] - (0.5+0.5*hl)/(0.8+0.2*hl)),2)*100} %')
##                     )
##   tab[[length(tab) + 1]] <- row
## }
## tab <- rbindlist(tab)


## library(kableExtra)
## tbl <- kableExtra::kable(tab, "html",
##                          col.names = c("Relativ smittsomhet",
##                                        "Høy risiko kontakter",
##                                        "Assortativ miksing",
##                                        "Andel smitted norske",
##                                        "Andel smitted innvandrere",
##                                        "Andel smitted høyrisiko",
##                                        "Andel smitted lavrisiko",
##                                        "Relative riskio høyrisiko",
##                                        "Økt forekomst innvandrere",
##                                         "Forklart økt forekomst",
##                                        "Uforklart økt forekomst"
##                                        )
##                                        ) %>%

##   kable_styling() 
## save_kable(tbl, file="table.html"
##            )


## for(r in results){
##   print(summary(r$fit))
## }




## png("hoylavriskI.png", width = 4.7, height = 3.55, units = "in", res = 800)
## par(mfrow = c(2,2), xaxs = "i", ps = 8.5, cex = 1, cex.main = 1, mgp = c(0.9,0.4,0), mar = c(2.0, 1.5, 0.55, 0.6))

## x2 = (1:length(results[[2]]$full_results$It_Nl))/10 
## x3 = (1:length(results[[3]]$full_results$It_Nl))/10
## x4 = (1:length(results[[4]]$full_results$It_Nl))/10 
## x5 = (1:length(results[[5]]$full_results$It_Nl))/10 
## plot(x = x2, y = results[[2]]$full_results$It_Nh, xlab = "", ylab = "Andel smittede", lwd = 2, type = "l", ylim = c(0, 0.052))
## title("a)", adj = 0, font.main = 1)
## lines(x = x3, y = results[[3]]$full_results$It_Nh, lwd = 2, col = "skyblue")
## lines(x = x4, y = results[[4]]$full_results$It_Nh, lwd = 2, col = "darkorange")
## lines(x = x5, y = results[[5]]$full_results$It_Nh, lwd = 2, col = "pink")


## plot(x = x2, y = results[[2]]$full_results$It_Nl, xlab = "", ylab = "", lwd = 2, type = "l", ylim = c(0, 0.052))
## title("b)", adj = 0, font.main = 1)
## lines(x = x3, y = results[[3]]$full_results$It_Nl, lwd = 2, col = "skyblue")
## lines(x = x4, y = results[[4]]$full_results$It_Nl, lwd = 2, col = "darkorange")
## lines(x = x5, y = results[[5]]$full_results$It_Nl, lwd = 2, col = "pink")
## legend("top", ncol = 2, c("Tilfeldig miksing", "25% assortativ", "50% assortativ", "75% assortativ"), col = c("black", "skyblue", "darkorange", "pink"), lwd = 2, bty = "n", x.intersp=0.1, y.intersp=0.8, seg.len = 0.5,)


## plot(x = x2, y = results[[2]]$full_results$It_Ih, xlab = "Dager", ylab = "Andel smittede", lwd = 2, type = "l", ylim = c(0, 0.052))
## title("c)", adj = 0, font.main = 1)
## lines(x = x3, y = results[[3]]$full_results$It_Ih, lwd = 2, col = "skyblue")
## lines(x = x4, y = results[[4]]$full_results$It_Ih, lwd = 2, col = "darkorange")
## lines(x = x5, y = results[[5]]$full_results$It_Ih, lwd = 2, col = "pink")

## plot(x = x2, y =results[[2]]$full_results$ It_Il, xlab = "Dager", ylab = "", lwd = 2, type = "l", ylim = c(0, 0.052))
## title("d)", adj = 0, font.main = 1)
## lines(x = x3, y = results[[3]]$full_results$It_Il, lwd = 2, col = "skyblue")
## lines(x = x4, y = results[[4]]$full_results$It_Il, lwd = 2, col = "darkorange")
## lines(x = x5, y = results[[5]]$full_results$It_Il, lwd = 2, col = "pink")
## dev.off()


## # Investigations
## res <- list()
## for(ci in c(1, 1.25, 1.5, 1.75)){
##   for(rel_beta_h in c(1, 2, 5)){
##    # for(inc_cont in c(1, 1.5, 2)){
##     #  for(assort in c(1, 1.5, 2)){
##         for(R in c(1.1, 1.15, 1.2)){
##           cimat <- matrix(c(ci, 1, 1, ci), nrow = 2) 
##           cimat <- cimat / rowSums(cimat)
##           crmat <- matrix(c(assort, 1, 1,assort), nrow = 2) 
##           crmat <- crmat / rowSums(crmat) * inc_cont
##           crmat <- matrix(c(1,1,1,1), nrow=2)
##           res[[length(res) +1]] <- run_model(cimat,
##                                              crmat,
##                                              R=R/(0.7+ 0.3*rel_beta_h)/inc_cont,
##                                              rel_beta_h=rel_beta_h
##                                              ) %>% mutate(
##                                                      ci=ci,
##                                                      rel_beta_h=rel_beta_h,
##                                                      inc_cont=1, #, inc_cont,
##                                                      assort=1, #assort,
##                                                      R=R)$fractions
##         }
##      # }
##    # }
##   }
## }



## res <- data.table::rbindlist(res)

## res[, frac:=I/N]
## res[, exp_frac:=(0.5+0.5*rel_beta_h) / (0.8+0.2*rel_beta_h)]

## res[, unexp_frac:=frac/exp_frac]

## library(ggplot2)
## ggplot(res %>% filter(assort==1 & inc_cont==1)) + geom_line(aes(x=ci, y=frac, color=factor(R))) + facet_wrap(rel_beta_h~.)

## ggplot(res) + geom_line(aes(x=ci, y=frac, color=factor(R))) + facet_wrap(rel_beta_h~.)
