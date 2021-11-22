##sum the numbers of effect alleles weighted by beta
score <- kfun_ukb[, .(eid)]
for (i in 1:nrow(score)) {
  tmp <- data.table(SNP = colnames(kfun_ukb)[-1],
                    num = as.numeric(t(kfun_ukb[i, -1])))
  tmp <- kfun_snp[tmp, on = 'SNP'][, score := num*beta.exposure]
  score$score_kfun[i] <- with(tmp, sum(score, na.rm = T)/sum(!is.na(num)))
}

##z-standardize
m <- mean(score$score_kfun); sd <- sd(score$score_kfun)
score[, score_kfun := (score_kfun - m)/sd]

##examine the variance explained by GRS
biomarker <- fread('D:/tanglei/3.Biological samples/Assay results_2_0.csv')
biomarker <- biomarker[, .(eid, 
                           creatinine_f30700_0_0,
                           creatinine_enzymatic_in_urine_f30510_0_0,
                           microalbumin_in_urine_f30500_0_0)]
colnames(biomarker) <- c('eid', 'scr', 'ucr', 'alb')
biomarker[, `:=` (age = geno$age,
                  sex = geno$sex,
                  ethnicity = rep(0, nrow(biomarker)),
                  creatinine = scr*0.01131,
                  uacr = alb/(ucr*0.0001131))]#scr in umol/L to mg/dl; ucr in mmol/L to g/L
biomarker[, `:=` (egfr = nephro::CKDEpi.creat(creatinine, sex, age, ethnicity),
                  ma = as.factor(ifelse(uacr >= 30, 1, 0)))]

biomarker <- cbind(biomarker, geno[, c(3, 5:14, 19)])[score, on = 'eid', nomatch = 0]

#kidney function
mod <- lm(log(egfr) ~ age + sex + rel3 + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, biomarker)
r0 <- summary(mod)$r.squared
mod <- lm(log(egfr) ~ score_kfun + age + sex + rel3 + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, biomarker)
r <- summary(mod)$r.squared
rsq <- r - r0

#ma
biomarker[, score_ma1 := cut(score_ma,
                             breaks = quantile(score_ma, 
                                               probs = c(0, 0.2, 0.8, 1)),
                             include.lowest = T,
                             labels = 1:3)]
mod <- glm(ma ~ score_ma1 + age + sex + rel3 + batch + pc1 + pc2 + pc3+ pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, biomarker, family = binomial())
rsq <- 1 - with(mod, deviance/null.deviance)
coeff <- summary(mod)$coeff[2:3, ]
or <- round(exp(coeff[, 1]), 2)
ci <- paste(round(exp(coeff[, 1] - 1.96*coeff[, 2]), 2),
            round(exp(coeff[, 1] + 1.96*coeff[, 2]), 2),
            '-')
p <- round(coeff[, 4], 3)

###individual-level MR
#scale effect size to per 10% decrement in eGFR
#a SD increment in GRS corresponds to ~3% increment in eGFR
nsd <- log(0.9)/summary(mod)$coeff[2, 1]

ca_all <- ca_all[score, on = 'eid']
ca_f <- ca_f[score, on = 'eid']
ca_m <- ca_m[score, on = 'eid']

grs_results <- data.table()
#all
for (i in 2:16) {
  ca <- as.factor(ca_all[[i]])
  mod <- glm(ca ~ score_kfun/nsd + age + sex + rel3 + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, ca_all, family = binomial())
  coeff <- summary(mod)$coeff[2, ]
  or1 <- round(exp(coeff[1]), 2)
  ci1 <- paste(round(exp(coeff[1] - 1.96*coeff[2]), 2),
              round(exp(coeff[1] + 1.96*coeff[2]), 2),
              '-')
  p1 <- round(coeff[4], 3)
  
  mod <- glm(ca ~ score_ma1 + age + sex + rel3 + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, ca_all, family = binomial())
  coeff <- summary(mod)$coeff[3, ]
  or2 <- round(exp(coeff[1]), 2)
  ci2 <- paste(round(exp(coeff[1] - 1.96*coeff[2]), 2),
              round(exp(coeff[1] + 1.96*coeff[2]), 2),
              '-')
  p2 <- round(coeff[4], 3)
  
  grs_results <- rbind(grs_results,
                       data.table(ca = c(colnames(ca_all)[i], NA),
                                  or = c(or1, or2), 
                                  ci = c(ci1, ci2), 
                                  p = c(p1, p2),
                                  exp = c('kidney function', 'ma')))
}
#female
for (i in 2:5) {
  ca <- as.factor(ca_f[[i]])
  mod <- glm(ca ~ score_kfun/nsd + age + rel3 + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, ca_f, family = binomial())
  coeff <- summary(mod)$coeff[2, ]
  or1 <- round(exp(coeff[1]), 2)
  ci1 <- paste(round(exp(coeff[1] - 1.96*coeff[2]), 2),
               round(exp(coeff[1] + 1.96*coeff[2]), 2),
               '-')
  p1 <- round(coeff[4], 3)
  
  mod <- glm(ca ~ score_ma1 + age + rel3 + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, ca_f, family = binomial())
  coeff <- summary(mod)$coeff[3, ]
  or2 <- round(exp(coeff[1]), 2)
  ci2 <- paste(round(exp(coeff[1] - 1.96*coeff[2]), 2),
               round(exp(coeff[1] + 1.96*coeff[2]), 2),
               '-')
  p2 <- round(coeff[4], 3)
  
  grs_results <- rbind(grs_results,
                       data.table(ca = c(colnames(ca_f)[i], NA),
                                  or = c(or1, or2), 
                                  ci = c(ci1, ci2), 
                                  p = c(p1, p2),
                                  exp = c('kidney function', 'ma')))
}
#male
mod <- glm(prostate ~ score_kfun/nsd + age + rel3 + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, ca_m, family = binomial())
coeff <- summary(mod)$coeff[2, ]
or1 <- round(exp(coeff[1]), 2)
ci1 <- paste(round(exp(coeff[1] - 1.96*coeff[2]), 2),
             round(exp(coeff[1] + 1.96*coeff[2]), 2),
             '-')
p1 <- round(coeff[4], 3)

mod <- glm(prostate ~ score_ma1 + age + rel3 + batch + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, ca_m, family = binomial())
coeff <- summary(mod)$coeff[3, ]
or2 <- round(exp(coeff[1]), 2)
ci2 <- paste(round(exp(coeff[1] - 1.96*coeff[2]), 2),
             round(exp(coeff[1] + 1.96*coeff[2]), 2),
             '-')
p2 <- round(coeff[4], 3)

grs_results <- rbind(grs_results,
                     data.table(ca = c('prostate', NA),
                                or = c(or1, or2), 
                                ci = c(ci1, ci2), 
                                p = c(p1, p2),
                                exp = c('kidney function', 'ma')))



