library(TwoSampleMR)
library(MendelianRandomization)
##kidney function/microalbuminuria
ca <- c('lu', 'uri', 'ki', 'thy', 'eso', 'stom', 'colo', 'liv', 'panc', 'mela', 'nmela', 'nhl', 'leu', 'mm', 'all', 'brea', 'endo', 'ova', 'cerv', 'pros', 'ki_f', 'ki_m')
mr_results <- data.table()
for (i in ca) {
  tmp <- fread(paste0('~/ckd_ca/glm/plink2.', i, '.glm.logistic.hybrid'))
  mr_out <- tmp[ID %in% kfun_snp$SNP][, beta := log(OR)]
  setnames(mr_out, 'ID', 'SNP')
  mr_out <- kfun_snp[, .(SNP, eaf.exposure)][mr_out, on = 'SNP']
  mr_out <- format_data(mr_out,
                        type = 'outcome',
                        snp_col = 'SNP',
                        chr_col = '#CHROM',
                        effect_allele_col = 'ALT',
                        other_allele_col = 'REF',
                        beta_col = 'beta',
                        se_col = 'LOG(OR)_SE',
                        eaf_col = 'eaf.exposure',
                        pval_col = 'P')
  
  mr_input <- harmonise_data(kfun_snp, mr_out, action = 1)
  mr_input <- mr_input %>% filter(ambiguous == FALSE)
  mr_output <- mr(mr_input, method_list = c('mr_ivw', 
                                            'mr_two_sample_ml',
                                            'mr_egger_regression',
                                            'mr_weighted_median',
                                            'mr_weighted_mode'))
  
  mr_i <- with(mr_input, MendelianRandomization::mr_ivw(mr_input(bx = beta.exposure,
                                                                 by = beta.outcome,
                                                                 bxse = se.exposure,
                                                                 byse = se.outcome)))
  print(mr_i)
  
  mr_e <- with(mr_input, mr_egger_regression(b_exp = beta.exposure,
                                             b_out = beta.outcome,
                                             se_exp = se.exposure,
                                             se_out = se.outcome))
  
  or = round(0.9^mr_output$b, 2)
  lb = with(mr_output, round(0.9^(b + 1.96*se), 2))
  ub = with(mr_output, round(0.9^(b - 1.96*se), 2))
  
  #or = round(exp(0.693*mr_output$b), 2)
  #lb = with(mr_output, round(exp(0.693*(b - 1.96*se)), 2))
  #ub = with(mr_output, round(exp(0.693*(b + 1.96*se)), 2))
  ci = paste(lb, ub, sep = '-')
  
  cochranq = round(mr_i$Heter.Stat, 1)[1]
  het = round(mr_i$Heter.Stat, 3)[2]
  cochranq_p = paste0(cochranq, ' (', 
                      ifelse(het == 0, '< 0.001', het), ')')
            
  mr_results <- rbind(mr_results,
                      data.table(ca = c(i, rep(NA, 4)),
                                 method = mr_output$method,
                                 or_95ci = paste0(or, ' (', ci, ')'),
                                 or, lb, ub,
                                 p = round(mr_output$pval, 3),
                                 plei = c(round(mr_e$pval_i, 2), rep(NA, 4)),
                                 cochranq = c(cochranq_p, rep(NA, 4))))
}
fwrite(mr_results, '~/ckd_ca/results/eGFR_mr_ukb.csv')
fwrite(mr_results, '~/ckd_ca/results/ma_mr_ukb.csv')

##MRPRESSO
mrpresso <- data.table()
for (i in ca) {
  tmp <- fread(paste0('~/ckd_ca/glm/plink2.', i, '.glm.logistic.hybrid'))
  mr_out <- tmp[ID %in% ma_snp$SNP][, beta := log(OR)]
  setnames(mr_out, 'ID', 'SNP')
  mr_out <- ma_snp[, .(SNP, eaf.exposure)][mr_out, on = 'SNP']
  mr_out <- format_data(mr_out,
                        type = 'outcome',
                        snp_col = 'SNP',
                        chr_col = '#CHROM',
                        effect_allele_col = 'ALT',
                        other_allele_col = 'REF',
                        beta_col = 'beta',
                        se_col = 'LOG(OR)_SE',
                        eaf_col = 'eaf.exposure',
                        pval_col = 'P')
  
  mr_input <- harmonise_data(ma_snp, mr_out, action = 1)
  mr_input <- mr_input %>% filter(ambiguous == FALSE)
  
  mr_p <- MRPRESSO::mr_presso(BetaOutcome = 'beta.outcome',
                              BetaExposure = 'beta.exposure',
                              SdOutcome = 'se.outcome',
                              SdExposure = 'se.exposure',
                              data = mr_input,
                              OUTLIERtest = T, 
                              DISTORTIONtest = T, 
                              NbDistribution = 5000, 
                              SignifThreshold = 0.05)
  
  or = with(mr_p$`Main MR results`, round(ifelse(is.na(`Causal Estimate`[2]),
                                                 exp(0.693*`Causal Estimate`)[1],
                                                 exp(0.693*`Causal Estimate`)[2]), 2))
  lb = with(mr_p$`Main MR results`, round(ifelse(is.na(Sd[2]), 
                                                 exp(0.693*(`Causal Estimate` - 1.96*Sd))[1], 
                                                 exp(0.693*(`Causal Estimate` - 1.96*Sd))[2]), 2))
  ub = with(mr_p$`Main MR results`, round(ifelse(is.na(Sd[2]), 
                                                 exp(0.693*(`Causal Estimate` + 1.96*Sd))[1], 
                                                 exp(0.693*(`Causal Estimate` + 1.96*Sd))[2]), 2))
  ci = paste(lb, ub, sep = '-')
  outlier = paste(mr_input$SNP[mr_p$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`], collapse = ', ')
  
  mrpresso <- rbind(mrpresso, data.table(ca = i,
                                         or_95ci = paste0(or, ' (', ci, ')'),
                                         or, lb, ub,
                                         outlier))
}
fwrite(mrpresso, '~/ckd_ca/results/eGFR_mrpresso_ukb.csv')
fwrite(mrpresso, '~/ckd_ca/results/ma_mrpresso_ukb.csv')

##multivariable mr of eGFR and microalbuminuria
library(MVMR)
mv_snp <- fread('~/ckd_ca/mv_ma_snp.csv')
mvmr_results <- data.table()
for (i in ca) {
  tmp <- fread(paste0('~/ckd_ca/glm/plink2.', i, '.glm.logistic.hybrid'))
  mv_out <- tmp[ID %in% mv_snp$SNP][, `:=` (beta = log(OR),
                                                 se = `LOG(OR)_SE`)]
  setnames(mv_out, 'ID', 'SNP')
  mv_out <- mv_snp[, .(SNP, eaf.exposure)][mv_out, on = 'SNP']
  mv_out <- format_data(mv_out,
                        type = 'outcome',
                        snp_col = 'SNP',
                        chr_col = '#CHROM',
                        effect_allele_col = 'ALT',
                        other_allele_col = 'REF',
                        beta_col = 'beta',
                        se_col = 'LOG(OR)_SE',
                        eaf_col = 'eaf.exposure',
                        pval_col = 'P')
  mv <- as.data.table(harmonise_data(mv_snp, mv_out, action = 1))
  mv_input <- format_mvmr(
    BXGs = mv[, .(beta_kfun,
                  beta_ma)],
    BYG = mv$beta.outcome,
    seBXGs = mv[, .(se_kfun,
                    se_ma)],
    seBYG = mv$se.outcome,
    RSID = mv$SNP
  )
  
  mv_output <- as.data.table(ivw_mvmr(mv_input))
  or = with(mv_output, c(round(0.9^Estimate[1], 2),
                         round(exp(0.693*Estimate[2]), 2)))
  ci = with(mv_output, paste(c(round(0.9^(Estimate[1] + 1.96*`Std. Error`[1]), 2),
                               round(exp(0.693*(Estimate[2] - 1.96*`Std. Error`[2])), 2)),
                             c(round(0.9^(Estimate[1] - 1.96*`Std. Error`[1]), 2),
                               round(exp(0.693*(Estimate[2] + 1.96*`Std. Error`[2])), 2)),
                             sep = '-'))
  
  #confitional F statistics
  con_fstat <- strength_mvmr(mv_input)
  
  mvmr_results <- rbind(mvmr_results,
                      data.table(ca = c(i, NA),
                                 exp = c('kfun', 'ma'),
                                 or_95ci = paste0(or, ' (', ci, ')'),
                                 p = round(mv_output$`Pr(>|t|)`, 3),
                                 con_fstat = t(con_fstat),
                                 nsnp = summary(as.factor(mv$exp))))
}
fwrite(mvmr_results, '~/ckd_ca/results/mv_ma_ukb.csv')



##test of difference in glm between plink and R
#get genotypes in hpc
#plink \
#--bfile ukb_imp \
#--extract kfun.snplist \
#--keep ukb.qc.id \
#--recode A \
#--recode-allele kfun_snp.a1 \
#--out kfun
kfun_ukb <- fread('~/ckd_ca/test/kfun.raw') %>%
  select('FID' |
           starts_with('rs'))

for (i in 1:213) {
  colnames(kfun_ukb)[i+1] <- unlist(
    strsplit(
      colnames(kfun_ukb)[i+1], '_'))[1]
}

kfun_beta <- data.table(snp = colnames(kfun_ukb)[-1])
kfun_se <- data.table(snp = colnames(kfun_ukb)[-1])
for (i in 3:17) {
  tmp <- kfun_ukb[ca_all[, c(1, 17:31)], on = 'FID'][, ca := ca_all[[i]]]
  for (j in 1:213) {
    snp <- tmp[[j+1]]
    mod <- glm(as.factor(ca) ~ snp + age + sex + rel3 + batch + pc1 + pc2 + pc3+ pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, tmp, family = binomial())
    kfun_beta$ca[j] <- summary(mod)$coef[2, 1]
    kfun_se$ca[j] <- summary(mod)$coef[2, 2]
  }
  setnames(kfun_beta, 'ca', colnames(ca_all)[i])
  setnames(kfun_se, 'ca', colnames(ca_all)[i])
}
kfun_snp <- fread('~/ckd_ca/test/kidney_function_snp.csv')
mr_out <- data.frame(SNP = kfun_beta$snp,
                         beta_ca = kfun_beta$ca,
                         se_ca = kfun_se$ca) 
mr_out <- kfun_snp[mr_out, on = 'SNP']
mr_out<- format_data(mr_out,
                          type = 'outcome',
                          snp_col = 'SNP',
                          effect_allele_col = 'effect_allele.exposure',
                          other_allele_col = 'other_allele.exposure',
                          eaf_col = 'eaf.exposure',
                          beta_col = 'beta_ca',
                          se_col = 'se_ca')
mr_input <- harmonise_data(exposure_dat = kfun_snp,
                           outcome_dat = mr_out)
