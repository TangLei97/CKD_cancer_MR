library(TwoSampleMR)
library(MendelianRandomization)
mr_results <- data.table()
  for (i in out_dat) {
    ##in mrcieu
    #mr_out <- extract_outcome_data(snps = ma_snp$SNP, 
    #                               outcomes = i, 
    #                               access_token = NULL)
    
    ##rcc
    mr_out <- format_data(rcc_m[rsid %in% ma_snp$SNP], type = 'outcome',
                          snp_col = 'rsid',
                          chr_col = 'chr',
                          effect_allele_col = 'ea',
                          other_allele_col = 'nea',
                          eaf_col = 'eaf',
                          beta_col = 'beta',
                          se_col = 'se',
                          pval_col = 'p')
    
    #effect of a snp on an outcome and exposure must be harmonised to be relative to the same allele
    mr_input <- harmonise_data(ma_snp, mr_out, action = 1)
    #removing palindromic SNPs with intermediate allele frequency
    mr_input <- mr_input %>% filter(ambiguous == FALSE)
    
    ##obtaining estimates from TwoSampleMR
    mr_output <- mr(mr_input, method_list = c('mr_ivw',
                                              'mr_two_sample_ml',
                                              'mr_egger_regression',
                                              'mr_weighted_median',
                                              'mr_weighted_mode'))
    
    ##obtaining pleiotropy
    mr_e <- with(mr_input, mr_egger_regression(b_exp = beta.exposure,
                                               b_out = beta.outcome,
                                               se_exp = se.exposure,
                                               se_out = se.outcome))
    
    ##obtaining i2/(cochran's 1 statistic) and associated p-value from MendelianRandomization::
    mr_i <- with(mr_input, MendelianRandomization::mr_ivw(mr_input(bx = beta.exposure,
                                                                   by = beta.outcome,
                                                                   bxse = se.exposure,
                                                                   byse = se.outcome)))
    print(mr_i)
    
    ##calculate F-statistics
    tstat <- with(mr_input, beta.exposure/se.exposure) 
    fstat <- mean((tstat)^2)

    ##scale results (per 10% decrement in eGFR/per doubling odds of ma)
    #or = round(0.9^mr_output$b, 2)
    #lb = with(mr_output, round(0.9^(b + 1.96*se), 2))
    #ub = with(mr_output, round(0.9^(b - 1.96*se), 2))
    
    or = round(exp(0.693*mr_output$b), 2)
    lb = with(mr_output, round(exp(0.693*(b - 1.96*se)), 2))
    ub = with(mr_output, round(exp(0.693*(b + 1.96*se)), 2))
    ci = paste(lb, ub, sep = '-')
   
    cochranq = round(mr_i$Heter.Stat, 1)[1]
    het = round(mr_i$Heter.Stat, 3)[2]
    cochranq_p = paste0(cochranq, ' (',
                        ifelse(het == 0, '< 0.001', het), ')')
    
    mr_results <- rbind(mr_results,
                        data.table(ca = c('rccm', rep(NA, 4)),
                                   method = mr_output$method,
                                   or_95ci = paste0(or, ' (', ci, ')'),
                                   or, lb, ub,
                                   p = round(mr_output$pval, 3),
                                   plei = c(round(mr_e$pval_i, 2), rep(NA, 4)),
                                   cochranq = c(cochranq_p, rep(NA, 4)),
                                   nsnp = c(mr_output$nsnp[1], rep(NA, 4)),
                                   f = c(round(fstat, 1), rep(NA, 4))))
  }
fwrite(mr_results, '~/desktop/ckd cancer/results/eGFR_mr.csv')
fwrite(mr_results, '~/desktop/ckd cancer/results/ma_mr.csv')

##obtaining outlier(s) & distortion test from MRPRESSO
mrpresso <- data.table()
for (i in out_dat) {
  mr_out <- extract_outcome_data(snps = ma_snp$SNP, 
                                 outcomes = i, 
                                 access_token = NULL)
  
  #mr_out <- format_data(rcc_m[rsid %in% ma_snp$SNP], type = 'outcome',
  #                      snp_col = 'rsid',
  #                      chr_col = 'chr',
  #                      effect_allele_col = 'ea',
  #                      other_allele_col = 'nea',
  #                      eaf_col = 'eaf',
  #                      beta_col = 'beta',
  #                      se_col = 'se',
  #                      pval_col = 'p')
  
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
fwrite(mrpresso, '~/desktop/ckd cancer/results/eGFR_mrpresso.csv')
fwrite(mrpresso, '~/desktop/ckd cancer/results/ma_mrpresso.csv')

##mvmr
library(MVMR)
mv_snp <- fread('~/desktop/ckd cancer/ivs/mvmr/mv_kfun_snp.csv')
mvmr_results <- data.table()
for (i in out_dat) {
  #mv_out <- extract_outcome_data(snps = mv_snp$SNP, 
  #                                      outcomes = i, 
  #                                      proxies = FALSE,
  #                                      access_token = NULL)
  
  mv_out <- format_data(rcc_f[rsid %in% mv_snp$SNP], type = 'outcome',
                        snp_col = 'rsid',
                        chr_col = 'chr',
                        effect_allele_col = 'ea',
                        other_allele_col = 'nea',
                        eaf_col = 'eaf',
                        beta_col = 'beta',
                        se_col = 'se',
                        pval_col = 'p')
  
  mv <- as.data.table(harmonise_data(mv_snp, mv_out, action = 1))
  
  ###read in MVMR::
  mv_input <- format_mvmr(BXGs = mv[, .(beta_kfun,
                                        beta_ma)],
                          BYG = mv$beta.outcome,
                          seBXGs = mv[, .(se_kfun,
                                          se_ma)],
                          seBYG = mv$se.outcome,
                          RSID = mv$SNP)
  
  ##conditional F statistics (ref: https://doi.org/10.1101/2020.04.02.021980, "Testing and Correcting for Weak and Pleiotropic Instruments in Two-Sample Multivariable Mendelian Randomisation")
  con_fstat <- strength_mvmr(mv_input)
  
  ##‘calculating the weak instrument test’
  #delta1 <- summary(lm(mvmrdata$x1.beta ~ -1 + mvmrdata$x2.beta))$coefficients[1,1]
  #v1 <- mvmrdata$x1.se^2 + delta1^2*(mvmrdata$x2.se^2)
  #F.x1 <- sum((1/v1)*(mvmrdata$x1.beta - (delta1*mvmrdata$x2.beta))^2)/(70-1)
  
  mv_output <- as.data.table(ivw_mvmr(mv_input))
  or = with(mv_output, c(round(0.9^Estimate[1], 2),
                         round(exp(0.693*Estimate[2]), 2)))
  ci = with(mv_output, paste(c(round(0.9^(Estimate[1] + 1.96*`Std. Error`[1]), 2),
                               round(exp(0.693*(Estimate[2] - 1.96*`Std. Error`[2])), 2)),
                             c(round(0.9^(Estimate[1] - 1.96*`Std. Error`[1]), 2),
                               round(exp(0.693*(Estimate[2] + 1.96*`Std. Error`[2])), 2)),
                             sep = '-'))
  
  ##or of cancer per 10% decrease in eGFR/sd (4.8) increase in bmi/doubling in odds of smkinit
  mvmr_results <- rbind(mvmr_results,
                        data.table(ca = c('rccf', NA),
                                   exp = c('kfun', 'ma'),
                                   or_95ci = paste0(or, ' (', ci, ')'),
                                   p = round(mv_output$`Pr(>|t|)`, 3),
                                   con_fstat = t(con_fstat),
                                   nsnp = summary(as.factor(mv$exp))))
}

fwrite(mvmr_results, '~/desktop/ckd cancer/results/mv_kfun.csv')
