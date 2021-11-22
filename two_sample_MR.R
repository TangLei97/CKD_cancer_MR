exp <- c('kidney_function', 'ma')
mr_results <- data.table()
for (i in exp) {
  mr_exp <- fread(paste0(i,'_snp.csv'))
  for (j in out_dat) {
    mr_out <- extract_outcome_data(snps = mr_exp$SNP, 
                                   outcomes = j, 
                                   access_token = NULL)
    
    mr_input <- harmonise_data(exposure_dat = mr_exp, 
                               outcome_dat = mr_out)#default action = 2
    
    ##obtain estimates from TwoSampleMR
    mr_output <- mr(mr_input,
                    method_list = c('mr_ivw',
                                    'mr_two_sample_ml',
                                    'mr_egger_regression',
                                    'mr_weighted_median',
                                    'mr_weighted_mode'))
    
    ##obtain pleiotropy
    mr_e <- mr_egger_regression(b_exp = mr_input$beta.exposure,
                                b_out = mr_input$beta.outcome,
                                se_exp = mr_input$se.exposure,
                                se_out = mr_input$se.outcome)
    
    ##obtain I2/(Cochran's Q statistic) and associated p-value from MendelianRandomization::
    mr_input <- mr_input[mr_keep == 'TRUE', ]#remove palindromic SNPs with intermediate allele frequency
    
    mr_i <- MendelianRandomization::mr_ivw(mr_input(bx = mr_input$beta.exposure,
                                                    by = mr_input$beta.outcome,
                                                    bxse = mr_input$se.exposure,
                                                    byse = mr_input$se.outcome))
    print(mr_i)
    
    ##calculate F-statistics
    tstat <- with(mr_input, beta.exposure/se.exposure) 
    fstat <- mean((tstat)^2)
    
    ##obtain outlier & distortion test from MRPRESSO
    mr_p <- MRPRESSO::mr_presso(BetaOutcome = 'beta.outcome',
                                BetaExposure = 'beta.exposure',
                                SdOutcome = 'se.outcome',
                                SdExposure = 'se.exposure',
                                data = mr_input,
                                OUTLIERtest = T, 
                                DISTORTIONtest = T, 
                                NbDistribution = 5000, 
                                SignifThreshold = 0.05)
    
    ##scale results (per 10% decrement in eGFR/per doubling odds of ma)
    if (i == 'kidney_function') {
      or = c(round(0.9^mr_output$b, 2),
             round(0.9^mr_p$`main MR results`$`Causal Estimate`,
                   2)[2])
      ci = c(paste(round(0.9^(mr_output$b + 1.96*mr_output$se), 2),
                   round(0.9^(mr_output$b - 1.96*mr_output$se), 2),
                   sep = '-'),
             paste(round(0.9^(mr_p$`main MR results`$`Causal Estimate` -
                                1.96*mr_p$`main MR results`$Sd), 2)[2],
                   round(0.9^(mr_p$`main MR results`$`Causal Estimate` -
                                1.96*mr_p$`main MR results`$Sd), 2)[2],
                   sep = '-'))
    } else {
      or = c(round(exp(0.693*mr_output$b), 2),
             round(exp(0.693*mr_p$`Main MR results`$`Causal Estimate`),
                   2)[2])
      ci = c(paste(round(exp(0.693*(mr_output$b - 1.96*mr_output$se)), 2),
                   round(exp(0.693*(mr_output$b + 1.96*mr_output$se)), 2),
                   sep = '-'),
             paste(round(exp(0.693*(mr_p$`Main MR results`$`Causal Estimate` -
                                      1.96*mr_p$`Main MR results`$Sd)), 2)[2],
                   round(exp(0.693*(mr_p$`Main MR results`$`Causal Estimate` +
                                      1.96*mr_p$`Main MR results`$Sd)), 2)[2],
                   sep = '-'))
    }
      mr_results <- rbind(mr_results,
                          data.table(method = c(mr_output$method, 'MR PRESSO'),
                                     or,
                                     ci,
                                     p_val = c(round(mr_output$pval, 3), NA),
                                     p_plei = c(round(mr_e$pval_i, 2), rep(NA, 5)),
                                     outlier = c(length(mr_p$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`),
                                                 rep(NA, 5)),
                                     exp = c(i, rep(NA, 5)),
                                     ca = c(j, rep(NA, 5)),
                                     fstat = c(fstat, rep(NA, 5)),
                                     nsnp = c(unique(mr_output$nsnp), rep(NA, 5)),
                                     cochranq_stat = c(round(mr_i$Heter.Stat, 1)[1], rep(NA, 5)),
                                     p_het = c(round(mr_i$Heter.Stat, 3)[2], rep(NA, 5))))
  }
}

