bi_ma_results <- data.table()
bi_kfun_results <- data.table()
###rcc
##Supplementary table S2-3 "Following MAF and LD based QC, a total of 17 sex-specific SNPs (6 significant and 11 suggestive) were selected for follow-up."
##7/15 associations were male-specific, whereas, 8/15 SNPs were female-specific
library(readxl)
library(purrr)
path = '~/Desktop/ckd cancer/ivs/ca/rcc.xlsx'
walk(
  .x = excel_sheets(path),
  .f = ~ assign(., read_excel(path, sheet = .), envir = .GlobalEnv)
)                     

f <- f %>%
  mutate(eaf = (eaf_cas + eaf_con)/2,
         beta = log(or),
         se = (log(ub) - log(or))/1.96)
rcc_snp <- format_data(f, type = 'exposure',
                       snp_col = 'snp',
                       chr_col = 'chr',
                       effect_allele_col = 'ea',
                       other_allele_col = 'nea',
                       eaf_col = 'eaf',
                       beta_col = 'beta',
                       se_col = 'se',
                       pval_col = 'p')

rcc_snp_f <- rcc_snp %>% mutate(ca = 'rcc_f')

###crc, leu, cerv, nhl
colo <- fread('~/ckd_ca/gwas/gwas.colo.glm.logistic.hybrid') %>%
  filter(P < 5e-08)
colo[, beta := log(OR)]
colo_snp <- clump_data(format_data(colo, type = 'exposure',
                                  snp_col = 'ID',
                                  chr_col = '#CHROM',
                                  effect_allele_col = 'ALT',
                                  other_allele_col = 'REF',
                                  beta_col = 'beta',
                                  se_col = 'LOG(OR)_SE',
                                  pval_col = 'P'), clump_r2 = 0.01, pop = "EUR")

nhl_snp <- nhl_snp %>% mutate(ca = 'nhl')

mr_out <- format_data(teumer_2019_sum[RSID %in% colo_snp$SNP],
                      type = 'outcome',
                      snp_col = 'RSID',
                      chr_col = 'Chr',
                      effect_allele_col = 'Allele1',
                      other_allele_col = 'Allele2',
                      eaf_col = 'Freq1',
                      beta_col = 'Effect',
                      se_col = 'StdErr',
                      pval_col = 'P-value')

mr_input <- harmonise_data(colo_snp, mr_out, action = 1)
mr_output <- mr(mr_input, method_list = c('mr_ivw',
                                          'mr_two_sample_ml',
                                          'mr_egger_regression',
                                          'mr_weighted_median',
                                          'mr_weighted_mode'))
mr_output <- mr(mr_input)

mr_e <- with(mr_input, mr_egger_regression(b_exp = beta.exposure,
                                           b_out = beta.outcome,
                                           se_exp = se.exposure,
                                           se_out = se.outcome))

mr_i <- with(mr_input, MendelianRandomization::mr_ivw(mr_input(bx = beta.exposure,
                                                               by = beta.outcome,
                                                               bxse = se.exposure,
                                                               byse = se.outcome)))
print(mr_i)

mr_p <- MRPRESSO::mr_presso(BetaOutcome = 'beta.outcome',
                            BetaExposure = 'beta.exposure',
                            SdOutcome = 'se.outcome',
                            SdExposure = 'se.exposure',
                            data = mr_input,
                            OUTLIERtest = T, 
                            DISTORTIONtest = T, 
                            NbDistribution = 5000, 
                            SignifThreshold = 0.05)

tstat <- with(mr_input, beta.exposure/se.exposure) 
fstat <- mean((tstat)^2)

or = round(exp(0.693*mr_output$b), 2)
lb = with(mr_output, round(exp(0.693*(b - 1.96*se)), 2))
ub = with(mr_output, round(exp(0.693*(b + 1.96*se)), 2))
ci = paste(lb, ub, sep = '-')

beta = with(mr_p$`Main MR results`, ifelse(is.na(`Causal Estimate`[2]), `Causal Estimate`[1], `Causal Estimate`[2]))
p = with(mr_p$`Main MR results`, ifelse(is.na(`Causal Estimate`[2]), `P-value`[1], `P-value`[2]))
outlier = paste(mr_input$SNP[mr_p$`MR-PRESSO results`$`Distortion Test`$`Outliers Indices`], collapse = ', ')
se = sqrt(((beta)^2)/qchisq(p, 1, lower.tail = F))
or = round(exp(0.693*beta), 2)
lb = round(exp(0.693*(beta - 1.96*se)), 2)
ub = round(exp(0.693*(beta + 1.96*se)), 2)
ci = paste(lb, ub, sep = '-')
p = round(p, 3)

cochranq = round(mr_i$Heter.Stat, 1)[1]
het = round(mr_i$Heter.Stat, 3)[2]
cochranq_p = paste0(cochranq, ' (',
                    ifelse(het == 0, '< 0.001', het), ')')

bi_ma_results <- rbind(bi_ma_results,
                    data.table(ca = c('crc', rep(NA, 5)),
                               method = c(mr_output$method, 'MRPRESSO'),
                               or_95ci = c(paste0(or, ' (', ci, ')'), paste0(or_p, ' (', ci_p, ')')),
                               p = c(round(mr_output$pval, 3), p),
                               plei = c(round(mr_e$pval_i, 2), rep(NA, 5)),
                               cochranq = c(cochranq_p, rep(NA, 5)),
                               outlier = c(outlier, rep(NA, 5)),
                               nsnp = c(mr_output$nsnp[1], rep(NA, 5)),
                               f = c(round(fstat, 1), rep(NA, 5))))

bi_ma_results <- rbind(bi_ma_results,
                       data.table(ca = 'nhl',
                                  method = 'Wald ratio',
                                  or_95ci = paste0(or, ' (', ci, ')'),
                                  p = round(mr_output$pval, 3),
                                  plei = NA,
                                  cochranq = NA,
                                  outlier = NA,
                                  nsnp = 1,
                                  f = round(fstat, 1)))
fwrite(rbind(rcc_snp_f, rcc_snp_m, leu_snp, cerv_snp, colo_snp, nhl_snp),
       '~/desktop/ckd cancer/ivs/ca/ca_snp.csv')
fwrite(bi_kfun_results, '~/desktop/ckd cancer/results/bi_kfun.csv')
fwrite(bi_ma_results, '~/desktop/ckd cancer/results/bi_ma.csv')
                       