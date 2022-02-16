dir <- '~/desktop/ckd cancer/finngen/'
dir <- '~/desktop/ckd cancer/yanglab/'
file <- dir(dir)

mr_results <- data.table()
for (f in file[12]) {
  summary <- fread(paste0(dir, f))
  phe <- unlist(strsplit(f, '_'))[1]
  
  for (i in c('kidney_function', 'ma')) {
    mr_exp <- fread(paste0('~/desktop/ckd cancer/', i, '_snp.csv'))
    mr_out <- format_data(summary[SNP %in% mr_exp$SNP], type = 'outcome',
                          snp_col = 'SNP',
                          chr_col = 'CHR',
                          effect_allele_col = 'A1',
                          other_allele_col = 'A2',
                          eaf_col = 'AF1',
                          beta_col = 'BETA',
                          se_col = 'SE',
                          pval_col = 'P')
    
    mr_input <- harmonise_data(mr_exp, mr_out, action = 1) %>% filter(ambiguous == FALSE)
    
    mr_output <- mr(mr_input, method_list = c('mr_ivw',
                                              'mr_two_sample_ml',
                                              'mr_egger_regression',
                                              'mr_weighted_median',
                                              'mr_weighted_mode'))
    mr_results <- rbind(mr_results, mr_output %>% 
                          select(method, nsnp, b, se, pval) %>%
                          mutate(phecode = phe, exp = i))
  }
}

kfun_mr <- mr_results[exp == 'kidney_function']
kfun_mr[, `:=` (or = round(0.9^b, 2),
                ub = round(0.9^(b + 1.96*se), 2),
                lb = round(0.9^(b - 1.96*se), 2),
                pval = round(pval, 3))]
kfun_mr[, or_95ci := paste0(or, ' (', ub, '-', lb, ')')]

ma_mr <- mr_results[exp == 'ma']
ma_mr[, `:=` (or = round(exp(0.693*b), 2),
              ub = round(exp(0.693*(b - 1.96*se)), 2),
              lb = round(exp(0.693*(b + 1.96*se)), 2),
              pval = round(pval, 3))]
ma_mr[, or_95ci := paste0(or, ' (', ub, '-', lb, ')')]

xlsx <- rbind(kfun_mr, ma_mr) %>% 
  group_split(exp) %>%
  set_names(levels(as.factor(rbind(kfun_mr, ma_mr)$exp)))
write_xlsx(xlsx, "~/desktop/ckd cancer/yanglab.xlsx")
