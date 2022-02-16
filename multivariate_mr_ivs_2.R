###mvmr incorporating kfun, bmi and smk init
mv_snp <- rbind(kidney_function_snp, bmi_snp, smkinit_snp)
mv_snp <- mv_snp %>% distinct(SNP, .keep_all = T) 

cl <- clump_data(
  format_data(subset(wuttke_2019_sum, RSID %in% mv_snp$SNP),
              type = 'exposure',
              snp_col = 'RSID',
              chr_col = 'Chr',
              effect_allele_col = 'Allele1',
              other_allele_col = 'Allele2',
              eaf_col = 'Freq1',
              beta_col = 'Effect',
              se_col = 'StdErr',
              pval_col = 'P-value'),
  clump_r2 = 0.01, pop = "EUR"
) 

p <- subset(
  mv_snp,
  ((effect_allele.exposure == 'A' & other_allele.exposure == 'T') | 
     (effect_allele.exposure == 'T' & other_allele.exposure == 'A') |
     (effect_allele.exposure == 'G' & other_allele.exposure == 'C') |
     (effect_allele.exposure == 'C' & other_allele.exposure == 'G')) &
    (eaf.exposure > 0.42 & eaf.exposure < 0.58),
  select = c(
    'SNP',
    'effect_allele.exposure',
    'other_allele.exposure',
    'eaf.exposure'
  )
)
mv_snp <- mv_snp %>% filter(!SNP %in% p$SNP & SNP %in% cl$SNP)

##kfun
kfun_mv <- format_data(subset(wuttke_2019_sum, RSID %in% mv_snp$SNP),
                       type = 'outcome',
                       snp_col = 'RSID',
                       chr_col = 'Chr',
                       effect_allele_col = 'Allele1',
                       other_allele_col = 'Allele2',
                       eaf_col = 'Freq1',
                       beta_col = 'Effect',
                       se_col = 'StdErr',
                       pval_col = 'P-value')

tmp <- harmonise_data(mv_snp, kfun_mv, action = 1) %>% 
  select(
    'SNP' | 
      ends_with('exposure') | 
      ends_with('stat') |
      starts_with('beta') |
      starts_with('se') |
      starts_with('pval.') |
      'exp'
  )
colnames(tmp)[13:15] <- c('beta_kfun', 'se_kfun', 'pval_kfun')

##bmi
bmi_mv <- extract_outcome_data(snps = tmp$SNP, 
                               outcomes = 'ieu-a-835',
                               proxies = FALSE,
                               access_token = NULL)

tmp <- harmonise_data(tmp, bmi_mv, action = 1) %>% 
  select(
    'SNP' | 
      ends_with('exposure') | 
      ends_with('stat') |
      ends_with('_kfun') |
      starts_with('beta.') |
      starts_with('se.') |
      starts_with('pval.') |
      'exp'
  )
colnames(tmp)[16:18] <- c('beta_bmi', 'se_bmi', 'pval_bmi')

#smk
liu_2019_sum <- fread('~/Desktop/ckd cancer/summary data/smk_init.gz')

smkinit_mv <- format_data(
  subset(liu_2019_sum, RSID %in% tmp$SNP),
  type = 'outcome',
  snp_col = 'RSID',
  effect_allele_col = 'ALT',
  other_allele_col = 'REF',
  beta_col = 'BETA',
  se_col = 'SE',
  pval_col = 'PVALUE',
)

tmp <- harmonise_data(tmp, smkinit_mv, action = 1) %>% 
  select(
    'SNP' | 
      ends_with('exposure') | 
      ends_with('stat') |
      ends_with('_kfun') |
      ends_with('_bmi') |
      starts_with('beta.') |
      starts_with('se.') |
      starts_with('pval.') |
      'exp'
  )
colnames(tmp)[19:21] <- c('beta_smkinit', 'se_smkinit', 'pval_smkinit')

mv_snp <- rbind(
  subset(tmp, exp == 's') %>%
    mutate(
      beta_smkinit = beta.exposure,
      se_smkinit = se.exposure,
      pval_smkinit = pval.exposure),
  subset(tmp, exp != 's')
)

fwrite(mv_snp, '~/desktop/ckd cancer/ivs/mvmr/mv_kbs_snp.csv')
