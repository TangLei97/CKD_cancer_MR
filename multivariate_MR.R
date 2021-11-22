###SNPs for BMI
##GWAS(Locke et al 2015 GIANT) identified xxx SNPs significant in European ancenstry: 
#‘Our primary meta-analysis of [European-descent individuals from GWAS and Metabochip studies] (n = 322,154) identified 77* loci reaching genome-wide significance (GWS) and separated by at least 500 kilo-bases (kb)’
#Supplementary Table 4 Significant loci for BMI at P < 1 × 10-5* for European sex-combined analysis
locke_2015_snp <- readxl::read_xlsx('~/Desktop/ckd cancer/Locke2015-supplement-2.xlsx', sheet = 4) %>% 
  filter(`European sex-combined GWAS I+II + Metabochip P value` < 5e-08)#77 SNPs

##format & further clump data using a window of 10,000 kb and maximal linkage disequilibrium of r2 < 0.01 as threshold
bmi_snp <- clump_data(
  format_data(
    locke_2015_snp,
    type = 'exposure',
    snp_col = 'SNP',
    chr_col = 'Chr.',
    effect_allele_col = 'Effect alleles',
    other_allele_col = 'Other alleles',
    eaf_col = 'Effect allele frequency',
    beta_col = 'BETA',
    se_col = 'SE',
    pval_col = 'European sex-combined GWAS I+II + Metabochip P value',
  ),
  clump_r2 = 0.01,pop = "EUR"
)

###SNPs for smoking
##GWAS(Liu et al 2019) identified xxx SNPs significant in European ancenstry: 
#‘Smoking initiation phenotypes included a binary phenotype indicating whether an individual had ever smoked regularly (SmkInit, N=1,232,091; 378 associated variants)’
#Supplementary Table 4 Genome-Wide Significant Conditionally Independent Association Results for Smoking Initiation
liu_2019_snp <- readxl::read_xlsx('~/Desktop/ckd cancer/Liu2019-supplement-2.xlsx',sheet = 4) %>% 
  filter(`Pvalue` < 5e-08)#378 SNPs

##format & further clump data using a window of 10,000 kb and maximal linkage disequilibrium of r2 < 0.01 as threshold
smkinit_snp <- clump_data(
  format_data(
    liu_2019_snp,
    type = 'exposure',
    snp_col = 'rsID',
    chr_col = 'Chr',
    effect_allele_col = 'Alternate Allele',
    other_allele_col = 'Reference Allele',
    eaf_col = 'Alternate Allele Frequency',
    beta_col = 'Beta',
    se_col = 'SE',
    pval_col = 'Pvalue',
  ),
  clump_r2 = 0.01,pop = "EUR"
)

kidney_function_snp[, exp := 'k']
bmi_snp[, exp := 'b']
smkinit_snp[, exp := 's']

multi_mr_snp <- rbind(kidney_function_snp,
                      bmi_snp,
                      smkinit_snp)

#check for duplicated SNPs
length(unique(multi_mr_snp$SNP))

##Further identify independent SNPs (ref: doi: 10.7554/eLife.43990)
#Combined results clumped by selecting the SNP within a region most strongly associated with kidney function in the Wuttke et al GWAS (i.e. lowest p-value). 
multi_mr_snp_cl <- clump_data(
  format_data(
    subset(wuttke_2019_sum,
           RSID %in% multi_mr_snp$SNP,
           select = c(
             'RSID',
             'Chr',
             'Allele1',
             'Allele2',
             'Freq1',
             'Effect',
             'StdErr',
             'P-value'
           )),
    type = 'exposure',
    snp_col = 'RSID',
    chr_col = 'Chr',
    effect_allele_col = 'Allele1',
    other_allele_col = 'Allele2',
    eaf_col = 'Freq1',
    beta_col = 'Effect',
    se_col = 'StdErr',
    pval_col = 'P-value',
  ),
  clump_r2 = 0.01,pop = "EUR"
)

##Remove palindromic SNPs (i.e., SNPs with A/T or C/G) with intermediate allele frequency (0.42-0.58, inferred from TwoSampleMR::harmonise_data(action = 2))
multi_mr_snp_pal <- subset(
  multi_mr_snp,
  ((effect_allele.exposure == 'A' & other_allele.exposure == 'T') | 
     (effect_allele.exposure == 'T' & other_allele.exposure == 'A') |
     (effect_allele.exposure == 'G' & other_allele.exposure == 'C') |
     (effect_allele.exposure == 'C' & other_allele.exposure == 'G')) &
    (eaf.exposure >= 0.42 & eaf.exposure <= 0.58),
  select = c(
    'SNP',
    'effect_allele.exposure',
    'other_allele.exposure',
    'eaf.exposure'
  )
)

###Extract beta, se of each exposure from GWAS (effect of a SNP on an outcome and exposure must be harmonised to be relative to the same allele)
multi_mr_dat <- multi_mr_snp %>% 
  filter((SNP %in% multi_mr_snp_cl$SNP) &
           (!SNP %in% multi_mr_snp_pal$SNP))

##kidney function
kidney_function_mv_dat <- format_data(
  subset(wuttke_2019_sum,
         RSID %in% multi_mr_dat$SNP,
         select = c(
           'RSID',
           'Chr',
           'Allele1',
           'Allele2',
           'Freq1',
           'Effect',
           'StdErr',
           'P.value'
         )),
  type = 'outcome',
  snp_col = 'RSID',
  chr_col = 'Chr',
  effect_allele_col = 'Allele1',
  other_allele_col = 'Allele2',
  eaf_col = 'Freq1',
  beta_col = 'Effect',
  se_col = 'StdErr',
  pval_col = 'P.value',
)

multi_mr_dat_tmp <- harmonise_data(exposure_dat = multi_mr_dat,
                                   outcome_dat = kidney_function_mv_dat,
                                   action = 1) %>% 
  select(
    'SNP' | 
      ends_with('exposure') | 
      ends_with('stat') |
      starts_with('beta') |
      starts_with('se') |
      starts_with('pval.') |
      'trait'
  )
colnames(multi_mr_dat_tmp)[13:15] <- c('beta_kidney_function',
                                       'se_kidney_function',
                                       'pval_kidney_function')

##BMI
bmi_mv_dat <- extract_outcome_data(snps = multi_mr_dat_tmp$SNP, 
                                   outcomes = 'ieu-a-835',
                                   proxies = FALSE,
                                   access_token = NULL)

multi_mr_dat_tmp <- harmonise_data(exposure_dat = multi_mr_dat_tmp,
                                   outcome_dat = bmi_mv_dat,
                                   action = 1) %>% 
  select(
    'SNP' | 
      ends_with('exposure') | 
      ends_with('stat') |
      ends_with('_kidney_function') |
      starts_with('beta.') |
      starts_with('se.') |
      starts_with('pval.') |
      'trait'
  )
colnames(multi_mr_dat_tmp)[16:18] <- c('beta_bmi',
                                       'se_bmi',
                                       'pval_bmi')

#smk
liu_2019_sum <- fread('~/Desktop/ckd cancer/summary data/smk_init.gz')

smkinit_mv_dat <- format_data(
  subset(liu_2019_sum,
         RSID %in% multi_mr_dat_tmp$SNP),
  type = 'outcome',
  snp_col = 'RSID',
  effect_allele_col = 'ALT',
  other_allele_col = 'REF',
  beta_col = 'BETA',
  se_col = 'SE',
  pval_col = 'PVALUE',
)

multi_mr_dat_tmp <- harmonise_data(exposure_dat = multi_mr_dat_tmp,
                                   outcome_dat = smkinit_mv_dat,
                                   action = 1) %>% 
  select(
    'SNP' | 
      ends_with('exposure') | 
      ends_with('stat') |
      ends_with('_kidney_function') |
      ends_with('_bmi') |
      starts_with('beta.') |
      starts_with('se.') |
      starts_with('pval.') |
      'trait'
  )
colnames(multi_mr_dat_tmp)[19:21] <- c('beta_smkinit',
                                       'se_smkinit',
                                       'pval_smkinit')

multi_mr_exp_dat <- rbind(
  subset(multi_mr_dat_tmp,trait == 's') %>%
    mutate(
      beta_smkinit = beta.exposure,
      se_smkinit = se.exposure,
      pval_smkinit = pval.exposure),
  subset(multi_mr_dat_tmp,trait != 's')
)

for (i in out_dat) {
  multi_mr_out_dat <- extract_outcome_data(snps = multi_mr_exp_dat$SNP, 
                                           outcomes = i, 
                                           proxies = FALSE,
                                           access_token = NULL)
  
  multi_mr_comb_dat <- harmonise_data(exposure_dat = multi_mr_exp_dat,
                                      outcome_dat = multi_mr_out_dat,
                                      action = 1)
  
  ###read in MVMR::
  #library(MVMR)
  mvmr_input <- format_mvmr(
    BXGs = multi_mr_comb_dat[, c('beta_kidney_function',
                                'beta_bmi',
                                'beta_smkinit')],
    BYG = multi_mr_comb_dat$beta.outcome,
    seBXGs = multi_mr_comb_dat[, c('se_kidney_function',
                                  'se_bmi',
                                  'se_smkinit')],
    seBYG = multi_mr_comb_dat$se.outcome,
    RSID = multi_mr_comb_dat$SNP
  )
  
  ##conditional F statistics (ref: https://doi.org/10.1101/2020.04.02.021980, "Testing and Correcting for Weak and Pleiotropic Instruments in Two-Sample Multivariable Mendelian Randomisation")
  con_fstat <- strength_mvmr(mvmr_input)#26/14/22
  colnames(con_fstat) <- c('kidney_function',
                           'bmi',
                           'smkinit')
  ##‘calculating the weak instrument test’
  #delta1 <- summary(lm(mvmrdata$x1.beta ~ -1 + mvmrdata$x2.beta))$coefficients[1,1]
  #v1 <- mvmrdata$x1.se^2 + delta1^2*(mvmrdata$x2.se^2)
  #F.x1 <- sum((1/v1)*(mvmrdata$x1.beta - (delta1*mvmrdata$x2.beta))^2)/(70-1)
  
  ##OR of cancer per 10% decrease in eGFR/SD (4.8) increase in BMI/doubling in odds of smkinit
  mvmr_results <- rbind(
    mvmr_results,
    mutate(
      as.data.table(ivw_mvmr(mvmr_input)),
      or = c(round(0.9^Estimate[1], 2),
             round(exp(Estimate[2]), 2),
             round(exp(0.693*Estimate[3]), 2)),
      ci = paste(
        c(round(0.9^(Estimate[1] + 1.96*`Std. Error`[1]), 2),
          round(exp(Estimate[2] - 1.96*`Std. Error`[2]), 2),
          round(exp(0.693*(Estimate[3] - 1.96*`Std. Error`[3])), 2)),
        c(round(0.9^(Estimate[1] - 1.96*`Std. Error`[1]), 2),
          round(exp(Estimate[2] + 1.96*`Std. Error`[2]), 2),
          round(exp(0.693*(Estimate[3] + 1.96*`Std. Error`[3])), 2)),
        sep = '-'
      ),
      p_val = round(`Pr(>|t|)`, 3),
      ca = c(i, NA, NA),
      exp = c('kidney function',
              'BMI',
              'smoking initiation'),
      con_fstat = t(con_fstat),
      nsnp = summary(as.factor(multi_mr_comb_dat$exp))
    )
  )
}
