###snps for BMI
##gwas (Locke et al 2015 GIANT) identified xxx snps significant in European ancenstry: 
#‘Our primary meta-analysis of [European-descent individuals from GWAS and Metabochip studies] (n = 322,154) identified 77* loci reaching genome-wide significance (GWS) and separated by at least 500 kilo-bases (kb)’
#Supplementary Table 4 Significant loci for BMI at P < 1 × 10-5* for European sex-combined analysis
locke_2015_snp <- readxl::read_xlsx('~/Desktop/ckd cancer/ivs/mvmr/Locke2015-supplement-2.xlsx', sheet = 4) %>% 
  filter(`European sex-combined GWAS I+II + Metabochip P value` < 5e-08) #77 SNPs

##formatting & further clumping data using a window of 10,000 kb and maximal linkage disequilibrium of r2 < 0.01 as threshold
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
  clump_r2 = 0.01, pop = "EUR"
)

###snps for smoking
##gwas (Liu et al 2019) identified xxx snps significant in European ancenstry: 
#‘Smoking initiation phenotypes included a binary phenotype indicating whether an individual had ever smoked regularly (SmkInit, N=1,232,091; 378 associated variants)’
#Supplementary Table 4 Genome-Wide Significant Conditionally Independent Association Results for Smoking Initiation
liu_2019_snp <- readxl::read_xlsx('~/Desktop/ckd cancer/ivs/mvmr/Liu2019-supplement-2.xlsx',sheet = 4) %>% 
  filter(`Pvalue` < 5e-08) #378 SNPs

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

kidney_function_snp <- kidney_function_snp %>%
  mutate(exp = 'k')
ma_snp <- ma_snp %>%
  mutate(exp = 'a')
bmi_snp <- bmi_snp %>%
  mutate(exp = 'b')
smkinit_snp <- smkinit_snp %>%
  mutate(exp = 's')

###mvmr of eGFR and ma
mv_snp <- rbind(kidney_function_snp, ma_snp)

## 1.checking for overlapping snp(s)
mv_snp <- mv_snp %>% distinct(SNP, .keep_all = T) 
intersect(kidney_function_snp$SNP, ma_snp$SNP) #one overlapping locus (rs4410790)

## 2.identifying independent snps 
#(ref: doi: 10.7554/eLife.43990: combined results clumped by selecting the SNP within a region most strongly associated with kidney function in the Wuttke et al GWAS.) 
cl_kfun <- clump_data(
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
# in the sensitivity analyses, clumping according to microalbuminuria
cl_ma <- clump_data(
  format_data(subset(teumer_2019_sum, RSID %in% mv_snp$SNP),
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

## 3.removing palindromic snps (i.e., snps with A/T or C/G) with intermediate allele frequency (0.42-0.58, inferred from TwoSampleMR::harmonise_data(action = 2))
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
mv_kfun_snp <- mv_snp %>% 
  filter(!SNP %in% p$SNP & SNP %in% cl_kfun$SNP)
mv_ma_snp <- mv_snp %>% 
  filter(!SNP %in% p$SNP & SNP %in% cl_ma$SNP)

## 4.extracting beta, se of each exposure from gwas (effect of a snp on an outcome and exposure must be harmonised to be relative to the same allele)
#kfun
kfun_mv <- format_data(subset(wuttke_2019_sum, RSID %in% mv_kfun_snp$SNP),
                       type = 'outcome',
                       snp_col = 'RSID',
                       chr_col = 'Chr',
                       effect_allele_col = 'Allele1',
                       other_allele_col = 'Allele2',
                       eaf_col = 'Freq1',
                       beta_col = 'Effect',
                       se_col = 'StdErr',
                       pval_col = 'P-value')

tmp <- harmonise_data(mv_kfun_snp, kfun_mv, action = 1) %>% 
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

#ma
ma_mv <- format_data(subset(teumer_2019_sum, RSID %in% tmp$SNP),
                       type = 'outcome',
                       snp_col = 'RSID',
                       chr_col = 'Chr',
                       effect_allele_col = 'Allele1',
                       other_allele_col = 'Allele2',
                       eaf_col = 'Freq1',
                       beta_col = 'Effect',
                       se_col = 'StdErr',
                       pval_col = 'P-value')

tmp <- harmonise_data(tmp, ma_mv, action = 1) %>% 
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
colnames(tmp)[16:18] <- c('beta_ma', 'se_ma', 'pval_ma')

fwrite(tmp, '~/desktop/ckd cancer/ivs/mvmr/mv_kfun_snp.csv')
fwrite(tmp, '~/desktop/ckd cancer/ivs/mvmr/mv_ma_snp.csv')