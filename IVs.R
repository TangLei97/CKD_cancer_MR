###SNPs for eGFR
##GWAS(Wuttke et al 2019) identified 243* SNPs significant in European ancenstry (marginal(p_value) & conditional(Joint p_value))
#Supplementary Table 11	Characteristics of 99% credible sets for the 253* independent SNPs identified in the EA sample.
library(data.table)
library(dplyr)
wuttke_2019_snp <- readxl::read_xlsx('~/Desktop/ckd cancer/Wuttke2019-supplement-Tables.xlsx',sheet = 12) %>% 
  filter(p_value < 5e-08)
#wuttke_2019_snp <- mutate(wuttke_2019_snp,chr=matrix(unlist(strsplit(`Chr/Pos (b37)`,split = ':')),byrow = T,ncol = 2)[,1])
#Summary data
wuttke_2019_sum <- fread('~/Desktop/ckd cancer/summary data/20171017_MW_eGFR_overall_EA_nstud42.dbgap.txt')

##format & clump data using a window of 10,000 kb and maximal linkage disequilibrium of r2 < 0.01 as threshold
kidney_function_snp <- wuttke_2019_sum[RSID %in% wuttke_2019_snp$`RS number`,
                                       .(RSID, Chr, Allele1, Allele2, Freq1, Effect, StdErr, `P-value`, n_total_sum)]
library(TwoSampleMR)
kidney_function_snp <- clump_data(
  format_data(
    kidney_function_snp,
    type = 'exposure',
    snp_col = 'RSID',
    chr_col = 'Chr',
    effect_allele_col = 'Allele1',
    other_allele_col = 'Allele2',
    eaf_col = 'Freq1',
    beta_col = 'Effect',
    se_col = 'StdErr',
    pval_col = 'P-value',
    samplesize_col = 'n_total_sum'
  ),
  clump_r2 = 0.01,pop = "EUR"
)#214 SNPs
fwrite(kidney_function_snp, 'kidney_function_snp.csv')

###SNPs for microalbuminuria (ma)
##GWAS(Teumer et al 2019) identified 17 SNPs significant in trans-ethnic population
#Supplementary Table 3	Index SNP characteristics at 59 loci associated with UACR in trans-ethnic meta-analysis
library(dplyr)
teumer_2019_snp <- readxl::read_xlsx('~/Desktop/ckd cancer/Teumer2019-supplement-3.xlsx') %>% 
  filter(`P-value MA` < 5e-08)
#Summary data
teumer_2019_sum <- fread('~/Desktop/ckd cancer/summary data/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.gz', 
                         fill = T)

##format data
ma_snp <- format_data(
  teumer_2019_snp,
  type = 'exposure',
  snp_col = 'SNP',
  chr_col = 'chr',
  effect_allele_col = 'EA',
  other_allele_col = 'NEA',
  eaf_col = 'FreqEA',
  beta_col = 'Effect MA',
  se_col = 'SE MA',
  pval_col = 'P-value MA',
  samplesize_col = 'N MA'
)#17 SNPs
write.csv(ma_snp, 'ma_snp.csv')

###calculate F statistics
##1. ref: https://github.com/cb12104/adiposity_metabolites_crc
#the F statistic can be calculated as the average of (beta/se), (t statistic) squared for each snp 
kidney_function_snp[, t_stat := beta.exposure/se.exposure]
kidney_function_snp[, f_stat := t_stat^2]
F = mean(kidney_function_snp$f_stat)

ma_snp[, t_stat := beta.exposure/se.exposure]
ma_snp[, f_stat := t_stat^2]
F = mean(ma_snp$f_stat)

##2. 
rsq <- sum(with(kidney_function_snp, 
            (get_r_from_pn(pval.exposure,
                          samplesize.exposure))^2))
k <- nrow(kidney_function_snp); N <- 480698
F <- (rsq*(N - k - 1))/((1 - rsq)*k)

#Use get_r_from_lor for binary traits.
#no available ncase/ncontrol...


