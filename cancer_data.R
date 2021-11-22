####preparation of cancer data
###cancers with available summary statistics from MRCIEU
##lung cancer: ad ieu-a-965; sc 967; all 966
##breast cancer: all ieu-a-1126; ER+ 1127; ER- 1128
##ovary cancer: ieu-a-1120
##endometrium cancer: 
##prostate cancer: ieu-b-85
out_dat <- c('ieu-a-966',
             'ieu-a-965',
             'ieu-a-967',
             'ieu-a-1126',
             'ieu-a-1127',
             'ieu-a-1128',
             'ieu-a-1120',
             'ieu-b-85',
             'ebi-a-GCST006464')

##renal carcinoma
rcc_f <- fread('~/Downloads/Laskar_31231134_Females.gz')
rcc_f[, `:=` (Beta = log(Odds_ratio), 
              Other_allele = ifelse(Allele_1 == Effect_allele, Allele_2, Allele_1))]
rcc_m <- fread('~/Downloads/Laskar_31231134_Males.gz') 
rcc_m[, `:=` (Beta = log(Odds_ratio), 
              Other_allele = ifelse(Allele_1 == Effect_allele, Allele_2, Allele_1))]


###cancers in ukbb
##ovarall diagnosis
#1. from cancer register (NHS Digital for England & Wales, and NHS Central Register for Scotland)
ca_reg <- data.table()
for (i in 0:16) {
  tmp <- fread(paste0('D:/tanglei/7.Health-related outcome/Cancer register_9000002_',
                      i, '.csv'))
  icd10 <- tmp %>% 
    select('eid' |
             starts_with('date_of_cancer_diagnosis_f40005_') |
             starts_with('type_of_cancer_icd10_f40006_'))
  colnames(icd10) <- c('eid', 'dt', 'ca')
  #icd 9 is missing in file 15 & 16
  if (i %in% 0:14) {
    icd9 <- tmp %>% 
      select('eid' |
               starts_with('date_of_cancer_diagnosis_f40005_') |
               starts_with('type_of_cancer_icd9_f40013_'))
    colnames(icd9) <- c('eid', 'dt', 'ca')
    tmp <- rbind(icd9, icd10)
  } else {
    tmp <- icd10
  }
  ca_reg <- rbind(ca_reg,
                  tmp[!is.na(ca)])
}
#Dxx & 2(1-3) refer to benige or uncertain neoplasma
ca_reg <- ca_reg[substr(ca, 1, 1) == 'C' |
                   substr(ca, 1, 2) == '14' |
                   substr(ca, 1, 2) == '15' |
                   substr(ca, 1, 2) == '16' |
                   substr(ca, 1, 2) == '17' |
                   substr(ca, 1, 2) == '18' |
                   substr(ca, 1, 2) == '19' |
                   substr(ca, 1, 2) == '20']

#2. from hospital inpatient (Hospital Episode Statistics for England, Scottish Morbidity Record (SMR), and Patient Episode Database for Wales (PEDW))
ip <- fread('D:/tanglei/7.Health-related outcome/Hospital inpatient_0_0.csv') %>%
  select('eid' |
           starts_with('diagnoses_icd10_f41270_0_') |
           starts_with('date_of_first_inpatient_diagnosis_icd10_f41280_0_') |
           starts_with('diagnoses_icd9_f41271_0_') |
           starts_with('date_of_first_inpatient_diagnosis_icd9_f41281_0_'))

ca_ip <- data.table()
for (i in 2:214) {
  ca <- grep('C|^Z85',
             ip[[i]])
  ca_ip <- rbind(ca_ip,
                 data.table(eid = ip$eid[ca],
                            dt = ip[[i+213]][ca],
                            ca = ip[[i]][ca]))
}
for (i in 428:474) {
  ca <- grep('^1[4-9]|^20|^V10',
             ip[[i]])
  ca_ip <- rbind(ca_ip,
                 data.table(eid = ip$eid[ca],
                            dt = ip[[i+47]][ca],
                            ca = ip[[i]][ca]))
}

#3. from death register
ca_death <- data.table()
for (i in 0:1) {
  death <- fread(paste0('D:/tanglei/7.Health-related outcome/Death register_9000001_',
                        i, '.csv')) %>%
    select('eid' |
             starts_with('date_of_death_f40000_') |
             starts_with('underlying_primary_cause_of_death_icd10_f40001_') |
             starts_with('contributory_secondary_ca_f40002_'))
  for (j in 3:17) {
    #only icd 10
    ca <- grep('C|^Z85',
               death[[j]])
    ca_death <- rbind(ca_death,
                      data.table(eid = death$eid[ca],
                                 dt = death[[2]][ca],
                                 ca = death[[j]][ca]))
  }
}

#4. from self report
ca_self <- data.table()
for (i in 0:3) {
  self <- fread(paste0('D:/tanglei/2.UK Biobank Assessment Centre/Verbal interview_2_',
                       i, '.csv')) %>%
    select('eid' |
             starts_with('cancer_yearage_first_occurred_f84_') |
             starts_with('cancer_code_selfreported_f20001_'))
  for (j in 8:13) {
    ca <- which(!is.na(self[[j]]))
    ca_self <- rbind(ca_self,
                     data.table(eid = self$eid[ca],
                                dt = self[[j-6]][ca],
                                ca = self[[j]][ca]))
  }
}

ca <- rbind(ca_reg[, .(eid, ca)],
            ca_ip[, .(eid, ca)],
            ca_death[, .(eid, ca)],
            ca_self[, .(eid, ca)])
ca_id <- unique(ca$eid)


##site-specific cancers
ca_code <- fread('/Users/tl/desktop/ca_code.csv') [, code := paste0(icd9, '|',
                                                                    icd10, '|',
                                                                    self)]

ca_spec <- ip[, .(eid)]
for (i in 1:19) {
  id <- ca$eid[grep(ca_code$code[i], ca$ca)]
  ca_spec[, ca := ifelse(eid %in% id, 1, 0)]
  setnames(ca_spec, 'ca', ca_code$site[i])
}
#sex-combined
ca_all <- ca_spec[, 1:15][, all := ifelse(eid %in% ca_id, 1, 0)]
#sex-specific
ca_f <- ca_spec[, .(eid, breast, endometrium, ovary, cervix)]
ca_m <- ca_spec[, .(eid, prostate)]

##derive the genetic analysis cohort
geno <- fread('D:/tanglei/4.Genomics/Genotypes_0_0.csv') %>%
  select('eid' |
           'genetic_ethnic_grouping_f22006_0_0',
           'genotype_measurement_batch_f22000_0_0' |
           'genetic_sex_f22001_0_0' |
           starts_with('genetic_principal_components_f22009_0_')) 
geno <- geno[, 1:14]
colnames(geno) <- c('eid', 'caucasian', 'batch', 'sex', paste0('pc', 1:10))

baseline <- fread('D:/tanglei/1.Population characteristics/Baseline characteristics_0_0.csv')
geno[, `:=` (age = baseline$age_at_recruitment_f21022_0_0,
             batch = ifelse(batch > 0, 1, 2))]

#1. Caucasian
#2. matched genetic and self-reported sex
geno[, sex_match := ifelse(sex == baseline$sex_f31_0_0, 1, 0)]
#3. pass QC of genetic data in hpc (missing in a high fraction of subjects/high heterozygosity)
qc <- fread('/Users/tl/desktop/ukb.valid.sample')
geno[, qc := ifelse(eid %in% qc$FID, 1, 0)]
#4. not in 1st/2nd degree of relationships (kinship coefficient > 0.0884)
rel <- fread('/Users/tl/desktop/ukb.rel.id')
#adjust for 3rd relationship
rel3 <- fread('/Users/tl/desktop/ukb.rel3.id')
geno[, `:=` (rel = ifelse(eid %in% rel$ID1, 1, 0),
             rel3 = ifelse(eid %in% rel3$ID1, 1, 0))]
#5. did not withdraw consent
withdraw <- fread('/Users/tl/desktop/w54803_20210809.csv')

geno_qc <- geno[caucasian == 1 & sex_match == 1 & qc == 1 & rel == 0 & !eid %in% withdraw$eid, c(1:15, 19)]

ca_all <- ca_all[geno_qc, on = 'eid', nomatch = 0]
ca_f <- ca_f[geno_qc, on = 'eid', nomatch = 0][sex == 0]
ca_m <- ca_m[geno_qc, on = 'eid', nomatch = 0][sex == 1]

apply(ca_all[, 2:16], 2, sum)
apply(ca_f[, 2:5], 2, sum)
apply(ca_m[, 2], 2, sum)

###calculate snp-ca effect size in ukb
kfun_snp <- fread('/Users/tl/desktop/kidney_function_snp.csv')
fwrite(kfun_snp[, .(SNP, effect_allele.exposure)],
       'kfun_snp.a1', 
       quote = F, col.names = F)
fwrite(kfun_snp[, .(SNP)],
       'kfun.snplist',
       quote = F, col.names = F)
ukb.qc.id <- data.table(FID = geno_qc$eid,
                        IID = geno_qc$eid)
fwrite(ukb.qc.id, 'ukb.qc.id', sep = '\t')

##get genotypes in hpc
##plink \
##--bfile ukb_imp \
##--extract kfun.snplist \
##--keep ukb.qc.id \
##--recode A \
##--recode-allele kfun_snp.a1 \
##--out kfun_snp

kfun_ukb <- fread('kfun_snp.raw') %>%
  select('FID' |
           starts_with('rs'))
setnames(kfun_ukb, 'FID', 'eid')

for (i in 1:213) {
  colnames(kfun_ukb)[i+1] <- unlist(
    strsplit(
      colnames(kfun_ukb)[i+1], '_'))[1]
}

kfun_beta <- data.table(SNP = colnames(kfun_ukb)[-1])
kfun_se <- data.table(SNP = colnames(kfun_ukb)[-1])
#all
for (i in 2:16) {
  tmp <- kfun_ukb[ca_all[, c(1, 17:31)], on = 'eid'][, ca := ca_all[[i]]]
  for (j in 1:213) {
    snp <- tmp[[j+1]]
    mod <- glm(as.factor(ca) ~ snp + age + sex + rel3 + batch + pc1 + pc2 + pc3+ pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, tmp, family = binomial())
    kfun_beta$ca[j] <- summary(mod)$coef[2, 1]
    kfun_se$ca[j] <- summary(mod)$coef[2, 2]
  }
  setnames(kfun_beta, 'ca', colnames(ca_all)[i])
  setnames(kfun_se, 'ca', colnames(ca_all)[i])
}
#female
for (i in 2:5) {
  tmp <- kfun_ukb[ca_f[, c(1, 6:20)], on = 'eid'][, ca := ca_f[[i]]]
  for (j in 1:213) {
    snp <- tmp[[j+1]]
    mod <- glm(as.factor(ca) ~ snp + age + rel3 + batch + pc1 + pc2 + pc3+ pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, tmp, family = binomial())
    kfun_beta$ca[j] <- summary(mod)$coef[2, 1]
    kfun_se$ca[j] <- summary(mod)$coef[2, 2]
  }
  setnames(kfun_beta, 'ca', colnames(ca_f)[i])
  setnames(kfun_se, 'ca', colnames(ca_f)[i])
}
#male
tmp <- kfun_ukb[ca_m, on = 'eid']
for (j in 1:213) {
  snp <- tmp[[j+1]]
  mod <- glm(as.factor(prostate) ~ snp + age + rel3 + batch + pc1 + pc2 + pc3+ pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10, tmp, family = binomial())
  kfun_beta$prostate[j] <- summary(mod)$coef[2, 1]
  kfun_se$prostate[j] <- summary(mod)$coef[2, 2]
}


###same with ma
#ma_beta; ma_se
