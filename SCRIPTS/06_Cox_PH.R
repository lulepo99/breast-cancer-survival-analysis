library(survival)
library(survminer)
library(ggplot2)
library(tidyverse)


METABRIC_SUBSET$tp53_mut <- as.factor(ifelse(METABRIC_SUBSET$tp53_mut=='0', 0, 1)) 
METABRIC_SUBSET$brca1_mut <- as.factor(ifelse(METABRIC_SUBSET$brca1_mut=='0', 0, 1)) 
METABRIC_SUBSET$brca2_mut <- as.factor(ifelse(METABRIC_SUBSET$brca2_mut=='0', 0, 1)) 
METABRIC_SUBSET$pten_mut <- as.factor(ifelse(METABRIC_SUBSET$pten_mut=='0', 0, 1)) 
METABRIC_SUBSET$pik3ca_mut <- as.factor(ifelse(METABRIC_SUBSET$pik3ca_mut=='0', 0, 1)) 
METABRIC_SUBSET$cdh1_mut <- as.factor(ifelse(METABRIC_SUBSET$cdh1_mut=='0', 0, 1)) 
METABRIC_SUBSET$akt1_mut <- as.factor(ifelse(METABRIC_SUBSET$akt1_mut=='0', 0, 1)) 
METABRIC_SUBSET$akt2_mut <- as.factor(ifelse(METABRIC_SUBSET$akt2_mut=='0', 0, 1)) 

METABRIC_SUBSET$brca_mut <- as.factor(ifelse(METABRIC_SUBSET$brca1_mut == 1 | METABRIC_SUBSET$brca2_mut==1, 1, 0))
METABRIC_SUBSET$akt_mut <- as.factor(ifelse(METABRIC_SUBSET$akt1_mut == 1 | METABRIC_SUBSET$akt2_mut==1, 1, 0))


#METABRIC_wo_NC <- METABRIC_SUBSET[METABRIC_SUBSET$receptor_subtype!='NC',]

METABRIC_COX_SUBSET <- METABRIC_SUBSET[,-c(1,4,8,11,13,22,28,35,36)]
METABRIC_COX_SUBSET <- METABRIC_COX_SUBSET[,c(1:29, 201, 202, 59, 38)]


cox_model <- coxph(Surv(overall_survival_years, overall_survival) ~  age_at_diagnosis + tumor_size + tp53_mut, 
                   data = METABRIC_COX_SUBSET)
summary(cox_model)

list_cox_ph <- list( NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                     NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)


for (i in 1:ncol(METABRIC_COX_SUBSET)) {
  
  if (i == 18 | i == 19) {
    next
  }
  curr_dataset <- METABRIC_COX_SUBSET[,c(18,19,i)]
  print(i)
  curr_model <- coxph(Surv(overall_survival_years, overall_survival) ~ .,
                      data = curr_dataset)
  ph_table <- cox.zph(curr_model)
  list_cox_ph[[i]] <- ph_table
  
}
  
x <- coxph(Surv(overall_survival_years, overall_survival) ~ brca_mut + strata(tp53_mut),
      data = METABRIC_COX_SUBSET[,c(18,19,29,30)])
summary(x)
h_ph <- cox.zph(x)
h_ph
  
  
  
  