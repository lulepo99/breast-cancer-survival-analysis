chi.test <- chisq.test(METABRIC_SUBSET$overall_survival, METABRIC_SUBSET$RFS_STATUS)
chi.test$observed
chi.test$expected

cor.test(METABRIC_SUBSET$overall_survival_years[METABRIC_SUBSET$overall_survival==1 & METABRIC_SUBSET$RFS_STATUS==1 ], METABRIC_SUBSET$RFS_years[METABRIC_SUBSET$overall_survival==1 & METABRIC_SUBSET$RFS_STATUS==1], method = "pearson")
cor.test(METABRIC_SUBSET$overall_survival_years[METABRIC_SUBSET$overall_survival==1 & METABRIC_SUBSET$RFS_STATUS==1 ], METABRIC_SUBSET$RFS_years[METABRIC_SUBSET$overall_survival==1 & METABRIC_SUBSET$RFS_STATUS==1], method = "spearman")
cor.test(METABRIC_SUBSET$overall_survival_years[METABRIC_SUBSET$overall_survival==1 & METABRIC_SUBSET$RFS_STATUS==1 ], METABRIC_SUBSET$RFS_years[METABRIC_SUBSET$overall_survival==1 & METABRIC_SUBSET$RFS_STATUS==1], method = "kendall")

# Frequency of remission among dead patients
nrow(METABRIC_NEW[METABRIC_NEW$overall_survival==1 & METABRIC_NEW$RFS_STATUS==1,]) / nrow(METABRIC_NEW[METABRIC_NEW$overall_survival==1,])

# Threshold for RFS_years = 10
METABRIC_SUBSET_for_relapse <- METABRIC_SUBSET[!(METABRIC_SUBSET$overall_survival_years <= 20 & 
                                                 METABRIC_SUBSET$overall_survival==0) & 
                                                 METABRIC_SUBSET$RFS_years<=20,]
table(METABRIC_SUBSET_for_relapse$overall_survival, METABRIC_SUBSET_for_relapse$RFS_STATUS)

fisher.test <- fisher.test(METABRIC_SUBSET_for_relapse$overall_survival, 
                       METABRIC_SUBSET_for_relapse$RFS_STATUS)
fisher.test$p.value


library(survival)
library(survminer)

fit_KM_by_RFS_status <- survfit(Surv(overall_survival_years, overall_survival) ~ RFS_STATUS,
                                data=METABRIC_SUBSET)


ggsurvplot(fit_KM_by_RFS_status)









