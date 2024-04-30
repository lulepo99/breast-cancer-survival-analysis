# KAPLAN MEIER CURVES AND LOG RANK TESTS

library(survival)
library(survminer)
library(dplyr)

METABRIC_SUBSET$RFS_STATUS <- as.numeric(METABRIC_SUBSET$RFS_STATUS)

# KM curve based on receptor subtype and overall survival (p= 9e-14). From the curve, we can see
# lower survival rates for HER2+ and LumB. This observation aligns with the absence of Trastuzumab 
# in the therapeutic regimens for patients testing positive for HER2

fit_KM_by_receptor <- survfit(Surv(overall_survival_years, overall_survival) ~ 
                                receptor_subtype, data= METABRIC_SUBSET,
                              subset = receptor_subtype %in% c("HER2+", "LumA", "LumB", "TNBC"))

ggsurvplot(fit_KM_by_receptor, risk.table.col = "strata", surv.median.line = "hv", conf.int = FALSE) 

fit_log_rank_receptor <- survdiff(Surv(overall_survival_years, overall_survival) ~ 
                           receptor_subtype, data= METABRIC_SUBSET,
                         subset =  receptor_subtype %in% c("HER2+", "LumA", "LumB", "TNBC"))

fit_log_rank_receptor   



# KM curve based on receptor subtype and RFS (p= 1e-08). The pattern is very similar
# to the previous one

fit_KM_by_RFS_receptor <- survfit(Surv(RFS_years, RFS_STATUS) ~ 
                                receptor_subtype, data= METABRIC_SUBSET,
                              subset = receptor_subtype %in% c("HER2+", "LumA", "LumB", "TNBC"))

ggsurvplot(fit_KM_by_RFS_receptor, risk.table.col = "strata", surv.median.line = "hv", conf.int = FALSE) 

fit_log_rank_RFS_receptor <- survdiff(Surv(RFS_years, RFS_STATUS) ~ 
                           receptor_subtype, data= METABRIC_SUBSET,
                         subset =  receptor_subtype %in% c("HER2+", "LumA", "LumB", "TNBC"))

fit_log_rank_RFS_receptor



# KM curve based on neoplasm histologic grade and overall survival (p= 8e-10). 
# The test is performed considering all the patients. This difference is surely related to
# the distribution of the grades among the receptor subtypes (see the chi-squared test
# previously performed)

fit_KM_by_grade <- survfit(Surv(overall_survival_years, overall_survival) ~ 
                                    neoplasm_histologic_grade, data= METABRIC_SUBSET)

ggsurvplot(fit_KM_by_grade, risk.table.col = "strata", surv.median.line = "hv", conf.int = FALSE) 

fit_log_rank_grade <- survdiff(Surv(overall_survival_years, overall_survival) ~ 
                                        neoplasm_histologic_grade, data= METABRIC_SUBSET)

fit_log_rank_grade



# KM curve based on neoplasm histologic grade and RFS (p= 2e-06). The pattern is very similar
# to the previous one

fit_KM_by_RFS_grade <- survfit(Surv(RFS_years, RFS_STATUS) ~ 
                             neoplasm_histologic_grade, data= METABRIC_SUBSET)

ggsurvplot(fit_KM_by_RFS_grade, risk.table.col = "strata", surv.median.line = "hv", conf.int = FALSE) 

fit_log_rank_RFS_grade <- survdiff(Surv(RFS_years, RFS_STATUS) ~ 
                                 neoplasm_histologic_grade, data= METABRIC_SUBSET)

fit_log_rank_RFS_grade


# The same tests were performed also for cellularity, but in that case we do not have 
# evidence to reject H0

# Laterality is non significant, as we already know from previous from the last semester.
# Same destiny for the histological subtypes (we never assessed it, so I tried just to be sure)




# Despite the number of NA, I also tried for tumor stage, considering its importance
# from a prognostic point of view. We have 975 patients out of 1285 with an available tumor stage
# (1 = 320, 2 = 559, 3= 88, 4 = 8). We can ask to Alessia or Cappozzo if we can proceed
# but I don't think there will be any problems (maybe if we want to introduce it 
# in the Cox model). 

# KM curve based on tumor stage and overall survival (p= <2e-16). It was clearly expected this result

fit_KM_by_stage <- survfit(Surv(overall_survival_years, overall_survival) ~ 
                             tumor_stage, data= METABRIC_SUBSET,
                           subset = tumor_stage %in% c("1", "2", "3", "4"))

ggsurvplot(fit_KM_by_stage, risk.table.col = "strata", surv.median.line = "hv", conf.int = FALSE) 

fit_log_rank_stage <- survdiff(Surv(overall_survival_years, overall_survival) ~ 
                                 tumor_stage, data= METABRIC_SUBSET,
                               subset = tumor_stage %in% c("1", "2", "3", "4"))

fit_log_rank_stage


# KM curve based on tumor stage and RFS (p= <2e-16)

fit_KM_by_RFS_stage <- survfit(Surv(RFS_years, RFS_STATUS) ~ 
                             tumor_stage, data= METABRIC_SUBSET,
                           subset = tumor_stage %in% c("1", "2", "3", "4"))

ggsurvplot(fit_KM_by_RFS_stage, risk.table.col = "strata", surv.median.line = "hv", conf.int = FALSE) 

fit_log_rank_RFS_stage <- survdiff(Surv(RFS_years, RFS_STATUS) ~ 
                                 tumor_stage, data= METABRIC_SUBSET,
                               subset = tumor_stage %in% c("1", "2", "3", "4"))

fit_log_rank_RFS_stage

# Considering the results obtained regarding tumor stage, we can surely stratify 
# within tumor size and pos lymph nodes. Besides, we can also evaluate a KM curve
# based on NPI



# Creation of the column for the categorization of tumor size based on TNM classification 

METABRIC_SUBSET <- METABRIC_SUBSET %>%
  mutate(tumor_size_stage = case_when (
     tumor_size <= 20 ~ "T1",
     tumor_size > 20 & tumor_size <= 50 ~ "T2",
     tumor_size > 50 ~ "T3"
  ))

METABRIC_SUBSET$tumor_size_stage <- factor(METABRIC_SUBSET$tumor_size_stage, level= c("T1", "T2", "T3"))



# KM curve based on tumor_size and overall survival (p= <2e-16)

fit_KM_by_size <- survfit(Surv(overall_survival_years, overall_survival) ~ 
                             tumor_size_stage, data= METABRIC_SUBSET)

ggsurvplot(fit_KM_by_size, risk.table.col = "strata", surv.median.line = "hv", conf.int = FALSE) 

fit_log_rank_size <- survdiff(Surv(overall_survival_years, overall_survival) ~ 
                                tumor_size_stage, data= METABRIC_SUBSET)

fit_log_rank_size