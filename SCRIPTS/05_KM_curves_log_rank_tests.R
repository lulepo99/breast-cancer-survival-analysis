# KAPLAN MEYER CURVES AND LOG RANK TESTS

library(survival)
library(survminer)


# KM curve based on receptor subtype and overall survival (p= 9e-14). From the curve, we can see
# lower survival rates for HER2+ and LumB. This observation aligns with the absence of Trastuzumab 
# in the therapeutic regimens for patients testing positive for HER2

fit_KM_by_receptor <- survfit(Surv(overall_survival_years, overall_survival) ~ 
                                receptor_subtype, data= METABRIC_SUBSET,
                              subset = receptor_subtype %in% c("HER2+", "LumA", "LumB", "TNBC"))

ggsurvplot(fit_KM_by_receptor, risk.table.col = "strata", surv.median.line = "hv", conf.int = FALSE) 

fit_log_rank <- survdiff(Surv(overall_survival_years, overall_survival) ~ 
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


