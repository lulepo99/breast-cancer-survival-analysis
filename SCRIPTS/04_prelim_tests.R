
# CORRELATION TESTS

# Correlation between positive_lymph_nodes and overall_survival_years in breast cancer deaths
# result: weak negative correlation (-0.2177), but statistically significant (p-value = 1.4*10^-7)

qqnorm(METABRIC_SUBSET$lymph_nodes_examined_positive[METABRIC_SUBSET$overall_survival==1])
qqline(METABRIC_SUBSET$lymph_nodes_examined_positive[METABRIC_SUBSET$overall_survival==1], col="red")

hist(METABRIC_SUBSET$lymph_nodes_examined_positive[METABRIC_SUBSET$overall_survival==1], prob=TRUE, main="Hystogram of Pos Lymph Nodes")

qqnorm(METABRIC_SUBSET$overall_survival_years[METABRIC_SUBSET$overall_survival==1])
qqline(METABRIC_SUBSET$overall_survival_years[METABRIC_SUBSET$overall_survival==1], col="red")

hist(METABRIC_SUBSET$overall_survival_years[METABRIC_SUBSET$overall_survival==1], prob=TRUE, main="Hystogram of OS")

cor.test(METABRIC_SUBSET$lymph_nodes_examined_positive[METABRIC_SUBSET$overall_survival==1], METABRIC_SUBSET$overall_survival_years[METABRIC_SUBSET$overall_survival==1], 
         method="spearman")


# Correlation between tumor_size and overall_survival_years in breast cancer deaths
# result: weak negative correlation (rho= -0.189), but statistically significant (p-value = 5.28*10^-6)

qqnorm(METABRIC_SUBSET$tumor_size[METABRIC_SUBSET$overall_survival==1])
qqline(METABRIC_SUBSET$tumor_size[METABRIC_SUBSET$overall_survival==1], col="red")

hist(METABRIC_SUBSET$tumor_size[METABRIC_SUBSET$overall_survival==1], prob=TRUE, main="Hystogram of Tumor size")

cor.test(METABRIC_SUBSET$tumor_size[METABRIC_SUBSET$overall_survival==1], METABRIC_SUBSET$overall_survival_years[METABRIC_SUBSET$overall_survival==1], 
         method="spearman")


# Correlation between age_at_diagnosis and tumor_size in breast cancer deaths
# result: weak positive correlation (rho= 0.083), but statistically significant (p-value = 0.045)

qqnorm(METABRIC_SUBSET$age_at_diagnosis[METABRIC_SUBSET$overall_survival==1])
qqline(METABRIC_SUBSET$age_at_diagnosis[METABRIC_SUBSET$overall_survival==1], col="red")

hist(METABRIC_SUBSET$age_at_diagnosis[METABRIC_SUBSET$overall_survival==1], prob="TRUE", main="Hystogram of Age at diagnosis")

cor.test(METABRIC_SUBSET$age_at_diagnosis[METABRIC_SUBSET$overall_survival==1], METABRIC_SUBSET$tumor_size[METABRIC_SUBSET$overall_survival==1], 
         method="spearman")


# Correlation between age_at_diagnosis and overall_survival_years in breast cancer deaths
# result: non-significant

cor.test(METABRIC_SUBSET$age_at_diagnosis[METABRIC_SUBSET$overall_survival==1], METABRIC_SUBSET$overall_survival_years[METABRIC_SUBSET$overall_survival==1], 
         method="spearman")


# Correlation between age_at_diagnosis and RFS_years in breast cancer deaths
# result: weak positive correlation (rho= 0.124), but statistically significant (p-value = 0.003)

qqnorm(METABRIC_SUBSET$RFS_years[METABRIC_SUBSET$overall_survival==1])
qqline(METABRIC_SUBSET$RFS_years[METABRIC_SUBSET$overall_survival==1], col="red")

hist(METABRIC_SUBSET$RFS_years[METABRIC_SUBSET$overall_survival==1], prob="TRUE", main="Hystogram of RFS years")

cor.test(METABRIC_SUBSET$age_at_diagnosis[METABRIC_SUBSET$overall_survival==1], METABRIC_SUBSET$RFS_years[METABRIC_SUBSET$overall_survival==1], 
         method="spearman")


# Correlation between tumor size and pos lymph nodes in breast cancer deaths 
# weak positive correlation (rho = 0.267), but statistically significant (p-value = 7.742*10^-11)

cor.test(METABRIC_SUBSET$tumor_size[METABRIC_SUBSET$overall_survival==1], METABRIC_SUBSET$lymph_nodes_examined_positive[METABRIC_SUBSET$overall_survival==1], 
         method="spearman")


# Correlation between tumor size and pos lymph nodes considering all the subset dataset
# weak positive correlation (rho = 0.34), but statistically significant (p-value < 2.2*10^-16).
# It makes more sense using all the patients in this case, considering the absense of variables related
# to censoring.

cor.test(METABRIC_SUBSET$tumor_size, METABRIC_SUBSET$lymph_nodes_examined_positive, 
         method="spearman")


# Correlation between age and pos lymph nodes is non-significant

# Correlation using mutation count and TMB are never significant. Specifically, 
# TMB is computed from mutation count so beware to correlate with each other

# Next step: transformation using logarithms etc to see if we are able to obtain 
# stronger correlations





# ASSOCIATION TESTS

# Association between the receptor subtype and the histologic grade (p-value < 2.2e-16). All the subset patients are considered.

table(METABRIC_SUBSET$receptor_subtype[METABRIC_SUBSET$receptor_subtype!="NC"], METABRIC_SUBSET$neoplasm_histologic_grade[METABRIC_SUBSET$receptor_subtype!="NC"])

chisq.test(METABRIC_SUBSET$receptor_subtype[METABRIC_SUBSET$receptor_subtype!="NC"], METABRIC_SUBSET$neoplasm_histologic_grade[METABRIC_SUBSET$receptor_subtype!="NC"], correct=FALSE)

grade=table(METABRIC_SUBSET$neoplasm_histologic_grade[METABRIC_SUBSET$receptor_subtype!="NC"], METABRIC_SUBSET$receptor_subtype[METABRIC_SUBSET$receptor_subtype!="NC"])

barplot(grade, beside=TRUE, legend=c("Grade 1","Grade 2", "Grade 3"), args.legend=list(x="topleft"), main="Neoplasm Histologic Grade")


# Association between the receptor subtype and cellularity (p-value = 10^-12). All the subset patients are considered.

METABRIC_SUBSET$cellularity <- factor(METABRIC_SUBSET$cellularity, level=c("Low", "Moderate", "High"))

table(METABRIC_SUBSET$receptor_subtype[METABRIC_SUBSET$receptor_subtype!="NC"], METABRIC_SUBSET$cellularity[METABRIC_SUBSET$receptor_subtype!="NC"])

chisq.test(METABRIC_SUBSET$receptor_subtype[METABRIC_SUBSET$receptor_subtype!="NC"], METABRIC_SUBSET$cellularity[METABRIC_SUBSET$receptor_subtype!="NC"], correct=FALSE)

cellularity=table(METABRIC_SUBSET$cellularity[METABRIC_SUBSET$receptor_subtype!="NC"], METABRIC_SUBSET$receptor_subtype[METABRIC_SUBSET$receptor_subtype!="NC"])

barplot(cellularity, beside=TRUE, legend=c("Low", "Moderate", "High"), args.legend=list(x="topleft"), main="Cellularity")


