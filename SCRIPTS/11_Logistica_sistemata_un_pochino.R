library(pROC)
library(caret)
library(glmnet)
library(car)

# Creation of the variable "relapse_within_5yr"
METABRIC_SUBSET$relapse_within_5yr <- with(METABRIC_SUBSET, ifelse(overall_survival_years <= 5 & overall_survival == 1, 1, ifelse(overall_survival_years > 5, 0, NA)))

# Creation of the relapse dataframe and removal of patients with RFS status = 0 and RFS years < 5 years
METABRIC_SUBSET_relapse <- METABRIC_SUBSET[,c(1:39, 215, 40:45)]

METABRIC_SUBSET_relapse <- METABRIC_SUBSET_relapse[,-c(1,4,5,9,11,13,15,24,25,28,30,31,33,34,35,36)]

METABRIC_SUBSET_relapse <- METABRIC_SUBSET_relapse[!is.na(METABRIC_SUBSET$relapse_within_5yr) & METABRIC_SUBSET$primary_tumor_laterality!=""
                                                   & METABRIC_SUBSET$type_of_breast_surgery!="",]

METABRIC_SUBSET_relapse <- na.omit(METABRIC_SUBSET_relapse)


METABRIC_SUBSET_relapse$cellularity <- factor(METABRIC_SUBSET_relapse$cellularity, levels = c("Low", "Moderate", "High"))
METABRIC_SUBSET_relapse$type_of_breast_surgery <- as.factor(METABRIC_SUBSET_relapse$type_of_breast_surgery)
#METABRIC_SUBSET_relapse$cancer_type_detailed <- as.factor(METABRIC_SUBSET_relapse$cancer_type_detailed)
METABRIC_SUBSET_relapse$chemotherapy <- as.factor(METABRIC_SUBSET_relapse$chemotherapy)
#METABRIC_SUBSET_relapse$cohort <- as.factor(METABRIC_SUBSET_relapse$cohort)
METABRIC_SUBSET_relapse$er_status_measured_by_ihc <- as.factor(METABRIC_SUBSET_relapse$er_status_measured_by_ihc)
METABRIC_SUBSET_relapse$neoplasm_histologic_grade <- factor(METABRIC_SUBSET_relapse$neoplasm_histologic_grade, levels = c("1", "2", "3"))
METABRIC_SUBSET_relapse$her2_status <- as.factor(METABRIC_SUBSET_relapse$her2_status)
#METABRIC_SUBSET_relapse$tumor_other_histologic_subtype <- as.factor(METABRIC_SUBSET_relapse$tumor_other_histologic_subtype)
METABRIC_SUBSET_relapse$hormone_therapy <- as.factor(METABRIC_SUBSET_relapse$hormone_therapy)
METABRIC_SUBSET_relapse$inferred_menopausal_state <- as.factor(METABRIC_SUBSET_relapse$inferred_menopausal_state)
METABRIC_SUBSET_relapse$integrative_cluster <- as.factor(METABRIC_SUBSET_relapse$integrative_cluster)
METABRIC_SUBSET_relapse$primary_tumor_laterality <- as.factor(METABRIC_SUBSET_relapse$primary_tumor_laterality)
#METABRIC_SUBSET_relapse$oncotree_code <- as.factor(METABRIC_SUBSET_relapse$oncotree_code)
METABRIC_SUBSET_relapse$pr_status <- as.factor(METABRIC_SUBSET_relapse$pr_status)
METABRIC_SUBSET_relapse$radio_therapy <- as.factor(METABRIC_SUBSET_relapse$radio_therapy)
#METABRIC_SUBSET_relapse$RFS_STATUS <- as.factor(METABRIC_SUBSET_relapse$RFS_STATUS)
METABRIC_SUBSET_relapse$pik3ca_mut <- as.factor(METABRIC_SUBSET_relapse$pik3ca_mut)
METABRIC_SUBSET_relapse$tp53_mut <- as.factor(METABRIC_SUBSET_relapse$tp53_mut)
METABRIC_SUBSET_relapse$brca_mut <- as.factor(METABRIC_SUBSET_relapse$brca_mut)
METABRIC_SUBSET_relapse$akt_mut <- as.factor(METABRIC_SUBSET_relapse$akt_mut)
METABRIC_SUBSET_relapse$pten_mut <- as.factor(METABRIC_SUBSET_relapse$pten_mut)
METABRIC_SUBSET_relapse$cdh1_mut <- as.factor(METABRIC_SUBSET_relapse$cdh1_mut)

METABRIC_SUBSET_relapse$pam50_._claudin.low_subtype <- as.factor(METABRIC_SUBSET_relapse$pam50_._claudin.low_subtype)
METABRIC_SUBSET_relapse$pam50_._claudin.low_subtype <- relevel(METABRIC_SUBSET_relapse$pam50_._claudin.low_subtype, ref = "LumA")

METABRIC_SUBSET_relapse$NPI_stage <- as.character(METABRIC_SUBSET_relapse$NPI_stage)
METABRIC_SUBSET_relapse$NPI_stage[METABRIC_SUBSET_relapse$NPI_stage %in% c("Excellent", "Good")] <- "Excellent_Good"
METABRIC_SUBSET_relapse$NPI_stage <- factor(METABRIC_SUBSET_relapse$NPI_stage, levels = c("Excellent_Good", "Moderate", "Poor"))
METABRIC_SUBSET_relapse$NPI_stage <- relevel(METABRIC_SUBSET_relapse$NPI_stage, ref = "Excellent_Good")

METABRIC_SUBSET_relapse$tumor_size_stage <- as.character(METABRIC_SUBSET_relapse$tumor_size_stage)
METABRIC_SUBSET_relapse$tumor_size_stage[METABRIC_SUBSET_relapse$tumor_size_stage %in% c("T2", "T3")] <- "T2_T3"
METABRIC_SUBSET_relapse$tumor_size_stage <- factor(METABRIC_SUBSET_relapse$tumor_size_stage, levels = c("T1", "T2_T3"))


covariates <- colnames(METABRIC_SUBSET_relapse)[!colnames(METABRIC_SUBSET_relapse) %in% c("relapse_within_5yr")]
models <- list()
pvalues <- list()

for (covariate in covariates) {
  # glm model
  formula <- as.formula(paste("relapse_within_5yr ~", covariate))
  glm_model <- glm(formula, family=binomial(link=logit), data = METABRIC_SUBSET_relapse)
  
  # save the model
  models[[covariate]] <- glm_model
  
  # take the p-value
  pvalue <- summary(glm_model)$coefficients[2,4]
  
  # save the p-value p-value
  pvalues[[covariate]] <- pvalue
}


test <- p.adjust(pvalues, method = "BH")
p_selection <- ifelse(test <= 0.05, test, NA)
p_significant <- p_selection[!is.na(p_selection)]
significant_variables <- names(p_significant)


selected_variables <- c("age_at_diagnosis", significant_variables, "relapse_within_5yr")
METABRIC_SUBSET_new_relapse <- METABRIC_SUBSET_relapse[, selected_variables, drop = FALSE]





# Riproducibility seed 
set.seed(123)

# Index for the training and test groups
train_index <- createDataPartition(METABRIC_SUBSET_new_relapse$relapse_within_5yr, p = 0.8, list = FALSE, times = 1)

# Creation of the training and test groups
train_set <- METABRIC_SUBSET_new_relapse[train_index, ]
test_set <- METABRIC_SUBSET_new_relapse[-train_index, ]

# Costruire il modello di regressione logistica su train_set solo come prova per valutare co-linearità
glm_model <- glm(relapse_within_5yr ~ . , family = binomial(link = logit), data = train_set)

# Calcolare il VIF
print(vif(glm_model))


# Dato che ci sono variabili collineari, le rimuovo.
train_set <- train_set[,-c(2,10,3,6,8,9,11,14)]
test_set <- test_set[,-c(2,10,3,6,8,9,11,14)]



# Convertire i dati di training in formato matrice
x_train <- model.matrix(relapse_within_5yr ~ ., train_set)[,-1] # rimuovi l'intercept
y_train <- train_set$relapse_within_5yr

# Trovare il valore ottimale di lambda tramite cross-validation
cv.out <- cv.glmnet(x_train, y_train, alpha = 1, family = "binomial", type.measure = "mse")

# Lambda ottimale
lambda_min <- cv.out$lambda.min
cat("Lambda min:", lambda_min, "\n")

# Coefficienti del modello di regressione
print(coef(cv.out, s = lambda_min))





# Fare predizioni sul training set
train_prob <- predict(cv.out, newx = x_train, s = lambda_min, type = "response")

# Assicurati che train_prob sia un vettore
train_prob <- as.vector(train_prob)

# Calcolare la curva ROC e l'AUC sul training set
roc_curve_train <- roc(y_train, train_prob)
auc_train <- auc(roc_curve_train)
print(paste("AUC:", auc_train))
plot(roc_curve_train, main = "ROC Curve for Lasso Logistic Regression on Test Set")
abline(0, 1, col = "red", lwd = 2, lty = 2)

# Ottenere le coordinate del best threshold sul training set
best_coords_train <- coords(roc_curve_train, x = "best", transpose = TRUE)
best_threshold <- best_coords_train["threshold"]



# Convertire i dati di test in formato matrice
x_test <- model.matrix(relapse_within_5yr ~ ., test_set)[,-1] # rimuovi l'intercept
y_test <- test_set$relapse_within_5yr

# Fare predizioni sul test set
lasso_prob <- predict(cv.out, newx = x_test, s = lambda_min, type = "response")

# Assicurati che lasso_prob sia un vettore
lasso_prob <- as.vector(lasso_prob)

# Applicare il miglior threshold selezionato sul training set
lasso_predict <- ifelse(lasso_prob > best_threshold, "1", "0")

# Confusion matrix
conf_matrix <- table(pred = lasso_predict, true = y_test)
print(conf_matrix)

# Calcolare la curva ROC e l'AUC sul test set
roc_curve_test <- roc(y_test, lasso_prob)
auc_test <- auc(roc_curve_test)
cat("AUC sul test set:", auc_test, "\n")

# Plot della curva ROC per il test set
plot(roc_curve_test, main = "ROC Curve for Lasso Logistic Regression on Test Set")
abline(0, 1, col = "red", lwd = 2, lty = 2)
coords(roc_curve_test, x="best", transpose = TRUE)

