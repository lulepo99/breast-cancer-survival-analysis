library(glmnet)

# cox2 è il dataset con tutte le solite osservazioni e le colonne più interessanti

# Dataset ottenuto togliendo tutte le righe con almeno un NA
cox2_wo_null <- na.omit(cox2)

df_x <- cox2_wo_null[,-c(18,19,25)]    # Tolto anche RFS_STATUS

y <- Surv(cox2_wo_null$overall_survival_years, cox2_wo_null$overall_survival)

x <- makeX(df_x) # x <- makeX(df_x, na.impute = TRUE) se si volessero lasciare anche i valori nulli (prenderanno come valore il valore medio della rispettiva colonna)

fit_glmnet <- glmnet(x, y, family = "cox")

plot(fit_glmnet)

coef(fit_glmnet, s = 0.05)

fit_cv.glmnet <- cv.glmnet(x, y, family = "cox", type.measure = 'C')
best_lambda <- fit_cv.glmnet$lambda.min
plot(fit_cv.glmnet)

coef(fit_cv.glmnet, s = best_lambda)
