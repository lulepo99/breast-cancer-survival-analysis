assign("df_with_tmb", a_00_data_clinical_sample)
rm(a_00_data_clinical_sample)
df_with_relapse <- a_00_data_clinical_patient
rm(a_00_data_clinical_patient)


df_with_tmb <- df_with_tmb[-c(1986:2509),]
df_with_tmb$patient_id <- as.integer(substr(df_with_tmb$PATIENT_ID, 4, 7))

total <- merge(METABRIC_RNA_Mutation, df_with_tmb, by ="patient_id")
total <- total[,-c(694:705)]

df_with_relapse <- df_with_relapse[-c(1986:2509),]
df_with_relapse$patient_id <- as.integer(substr(df_with_relapse$PATIENT_ID, 4, 7))

total2 <- merge(total, df_with_relapse, by ="patient_id")
total2 <- total2[,-c(695:716)]

total2$RFS_STATUS <- factor(substr(total2$RFS_STATUS, 1, 1), levels = c("0","1"))

METABRIC_NEW <- total2[,c(1:31,694:696,32:693)]

rm(total,total2,total3,df_with_relapse,df_with_tmb)

METABRIC_NEW <- METABRIC_NEW[,-c(35:523)] 

expr_transposed <- t(a_00_data_mrna_illumina_microarray_zscores_ref_diploid_samples)
colnames(expr_transposed) <- expr_transposed[1,]
expr_transposed <- expr_transposed[-c(1,2),]
expr_transposed <- as.data.frame(expr_transposed)
expr_transposed$patient_id <- as.integer(substr(rownames(expr_transposed), 4, 7))

total3 <- merge(METABRIC_NEW, expr_transposed, by ="patient_id")
total3 <- total3[,c(1:34, 208:20810,35:207)]
METABRIC_NEW_GENES <- total3

rm(expr_transposed, total3, a_00_data_mrna_illumina_microarray_zscores_ref_diploid_samples, a_00_data_mutations)





