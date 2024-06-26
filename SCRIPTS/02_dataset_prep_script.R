library(dplyr)

# I used "NC" instead of "NA" or "" to create the new columns because these values are not 
# properly NA (they cannot simply be classified)


# Inversion of 0 and 1 in overall_survival

METABRIC_NEW$overall_survival <- ifelse(METABRIC_NEW$overall_survival==0, 1, 0)
METABRIC_NEW_GENES$overall_survival <- ifelse(METABRIC_NEW$overall_survival==0, 1, 0)


# creation of the column for the ki67 proliferation index

METABRIC_NEW <- METABRIC_NEW %>%
  mutate(ki67_proliferation_index = case_when(
    X3.gene_classifier_subtype == "ER+/HER2- High Prolif" ~ "High",
    X3.gene_classifier_subtype == "ER+/HER2- Low Prolif" ~ "Low",
    TRUE ~ "NC"
  ))

METABRIC_NEW_GENES <- METABRIC_NEW_GENES %>%
  mutate(ki67_proliferation_index = case_when(
    X3.gene_classifier_subtype == "ER+/HER2- High Prolif" ~ "High",
    X3.gene_classifier_subtype == "ER+/HER2- Low Prolif" ~ "Low",
    TRUE ~ "NC"
  ))


# Change the value "Positve" in "Positive" in the er_status_measured_by_ihc
# It took 1 hour to understand why the script did not work before. Hate them all .-.

METABRIC_NEW <- METABRIC_NEW %>%
  mutate(er_status_measured_by_ihc = ifelse(er_status_measured_by_ihc == "Positve", 
                                            "Positive", er_status_measured_by_ihc))

METABRIC_NEW_GENES <- METABRIC_NEW_GENES %>%
  mutate(er_status_measured_by_ihc = ifelse(er_status_measured_by_ihc == "Positve", 
                                            "Positive", er_status_measured_by_ihc))


# creation of the column for the receptor_subtype 

METABRIC_NEW <- METABRIC_NEW %>%
  mutate(receptor_subtype = case_when (
    er_status_measured_by_ihc == "Negative" & pr_status == "Negative" & her2_status== "Negative" ~ "TNBC",
    er_status_measured_by_ihc == "Negative" & pr_status == "Negative" & her2_status== "Positive" ~ "HER2+",
    er_status_measured_by_ihc == "Positive" & her2_status == "Negative" & ki67_proliferation_index == "Low" ~ "LumA",
    er_status_measured_by_ihc == "Positive" & her2_status == "Positive" ~ "LumB",
    er_status_measured_by_ihc == "Positive" & her2_status == "Negative" & ki67_proliferation_index == "High" ~ "LumB",
    TRUE ~ "NC"
  ))


METABRIC_NEW_GENES <- METABRIC_NEW_GENES %>%
  mutate(receptor_subtype = case_when (
    er_status_measured_by_ihc == "Negative" & pr_status == "Negative" & her2_status== "Negative" ~ "TNBC",
    er_status_measured_by_ihc == "Negative" & pr_status == "Negative" & her2_status== "Positive" ~ "HER2+",
    er_status_measured_by_ihc == "Positive" & her2_status == "Negative" & ki67_proliferation_index == "Low" ~ "LumA",
    er_status_measured_by_ihc == "Positive" & her2_status== "Positive" ~ "LumB",
    er_status_measured_by_ihc == "Positive" & her2_status == "Negative" & ki67_proliferation_index == "High" ~ "LumB",
    TRUE ~ "NC"
  ))


# moving the new two variables after the other clinical ones

METABRIC_NEW <- METABRIC_NEW[,c(1:34, 208:209,35:207)]

METABRIC_NEW_GENES <- METABRIC_NEW_GENES[,c(1:34, 20811:20812,35:20810)]



# Deleting patients based on NA and "Died of Other Causes" value 

METABRIC_SUBSET <- METABRIC_NEW %>%
  filter(death_from_cancer!="Died of Other Causes", death_from_cancer!="", 
         pam50_._claudin.low_subtype!="NC",
         er_status_measured_by_ihc!="", cellularity!="", !is.na(mutation_count), 
         !is.na(tumor_size), !is.na(neoplasm_histologic_grade))

METABRIC_SUBSET_GENES <- METABRIC_NEW_GENES %>%
  filter(death_from_cancer!="Died of Other Causes", death_from_cancer!="", 
         pam50_._claudin.low_subtype!="NC",
         er_status_measured_by_ihc!="", cellularity!="", !is.na(mutation_count), 
         !is.na(tumor_size), !is.na(neoplasm_histologic_grade))



# Converting age_at_diagnosis in a numeric variable

METABRIC_SUBSET$age_at_diagnosis=as.integer(METABRIC_SUBSET$age_at_diagnosis)

METABRIC_SUBSET_GENES$age_at_diagnosis=as.integer(METABRIC_SUBSET_GENES$age_at_diagnosis)



# Rename the column "overall_survival_months" to "overall_survival_years"
# and convert the survival time from months to year, to make the data easier to interpret and compare.
# The same is performed for the RFS columns

colnames(METABRIC_SUBSET)[colnames(METABRIC_SUBSET) == "overall_survival_months"] <- "overall_survival_years"
METABRIC_SUBSET$overall_survival_years <- METABRIC_SUBSET$overall_survival_years / 12

colnames(METABRIC_SUBSET_GENES)[colnames(METABRIC_SUBSET_GENES) == "overall_survival_months"] <- "overall_survival_years"
METABRIC_SUBSET_GENES$overall_survival_years <- METABRIC_SUBSET_GENES$overall_survival_years / 12


colnames(METABRIC_SUBSET)[colnames(METABRIC_SUBSET) == "RFS_MONTHS"] <- "RFS_years"
METABRIC_SUBSET$RFS_years <- METABRIC_SUBSET$RFS_years / 12

colnames(METABRIC_SUBSET_GENES)[colnames(METABRIC_SUBSET_GENES) == "RFS_MONTHS"] <- "RFS_years"
METABRIC_SUBSET_GENES$RFS_years <- METABRIC_SUBSET_GENES$RFS_years / 12



# Save the new dataframe

write.csv(METABRIC_SUBSET, file = "METABRIC_SUBSET.csv", row.names = FALSE)
