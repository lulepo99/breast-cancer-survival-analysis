library(ggplot2)

frequenza <- table(METABRIC_SUBSET$primary_tumor_laterality)

df <- data.frame(laterality = names(frequenza), frequency = as.numeric(frequenza))

colori <- c("Left" = "#9932CC", "Right" = "pink")

ggplot(df, aes(x = "", y = frequency, fill = laterality)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(frequency / sum(frequency) * 100), "%")), 
            position = position_stack(vjust = 0.5), size = 5) +  # Aggiungi le percentuali
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = colori) +
  labs(title = "Pie Chart: Primary Tumor Laterality", fill = "Laterality") +
  theme(legend.text = element_text(size = 13),  # Imposta la dimensione del testo della legenda
        legend.title = element_text(size = 15))  # Imposta la dimensione del titolo della legenda



# NPI
frequenza <- table(METABRIC_SUBSET$NPI_stage)

df <- data.frame(NPI_stage = names(frequenza), frequency = as.numeric(frequenza))

colori <- c("Poor" = "#9932CC", "Moderate" = "orchid", "Good" = "#ED98D1", "Excellent" = "pink")

ggplot(df, aes(x = "", y = frequency, fill = NPI_stage)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(frequency / sum(frequency) * 100), "%")), 
            position = position_stack(vjust = 0.5), size = 5) +  # Aggiungi le percentuali
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = colori) +
  labs(title = "Pie Chart: NPI", fill = "NPI stage") +
  theme(legend.text = element_text(size = 13),  # Imposta la dimensione del testo della legenda
        legend.title = element_text(size = 15))  # Imposta la dimensione del titolo della legenda


# Stage
frequenza <- table(METABRIC_SUBSET$tumor_stage)

df <- data.frame(tumor_stage = names(frequenza), frequency = as.numeric(frequenza))


colori <- c("1" = "#2c7a7b", "2" = "#38a89d", "3" = "#81c784", "4" = "#c8e6c9")

ggplot(df, aes(x = "", y = frequency, fill = tumor_stage)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(frequency / sum(frequency) * 100), "%")), 
            position = position_stack(vjust = 0.5), size = 5) +  # Aggiungi le percentuali
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = colori) +
  labs(title = "Pie Chart: stage", fill = "tumor stage") +
  theme(legend.text = element_text(size = 13),  # Imposta la dimensione del testo della legenda
        legend.title = element_text(size = 15))  # Imposta la dimensione del titolo della legenda


# BRCA
frequenza <- table(METABRIC_SUBSET$brca_mut)

df <- data.frame(brca_mut = names(frequenza), frequency = as.numeric(frequenza))


colori <- c("0" = "pink", "1" = "darkred")

ggplot(df, aes(x = "", y = frequency, fill = brca_mut)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(frequency / sum(frequency) * 100), "%")), 
            position = position_stack(vjust = 0.5), size = 5) +  # Aggiungi le percentuali
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = colori) +
  labs(title = "Pie Chart: brca", fill = "Brca mutated") +
  theme(legend.text = element_text(size = 13),  # Imposta la dimensione del testo della legenda
        legend.title = element_text(size = 15))  # Imposta la dimensione del titolo della legenda



# TP53
frequenza <- table(METABRIC_SUBSET$tp53_mut)

df <- data.frame(tp53_mut = names(frequenza), frequency = as.numeric(frequenza))
#c("#DDA0DD", "#DA70D6", "#BA55D3", "#9932CC")
#("#d4b4e2", "#c49ac6", "#b380b1", "#9d67a5")

colori <- c("0" = "#d4b4e2", "1" = "purple")

ggplot(df, aes(x = "", y = frequency, fill = tp53_mut)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(frequency / sum(frequency) * 100), "%")), 
            position = position_stack(vjust = 0.5), size = 5) +  # Aggiungi le percentuali
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = colori) +
  labs(title = "Pie Chart: tp53", fill = "tp53 mutated") +
  theme(legend.text = element_text(size = 13),  # Imposta la dimensione del testo della legenda
        legend.title = element_text(size = 15))  # Imposta la dimensione del titolo della legenda





#receptors

frequenza <- table(METABRIC_SUBSET$receptor_subtype)

df <- data.frame(receptor_subtype = names(frequenza), frequency = as.numeric(frequenza))

colori <- c("TNBC" = "#0047AB", "LumA" = "#1E90FF", "LumB" = "#87CEEB", "HER2+" = "#40E0D0")

ggplot(df, aes(x = "", y = frequency, fill = receptor_subtype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(round(frequency / sum(frequency) * 100), "%")), 
            position = position_stack(vjust = 0.5), size = 5) +  # Aggiungi le percentuali
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_manual(values = colori) +
  labs(title = "Pie Chart: receptors", fill = "receptor subtype") +
  theme(legend.text = element_text(size = 13),  # Imposta la dimensione del testo della legenda
        legend.title = element_text(size = 15))  # Imposta la dimensione del titolo della legenda


#age

# Istogramma con stratificazione per la variabile Menopausa
colori_personalizzati <- c("Pre" = "lightgreen", "Post" = "hotpink")

# Definizione delle età di riferimento per l'istogramma
age_ranges <- c(30, 40, 50, 60, 70, 80, 90, 100)

# Creazione del grafico ggplot
ggplot(METABRIC_RNA_Mutation, aes(x = age_at_diagnosis, fill = inferred_menopausal_state)) + 
  geom_histogram(breaks = c(20, age_ranges), color = "black") + 
  geom_vline(aes(xintercept = mean(age_at_diagnosis)), color = "white", linetype = "dashed", size = 1.5) + 
  scale_x_continuous(breaks = seq(20, 100, by = 10)) + 
  scale_y_continuous(breaks = seq(0, 600, by = 100)) + 
  labs(title = "Histogram of Age at Diagnosis", x = "Age at Diagnosis", y = "No. Patients",
       fill = "Inferred Menopausal State") + 
  guides(fill = guide_legend(title = "Inferred Menopausal State")) +
  theme_light() +
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5, size = 20),  # Ingrandisce il titolo
        legend.title = element_text(size = 15),              # Ingrandisce il titolo della legenda
        legend.text = element_text(size = 12),               # Ingrandisce il testo della legenda
        axis.title.x = element_text(size = 15),              # Ingrandisce il testo dell'asse x
        axis.title.y = element_text(size = 15)) +            # Ingrandisce il testo dell'asse y
  
  # Personalizzazione dei colori
  scale_fill_manual(values = colori_personalizzati)


####################################à
ggsurvplot(fit_KM_by_receptor, 
           risk.table.col = "strata", 
           surv.median.line = "hv", 
           conf.int = FALSE,
           xlab = "Time(years)",                  # Modifica il titolo dell'asse x
           ylab = "Survival probability",  # Modifica il titolo dell'asse y
           title = "Curve di sopravvivenza per i sottotipi di recettori",  # Modifica il titolo del grafico
           legend.title = "Strata",  # Modifica il titolo della legenda
           title.size = 16,  # Modifica le dimensioni del testo del titolo
           legend.title.size = 14,  # Modifica le dimensioni del testo del titolo della legenda
           axis.title.size = 12,   # Modifica le dimensioni del testo degli assi
           axis.text.size = 10,    # Modifica le dimensioni del testo sugli assi
           legend.text.size = 12   # Modifica le dimensioni del testo nella legenda
)
