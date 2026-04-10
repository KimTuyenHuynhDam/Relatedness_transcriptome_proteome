# ==============================================================================
# INTEGRATIVE OMICS: BRAIN RNA-SEQ & PLASMA PROTEOMICS 
# TARGET: Systemic Mean Drift and Paired Statistical Comparison
# ==============================================================================

library(tidyverse)
library(openxlsx)
library(ggpubr)
library(gprofiler2)

# 1. SETUP & DATA CLEANING 
# ------------------------------------------------------------------------------

out_v <- "Integrative omics (rm outliers)/Visualizations"
out_t <- "Integrative omics (rm outliers)/Data_Tables"
if(!dir.exists(out_v)) dir.create(out_v, recursive = TRUE)
if(!dir.exists(out_t)) dir.create(out_t, recursive = TRUE)

clean_id <- function(x) {
  x <- gsub("^X", "", x); x <- gsub("[[:punct:] ]", "", x); return(toupper(x))
}

# 2. RNA-SEQ: ALIGNMENT, FILTERING & RESIDUALIZATION
# ------------------------------------------------------------------------------
mice_info <- read.xlsx("RNAseq-brain-F1,F2/Rawdata/mice information.xlsx") %>%
  filter(Sex == "F", !Analysis.ID %in% c('F2Fc_Rep6', 'F2Mm_Rep3')) %>%
  mutate(Match_ID = clean_id(Analysis.ID))

RNA_seq <- read.csv("RNAseq-brain-F1,F2/Rawdata/20260107tpm_All.csv")
rownames(RNA_seq) <- toupper(RNA_seq$X); RNA_seq$X <- NULL
colnames(RNA_seq) <- clean_id(colnames(RNA_seq))

valid_ids <- intersect(mice_info$Match_ID, colnames(RNA_seq))
rna_raw <- RNA_seq[, valid_ids]

# Define groups based on kinship metadata
rna_f1_cols <- mice_info$Match_ID[mice_info$Relatedness == "0"] %>% intersect(valid_ids)
rna_f2_cols <- mice_info$Match_ID[mice_info$Relatedness == "50"] %>% intersect(valid_ids)

# Apply the 70% Group-Wise Filter to Raw RNA Data
keep_rna_groupwise <- apply(rna_raw, 1, function(row) {
  f1_presence <- sum(row[colnames(rna_raw) %in% rna_f1_cols] > 0) / length(rna_f1_cols)
  f2_presence <- sum(row[colnames(rna_raw) %in% rna_f2_cols] > 0) / length(rna_f2_cols)
  return(f1_presence >= 0.7 | f2_presence >= 0.7)
})

# Filter, Log Transform, and Residualize
rna_filtered_log <- log2(rna_raw[keep_rna_groupwise, ] + 1)
rna_meta_subset <- mice_info[match(colnames(rna_filtered_log), mice_info$Match_ID),]

rna_resid_mat <- t(apply(rna_filtered_log, 1, function(gene_row) {
  fit <- lm(as.numeric(gene_row) ~ as.factor(rna_meta_subset$Tissue))
  return(residuals(fit) + mean(gene_row))
}))
rna_resid <- as.data.frame(rna_resid_mat)
colnames(rna_resid) <- colnames(rna_filtered_log)
rownames(rna_resid) <- rownames(rna_filtered_log)

# 3. PROTEOMICS: FILTERING, IMPUTATION & NORMALIZATION
# ------------------------------------------------------------------------------
pro_exp <- read.xlsx("proteomics/All_expression.xlsx")

cols_to_drop <- grep("134_RA1", colnames(pro_exp))
if(length(cols_to_drop) > 0) {
  pro_exp <- pro_exp[, -cols_to_drop]
  message("Technical outlier (Sample 134) successfully removed.")
}

# Condense to unique genes (Mean of raw intensities, NO log transform yet)
pro_raw_means <- pro_exp %>% 
  mutate(Gene = toupper(PG.Genes)) %>% 
  dplyr::select(Gene, contains("PG.Quantity_F")) %>% 
  group_by(Gene) %>% 
  summarise(across(everything(), mean, na.rm = TRUE)) %>% 
  filter(!is.na(Gene), Gene != "NA") %>% 
  column_to_rownames("Gene")

# Define groups based on column names 
pro_f1_cols <- grep("F1$", colnames(pro_raw_means), value = TRUE)
pro_f2_cols <- grep("F2$", colnames(pro_raw_means), value = TRUE)

# Apply the 70% Group-Wise Filter to Raw Proteomic Data
keep_pro_groupwise <- apply(pro_raw_means, 1, function(row) {
  f1_presence <- sum(row[pro_f1_cols] > 0, na.rm = TRUE) / length(pro_f1_cols)
  f2_presence <- sum(row[pro_f2_cols] > 0, na.rm = TRUE) / length(pro_f2_cols)
  return(f1_presence >= 0.7 | f2_presence >= 0.7)
})

# Filter, Log Transform, and Impute
pro_filtered <- pro_raw_means[keep_pro_groupwise, ]
pro_log <- log2(pro_filtered + 1)

pro_log_imputed <- as.data.frame(apply(pro_log, 2, function(x) {
  temp <- x; temp[temp == 0] <- NA
  if(all(is.na(temp))) return(rep(0, length(x))) 
  x[is.na(temp)] <- min(temp, na.rm = TRUE) * 0.95; return(x)
}))
colnames(pro_log_imputed) <- clean_id(colnames(pro_log_imputed))

# 4. STATISTICAL CALCULATIONS (INTEGRATIVE GAP)
# ------------------------------------------------------------------------------
common_genes <- intersect(rownames(rna_resid), rownames(pro_log_imputed))

pro_f1_clean <- clean_id(pro_f1_cols)
pro_f2_clean <- clean_id(pro_f2_cols)

df_fidelity <- data.frame(
  Gene = common_genes,
  RNA_F1 = rowMeans(rna_resid[common_genes, rna_f1_cols], na.rm = TRUE), 
  RNA_F2 = rowMeans(rna_resid[common_genes, rna_f2_cols], na.rm = TRUE),
  PRO_F1 = rowMeans(pro_log_imputed[common_genes, pro_f1_clean], na.rm = TRUE), 
  PRO_F2 = rowMeans(pro_log_imputed[common_genes, pro_f2_clean], na.rm = TRUE)
) %>% mutate(
  Gap_F1 = abs(RNA_F1 - PRO_F1) / sqrt(2),
  Gap_F2 = abs(RNA_F2 - PRO_F2) / sqrt(2),
  Gap_Drift = Gap_F2 - Gap_F1,
  Rel_Change = (Gap_F2 + 0.01) / (Gap_F1 + 0.01)
)

# Threshold: Mean Drift + 1.5 SD + 20% relative change
drift_threshold <- mean(df_fidelity$Gap_Drift, na.rm = TRUE) + (1.5 * sd(df_fidelity$Gap_Drift, na.rm = TRUE))

df_fidelity <- df_fidelity %>%
  mutate(Status = ifelse(Gap_Drift > drift_threshold & Rel_Change > 1.2, 
                         "Decanalized Driver", "Stable"))

write.xlsx(df_fidelity, file.path(out_t, "Full_Integrative_Data_Final.xlsx"))

# 5. VISUALIZATION 1: SYSTEMIC MEAN DRIFT (DENSITY)
# ------------------------------------------------------------------------------
mean_f1 <- mean(df_fidelity$Gap_F1, na.rm = TRUE)
mean_f2 <- mean(df_fidelity$Gap_F2, na.rm = TRUE)
ks_test <- ks.test(df_fidelity$Gap_F1, df_fidelity$Gap_F2)

p1 <- ggplot(df_fidelity) +
  geom_density(aes(x = Gap_F1, fill = "F1 (0%)"), alpha = 0.4) +
  geom_density(aes(x = Gap_F2, fill = "F2 (50%)"), alpha = 0.4) +
  geom_vline(xintercept = mean_f1, color = "#377eb8", linetype = "dashed", size = 1) +
  geom_vline(xintercept = mean_f2, color = "#e41a1c", linetype = "dashed", size = 1) +
  annotate("text", x = mean_f2 * 1.5, y = 0.5, label = paste0("KS p-value: ", format.pval(ks_test$p.value)), fontface = "bold") +
  annotate("text", x = mean_f2 * 1.5, y = 0.4, label = paste0("Mean Shift: ", round(mean_f2 - mean_f1, 3))) +
  scale_fill_manual(values = c("F1 (0%)" = "#377eb8", "F2 (50%)" = "#e41a1c"), name = "Kinship") +
  theme_pubr() + labs(title = "Regulatory Fidelity Density", x = "mRNA-Protein Gap (Decoupling)", y = "Density")

ggsave(file.path(out_v, "Regulatory Fidelity Density.jpg"), p1, width = 7, height = 5)

# 6. VISUALIZATION 2: SYSTEMIC REGULATORY DRIFT (VIOLIN) 
# ------------------------------------------------------------------------------
df_long <- df_fidelity %>%
  dplyr::select(Gene, Gap_F1, Gap_F2) %>%
  pivot_longer(cols = c(Gap_F1, Gap_F2), names_to = "Generation", values_to = "Gap_Value") %>%
  mutate(Generation = factor(Generation, levels = c("Gap_F1", "Gap_F2"), 
                             labels = c("F1 (Canalized)", "F2 (Decanalized)")))

p1b_violin <- ggplot(df_long, aes(x = Generation, y = Gap_Value, fill = Generation)) +
  geom_violin(alpha = 0.5, trim = FALSE, color = "black") +
  geom_boxplot(width = 0.1, color = "black", outlier.shape = NA, alpha = 0.8) +
  geom_jitter(shape = 16, position = position_jitter(0.15), alpha = 0.05, size = 0.1) +
  stat_compare_means(method = "t.test", paired = TRUE, label = "p.format", 
                     label.x = 1.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("F1 (Canalized)" = "#377eb8", "F2 (Decanalized)" = "#e41a1c")) +
  theme_pubr() +
  labs(title = "Systemic Regulatory Drift (mRNA-Protein)", 
       subtitle = "Paired comparison showing distribution and individual gene decoupling",
       y = "Regulatory Gap (|RNA-PRO|/sqrt(2))", x = "")

ggsave(file.path(out_v, "Systemic Regulatory Drift (mRNA-Protein).jpg"), p1b_violin, width = 7, height = 8, dpi = 600)

# 7. DRIVER IDENTIFICATION & PATHWAY ENRICHMENT (BUBBLE PLOT)
# ------------------------------------------------------------------------------
shared_drivers <- df_fidelity %>% filter(Status == "Decanalized Driver")

gostres <- gost(query = shared_drivers$Gene, 
                organism = "pmbairdii", 
                significant = TRUE, 
                user_threshold = 0.05,
                sources = c("GO:BP", "REAC", "KEGG"))

if (!is.null(gostres$result)) {
  
  go_results <- as.data.frame(gostres$result) %>%
    mutate(Fold_Enrichment = (intersection_size / query_size) / (term_size / effective_domain_size)) %>%
    filter(source == "GO:BP") %>%
    arrange(desc(Fold_Enrichment))
  
  write.xlsx(go_results, file.path(out_t, "Pathway_Enrichment.xlsx"))
  
  p_bubble <- ggplot(head(go_results, 15), 
                     aes(x = Fold_Enrichment, y = reorder(term_name, Fold_Enrichment))) +
    geom_point(aes(size = intersection_size, color = p_value)) +
    scale_color_gradient(low = "red", high = "blue", name = "adj. p-value") +
    scale_size_continuous(name = "Gene Count") +
    theme_bw() +
    labs(title = "Functional Impact of Regulatory Drift",
         subtitle = paste0("Enriched Pathways in Decanalized Drivers (n = ", nrow(shared_drivers), ")"),
         x = "Fold Enrichment (Effect Magnitude)", y = "")
  
  ggsave(file.path(out_v, "Functional Impact of Systemic Regulatory Drift.jpg"), 
         p_bubble, width = 11, height = 8, dpi = 600)
  
  cat("Bubble plot generated. Significance encoded by color, Impact by size.")
}