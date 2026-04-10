# ==============================================================================
# PEROMYSCUS PROTEOMICS: DECANALIZATION & INSTABILITY PIPELINE
# Upgraded with Group-Wise Filtering & Technical Outlier Removal
# ==============================================================================

# 1. LOAD REQUIRED PACKAGES
# ------------------------------------------------------------------------------
library(tidyverse)
library(openxlsx)
library(vegan)
library(car)
library(ggrepel)
library(ggpubr) 

# 2. SETUP & DATA PREP
# ------------------------------------------------------------------------------

# Read Raw Expression Data
pro_exp <- read.xlsx("All_expression.xlsx")

# ==============================================================================
# REMOVE TECHNICAL OUTLIERS
# Removing F1 sample 134 (Sample [10]) due to massive mass-spec dropout
# ==============================================================================
cols_to_drop <- grep("134_RA1", colnames(pro_exp))
if(length(cols_to_drop) > 0) {
  pro_exp <- pro_exp[, -cols_to_drop]
  message("Technical outlier (Sample 134) successfully removed.")
}

# Define the Thresholds for sensitivity analysis
thresholds <- c(0, 0.3, 0.7, 0.9)
names(thresholds) <- c("No_Filter", "Filter_30", "Filter_70", "Filter_90")

# Publication Colors (Nature Publishing Group Hex Codes)
pub_colors <- c("F1" = "#E64B35FF", "F2" = "#4DBBD5FF")

# 3. MASTER LOOP: RUN ANALYSIS ACROSS ALL THRESHOLDS
# ------------------------------------------------------------------------------
for (i in seq_along(thresholds)) {
  
  thresh_val <- thresholds[i]
  thresh_name <- names(thresholds)[i]
  
  main_dir <- paste0("Decanalization_Analysis_rm_134_", thresh_name)
  dir.create(paste0(main_dir, "/Figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(paste0(main_dir, "/Data_Tables"), recursive = TRUE, showWarnings = FALSE)
  
  message(paste("Processing:", thresh_name, "..."))
  
  # A. ADVANCED FILTERING: "Present in X% of AT LEAST ONE group"
  exp_cols <- grep("Quantity_F", colnames(pro_exp), value = TRUE)
  df_quant <- pro_exp[, exp_cols]
  
  f1_cols <- grep("F1", colnames(df_quant))
  f2_cols <- grep("F2", colnames(df_quant))
  
  # Calculate presence fraction within each specific group
  f1_presence <- rowSums(df_quant[, f1_cols] > 0) / length(f1_cols)
  f2_presence <- rowSums(df_quant[, f2_cols] > 0) / length(f2_cols)
  
  # Keep protein if it meets the threshold in F1 OR F2
  keep_idx <- (f1_presence >= thresh_val) | (f2_presence >= thresh_val)
  df_filtered <- df_quant[keep_idx, ]
  
  # B. Log2 Transformation
  df_log <- log2(df_filtered + 1)
  
  # C. Imputation: Down-shifted Normal Distribution
  df_imputed <- apply(df_log, 2, function(x) {
    temp <- x
    temp[temp == 0] <- NA
    min_val <- min(temp, na.rm = TRUE)
    x[is.na(temp)] <- min_val * 0.95 
    return(x)
  })
  
  # D. Extract Metadata
  sample_names <- colnames(df_imputed)
  mice_ids <- str_extract(sample_names, "(?<=\\.)\\d+(?=_)")
  groups <- factor(ifelse(grepl("F1", sample_names), "F1", "F2"), levels = c("F1", "F2"))
  
  id_map <- data.frame(Sample = sample_names, Mice_ID = mice_ids, Group = groups)
  
  # ------------------------------------------------------------------------------
  # 4. GLOBAL METRICS: SHANNON ENTROPY & PCA
  # ------------------------------------------------------------------------------
  calc_entropy <- function(col_vec) {
    p <- (2^col_vec) / sum(2^col_vec)
    return(-sum(p * log2(p + 1e-12))) 
  }
  ent_scores <- apply(df_imputed, 2, calc_entropy)
  ent_df <- data.frame(Sample = names(ent_scores), Score = ent_scores) %>%
    left_join(id_map, by = "Sample")
  
  data_t <- t(df_imputed)
  perm_res <- adonis2(data_t ~ groups, method = "euclidean")
  perm_p_val <- perm_res$`Pr(>F)`[1]
  
  pca_res <- prcomp(data_t, scale. = TRUE)
  pca_df <- as.data.frame(pca_res$x) %>% 
    mutate(Sample = rownames(.)) %>%
    left_join(id_map, by = "Sample")
  
  # ------------------------------------------------------------------------------
  # 5. DECANALIZATION METRIC: UNBIASED INSTABILITY SCORE
  # ------------------------------------------------------------------------------
  f1_centroid <- rowMeans(df_imputed[, groups == "F1"])
  f2_centroid <- rowMeans(df_imputed[, groups == "F2"])
  
  calc_intra_dist <- function(x, group_label) {
    if(group_label == "F1") {
      return(sqrt(sum((x - f1_centroid)^2)))
    } else {
      return(sqrt(sum((x - f2_centroid)^2)))
    }
  }
  
  instab_scores <- sapply(1:ncol(df_imputed), function(j) {
    calc_intra_dist(df_imputed[, j], groups[j])
  })
  
  instab_df <- data.frame(Sample = colnames(df_imputed), Score = instab_scores) %>%
    left_join(id_map, by = "Sample")
  
  # Kolmogorov-Smirnov test for Density distributions
  f1_scores <- instab_df %>% filter(Group == "F1") %>% pull(Score)
  f2_scores <- instab_df %>% filter(Group == "F2") %>% pull(Score)
  ks_res <- ks.test(f1_scores, f2_scores)
  
  # ------------------------------------------------------------------------------
  # 6. PER-PROTEIN VARIABILITY: LEVENE'S VOLCANO
  # ------------------------------------------------------------------------------
  v_pvals <- sapply(1:nrow(df_imputed), function(j) {
    leveneTest(as.numeric(df_imputed[j, ]) ~ groups)$`Pr(>F)`[1]
  })
  
  var_ratio <- apply(df_imputed[, groups == "F2"], 1, var) / apply(df_imputed[, groups == "F1"], 1, var)
  
  volcano_df <- data.frame(
    Gene = pro_exp$PG.Genes[keep_idx],
    ProteinID = pro_exp$PG.ProteinGroups[keep_idx],
    Log2VarRatio = log2(var_ratio),
    Pval = v_pvals
  ) %>% 
    filter(!is.na(Gene), Gene != "NA") %>%
    mutate(Status = case_when(Pval < 0.05 ~ "Significantly Unstable (p < 0.05)", TRUE ~ "Stable"))
  
  # ------------------------------------------------------------------------------
  # 7. EXPORT PUBLICATION-READY FIGURES (300 DPI TIFF)
  # ------------------------------------------------------------------------------
  
  # Fig 1: Entropy Boxplot
  p1 <- ggboxplot(ent_df, x = "Group", y = "Score", fill = "Group", 
                  add = "jitter", add.params = list(size = 2, alpha = 0.6)) +
    scale_fill_manual(values = pub_colors) +
    stat_compare_means(method = "t.test", label = "p.format", label.x = 1.5, size = 5) +
    labs(title = "Global Proteomic Complexity", y = "Shannon Entropy (H)") +
    theme_pubr(base_size = 14) + theme(legend.position = "none")
  ggsave(paste0(main_dir, "/Figures/Fig1_Shannon_Entropy.tiff"), p1, width = 5, height = 6, dpi = 300)
  
  # Fig 2: PCA
  p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 4, alpha = 0.8) + 
    stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.15) +
    scale_color_manual(values = pub_colors) + scale_fill_manual(values = pub_colors) +
    theme_pubr(base_size = 14) + 
    labs(title = "Proteomic State Drift (PCA)", 
         subtitle = paste("PERMANOVA p =", format.pval(perm_p_val, digits=3)))
  ggsave(paste0(main_dir, "/Figures/Fig2_PCA_Drift.tiff"), p2, width = 7, height = 5, dpi = 300)
  
  # Fig 3: Instability Boxplot
  p3 <- ggboxplot(instab_df, x = "Group", y = "Score", fill = "Group", 
                  add = "jitter", add.params = list(size = 2, alpha = 0.6)) +
    scale_fill_manual(values = pub_colors) +
    stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.5, size = 5) +
    labs(title = "Proteomic Instability Score", 
         subtitle = "Intra-group stochastic noise",
         y = "Distance from Respective Centroid") +
    theme_pubr(base_size = 14) + theme(legend.position = "none")
  ggsave(paste0(main_dir, "/Figures/Fig3_Instability_Boxplot.tiff"), p3, width = 5, height = 6, dpi = 300)
  
  # Fig 4: Instability Density
  p4 <- ggplot(instab_df, aes(x = Score, fill = Group)) +
    geom_density(alpha = 0.5, color = "black") +
    scale_fill_manual(values = pub_colors) +
    theme_pubr(base_size = 14) +
    labs(title = "Distribution of Instability", 
         subtitle = paste("Kolmogorov-Smirnov test p =", format.pval(ks_res$p.value, digits = 3)),
         x = "Instability Score (Distance)", y = "Density")
  ggsave(paste0(main_dir, "/Figures/Fig4_Instability_Density.tiff"), p4, width = 7, height = 5, dpi = 300)
  
  # Fig 5: Levene's Volcano
  p5 <- ggplot(volcano_df, aes(x = Log2VarRatio, y = -log10(Pval))) +
    geom_point(aes(color = Status, shape = Status), size = 3, alpha = 0.8) +
    scale_color_manual(values = c("Significantly Unstable (p < 0.05)" = "#E64B35FF", "Stable" = "grey70")) +
    scale_shape_manual(values = c(17, 16)) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black", alpha = 0.5) +
    geom_text_repel(data = subset(volcano_df, Pval < 0.005), aes(label = Gene), max.overlaps = 15, size = 4) +
    theme_pubr(base_size = 14) +
    labs(title = "Differential Variability: F2 vs F1", 
         x = expression(Log[2]~Variance~Ratio~(F2/F1)), 
         y = expression(-Log[10]~Levene~P-value)) +
    theme(legend.position = "bottom", legend.title = element_blank())
  ggsave(paste0(main_dir, "/Figures/Fig5_Levene_Volcano.tiff"), p5, width = 8, height = 7, dpi = 300)
  
  # ------------------------------------------------------------------------------
  # 8. EXPORT DATA TABLES
  # ------------------------------------------------------------------------------
  summary_report <- instab_df %>%
    left_join(ent_df %>% dplyr:: select(Sample, Entropy = Score), by = "Sample") %>%
    group_by(Group) %>%
    summarize(
      Average_Entropy = mean(Entropy), Entropy_SD = sd(Entropy),
      Mean_Instability_Score = mean(Score), Instability_SD = sd(Score),
      CV_Entropy = (Entropy_SD / Average_Entropy) * 100,
      Group_Size = n()
    )
  
  write.xlsx(list(
    Summary_Report = summary_report,
    Instability_Scores = instab_df,
    Entropy_Scores = ent_df,
    PCA_Coordinates = pca_df,
    Protein_Variability_Stats = volcano_df
  ), file = paste0(main_dir, "/Data_Tables/Decanalization_Master_Results.xlsx"))
  
  unstable_genes <- volcano_df %>% filter(Pval < 0.05) %>% pull(Gene) %>% unique()
  writeLines(unstable_genes, con = paste0(main_dir, "/Data_Tables/Unstable_Genes_List.txt"))
  
}

message("Pipeline complete. Look for the folder 'Decanalization_Analysis_Filter_70' for your primary publication figures.")