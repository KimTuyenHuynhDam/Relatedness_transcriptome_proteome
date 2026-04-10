################################################################################
# 1. LIBRARIES & ENVIRONMENT SETUP
################################################################################
library(openxlsx)
library(tidyverse)
library(stringr)
library(ggplot2)
library(ggpubr)
library(patchwork)


# Define Male IDs for Sex Assignment
male_ids <- c(
  "273", "263", "272", "255", "276", "278", "259", "258", 
  "265", "274", "194", "202", "188", "209", "198", "204", 
  "196", "208", "192", "187", "207", "190"
)

################################################################################
# 2. DATA INGESTION & INITIAL CLEANING
################################################################################
qpcr <- read.xlsx('LUNA_RT-qPCR_SYBR_20260311_115052_796BR02056.xlsx')

# Handle technical replicates and remove water controls
qpcr_clean <- qpcr %>%
  group_by(Sample, Target) %>% 
  filter(!Sample == "water") %>% 
  summarise(Cq_Mean = mean(Cq, na.rm = TRUE), .groups = 'drop')

# Pivot for Delta Cq calculation
qpcr_wide <- qpcr_clean %>%
  pivot_wider(names_from = Target, values_from = Cq_Mean)

################################################################################
# 3. NORMALIZATION (Delta Cq & Delta Delta Cq)
################################################################################

# Calculate Delta Cq (Target - GAPDH)
analysis_results <- qpcr_wide %>%
  mutate(across(c(`IL-1B`, TNFa, IL6), 
                ~ . - GAPDH, 
                .names = "DeltaCq_{.col}"))

# Clean sample names and extract metadata (Treatment/Generation)
analysis_prepared <- analysis_results %>% 
  mutate(
    Treatment = case_when(
      str_detect(Sample, "Ctrl") ~ "Ctrl",
      str_detect(Sample, "LPS") ~ "LPS"
    ),
    Generation = case_when(
      str_detect(Sample, "F1") ~ "F1",
      str_detect(Sample, "F2") ~ "F2"
    )
  )

# Calculate Reference Mean (Global Control Delta Cq)
control_means <- analysis_prepared %>%
  filter(Treatment == "Ctrl") %>%
  summarise(across(starts_with("DeltaCq_"), mean, na.rm = TRUE))

# Calculate Fold Change (2^-DDCq)
final_results <- analysis_prepared %>%
  mutate(
    DDCq_IL1B = `DeltaCq_IL-1B` - control_means$`DeltaCq_IL-1B`,
    Fold_IL1B = 2^(-DDCq_IL1B),
    DDCq_TNFa = DeltaCq_TNFa - control_means$DeltaCq_TNFa,
    Fold_TNFa = 2^(-DDCq_TNFa),
    DDCq_IL6 = DeltaCq_IL6 - control_means$DeltaCq_IL6,
    Fold_IL6 = 2^(-DDCq_IL6)
  )

# Assign Sex based on IDs
final_results <- final_results %>%
  mutate(
    Mice_ID = str_extract(str_trim(Sample), "^[0-9]+"),
    Sex = ifelse(Mice_ID %in% male_ids, "Male", "Female")
  )

# Final Outlier Removal
final_results_no_outlier <- final_results %>%
  filter(Mice_ID != "205")

# Export Cleaned Data
write.xlsx(final_results, 'final_results.xlsx')
write.xlsx(final_results_no_outlier, 'final_results_remove_outlier(205-LPS-F2).xlsx')

################################################################################
# 4. DATA PREPARATION FOR PLOTTING (Long Format)
################################################################################

# Dataset 1: 2^-DeltaCq (Absolute Relative Expression)
baseline_2delta_data <- final_results_no_outlier %>%
  mutate(
    Group = paste(Sex, Generation, sep = "-"),
    Treatment = factor(Treatment, levels = c("Ctrl", "LPS"))
  ) %>%
  pivot_longer(cols = c(`DeltaCq_IL-1B`, DeltaCq_TNFa, DeltaCq_IL6), 
               names_to = "Gene", values_to = "DeltaCq") %>%
  mutate(
    Gene = gsub("DeltaCq_", "", Gene),
    Abs_Relative_Expr = 2^(-DeltaCq)
  )

# Dataset 2: Pairwise setup (Treatment-Generation)
plot_data_pairs <- baseline_2delta_data %>%
  mutate(
    Tx_Gen = factor(paste(Treatment, Generation, sep = "-"),
                    levels = c("Ctrl-F1", "LPS-F1", "Ctrl-F2", "LPS-F2"))
  )

# Dataset 3: Fold Change (LPS Induced Magnitude only)
plot_data_final <- final_results_no_outlier %>%
  filter(Treatment == "LPS") %>%
  mutate(Generation = factor(Generation, levels = c("F1", "F2"))) %>%
  pivot_longer(cols = c(Fold_IL1B, Fold_TNFa, Fold_IL6), 
               names_to = "Gene", values_to = "Expression") %>%
  mutate(Gene = case_when(
    Gene == "Fold_IL1B" ~ "IL-1B",
    Gene == "Fold_IL6"  ~ "IL6",
    Gene == "Fold_TNFa" ~ "TNFa"
  ))

################################################################################
# 5. CUSTOM PLOTTING FUNCTIONS
################################################################################

# Function 1: Multi-Generational Ctrl vs LPS Plot
make_gene_plot <- function(gene_name, data) {
  gene_data <- data %>% 
    filter(Gene == gene_name) %>%
    mutate(Tx_Gen = factor(paste(Treatment, Generation, sep = "-"),
                           levels = c("Ctrl-F1", "LPS-F1", "Ctrl-F2", "LPS-F2")))
  max_val <- max(gene_data$Abs_Relative_Expr, na.rm = TRUE)
  bracket_y <- max_val * 1.2
  anova_y <- max_val * 2.63   
  top_limit <- max_val * 2.8 
  my_comparisons <- list(c("Ctrl-F1", "Ctrl-F2"), c("LPS-F1", "LPS-F2"), 
                         c("Ctrl-F1", "LPS-F1"), c("Ctrl-F2", "LPS-F2"))
  
  p <- ggplot(gene_data, aes(x = Tx_Gen, y = Abs_Relative_Expr, fill = Tx_Gen)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.12, alpha = 0.4, size = 1.2) +
    facet_wrap(~Sex, scales = "free") + 
    stat_compare_means(method = "anova", label.y = anova_y, label.x = 2.5,
                       hjust = 0.5, color = "darkred", size = 5, fontface = "bold") +
    stat_compare_means(comparisons = my_comparisons, method = "t.test", 
                       label = "p.format", label.y = bracket_y,
                       step.increase = 0.35, tip.length = 0.05) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), limits = c(0, top_limit)) +
    scale_fill_manual(values = c("Ctrl-F1"="#fbb4ae", "LPS-F1"="#e41a1c", 
                                 "Ctrl-F2"="#ccebc5", "LPS-F2"="#4daf4a")) +
    theme_bw(base_size = 12) +
    labs(y = gene_name, x = "") + 
    theme(legend.position = "none",
          axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 10)),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
          strip.background = element_rect(fill = "gray95"),
          strip.text = element_text(face = "bold", size = 13),
          panel.grid.minor = element_blank(),
          panel.spacing.y = unit(1.5, "lines"))
  return(p)
}

# Function 2: Control-Only Comparison Plot
make_ctrl_plot <- function(gene_name, data) {
  ctrl_data <- data %>% 
    filter(Gene == gene_name, Treatment == "Ctrl") %>%
    mutate(Sex_Gen = factor(paste(Sex, Generation, sep = "-"),
                            levels = c("Female-F1", "Female-F2", "Male-F1", "Male-F2")))
  max_val <- max(ctrl_data$Abs_Relative_Expr, na.rm = TRUE)
  bracket_y <- max_val * 1.15
  anova_y <- max_val * 2.5
  top_limit <- max_val * 2.6
  ctrl_comparisons <- list(c("Female-F1", "Female-F2"), c("Male-F1", "Male-F2"),
                           c("Female-F1", "Male-F1"), c("Female-F2", "Male-F2"))
  
  p <- ggplot(ctrl_data, aes(x = Sex_Gen, y = Abs_Relative_Expr, fill = Sex_Gen)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    stat_compare_means(method = "anova", label.y = anova_y, label.x = 2.5,
                       hjust = 0.5, color = "darkred", size = 5, fontface = "bold") +
    stat_compare_means(comparisons = ctrl_comparisons, method = "t.test", 
                       label = "p.format", label.y = bracket_y,
                       step.increase = 0.3, tip.length = 0.01) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), limits = c(0, top_limit)) +
    scale_fill_manual(values = c("Female-F1"="#fbb4ae", "Female-F2"="#f781bf", 
                                 "Male-F1"="#d1e5f0", "Male-F2"="#a6cee3")) +
    theme_bw(base_size = 12) +
    labs(y = gene_name, x = "") + 
    theme(legend.position = "none",
          axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 10)),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
          panel.grid.minor = element_blank())
  return(p)
}

# Function 3: Fold Change Comparison Plot
make_fold_change_plot <- function(gene_name, data) {
  fold_data <- data %>% 
    filter(Gene == gene_name, Treatment == "LPS") %>%
    mutate(Sex_Gen = factor(paste(Sex, Generation, sep = "-"),
                            levels = c("Female-F1", "Female-F2", "Male-F1", "Male-F2")))
  max_val <- max(fold_data$Expression, na.rm = TRUE)
  bracket_y <- max_val * 1.15
  anova_y <- max_val * 2.3
  top_limit <- max_val * 2.6
  fold_comparisons <- list(c("Female-F1", "Female-F2"), c("Male-F1", "Male-F2"),
                           c("Female-F1", "Male-F1"), c("Female-F2", "Male-F2"))
  
  p <- ggplot(fold_data, aes(x = Sex_Gen, y = Expression, fill = Sex_Gen)) +
    geom_boxplot(alpha = 0.8, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.5, size = 1.5) +
    stat_compare_means(method = "anova", label.y = anova_y, label.x = 2.5,
                       hjust = 0.5, color = "darkred", size = 5, fontface = "bold") +
    stat_compare_means(comparisons = fold_comparisons, method = "t.test", 
                       label = "p.format", label.y = bracket_y,
                       step.increase = 0.25, tip.length = 0.01) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)), limits = c(0, top_limit)) +
    scale_fill_manual(values = c("Female-F1"="#f781bf", "Female-F2"="#e41a1c", 
                                 "Male-F1"="#a6cee3", "Male-F2"="#377eb8")) +
    theme_bw(base_size = 12) +
    labs(y = paste(gene_name, "\n(Fold Induction)"), x = "") + 
    theme(legend.position = "none",
          axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 10)),
          axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 12),
          panel.grid.minor = element_blank())
  return(p)
}

################################################################################
# 6. FIGURE GENERATION & STITCHING
################################################################################

# --- FIGURE 1: Ctrl vs LPS Stitched Plot ---
p_il1b <- make_gene_plot("IL-1B", plot_data_pairs)
p_il6  <- make_gene_plot("IL6", plot_data_pairs)
p_tnfa <- make_gene_plot("TNFa", plot_data_pairs)

fig1 <- (p_il1b / p_il6 / p_tnfa) + 
  plot_annotation(
    title = 'Multi-Generational Inflammatory Profile',
    subtitle = 'LPS Response vs. Baseline Shifts',
    theme = theme(plot.title = element_text(hjust = 0.5, face="bold", size=22),
                  plot.subtitle = element_text(hjust = 0.5, size=16, face="italic"))
  )
print(fig1)
ggsave("Multi-Generational_Inflammatory_Profile_Ctrl_vs_LPS.png", fig1, width = 8, height = 12, dpi = 300)

# --- FIGURE 2: Baseline Comparison (Ctrl Only) ---
p1 <- make_ctrl_plot("IL-1B", plot_data_pairs)
p2 <- make_ctrl_plot("IL6", plot_data_pairs)
p3 <- make_ctrl_plot("TNFa", plot_data_pairs)

fig2 <- (p1 / p2 / p3) + 
  plot_annotation(
    title = 'Baseline Inflammatory Profile',
    subtitle = 'Generational and Sexual Dimorphism in Resting Expression',
    theme = theme(plot.title = element_text(hjust = 0.5, face="bold", size=20),
                  plot.subtitle = element_text(hjust = 0.5, size=15))
  )
print(fig2)
ggsave("Figure_Ctrl_Baseline_Comparison.png", fig2, width = 8, height = 12, dpi = 300)

# --- FIGURE 3: Fold Change Comparison (LPS Only) ---
p1_fc <- make_fold_change_plot("IL-1B", plot_data_final)
p2_fc <- make_fold_change_plot("IL6", plot_data_final)
p3_fc <- make_fold_change_plot("TNFa", plot_data_final)

fig3 <- (p1_fc / p2_fc / p3_fc) + 
  plot_annotation(
    title = 'LPS-Induced Fold Change Comparison',
    subtitle = 'Magnitude of Inflammatory Response Relative to Internal Control',
    theme = theme(plot.title = element_text(hjust = 0.5, face="bold", size=20),
                  plot.subtitle = element_text(hjust = 0.5, size=15))
  )
print(fig3)
ggsave("Figure_LPS_FoldChange_Comparison.png", fig3, width = 8, height = 12, dpi = 300)