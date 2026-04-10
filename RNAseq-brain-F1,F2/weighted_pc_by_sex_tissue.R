library(tidyverse)
library(readxl)
library(glue)
library(ggpubr)
library(openxlsx)
library(ggplot2)
library(dplyr)

dir.create("Pc_analysis",  recursive = TRUE, showWarnings = FALSE)
dir.create("Output",  recursive = TRUE, showWarnings = FALSE)

# 1. Filtering Step: Keep genes expressed in at least 80% of samples

RNA_seq = read.csv("./Rawdata/20260107tpm_All.csv") %>% select(-c(66, 67)) 
# # Apply log2 transformation with a pseudocount
# mutate(across(-X, ~log2(.x + 1)))
mice_info <- read.xlsx("./Rawdata/mice information.xlsx") %>% 
  filter(!Analysis.ID %in% c('F2Fc_Rep6', 'F2Mm_Rep3'))

valid_sample_ids <- mice_info$Analysis.ID

# # Adjust the 0.8 threshold based on your smallest group size
# presence_threshold <- 0.8
# 
# min_tpm <- 0.1
# 
# # Identify genes that meet the criteria
# keep_genes <- RNA_seq %>%
#   pivot_longer(-X) %>%
#   group_by(X) %>%
#   summarize(prop_present = sum(value > min_tpm) / n()) %>%
#   filter(prop_present >= presence_threshold) %>%
#   pull(X)
# 
# # 2. Apply filter and then Log-transform
# RNA_seq_clean <- RNA_seq %>%
#   filter(X %in% keep_genes) %>%
#   mutate(across(-X, ~log2(.x + 1)))
# 
# message(glue("Filtered out {nrow(RNA_seq) - length(keep_genes)} genes due to low prevalence/sparsity."))



# 1. Strict Zero-Removal Filter
# This identifies genes where EVERY sample has a value > 0
genes_with_zero_presence <- RNA_seq %>%
  # Check across all columns except the Gene ID (X)
  filter(if_all(-X, ~ .x > 0)) %>%
  pull(X)

# 2. Apply the filter and log-transform


RNA_seq_strict_ori_clean <- RNA_seq %>%
  filter(X %in% genes_with_zero_presence) %>%
  select(X, all_of(valid_sample_ids))

# 3. SYNCHRONIZE AND EXPORT LOG2 TPM (Log Space)
RNA_seq_strict_clean <- RNA_seq %>%
  filter(X %in% genes_with_zero_presence) %>%
  select(X, all_of(valid_sample_ids)) %>%
  mutate(across(-X, ~log2(.x + 1)))

# 4. EXPORT TO EXCEL
write.xlsx(RNA_seq_strict_ori_clean, './Output/TPM_Cleaned_NoOutliers.xlsx')
write.xlsx(RNA_seq_strict_clean, './Output/Log2_TPM_Cleaned_NoOutliers.xlsx')

# 3. Report the loss to the console for your records
original_count <- nrow(RNA_seq)
final_count <- length(genes_with_zero_presence)
loss_percent <- round((1 - (final_count / original_count)) * 100, 2)




# Load expression and metadata
# Merge and preprocess


RNA_seq_long_all <- RNA_seq_strict_clean  %>%
  pivot_longer(cols = -X, names_to = "Analysis.ID", values_to = "Expression") %>%
  rename(Gene = X) %>%
  inner_join(mice_info, by = "Analysis.ID") %>%
  mutate(
    Relatedness = as.factor(Relatedness),
    Sex = as.factor(Sex),
    Tissue = as.factor(Tissue),
    Grouping = paste(Sex, Tissue, sep = "_")
  )


# --- Abundance-Weighted Pc Calculation Script ---

# 1. PRE-CALCULATE MEAN TPM PER GENE/GROUP/RELATEDNESS
# This serves as the 'Weighting Factor'
abundance_weights <- RNA_seq_long_all %>%
  group_by(Gene, Grouping, Relatedness) %>%
  summarize(Mean_TPM = mean(Expression, na.rm = TRUE), .groups = "drop")


# 2. REVISED FUNCTION: Weighted Pc Calculation
calculate_weighted_pc_by_grouping <- function(expr_data, weight_data, group_label, level_a = "0", level_b = "50") {
  message("→ Processing group (Weighted): ", group_label)
  
  # Subset and Filter
  df <- expr_data %>% filter(Grouping == group_label, Relatedness %in% c(level_a, level_b))
  if (n_distinct(df$Relatedness) < 2) return(NULL)
  
  # Create wide matrix
  wide <- df %>%
    select(Gene, Analysis.ID, Expression) %>%
    pivot_wider(names_from = Analysis.ID, values_from = Expression) %>%
    column_to_rownames("Gene") %>% as.matrix()
  
  samples_a <- df %>% filter(Relatedness == level_a) %>% pull(Analysis.ID) %>% unique()
  samples_b <- df %>% filter(Relatedness == level_b) %>% pull(Analysis.ID) %>% unique()
  
  # Step A: Standard Correlation Matrices
  cor_a <- cor(t(wide[, samples_a]), use = "p")
  cor_b <- cor(t(wide[, samples_b]), use = "p")
  
  # Step B: Extract Weights for the Partner Genes (columns of the matrix)
  # Ensure weights are in the same order as the matrix rows/cols
  gene_order <- rownames(cor_a)
  
  weights_a <- weight_data %>% 
    filter(Grouping == group_label, Relatedness == level_a) %>%
    slice(match(gene_order, Gene)) %>% pull(Mean_TPM)
  
  weights_b <- weight_data %>% 
    filter(Grouping == group_label, Relatedness == level_b) %>%
    slice(match(gene_order, Gene)) %>% pull(Mean_TPM)
  
  # Step C: Compute Weighted Pc
  # We use map_dbl to iterate through each gene (row)
  weighted_pcs <- map_dbl(1:nrow(cor_a), function(i) {
    # Multiply the correlation vector by the abundance weights
    v_a_weighted <- cor_a[i, ] * weights_a
    v_b_weighted <- cor_b[i, ] * weights_b
    
    # Calculate correlation between weighted profiles
    return(cor(v_a_weighted, v_b_weighted, use = "p"))
  })
  
  tibble(Gene = gene_order, Weighted_Pc = weighted_pcs, Grouping = group_label)
}

# 3. RUN FOR ALL GROUPS
groups <- unique(RNA_seq_long_all$Grouping)
pc_weighted_results <- map_dfr(groups, ~ calculate_weighted_pc_by_grouping(RNA_seq_long_all, abundance_weights, .x))

# 4. SAVE AND EXPORT
pc_weighted_save <- pc_weighted_results %>% 
  pivot_wider(names_from = Grouping, values_from = Weighted_Pc)

write.xlsx(pc_weighted_save, "./Pc_analysis/Abundance_Weighted_Pc_Results.xlsx")





#  VIOLIN PLOT: DISTRIBUTION OF WEIGHTED PC
ggplot(pc_weighted_results, aes(x = Grouping, y = Weighted_Pc, fill = Grouping)) +
  geom_violin(alpha = 0.6, trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.size = 0.5, fill = "white") +
  stat_summary(fun = median, geom = "point", shape = 21, size = 2, fill = "white") +
  stat_compare_means(
    comparisons = list(c("F_midbrain", "F_cortex"),
                       c("M_cortex", "M_midbrain"),
                       c("F_midbrain", "M_midbrain"),
                       c("F_cortex", "M_cortex")), 
                       
    method = "wilcox.test", label = "p.format"
    
  ) +
  theme_minimal() +
  labs(title = "Distribution of Abundance-Weighted Pc Scores",
       y = "Weighted Pc", x = "") +
  scale_fill_manual(values = c("F_cortex"="red", "F_midbrain"="blue", "M_cortex"="orange", "M_midbrain"="purple"))

ggsave("./Pc_analysis/Weighted_Pc_Violin.png", width = 8, height = 6)




# --- PLOTTING EXTREMES FOR WEIGHTED PC ---

# Define the list of thresholds requested
thresholds <- c(0.02, 0.05, 0.1, 0.15, 0.20, 0.25)

# 1. UPDATED FUNCTION for Weighted results
plot_ranked_weighted_extremes <- function(pc_df, top_percent, output_dir = "./Pc_analysis/") {
  
  # Ensure the directory exists
  if(!dir.exists(output_dir)) dir.create(output_dir)
  
  # Calculate ranks and identify extremes
  ranked <- pc_df %>%
    filter(!is.na(Weighted_Pc)) %>%
    group_by(Grouping) %>%
    arrange(Weighted_Pc) %>%
    mutate(
      Rank = row_number(),
      Total = n(),
      Percent = Rank / Total,
      is_extreme = Percent <= top_percent | Percent >= (1 - top_percent)
    ) %>%
    ungroup() %>%
    filter(is_extreme)
  
  # Create the plot
  p <- ggplot(ranked, aes(x = Rank, y = Weighted_Pc, color = Grouping)) +
    geom_line(size = 1.2) +
    # Note: Weighted Pc often has a wider range than standard Pc
    scale_y_continuous(limits = c(min(ranked$Weighted_Pc, na.rm=T), max(ranked$Weighted_Pc, na.rm=T))) +
    labs(
      x = "Ranked Transcripts (Extremes)",
      y = expression("Abundance-Weighted " * P[c]),
      title = glue::glue("Extremes of Functional Preservation ({top_percent * 100}%)"),
      subtitle = "Weighted by Mean TPM across Relatedness groups"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.title = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank() # Hide individual ranks for clarity
    ) +
    scale_color_manual(values = c("F_cortex" = "red", "F_midbrain" = "blue", 
                                  "M_cortex" = "orange", "M_midbrain" = "purple"))
  
  # Save the file
  file_name <- glue::glue("{output_dir}Weighted_Pc_Extremes_{top_percent*100}percent.png")
  ggsave(file_name, plot = p, width = 10, height = 8)
  
  return(p)
}


# 2. ITERATE THROUGH ALL REQUESTED THRESHOLDS
# This uses 'pc_weighted_results' from your previous calculation
map(thresholds, ~ plot_ranked_weighted_extremes(pc_weighted_results, .x))

message("All extreme-rank plots for Weighted Pc have been generated and saved.")