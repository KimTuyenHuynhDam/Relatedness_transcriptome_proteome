# ==============================================================================
# FINAL DECANALIZATION ANALYSIS & EXPORT PIPELINE
# ==============================================================================
library(tidyverse)
library(car)
library(ggpubr)
library(openxlsx)


# --- 1. SETUP & DATA PREP ---
dir.create("Entropy_Analysis_Results", showWarnings = FALSE)

# 1. DATA LOADING & PRE-PROCESSING ---------------------------------------------
RNA_seq <- read.csv("./Rawdata/20260107tpm_All.csv") %>% dplyr:: select(-c(66, 67)) 

# genes_with_zero_presence <- RNA_seq %>%
#   # Check across all columns except the Gene ID (X)
#   filter(if_all(-X, ~ .x > 0)) %>%
#   pull(X)

# 2. Apply the filter and log-transform

# 
# RNA_seq_strict_ori <- RNA_seq %>%
#   filter(X %in% genes_with_zero_presence)

# Load metadata and remove dendrogram-identified outliers
# F2Fc_Rep6 and F2Mm_Rep3 branched incorrectly, suggesting technical artifacts
mice_info <- read.xlsx("./Rawdata/mice information.xlsx") %>% 
  filter(!Analysis.ID %in% c('F2Fc_Rep6', 'F2Mm_Rep3'))

# Clean IDs and Merge: Extract Mouse_ID by removing 'M' or 'C'
RNA_seq_long <- RNA_seq %>%
  pivot_longer(cols = -X, names_to = "Analysis.ID", values_to = "Expression") %>%
  rename(Gene = X) %>%
  inner_join(mice_info, by = "Analysis.ID") %>%
  mutate(
    # Clean ID to link Cortex/Midbrain samples to the same individual
    Mouse_ID = str_remove(ID, "[MC]"),
    Relatedness = factor(Relatedness, levels = c("0", "50")),
    Sex = as.factor(Sex),
    Tissue = as.factor(Tissue),
    Mouse_ID = as.factor(Mouse_ID)
  )



# Calculate Entropy and Adjusted Metrics
entropy_data <- RNA_seq_long %>%
  group_by(Analysis.ID) %>%
  mutate(Probs = Expression / sum(Expression, na.rm = TRUE)) %>%
  summarize(Raw_Entropy = -sum(Probs * log2(Probs + 1e-12), na.rm = TRUE)) %>%
  inner_join(mice_info, by = "Analysis.ID") %>%
  mutate(Relatedness = factor(Relatedness, levels = c("0", "50")))



# Linear model to isolate noise (Residuals)
fit_global <- lm(Raw_Entropy ~ Sex + Tissue, data = entropy_data)

entropy_data <- entropy_data %>%
  mutate(
    Adj_Global = resid(fit_global),
    Instability_Score = abs(Adj_Global) # Absolute Deviation
  )


entropy_data$Adj_Tissue  <- resid(lm(Raw_Entropy ~ Tissue, data = entropy_data))
entropy_data$Adj_Sex     <- resid(lm(Raw_Entropy ~ Sex, data = entropy_data))

write.xlsx(entropy_data, 'Entropy_Analysis_Results/transcriptomic_Stability_Results.xlsx')

# --- 2. STATISTICAL CALCULATIONS ---

# Isolated Permutation Function
calc_emp_p <- function(df, sex_label) {
  sub_data <- df %>% filter(Sex == sex_label)
  real_f <- leveneTest(Adj_Global ~ Relatedness, data = sub_data)$`F value`[1]
  null_dist <- replicate(10000, {
    shuffled <- sub_data %>% mutate(Relatedness = sample(Relatedness))
    leveneTest(Adj_Global ~ Relatedness, data = shuffled)$`F value`[1]
  })
  return(sum(null_dist >= real_f) / 10000)
}

p_emp_f <- calc_emp_p(entropy_data, "F")
p_emp_m <- calc_emp_p(entropy_data, "M")

# --- 3. GRAPHICAL EXPORTS ---
# --- 1. GLOBAL STATE VS. GLOBAL INSTABILITY ---

# Plot A: Global Raw Entropy (The 'State' remains stable)
g_global_entropy <- ggboxplot(entropy_data, x = "Relatedness", y = "Raw_Entropy", 
                              fill = "Relatedness", palette = "npg", add = "jitter") +
  stat_compare_means(method = "t.test", label = "p.format", label.x = 1.5) +
  labs(title = "Global Transcriptomic State (Raw Entropy)",
       subtitle = "The absolute magnitude of complexity is preserved",
       y = "Shannon Entropy (H)", x = "Relatedness Group") +
  theme_minimal()
ggsave("Entropy_Analysis_Results/G1_Global_Raw_Entropy.png", g_global_entropy, width = 6, height = 6)

# Plot B: Global Residual Entropy (The 'Spread' centers at zero)
g_global_residual_entropy <- ggboxplot(entropy_data, x = "Relatedness", y = "Adj_Global", 
                                       fill = "Relatedness", palette = "npg", add = "jitter") +
  stat_compare_means(method = "t.test", label = "p.format", label.x = 1.5) +
  labs(title = "Global Residual Entropy",
       subtitle = "Means are centered at zero, but variance differs",
       y = "Residual Entropy (Adj for Sex/Tissue)", x = "Relatedness Group") +
  theme_minimal()
ggsave("Entropy_Analysis_Results/G2_Global_Residual_Entropy.png", g_global_residual_entropy, width = 6, height = 6)

# Plot C: Global Instability Score (The 'Decanalization Signal')
g_global_instability <- ggboxplot(entropy_data, x = "Relatedness", y = "Instability_Score", 
                                  fill = "Relatedness", palette = "npg", add = "jitter") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.5) +
  labs(title = "Global Transcriptomic Instability",
       subtitle = "Significant increase in stochastic noise (Levene p = 0.0168)",
       y = "Instability Score (|Residual|)", x = "Relatedness Group") +
  theme_minimal()
ggsave("Entropy_Analysis_Results/G3_Global_Instability.png", g_global_instability, width = 6, height = 6)

#GLOBAL LEVENE PLOT (THE SYSTEM-WIDE VIEW)

global_p <- leveneTest(Adj_Global ~ Relatedness, data = entropy_data)$`Pr(>F)`[1]

g_global <- ggboxplot(entropy_data, x = "Relatedness", y = "Instability_Score", 
                      fill = "Relatedness", palette = "npg", add = "jitter") +
  stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.5) +
  labs(title = "Global Transcriptomic Instability", 
       subtitle = paste("Global Levene's Test: p =", round(global_p, 4)),
       y = "Instability Score (|Residual Entropy|)",
       x = "Relatedness (0 vs 50)") +
  theme_minimal()


ggsave("Entropy_Analysis_Results/Global_Levene_Plot.png", g_global, width = 6, height = 6, dpi = 300)


# Raw Shannon Entropy (Transcriptomic State)

g1 <- ggboxplot(entropy_data, x = "Relatedness", y = "Raw_Entropy", 
                fill = "Relatedness", palette = "npg", facet.by = c("Sex", "Tissue")) +
  stat_compare_means(method = "t.test", label = "p.format") +
  labs(title = "Shannon Entropy by Relatedness", y = "Raw Entropy (H)")
ggsave("Entropy_Analysis_Results/Raw_Entropy_Boxplots.png", g1, width = 8, height = 8, dpi = 300)

# Density Distribution (Visualizing Decanalization)
g2 <- ggplot(entropy_data, aes(x = Adj_Global, fill = Relatedness)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Sex, labeller = as_labeller(c(
    F = paste0("Females (Empirical p = ", round(p_emp_f, 4), ")"),
    M = paste0("Males (Empirical p = ", round(p_emp_m, 4), ")")
  ))) +
  theme_minimal() +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  labs(title = "Residual Entropy Distribution", x = "Adjusted Entropy (Noise Component)")
ggsave("Entropy_Analysis_Results/Density_Decanalization.png", g2, width = 10, height = 5, dpi = 300)

# Instability Score (Tissue-Specific Noise)

g3 <- ggboxplot(entropy_data, x = "Relatedness", y = "Instability_Score", 
                fill = "Relatedness", palette = "npg", facet.by = c("Sex", "Tissue"), add = "jitter") +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  labs(title = "Transcriptional Instability Scores", y = "Instability Score (|Residual|)")
ggsave("Entropy_Analysis_Results/Instability_Scores.png", g3, width = 8, height = 8, dpi = 300)


g_instability_density <- ggplot(entropy_data, aes(x = Instability_Score, fill = Relatedness)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Sex, labeller = as_labeller(c(
    F = paste0("Females (Empirical p = ", round(p_emp_f, 4), ")"),
    M = paste0("Males (Empirical p = ", round(p_emp_m, 4), ")")
  ))) +
  theme_minimal() +
  scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
  labs(
    title = "Distribution of Transcriptional Instability Scores",
    subtitle = "Visualizing the 'long tail' of regulatory failure (Decanalization)",
    x = "Instability Score (|Residual Entropy|)",
    y = "Density"
  )

ggsave("Entropy_Analysis_Results/4_Instability_Density_Distribution.png", 
       g_instability_density, width = 10, height = 5, dpi = 300)

# --- 4. FINAL REPORT TABLE EXPORT ---
report_table <- entropy_data %>%
  group_by(Sex, Tissue, Relatedness) %>%
  summarize(
    Average_Entropy = mean(Raw_Entropy),
    Entropy_SD = sd(Raw_Entropy),
    Median_Instability = median(Instability_Score),
    Group_Size = n(),
    .groups = "drop"
  )

write.xlsx(report_table, "Entropy_Analysis_Results/Final_Report_Summary_Table.xlsx")


###############


# 1. Calculate the F-test p-values
f_test_results <- entropy_data %>%
  group_by(Sex, Tissue) %>%
  summarize(
    p_f_test = var.test(Raw_Entropy ~ Relatedness)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_label = paste0("F-test p = ", round(p_f_test, 4)))

# 2. Generate the plot with adjusted text positioning
g_cv <- report_table %>%
  mutate(CV = (Entropy_SD / Average_Entropy) * 100) %>%
  left_join(f_test_results, by = c("Sex", "Tissue")) %>%
  ggbarplot(x = "Relatedness", y = "CV", fill = "Relatedness", 
            palette = "npg", facet.by = c("Sex", "Tissue")) +
  
  # Adjusted geom_text parameters:
  # x = 1.1 moves the starting point slightly right of the first bar
  # hjust = 0 aligns the text to the left of that point
  # vjust = 1.5 pushes it down slightly from the top margin (Inf)
  geom_text(aes(label = p_label), 
            x = 1.1, y = Inf, 
            hjust = 0, vjust = 2, 
            size = 3.5, fontface = "italic",
            check_overlap = TRUE) +
  
  # Expand y-axis slightly to ensure text has "breathing room"
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  
  labs(title = "Coefficient of Variation (CV) & Variance Comparison",
       subtitle = "Variance stability analysis across biological groups",
       y = "CV (%) = (SD / Mean) * 100") +
  theme_bw()

ggsave("Entropy_Analysis_Results/5_CV_Comparison.png", g_cv, width = 8, height = 6)
