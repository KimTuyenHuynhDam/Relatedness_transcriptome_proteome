# ==============================================================================
# DEFINITIVE QC: TOTAL PROTEINS & STATISTICAL OUTLIER THRESHOLDS
# ==============================================================================
library(tidyverse)
library(openxlsx)
library(ggpubr)


pro_exp <- read.xlsx("All_expression.xlsx")

# Extract quantitative columns
exp_cols <- grep("Quantity_F", colnames(pro_exp), value = TRUE)
df_quant <- pro_exp[, exp_cols]

# 1. Calculate how many proteins have a value > 0 for each mouse
protein_counts <- colSums(df_quant > 0)

# 2. Create a clean data frame for plotting
count_df <- data.frame(
  Sample = colnames(df_quant),
  Proteins_Identified = protein_counts
) %>%
  mutate(
    Group = factor(ifelse(grepl("F1", Sample), "F1", "F2"), levels = c("F1", "F2")),
    Mouse_ID = stringr::str_extract(Sample, "(?<=\\.)\\d+(?=_)") 
  )

# 3. DYNAMICALLY CALCULATE IQR CUTOFFS FOR EACH GROUP
# Q1 - (1.5 * IQR)
iqr_thresholds <- count_df %>%
  group_by(Group) %>%
  summarize(
    Q1 = quantile(Proteins_Identified, 0.25),
    Q3 = quantile(Proteins_Identified, 0.75),
    IQR_Value = Q3 - Q1,
    Lower_Bound = Q1 - (1.5 * IQR_Value)
  )

# 4. Generate the Faceted QC Barplot
p_qc_stats <- ggplot(count_df, aes(x = Proteins_Identified, y = reorder(Mouse_ID, Proteins_Identified), fill = Group)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("F1" = "#E64B35FF", "F2" = "#4DBBD5FF")) +
  
  # Adds the numbers to the ends of the bars
  geom_text(aes(label = Proteins_Identified), hjust = -0.15, size = 4) +
  
  # Split the graph into F1 (top) and F2 (bottom)
  facet_wrap(~Group, scales = "free_y", ncol = 1) +
  
  # Draw the IQR threshold dashed lines dynamically for each group
  geom_vline(data = iqr_thresholds, aes(xintercept = Lower_Bound), 
             linetype = "dashed", color = "black", size = 1) +
  
  # Add text labels right next to the dashed lines to state the exact cutoff number
  geom_text(data = iqr_thresholds, 
            aes(x = Lower_Bound, y = 1.5, label = paste("Outlier Threshold:", round(Lower_Bound, 0))), 
            angle = 90, vjust = -1, size = 4.5, fontface = "italic", color = "black") +
  
  theme_pubr(base_size = 14) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "QC: Proteomic Detection & Outlier Analysis",
       subtitle = "Dashed lines represent the strict lower bound (Q1 - 1.5 * IQR) for each respective genotype.",
       x = "Number of Proteins Detected (>0)",
       y = "Mouse ID") +
  theme(legend.position = "none", # Legend not needed since we faceted
        strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold", size = 14))

# Save the plot
dir.create("Supplementary_Figures", showWarnings = FALSE)
ggsave("Supplementary_Figures/SuppFig_QC_ProteinCounts_Stats.tiff", p_qc_stats, width = 12, height = 7, dpi = 300)

print("Statistical QC Barplot generated!")