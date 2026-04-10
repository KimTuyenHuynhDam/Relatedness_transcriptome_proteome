
library(openxlsx)
library(tidyverse)
library(broom)
library(lubridate)
library(glue)
library(ggpubr)
library(purrr)
library(readxl)
library(dplyr)
library(ggplot2)




RNA_seq = read.csv("./Rawdata/20260107tpm_All.csv") %>% dplyr:: select(-c(66, 67))


mice_info = read.xlsx("./Rawdata/mice information.xlsx") %>% 
  filter(! Analysis.ID %in% c(  'F2Fc_Rep6', 'F2Mm_Rep3')) #remove outliers





# Merge RNA_seq and mice_info
# Assume the column "X" in RNA_seq corresponds to "Analysis.ID" in mice_info
RNA_seq_long <- RNA_seq %>%
  pivot_longer(cols = -X, names_to = "Analysis.ID", values_to = "Expression") %>%
  inner_join(mice_info, by = c("Analysis.ID" = "Analysis.ID")) %>% rename(Gene = X)

# Ensure metadata columns are in appropriate formats
RNA_seq_long <- RNA_seq_long %>%
  mutate(
    Relatedness = as.factor(Relatedness),
    Sex = as.factor(Sex),
    Tissue = as.factor(Tissue)
  )

# This function isolates a subgroup, determines the covariate, and runs the models
run_subgroup_analysis <- function(df, subgroup_col, subgroup_val) {
  
  message(glue("Processing Subgroup: {subgroup_col} = {subgroup_val}..."))
  
  # Determine which variable remains as a covariate
  # If we filter by Sex, Tissue is the remaining covariate, and vice versa.
  remaining_var <- ifelse(subgroup_col == "Sex", "Tissue", "Sex")
  
  # Define dynamic formulas
  reduced_form <- as.formula(glue("Expression ~ {remaining_var}"))
  full_form    <- as.formula(glue("Expression ~ Relatedness + {remaining_var}"))
  
  # Run Nested ANOVA
  results <- df %>%
    filter(!!sym(subgroup_col) == subgroup_val) %>%
    group_by(Gene) %>%
    nest() %>%
    mutate(
      reduced_model = map(data, ~ lm(reduced_form, data = .)),
      full_model    = map(data, ~ lm(full_form, data = .)),
      anova_result  = map2(full_model, reduced_model, anova),
      
      # Extract p-value for the 'Relatedness' addition
      p_value = map_dbl(anova_result, ~ .x$`Pr(>F)`[2]),
      
      # Extract Coefficient for directionality
      full_summary = map(full_model, tidy),
      rel_coef = map_dbl(full_summary, ~ {
        coef_row <- .x %>% filter(term == "Relatedness50")
        if(nrow(coef_row) > 0) coef_row$estimate else NA_real_
      })
    ) %>%
    ungroup() %>%
    mutate(
      FDR = p.adjust(p_value, method = "fdr"),
      Correlation_direction = case_when(
        rel_coef > 0 ~ "positive",
        rel_coef < 0 ~ "negative",
        TRUE ~ "NA"
      )
    ) %>%
    dplyr:: select(Gene, p_value, FDR, rel_coef, Correlation_direction)
  
  # Export to Excel
  file_name <- glue("./Output/Nested_ANOVA_Results_{subgroup_col}_{subgroup_val}_adj_by_{remaining_var}.xlsx")
  write.xlsx(results, file_name)
  
  return(results)
}

# --- 3. Execute for All Subgroups ---

# Define the levels we want to iterate over
sex_levels    <- levels(RNA_seq_long$Sex)    # e.g., "M", "F"
tissue_levels <- levels(RNA_seq_long$Tissue) # e.g., "cortex", "midbrain"

# Run for Sex subgroups (Adjusting for Tissue)
sex_results_list <- map(sex_levels, ~ run_subgroup_analysis(RNA_seq_long, "Sex", .x))

# Run for Tissue subgroups (Adjusting for Sex)
tissue_results_list <- map(tissue_levels, ~ run_subgroup_analysis(RNA_seq_long, "Tissue", .x))

message("Analysis complete. All files saved to your working directory.")