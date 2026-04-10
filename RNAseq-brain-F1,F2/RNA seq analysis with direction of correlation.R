
library(openxlsx)
library(tidyverse)
library(broom)
library(lubridate)
library(glue)
library(ggpubr)
library(purrr)
library(readxl)
library(dplyr)





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

# Extend nested model to extract coefficient and direction of correlation
nested_anova <- RNA_seq_long %>%
  group_by(Gene) %>%
  nest() %>%
  mutate(
    reduced_model = map(data, ~ lm(Expression ~ Sex + Tissue, data = .)),
    full_model = map(data, ~ lm(Expression ~ Relatedness + Sex + Tissue, data = .)),
    anova_result = map2(full_model, reduced_model, anova),
    p_value = map_dbl(anova_result, ~ .x$`Pr(>F)`[2]),  # Extract p-value from ANOVA
    
    # Extract full model summary
    full_model_summary = map(full_model, tidy),
    
    # Extract Relatedness coefficient & p-value from full model
    relatedness_coef = map_dbl(full_model_summary, ~ .x$estimate[.x$term == "Relatedness50"]),
    relatedness_p = map_dbl(full_model_summary, ~ .x$p.value[.x$term == "Relatedness50"])
  )




# Adjust p-values for multiple testing (False Discovery Rate)
filtered_nested_anova <- nested_anova %>%
  dplyr:: select(Gene, p_value, relatedness_coef ) %>%
  as.data.frame() %>%
  mutate(FDR = p.adjust(p_value, method = 'fdr'))  %>%
  mutate(Correlation_direction = ifelse(relatedness_coef > 0, "positive", ifelse(relatedness_coef < 0, "negative", "NA" )))




# Save results
write.xlsx(filtered_nested_anova ,
           "./Output/nested_anova_results with direction of correlation (remove 2 outliers).xlsx")


##################


# Get top significant genes (FDR < 0.05)
significant_genes <- filtered_nested_anova %>%
  filter(FDR < 0.05) %>%
  pull(Gene)  # Extract list of significant genes




# Compute adjusted Expression after removing Sex and Tissue effects
RNA_seq_long_adjusted <- RNA_seq_long %>%
  filter(Gene %in% significant_genes) %>%
  
  # Make sure Relatedness is not lost before nesting
  group_by(Gene) %>%
  nest() %>%
  
  # Fit models within the nested data
  mutate(
    # Fit full model (Expression ~ Sex + Tissue) for each gene
    adjusted_model = map(data, ~ lm(Expression ~ Sex + Tissue, data = .)),
    
    # Extract residuals (Expression adjusted for Sex & Tissue)
    residuals = map(adjusted_model, resid),
    
    # Extract mean Expression per gene
    mean_expression = map_dbl(data, ~ mean(.$Expression, na.rm = TRUE))
  ) %>%
  
  # Unnest the data while keeping Relatedness
  unnest(cols = c(data, residuals)) %>%
  
  # Compute adjusted Expression
  mutate(adjusted_expression = residuals + mean_expression) %>%
  
  # Explicitly keep Relatedness in the final dataset
  dplyr:: select(Gene, Analysis.ID, Relatedness, adjusted_expression)


# Save adjusted Expression values
write.csv(RNA_seq_long_adjusted, "./Output/Adjusted_Expression.csv", row.names = FALSE)

# View first few rows
print(head(RNA_seq_long_adjusted))
