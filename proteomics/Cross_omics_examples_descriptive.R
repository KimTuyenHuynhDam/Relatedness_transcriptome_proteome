# Descriptive RNA/plasma concordance for three named manuscript examples.
# Run from the proteomics directory; uses base R only.
# Plasma inputs are source cells from All_expression.xlsx.
# RNA inputs are saved model rows from the overall nested-model output in
# ../RNAseq-brain-F1,F2/Output/nested_anova_results with direction of correlation (remove 2 outliers).xlsx
# This is a targeted description, not an exhaustive overlap screen or a new RNA fit.

plasma <- read.csv("Cross_omics_examples_plasma_inputs.csv", stringsAsFactors = FALSE)
rna <- read.csv("Cross_omics_examples_RNA_results.csv", stringsAsFactors = FALSE)
genes <- c("Tgfb1", "Uap1l1", "Rpl12")
stopifnot(setequal(unique(plasma$Gene), genes), setequal(rna$Gene, genes),
          !anyDuplicated(plasma[c("Gene", "Generation", "SampleID")]),
          !anyDuplicated(rna$Gene), all(plasma$Generation %in% c("F1", "F2")))

retained <- plasma[!(plasma$Generation == "F1" & plasma$SampleID == 134), ]
results <- lapply(genes, function(gene) {
  x <- retained[retained$Gene == gene, ]
  observed <- is.finite(x$Raw_intensity) & x$Raw_intensity > 0
  f1 <- x$Raw_intensity[observed & x$Generation == "F1"]
  f2 <- x$Raw_intensity[observed & x$Generation == "F2"]
  r <- rna[rna$Gene == gene, ]
  stopifnot(length(f1) == 5L, length(f2) == if (gene == "Tgfb1") 6L else 5L,
            length(unique(x$Accession)) == 1L,
            r$RNA_P < 0.05, r$RNA_FDR < 0.10, r$RNA_generation_coefficient < 0)
  data.frame(Gene = gene, Accession = unique(x$Accession),
             nF1_observed = length(f1), nF2_observed = length(f2),
             meanF1_observed = mean(f1), meanF2_observed = mean(f2),
             F2_over_F1_observed = mean(f2) / mean(f1),
             RNA_generation_coefficient = r$RNA_generation_coefficient,
             RNA_P = r$RNA_P, RNA_FDR = r$RNA_FDR,
             Concordant_lower_F2 = mean(f2) < mean(f1) && r$RNA_generation_coefficient < 0)
})
results <- do.call(rbind, results)
write.csv(results, "Cross_omics_examples_descriptive_results.csv", row.names = FALSE)
print(results, row.names = FALSE)
