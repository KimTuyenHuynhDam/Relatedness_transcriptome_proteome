# Four named plasma proteins: observed positive raw intensities, without imputation.
# Two-sided, unpaired equal-variance Student t tests; nominal, unadjusted P values.
# Run from the proteomics directory. Source cells refer to All_expression.xlsx.

input <- read.csv("Targeted_abundance_inputs.csv", stringsAsFactors = FALSE)
genes <- c("Rassf2", "Akt1", "Sema3f", "Cpn1")
stopifnot(all(genes %in% input$Gene),
          !anyDuplicated(input[c("Gene", "Generation", "SampleID")]),
          all(input$Generation %in% c("F1", "F2")))

# Exclude the F1 specimen omitted from the main proteomic analyses.
retained <- input[!(input$Generation == "F1" & input$SampleID == 134), ]

results <- lapply(genes, function(gene) {
  x <- retained[retained$Gene == gene, ]
  measured <- is.finite(x$Raw_intensity) & x$Raw_intensity > 0
  f1 <- x$Raw_intensity[measured & x$Generation == "F1"]
  f2 <- x$Raw_intensity[measured & x$Generation == "F2"]
  stopifnot(length(f1) == 5L,
            length(f2) == if (gene == "Rassf2") 5L else 6L,
            length(unique(x$Accession)) == 1L)
  test <- t.test(f1, f2, alternative = "two.sided", paired = FALSE,
                 var.equal = TRUE)
  data.frame(
    Gene = gene, Accession = unique(x$Accession),
    nF1 = length(f1), nF2 = length(f2),
    meanF1 = mean(f1), meanF2 = mean(f2),
    F2_over_F1 = mean(f2) / mean(f1),
    F1_over_F2 = mean(f1) / mean(f2),
    Student_nominal_P = test$p.value,
    Non_detected_F1 = sum(!measured & x$Generation == "F1"),
    Non_detected_F2 = sum(!measured & x$Generation == "F2")
  )
})
results <- do.call(rbind, results)
write.csv(results, "Targeted_abundance_results.csv", row.names = FALSE)
print(results, row.names = FALSE)
