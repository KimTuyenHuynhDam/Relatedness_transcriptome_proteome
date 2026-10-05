# Relatedness, brain transcriptomes and plasma proteomes in Peromyscus

Analysis code and data for brain RNA sequencing, plasma proteomics and brain qPCR responses to LPS in *Peromyscus*.

## Study design

Crosses between the BW and SM2 stocks of *P. maniculatus* produced F1 hybrids. F1 sibling crosses produced F2 hybrids. The metadata codes `Relatedness = 0` and `Relatedness = 50` identify F1 and F2, respectively; they are breeding-group labels, not individual genomic measurements.

The brain dataset contains cortex and midbrain specimens from 32 donors, with eight donors in each generation-by-sex group. Of 64 collected specimens, 62 were retained after excluding `F2Fc_Rep6` and `F2Mm_Rep3`. Plasma came from a separate cohort of 12 female donors. Exclusion of F1-134 leaves five F1 and six F2 plasma specimens.

Brain and plasma donors were virgin. After weaning, animals were housed by generation and sex, with up to six animals per cage. No sampled F1 donor was a recorded parent of a sampled F2 donor. Donor ages, recorded parental pairs and founder information are in [Supplementary Table S2](metadata/Supplementary_Table_S2.xlsx). Library-level sequencing and mapping metrics are in [RNA-seq QC](metadata/RNAseq_QC.csv).

## Running the analyses

Run each script from the working directory below so its relative input paths resolve. Required R packages are declared in the scripts. The targeted abundance script uses base R. g:Profiler queries require an internet connection; annotation releases and unseeded permutations can change between executions.

| Analysis | Working directory | Script |
|---|---|---|
| RNA entropy and residual distributions | `RNAseq-brain-F1,F2` | `FINAL DECANALIZATION ANALYSIS and  EXPORT PIPELINE.R` |
| Abundance-weighted coordination (Pc) | `RNAseq-brain-F1,F2` | `weighted_pc_by_sex_tissue.R` |
| Overall expression models | `RNAseq-brain-F1,F2` | `RNA seq analysis with direction of correlation.R` |
| Expression models within sex and brain region | `RNAseq-brain-F1,F2` | `RNA seq analysis with direction of correlation - Sex, Tissue subgroups.R` |
| Coexpression networks | `RNAseq-brain-F1,F2` | `WGCNA analysis (updated).R` |
| Plasma entropy, centroid distances and variance | `proteomics` | `Decanalization_Analysis_rm_outlier.R` |
| Protein detection counts | `proteomics` | `QC Protein Count - detect outliers.R` |
| Four named plasma protein comparisons | `proteomics` | `Targeted_abundance_comparisons.R` |
| Three descriptive RNA/plasma examples | `proteomics` | `Cross_omics_examples_descriptive.R` |
| Brain RNA and plasma protein integration | Repository root | `integrative omics -RNA-SEQ (BRAIN) and PROTEOMICS (PLASMA) (Sensitivity to Coordination).R` |
| Brain qPCR under control and LPS conditions | `LPS treatment` | `qPCR results analysis.R` |

The scripts save numeric tables and plots in their specified output directories. The main plasma figures use `Decanalization_Analysis_rm_134_Filter_70`; the 0%, 30% and 90% filtering branches are separate sensitivity outputs.

## Measurements and statistical comparisons

- **RNA entropy:** Shannon entropy is calculated per specimen. Signed residuals come from `Raw_Entropy ~ Sex + Tissue`; absolute residuals define the instability score. Signed-residual Levene tests and comparisons of absolute residuals assess different quantities. The sex-specific permutation tests shuffle generation labels within sex, pooling brain regions.
- **Pc:** Genes must have positive TPM in all 64 initial specimens before specimen exclusions. Within each generation, sex and region, correlation profiles include self-correlations and are weighted by the partner genes' mean `log2(TPM + 1)` abundance. Pc measures agreement between the weighted F1 and F2 profiles. The Figure 2C gene-selection cutoff is Pc < −0.25; the ranked-tail plots use separate fractional-rank cutoffs. Gene-wise distribution comparisons share genes and correlation inputs and are exploratory; no sample-level randomization calibration of Pc is implemented here.
- **Expression screens:** Nested models use raw TPM. The overall 641-gene set and Figure 3 use nominal P < 0.05 and Benjamini–Hochberg FDR < 0.10. The separate `Adjusted_Expression.csv` export uses FDR < 0.05. Sex-specific models adjust for region; region-specific models adjust for sex. These are separate subset comparisons, not formal interaction tests.
- **WGCNA:** Module–trait correlations use nominal P < 0.05. ME0 is the unassigned gene set. R and P subtitles on eigengene plots describe module–trait correlations; stars denote nominal Welch comparisons between generations. The expression-screen FDR threshold does not apply to module selection.
- **Enrichment tools:** The R pathway branches use g:Profiler. The supplied GO Biological Process results in Supplementary Figure 3 were generated using ShinyGO 0.82 ([Ge et al., 2020](https://doi.org/10.1093/bioinformatics/btz931)). These are distinct enrichment outputs.
- **Plasma:** Proteins are retained if detected in at least 70% of specimens in either generation. Intensities are transformed as `log2(intensity + 1)`. Missing values and zeros are replaced by 0.95 times the minimum positive transformed value within each specimen. Distances are to each specimen's own generation centroid, including that specimen. PCA scales protein variables; PERMANOVA uses the unscaled transformed matrix and its default 999 permutations. Protein variance comparisons use nominal Levene P values and separately calculated log2 variance ratios.
- **Named protein abundances:** `Targeted_abundance_inputs.csv` contains the exact cells for Rassf2, Akt1, Sema3f and Cpn1 from `All_expression.xlsx`, including source-cell references. `Targeted_abundance_comparisons.R` excludes F1-134 and nonpositive or missing intensities, with no imputation or transformation. Two-sided, equal-variance Student tests produce unadjusted P values; fold changes are ratios of arithmetic group means. Rassf2 has five observations in each generation; the other proteins have five F1 and six F2 observations. These four comparisons describe detected values and are separate from the imputed variance analysis. Results are in `Targeted_abundance_results.csv`.
- **RNA/plasma examples:** `Cross_omics_examples_descriptive.R` summarizes Tgfb1, Uap1l1 and Rpl12 using positive observed raw plasma intensities after excluding F1-134, with stored zeros treated as missing and no imputation. Tgfb1 has five F1 and six F2 observations; Uap1l1 and Rpl12 each have five in each generation. Source-cell inputs, saved overall RNA model coefficients/P/FDR, and arithmetic plasma group means are supplied in the accompanying `Cross_omics_examples_*.csv` files. The separate brain and plasma cohorts show concordant lower RNA expression and mean plasma abundance for these examples. No new RNA model is fitted, no plasma significance test is performed by this descriptive script, and it does not establish an exhaustive set of shared targets.
- **Integration:** Gene-level scores compare generation means from separate female brain RNA and female plasma cohorts. RNA adjustment restores each gene's mean after removing the tissue effect. Paired comparisons match each gene's F1 and F2 score. The score depends on assay scales and tissue sources and does not measure post-transcriptional buffering or translation efficiency.
- **LPS qPCR:** Post-LPS relative expression uses a shared control calibrator pooled across sex and generation.

Generation contrasts also encompass segregation and shared ancestry. Two regions from one donor and donors sharing parental ancestry are not independent biological units; the ordinary models do not fully account for these dependencies. A nonsignificant comparison does not establish equivalence.

## Data availability

RNA-sequencing data are deposited in GEO under **GSE319000**, with public release scheduled upon manuscript acceptance. Mass spectrometry proteomics data are deposited through PRIDE under **PXD077328**.


## Manuscript figure sources

The saved files below identify the manuscript outputs. Additional exploratory plots and sensitivity branches are retained for provenance; a similar filename does not make an alternative plot a manuscript panel.

| Figure | Analysis and saved source |
|---|---|
| 1 | RNA entropy script; `RNAseq-brain-F1,F2/Entropy_Analysis_Results` (global, subgroup, density and individual-value displays) |
| 2 | `weighted_pc_by_sex_tissue.R`; `RNAseq-brain-F1,F2/Pc_analysis`; the Venn selection is Pc < -0.25 |
| 3 | Overall/subgroup nested-expression scripts; `RNAseq-brain-F1,F2/Output/Nested_ANOVA_Results_*` and the overall nested-model workbook |
| 4 | `WGCNA analysis (updated).R`; `WGCNA_Final_Results/p_wgcna_heatmap_600dpi.jpg`, `Pathway_Plots_Final/Pathway_Global_Pos.jpg`, `Pathway_Global_Neg.jpg`, and `p_boxplot_ME*.jpg` without computer-name suffixes |
| 5 | Plasma decanalization script, 70% detection branch: `proteomics/Decanalization_Analysis_rm_134_Filter_70`; variance enrichment derives from the saved unstable-protein list |
| 6 | Root integrative script; `integrative omics (rm outliers)/Data_Tables` and `Visualizations` |
| 7A/B | qPCR script and saved qPCR summaries in `LPS treatment`; B uses the pooled-control calibrator |
| Supplementary 4 | WGCNA script; `WGCNA_Final_Results/p_boxplot_ME77.jpg`, `p_boxplot_ME76.jpg` and `p_boxplot_ME6.jpg` |
| Supplementary 5 | Positive-intensity count QC; `proteomics/Supplementary_Figures/SuppFig_QC_ProteinCounts_Stats.tiff`; `plot_qc_detection_counts.R` reproduces the saved display from verified `Protein_detection_counts.csv` |
| Supplementary 6 | qPCR script; `LPS treatment/Multi-Generational_Inflammatory_Profile_Ctrl_vs_LPS.png`; values are `2^-DeltaCq` |

`WGCNA analysis (updated5).R` has the same upstream expression adjustment and module/eigengene calculations, with different network-export steps. It is retained as an earlier network-export implementation; the entry script above is the supported route for the figure-4 outputs. The additional network visualizations are not figure-4 panels. Files named `test_*`, containing a computer-name suffix, and the 0%, 30% and 90% plasma-filter outputs are auxiliary saved outputs, not substitutes for the manuscript sources identified above.

Other experimental assays and assembled figure panels are outside the scope of this repository. Plot layout and labels in the saved manuscript assets may differ from a fresh script export while using the same numerical results. Permutations, jitter, annotation versions and software versions can affect newly generated outputs.

The saved `Pc_analysis/Abundance_Weighted_Pc_Results.xlsx` export and the original Venn/Fisher contingency tables for Figure 2C/D were not present in the supplied analysis files. The Pc script exports that workbook when run, but it does not implement the Venn/Fisher panels. The subgroup expression workbooks underlying the Figure 3A counts are included; the original Fisher test tables/settings for that panel are not separately supplied. The assembled figures' Fisher annotations have therefore not been regenerated from original saved contingency tables in this repository.
