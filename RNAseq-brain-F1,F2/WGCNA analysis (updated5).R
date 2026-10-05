# WGCNA MASTER PIPELINE: PEROMYSCUS HYBRID RELATEDNESS

# Includes: Residualization, Cluster Dendrogram, Boxplots, and Native Pathways

library(WGCNA)
library(tidyverse)
library(openxlsx)
library(biomaRt)
library(ggpubr)
library(broom)
library(gprofiler2)

# Force WGCNA correlation to avoid namespace conflicts
cor <- WGCNA::cor

# 1. SETUP & DATA IMPORT --------------------------------------------------
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Set your paths
output_dir <- "./WGCNA_Final_Results"
if(!dir.exists(output_dir)) dir.create(output_dir)

# Data Loading
RNA_seq <- read.csv("./Rawdata/20260107tpm_All.csv") %>% dplyr::select(-c(66, 67))
mice_info <- read.xlsx("./Rawdata/mice information.xlsx") %>% 
  filter(!Analysis.ID %in% c('F2Fc_Rep6', 'F2Mm_Rep3'))

rownames(RNA_seq) <- RNA_seq$X
RNA_seq$X <- NULL
common_samples <- intersect(colnames(RNA_seq), mice_info$Analysis.ID)
RNA_seq <- RNA_seq[, common_samples]
mice_info <- mice_info[match(common_samples, mice_info$Analysis.ID), ]

# 2. COVARIATE ADJUSTMENT & DATA CLEANING ---------------------------------
print("Step 1: Residualizing Expression for Sex and Tissue...")
datExpr_log <- log2(as.matrix(RNA_seq) + 1)

adjust_expression <- function(gene_row, meta) {
  fit <- lm(gene_row ~ as.factor(meta$Sex) + as.factor(meta$Tissue))
  return(residuals(fit) + mean(gene_row)) 
}

datExpr_adjusted <- t(apply(datExpr_log, 1, adjust_expression, meta = mice_info))
datExpr <- as.data.frame(t(datExpr_adjusted))

# Remove genes with zero variance after residualization
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
  print(paste("Cleaned data. Remaining genes:", ncol(datExpr)))
}

# Define the Trait
datTraits <- mice_info %>%
  mutate(Relatedness_Num = ifelse(Relatedness == "50", 1, 0)) %>%
  column_to_rownames("Analysis.ID") %>%
  dplyr::select(Relatedness_Num)

# 3. NETWORK CONSTRUCTION & DENDROGRAM (The Requested Graph) --------------
print("Step 2: Building Network and Dendrogram...")
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

softPower <- ifelse(is.na(sft$powerEstimate), 12, sft$powerEstimate)

# This is the heavy calculation step
net <- blockwiseModules(datExpr, power = softPower, TOMType = "unsigned", 
                        minModuleSize = 30, numericLabels = TRUE, 
                        saveTOMs = TRUE, verbose = 3)

# SAVE THE CLUSTER DENDROGRAM (FIGURE FROM YOUR IMAGE)
jpeg(file.path(output_dir, "p_wgcna_cluster_dendrogram_600dpi.jpg"), 
     width = 10, height = 7, units = "in", res = 600)

# We plot the tree for the first block (standard WGCNA visual)
plotDendroAndColors(net$dendrograms[[1]], 
                    labels2colors(net$colors)[net$blockGenes[[1]]],
                    "Module Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Cluster Dendrogram: Adjusted for Sex/Tissue")
dev.off()

MEs <- net$MEs
moduleColors <- labels2colors(net$colors)
gene_info <- data.frame(Gene = colnames(datExpr), Module = moduleColors)

# 4. MODULE-TRAIT CORRELATION & HEATMAP -----------------------------------
print("Step 3: Calculating Correlations and Heatmap...")
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

rel_cor <- as.matrix(moduleTraitCor[, 1, drop = FALSE])
rel_p   <- as.matrix(moduleTraitPvalue[, 1, drop = FALSE])
sort_order <- order(rel_cor[,1], decreasing = TRUE)
rel_cor_sorted <- rel_cor[sort_order, , drop = FALSE]
rel_p_sorted <- rel_p[sort_order, , drop = FALSE]

jpeg(file.path(output_dir, "p_wgcna_heatmap_600dpi.jpg"), width = 8, height = 15, units = "in", res = 600)
textMatrix <- paste("R: ", round(rel_cor_sorted, 2), " (p: ", signif(rel_p_sorted, 2), ")", sep = "")
dim(textMatrix) <- dim(rel_cor_sorted)

labeledHeatmap(Matrix = rel_cor_sorted, xLabels = "Relatedness (F1 vs F2)", 
               yLabels = rownames(rel_cor_sorted), ySymbols = rownames(rel_cor_sorted), 
               colors = blueWhiteRed(50), textMatrix = textMatrix,
               xLabelsAngle = 0, xLabelsAdj = 0.5, cex.text = 0.6, zlim = c(-1,1))
dev.off()

# 5. MAPPING & SUMMARY EXPORT ---------------------------------------------
mapping_table <- data.frame(Number = as.character(net$colors),
                            Color = as.character(labels2colors(net$colors))) %>% distinct()

results_summary <- data.frame(ME_Name = rownames(rel_cor_sorted),
                              R_Value = rel_cor_sorted[,1],
                              P_Value = rel_p_sorted[,1]) %>%
  filter(P_Value < 0.05) %>%
  arrange(desc(R_Value)) %>%
  mutate(Mod_Num = as.character(gsub("ME", "", ME_Name))) %>%
  left_join(mapping_table, by = c("Mod_Num" = "Number")) %>%
  mutate(Rank_Pos = row_number(), Rank_Neg = rev(row_number()))

write.xlsx(results_summary, file.path(output_dir, "table_SIGNIFICANT_MODULE_SUMMARY.xlsx"))

# 6. MASTER PATHWAY LOOP (FDR < 0.05, Sort by Fold Enrichment) ------------
print("Step 4: Running Native Pathway Analysis...")
pathway_excel_list <- list()
fig_path <- file.path(output_dir, "Pathway_Plots_Final")
if(!dir.exists(fig_path)) dir.create(fig_path)

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.0001) return(formatC(p, format = "e", digits = 1))
  return(as.character(round(p, 4)))
}

# Group creation
groups_to_analyze <- list("Global_Pos" = (gene_info %>% filter(Module %in% results_summary$Color[results_summary$R_Value > 0]) %>% pull(Gene)),
                          "Global_Neg" = (gene_info %>% filter(Module %in% results_summary$Color[results_summary$R_Value < 0]) %>% pull(Gene)))

for(i in seq_len(nrow(results_summary))) {
  groups_to_analyze[[results_summary$ME_Name[i]]] <- (gene_info %>% filter(Module == results_summary$Color[i]) %>% pull(Gene))
}

for(group_name in names(groups_to_analyze)) {
  query_genes <- as.character(na.omit(groups_to_analyze[[group_name]]))
  if(length(query_genes) < 5) next 
  
  gostres <- gost(query = query_genes, organism = "pmbairdii", significant = TRUE, sources = c("GO:BP", "REAC", "KEGG"))
  if(is.null(gostres)) next
  
  res_df <- as.data.frame(gostres$result) %>%
    mutate(Fold_Enrichment = (intersection_size / query_size) / (term_size / 31878)) %>%
    dplyr::select(term_id, term_name, p_value, Fold_Enrichment, intersection_size, term_size, source) %>%
    arrange(desc(Fold_Enrichment)) 
  
  pathway_excel_list[[substr(group_name, 1, 31)]] <- res_df
  
  # Visualization
  if(grepl("Global", group_name)) {
    p_title <- group_name; p_sub <- paste(length(query_genes), "genes")
  } else {
    m_info <- results_summary %>% filter(ME_Name == group_name)
    r_lbl <- ifelse(m_info$Rank_Pos == 1, " (Top Pos)", ifelse(m_info$Rank_Neg == 1, " (Top Neg)", ""))
    p_title <- paste("Module:", group_name, r_lbl)
    p_sub <- paste0("R = ", round(m_info$R_Value, 2), " | p = ", format_p(m_info$P_Value))
  }
  
  p <- ggplot(head(res_df, 15), aes(x = Fold_Enrichment, y = reorder(term_name, Fold_Enrichment))) +
    geom_point(aes(size = intersection_size, color = p_value)) +
    scale_color_gradient(low = "red", high = "blue") + theme_bw() +
    labs(title = p_title, subtitle = p_sub, x = "Fold Enrichment", y = "Pathway")
  
  ggsave(filename = paste0("Pathway_", group_name, ".jpg"), plot = p, path = fig_path, width = 9, height = 7, dpi = 600)
}

write.xlsx(pathway_excel_list, file.path(output_dir, "MASTER_PATHWAY_RESULTS.xlsx"))

# 7. VALIDATION BOXPLOTS (p < 0.05) ----------------------------------------
print("Step 5: Generating Boxplots for all Significant Modules...")

for(i in seq_len(nrow(results_summary))) {
  mod_name <- results_summary$ME_Name[i]
  boxplot_data <- data.frame(Relatedness = factor(mice_info$Relatedness, levels = c("0", "50")),
                             Expression = MEs[[mod_name]])
  
  p <- ggplot(boxplot_data, aes(x = Relatedness, y = Expression, fill = Relatedness)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) + geom_jitter(width = 0.2, size = 2, alpha = 0.5) +
    theme_pubr() + scale_fill_manual(values = c("0" = "#66c2a5", "50" = "#fc8d62")) +
    labs(title = paste("Module Validation:", mod_name), 
         subtitle = paste0("R = ", round(results_summary$R_Value[i], 2), " | p = ", format_p(results_summary$P_Value[i])),
         y = "Adjusted Eigengene", x = "Relatedness (F1=0 vs F2=50)") + 
    stat_compare_means(method = "t.test", label = "p.signif")
  
  ggsave(filename = paste0("p_boxplot_", mod_name, ".jpg"), plot = p, path = output_dir, width = 6, height = 5, dpi = 600)
}

print("Workflow Complete. Dendrogram, Heatmap, Pathways, and Boxplots generated.")

###################

# 9. MODULE INTERNAL CONNECTIVITY (Internal Network Plots) ----------------

print("Step 7: Calculating Internal Module Networks...")

# Create a folder for the network data
net_dir <- file.path(output_dir, "Module_Networks")
if(!dir.exists(net_dir)) dir.create(net_dir)

# We must recalculate the TOM (Topological Overlap Matrix) for the network
# This takes the 'net' object from your previous step
TOM = TOMsimilarityFromExpr(datExpr, power = softPower)

for(i in seq_len(nrow(results_summary))) {
  
  mod_color <- results_summary$Color[i]
  mod_name <- results_summary$ME_Name[i]
  
  # Select genes belonging to this module
  probes = colnames(datExpr)
  inModule = (moduleColors == mod_color)
  modProbes = probes[inModule]
  modTOM = TOM[inModule, inModule]
  
  dimnames(modTOM) = list(modProbes, modProbes)
  
  # 1. Export for Cytoscape (The PhD standard for these graphs)
  # This creates two files: an 'Edges' file (the lines) and a 'Nodes' file (the dots)
  exportNetworkToCytoscape(modTOM,
                           edgeFile = file.path(net_dir, paste0("Edges-", mod_name, ".txt")),
                           nodeFile = file.path(net_dir, paste0("Nodes-", mod_name, ".txt")),
                           weighted = TRUE,
                           threshold = 0.1, # Adjust this higher (0.2) if the graph is too crowded
                           nodeNames = modProbes,
                           altNodeNames = modProbes,
                           nodeAttr = moduleColors[inModule])
  
  # 2. Quick Visual in R (Circular Plot)
  # We limit to the top 30 most connected genes so the graph is readable
  nTop = 30
  IMConn = softConnectivity(datExpr[, modProbes])
  top = (rank(-IMConn) <= nTop)
  
  jpeg(file.path(net_dir, paste0("p_internal_network_", mod_name, ".jpg")), 
       width = 7, height = 7, units = "in", res = 600)
  
  visLabels = modProbes[top]
  # Simplify LOC names for the graph display
  visLabels = gsub("LOC", "L-", visLabels) 
  
  plotNetworkHeatmap(datExpr[, modProbes[top]],
                     plotGenes = visLabels,
                     networkType = "unsigned",
                     
                     main = paste("Top 30 Hubs in", mod_name))
  dev.off()
}

print("Internal networks exported. Check 'Module_Networks' folder.")


################################
# 6. INTERNAL MODULE NETWORK (The "Circle and Lines" Graph) --------------
print("Step 6: Generating Internal Network Graphs in R...")
net_fig_path <- file.path(output_dir, "Module_Internal_Networks")
if(!dir.exists(net_fig_path)) dir.create(net_fig_path)

# Recalculate TOM specifically for network visualization
load(net$TOMFiles[1]) # Load the TOM calculated during blockwiseModules
TOM <- as.matrix(TOM)

for(i in seq_len(nrow(results_summary))) {
  mod_color <- results_summary$Color[i]; mod_name <- results_summary$ME_Name[i]
  inModule = (moduleColors == mod_color); modProbes = colnames(datExpr)[inModule]
  modTOM = TOM[inModule, inModule]; dimnames(modTOM) = list(modProbes, modProbes)
  
  # Select Top 30 Hubs to keep the graph readable
  nTop = 30; IMConn = softConnectivity(datExpr[, modProbes])
  top = (rank(-IMConn) <= nTop); visTOM = modTOM[top, top]
  
  # Create Igraph Object
  g <- graph_from_adjacency_matrix(visTOM, mode = "undirected", weighted = TRUE, diag = FALSE)
  # Remove weak edges for a cleaner visual
  g <- delete_edges(g, E(g)[weight < 0.15]) 
  
  jpeg(file.path(net_fig_path, paste0("p_network_", mod_name, ".jpg")), width = 8, height = 8, units = "in", res = 600)
  plot(g, layout = layout_with_fr(g), vertex.color = mod_color, 
       vertex.label.cex = 0.7, vertex.label.color = "black", vertex.frame.color = "white",
       edge.width = E(g)$weight * 5, edge.color = alpha("grey", 0.5),
       main = paste("Internal Network Hubs:", mod_name))
  dev.off()
}

