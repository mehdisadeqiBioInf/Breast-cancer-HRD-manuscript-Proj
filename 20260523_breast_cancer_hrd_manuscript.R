##### libraries #####
library(DESeq2)
library(MatrixGenerics)
library(WGCNA)
library(glmnet)
library(ranger)
library(caret)
library(pROC)
library(plotROC)
library(ggplot2)
library(cowplot)
library(ggplotify)
library(pheatmap)
library(survival)
library(broom)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)
library(reactome.db)
library(GO.db)
library(KEGGREST)
library(msigdbr)
library(limma)
library(readxl)
library(gridExtra)
library(janitor)

##### parameters #####
set.seed(123)
K_outer <- 5
K_inner <- 5
B_boot <- 500
rf_top_n <- 30
inner_freq_thresh <- 0.5
min_overlap_genes <- 4
freq_thresh_boot <- 0.6
primary_threshold <- 0.5
hrd_cutoffs <- c(42, 33)
tcga_count_directory <- "expression files"
tcga_hrd_file <- "TCGA BRCA with HRD scores.tsv"
tcga_metadata_file <- "data/tcga_sample_metadata_prepared.tsv"
tcga_cdr_file <- "data/tcga_clinical/TCGA-CDR-SupplementalTableS1.xlsx"
tcga_brca_status_file <- "data/tcga_brca_knijnenburg_brca_status.csv"
metabric_expression_file <- "data/data_mrna_illumina_microarray.txt"
metabric_hrd_file <- "data/METABRIC_HRD_scores.tsv"
metabric_metadata_file <- "data/metabric_sample_metadata_prepared.tsv"
metabric_clinical_file <- "data/metabric_clinical/brca_metabric_clinical_data.tsv"
celltype_fraction_file <- "data/celltype_fraction_matrix.tsv"
crispr_square_file <- "data/Breast cell lines-Olaparib treated/Breast cell lines-Olaparib treated/MAGeCKFlute_sump149pt_Day15.MLE/Selectionsump149pt_Day15.MLEsquareview_data.txt"
crispr_mle_file <- "data/Breast cell lines-Olaparib treated/Breast cell lines-Olaparib treated/sump149pt all samples MLE ready.txt"
jacobson_template_file <- "data/published_signatures/jacobson_2023_signature_templates.xlsx"
pan_hrd200_scores_file <- "data/published_signatures/METABRIC_Pan_HRD200_public_classifier_scores.tsv"
palette_main <- c("#868686FF", "#A73030FF", "#0073C2FF", "#EFC000FF", "#7AA6DCFF", "#8F7700FF", "#3B3B3BFF", "#7876B1FF", "#20854EFF", "#E18727FF")
pam50_cols <- c("LumA" = "#0073C2FF", "LumB" = "#EFC000FF", "Her2" = "#A73030FF", "Basal" = "#20854EFF", "Normal-like" = "#868686FF")

##### TCGA expression and WGCNA #####
datTraits <- read.table(tcga_hrd_file, stringsAsFactors = FALSE, header = TRUE, check.names = FALSE)
sampleFiles <- list.files(tcga_count_directory)
sampleFiles <- sampleFiles[sampleFiles %in% datTraits[, 1]]
sampleTable <- data.frame(sampleName = sampleFiles, fileName = sampleFiles, condition = factor("HRD"), stringsAsFactors = FALSE)
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = tcga_count_directory, design = ~ 1)
dds <- DESeq(dds)
vsd <- vst(dds)
counts1 <- counts(dds, normalized = FALSE)
norm.mat <- data.frame(assay(vsd), check.names = FALSE)
rv <- MatrixGenerics::rowVars(as.matrix(norm.mat))
norm.mat <- norm.mat[rv > as.numeric(quantile(rv, 0.75, na.rm = TRUE)), ]
datExpr <- data.frame(t(norm.mat), check.names = FALSE)
rownames(datTraits) <- datTraits[, 1]
datTraits <- datTraits[, -1, drop = FALSE]
datTraits <- datTraits[match(rownames(datExpr), rownames(datTraits)), , drop = FALSE]
colnames(datTraits)[1] <- "HRD_score"
y1 <- as.numeric(datTraits$HRD_score >= 42)
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
powers <- c(1:10, seq(12, 30, 2))
options(expressions = 5e5)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0, blockSize = 25000, networkType = "signed")
net <- blockwiseModules(datExpr, power = 14, TOMType = "signed", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = FALSE, pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 0, maxBlockSize = 25000)
MEs <- net$MEs
module.members <- data.frame(module = net$colors, stringsAsFactors = FALSE)
module.members$symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = substr(rownames(module.members), 1, 15), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
module.members$ENTREZ <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = substr(rownames(module.members), 1, 15), column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
moduleTraitCor <- as.data.frame(cor(MEs, datTraits, use = "p"))
moduleTraitPvalue <- as.data.frame(corPvalueStudent(as.matrix(moduleTraitCor), nSamples))
rownames(moduleTraitPvalue) <- rownames(moduleTraitCor)
module_order <- order(abs(moduleTraitCor[, 1]), decreasing = TRUE)
module_name_map <- setNames(paste0("Module ", seq_along(module_order)), rownames(moduleTraitCor)[module_order])
module.members$module_label <- module_name_map[paste0("ME", module.members$module)]
brown.mod <- module.members[module.members$module == "brown", , drop = FALSE]
candidate_genes <- substr(rownames(brown.mod), 1, 15)
rownames(counts1) <- substr(rownames(counts1), 1, 15)
counts1 <- counts1[!duplicated(rownames(counts1)), , drop = FALSE]
counts1 <- counts1[, rownames(datTraits), drop = FALSE]
genes_norm_train <- rownames(counts1)
candidate_genes_model <- intersect(candidate_genes, genes_norm_train)
supplementary_table_s1 <- split(module.members$symbol, module.members$module_label)
supplementary_table_s1 <- as.data.frame(lapply(supplementary_table_s1[paste0("Module ", seq_along(supplementary_table_s1))], "length<-", max(lengths(supplementary_table_s1))), stringsAsFactors = FALSE)

##### pathway annotations #####
reactome_t2g_raw <- AnnotationDbi::toTable(reactome.db::reactomeEXTID2PATHID)
reactome_t2g <- data.frame(source_id = as.character(reactome_t2g_raw[[2]]), gene = as.character(reactome_t2g_raw[[1]]), stringsAsFactors = FALSE)
reactome_t2g <- reactome_t2g[!is.na(reactome_t2g$source_id) & !is.na(reactome_t2g$gene), ]
reactome_t2n_raw <- AnnotationDbi::toTable(reactome.db::reactomePATHID2NAME)
reactome_meta <- data.frame(source_id = as.character(reactome_t2n_raw[[1]]), Database = "Reactome", name = as.character(reactome_t2n_raw[[2]]), stringsAsFactors = FALSE)
reactome_meta$term <- paste0("Reactome::", reactome_meta$source_id)
reactome_t2g$term <- paste0("Reactome::", reactome_t2g$source_id)
reactome_t2g <- reactome_t2g[, c("term", "gene")]
reactome_meta <- reactome_meta[, c("term", "source_id", "Database", "name")]
kegg_link <- KEGGREST::keggLink("hsa", "pathway")
kegg_t2g <- data.frame(source_id = sub("^path:", "", as.character(kegg_link)), gene = sub("^hsa:", "", names(kegg_link)), stringsAsFactors = FALSE)
kegg_names <- KEGGREST::keggList("pathway", "hsa")
kegg_meta <- data.frame(source_id = sub("^path:", "", names(kegg_names)), Database = "KEGG", name = sub(" - Homo sapiens \\(human\\)$", "", as.character(kegg_names)), stringsAsFactors = FALSE)
kegg_meta$term <- paste0("KEGG::", kegg_meta$source_id)
kegg_t2g$term <- paste0("KEGG::", kegg_t2g$source_id)
kegg_t2g <- kegg_t2g[, c("term", "gene")]
kegg_meta <- kegg_meta[, c("term", "source_id", "Database", "name")]
go_keys <- keys(org.Hs.eg.db, keytype = "ENTREZID")
go_map <- AnnotationDbi::select(org.Hs.eg.db, keys = go_keys, keytype = "ENTREZID", columns = c("ENTREZID", "GO", "ONTOLOGY"))
go_map <- go_map[!is.na(go_map$GO) & go_map$ONTOLOGY == "BP", ]
go_t2g <- unique(data.frame(term = paste0("GO_BP::", go_map$GO), gene = as.character(go_map$ENTREZID), stringsAsFactors = FALSE))
go_meta_raw <- AnnotationDbi::select(GO.db, keys = unique(go_map$GO), keytype = "GOID", columns = c("GOID", "TERM", "ONTOLOGY"))
go_meta <- unique(data.frame(term = paste0("GO_BP::", go_meta_raw$GOID), source_id = go_meta_raw$GOID, Database = "GO BP", name = go_meta_raw$TERM, stringsAsFactors = FALSE))
hallmark_t2g_raw <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
hallmark_t2g <- unique(data.frame(term = paste0("Hallmark::", hallmark_t2g_raw$gs_name), gene = as.character(hallmark_t2g_raw$entrez_gene), stringsAsFactors = FALSE))
hallmark_meta <- unique(data.frame(term = paste0("Hallmark::", hallmark_t2g_raw$gs_name), source_id = hallmark_t2g_raw$gs_name, Database = "MSigDB Hallmark", name = gsub("_", " ", sub("^HALLMARK_", "", hallmark_t2g_raw$gs_name)), stringsAsFactors = FALSE))
pooled_TERM2GENE <- unique(rbind(reactome_t2g, kegg_t2g, go_t2g, hallmark_t2g))
term_metadata <- unique(rbind(reactome_meta, kegg_meta, go_meta, hallmark_meta))
pooled_TERM2NAME <- unique(term_metadata[, c("term", "name")])

##### module enrichment #####
uni <- unique(na.omit(as.character(module.members$ENTREZ)))
sig.genes <- unique(na.omit(as.character(brown.mod$ENTREZ)))
pooled_TERM2GENE_wgcna <- pooled_TERM2GENE[pooled_TERM2GENE$gene %in% uni, ]
pooled_ora_module1 <- clusterProfiler::enricher(gene = sig.genes, universe = uni, TERM2GENE = pooled_TERM2GENE_wgcna, TERM2NAME = pooled_TERM2NAME, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
pooled_ora_module1_all <- as.data.frame(pooled_ora_module1)
pooled_ora_module1_all <- merge(pooled_ora_module1_all, term_metadata, by.x = "ID", by.y = "term", all.x = TRUE, sort = FALSE)
pooled_ora_module1_all$minus_log10_FDR <- -log10(pooled_ora_module1_all$p.adjust)
pooled_ora_module1_all <- pooled_ora_module1_all[order(pooled_ora_module1_all$p.adjust), ]
pooled_ora_module1_sig <- pooled_ora_module1_all[pooled_ora_module1_all$p.adjust <= 0.05, ]
module1_plot_terms <- c("Cell Cycle", "G2M Checkpoint", "DNA Repair", "DNA Damage Response", "Double-Strand Break Repair", "RHO GTPase Effectors", "MYC TARGETS V1", "Homologous Recombination")
module1_enrichment <- data.frame()
for (i in seq_along(module1_plot_terms)) {
  hit <- grep(module1_plot_terms[i], pooled_ora_module1_sig$Description, ignore.case = TRUE)
  if (length(hit) == 0) hit <- grep(module1_plot_terms[i], pooled_ora_module1_sig$name, ignore.case = TRUE)
  if (length(hit) > 0) module1_enrichment <- rbind(module1_enrichment, pooled_ora_module1_sig[hit[1], ])
}
module1_enrichment$plot_label <- module1_plot_terms[seq_len(nrow(module1_enrichment))]
module1_enrichment$Value <- -log10(module1_enrichment$p.adjust)

##### figure 1 objects #####
pca_res <- prcomp(datExpr, center = TRUE, scale. = FALSE)
pca_df <- data.frame(sample_id = rownames(datExpr), HRD_score = as.numeric(datTraits$HRD_score), PC1 = pca_res$x[, 1], PC2 = pca_res$x[, 2], stringsAsFactors = FALSE)
pc1_hrd_cor <- cor.test(pca_df$HRD_score, pca_df$PC1, method = "pearson")
figure1_panel_A <- ggplot(pca_df, aes(x = HRD_score, y = PC1)) + geom_point(color = "grey25", size = 0.8, alpha = 0.75) + geom_smooth(method = "lm", se = FALSE, color = "firebrick2", linewidth = 0.5) + annotate("text", x = quantile(pca_df$HRD_score, 0.82, na.rm = TRUE), y = quantile(pca_df$PC1, 0.12, na.rm = TRUE), label = paste0("R = ", round(unname(pc1_hrd_cor$estimate), 2)), fontface = "bold", size = 3.2) + theme_classic(base_size = 8) + labs(x = "HRD score", y = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"))
module_bar_df <- data.frame(module = module_name_map[rownames(moduleTraitCor)[module_order]], correlation = abs(moduleTraitCor[module_order, 1]), color = gsub("^ME", "", rownames(moduleTraitCor)[module_order]), stringsAsFactors = FALSE)
module_bar_df$module <- factor(module_bar_df$module, levels = module_bar_df$module)
figure1_panel_B <- ggplot(module_bar_df, aes(x = module, y = correlation, fill = module)) + geom_col(color = "black", linewidth = 0.25) + scale_fill_manual(values = c("#0073C2FF", "#EFC000FF", "#CD534CFF", "#8F7700FF", "#008080", "#A73030FF", "#868686FF", "#E18727FF", "#20854EFF", "#7876B1FF", "#003C67FF", "black")) + theme_classic(base_size = 8) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), axis.text.y = element_text(colour = "black")) + labs(x = NULL, y = "Absolute Correlation Coefficient")
figure1_panel_C <- ggplot(module1_enrichment, aes(x = reorder(plot_label, Value), y = Value)) + geom_col(fill = "#A73030FF", color = "black", linewidth = 0.25) + coord_flip() + theme_classic(base_size = 7) + theme(panel.border = element_rect(fill = NA, colour = "black"), axis.text = element_text(colour = "black")) + labs(title = "Module 1 enrichment results", x = NULL, y = expression(-log[10] * " " * italic(q)))
celltype_fraction <- read.table(celltype_fraction_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
celltype_fraction <- celltype_fraction[rownames(MEs), , drop = FALSE]
module_celltype_cor <- cor(MEs[, rownames(moduleTraitCor)[module_order], drop = FALSE], celltype_fraction, use = "pairwise.complete.obs")
rownames(module_celltype_cor) <- module_name_map[rownames(moduleTraitCor)[module_order]]
figure1_panel_D <- pheatmap::pheatmap(module_celltype_cor, cluster_rows = FALSE, cluster_cols = TRUE, display_numbers = TRUE, breaks = seq(-0.6, 0.6, length.out = 51), main = "", border_color = "black", silent = TRUE)
figure1 <- cowplot::plot_grid(figure1_panel_A, figure1_panel_B, figure1_panel_C, ggplotify::as.ggplot(figure1_panel_D$gtable), labels = c("A", "B", "C", "D"), ncol = 2, label_size = 16)

##### nested cross-validation #####
n_samples <- length(y1)
set.seed(123)
outer_folds <- caret::createFolds(y = factor(y1), k = K_outer, returnTrain = FALSE)
outer_results <- data.frame(fold = integer(0), auc = numeric(0), acc = numeric(0), sens = numeric(0), spec = numeric(0), stringsAsFactors = FALSE)
outer_gene_sets <- vector("list", K_outer)
outer_pred <- numeric(n_samples)
outer_true <- y1
for (k in seq_len(K_outer)) {
  test_idx <- outer_folds[[k]]
  train_idx <- setdiff(seq_len(n_samples), test_idx)
  counts_train <- counts1[, train_idx, drop = FALSE]
  counts_test <- counts1[, test_idx, drop = FALSE]
  y_outer_train <- y1[train_idx]
  y_outer_test <- y1[test_idx]
  coldata_train <- data.frame(row.names = colnames(counts_train), condition = factor(y_outer_train))
  dds_train <- DESeqDataSetFromMatrix(countData = counts_train, colData = coldata_train, design = ~ 1)
  geo_means <- exp(MatrixGenerics::rowMeans(log(counts(dds_train))))
  dds_train <- estimateSizeFactors(dds_train, geoMeans = geo_means)
  dds_train <- estimateDispersions(dds_train)
  vsd_train <- vst(dds_train, blind = TRUE)
  X_outer_train_all <- t(assay(vsd_train))
  disp_fun_outer <- dispersionFunction(dds_train)
  disp_outer <- dispersions(dds_train)
  names(disp_outer) <- rownames(dds_train)
  coldata_test <- data.frame(row.names = colnames(counts_test), condition = factor(y_outer_test))
  dds_test <- DESeqDataSetFromMatrix(countData = counts_test, colData = coldata_test, design = ~ 1)
  dds_test <- estimateSizeFactors(dds_test, geoMeans = geo_means[rownames(dds_test)])
  dispersions(dds_test) <- disp_outer[rownames(dds_test)]
  dispersionFunction(dds_test) <- disp_fun_outer
  vsd_test <- vst(dds_test, blind = TRUE)
  X_outer_test_all <- t(assay(vsd_test))
  X_outer_train <- X_outer_train_all[, candidate_genes_model, drop = FALSE]
  X_outer_test <- X_outer_test_all[, candidate_genes_model, drop = FALSE]
  set.seed(100 + k)
  inner_folds <- caret::createFolds(y = factor(y_outer_train), k = K_inner, returnTrain = FALSE)
  inner_gene_sets <- vector("list", K_inner)
  for (j in seq_len(K_inner)) {
    val_idx <- inner_folds[[j]]
    tr_idx_in <- setdiff(seq_len(length(y_outer_train)), val_idx)
    X_in_train <- X_outer_train[tr_idx_in, , drop = FALSE]
    y_in_train <- y_outer_train[tr_idx_in]
    rf_data <- data.frame(y = factor(y_in_train), X_in_train, check.names = FALSE)
    rf_fit <- ranger::ranger(y ~ ., data = rf_data, importance = "permutation", num.trees = 10000, probability = FALSE)
    rf_imp <- rf_fit$variable.importance
    rf_genes <- names(sort(rf_imp, decreasing = TRUE))[seq_len(min(rf_top_n, length(rf_imp)))]
    X_in_train_scaled <- scale(X_in_train)
    cv_lasso <- glmnet::cv.glmnet(x = as.matrix(X_in_train_scaled), y = y_in_train, family = "binomial", alpha = 1)
    coef_lasso <- coef(cv_lasso, s = "lambda.1se")
    lasso_genes <- setdiff(rownames(coef_lasso)[which(as.vector(coef_lasso) != 0)], "(Intercept)")
    inner_gene_sets[[j]] <- intersect(rf_genes, lasso_genes)
  }
  gene_counts <- table(unlist(inner_gene_sets))
  gene_freq <- gene_counts / length(inner_gene_sets)
  selected_genes <- names(gene_freq)[gene_freq >= inner_freq_thresh]
  if (length(selected_genes) < min_overlap_genes) stop("Few genes were selected in the inner loop.")
  outer_gene_sets[[k]] <- selected_genes
  X_train_sel <- X_outer_train[, selected_genes, drop = FALSE]
  X_test_sel <- X_outer_test[, selected_genes, drop = FALSE]
  X_train_mu <- colMeans(X_train_sel)
  X_train_sd <- apply(X_train_sel, 2, sd)
  X_train_sd[X_train_sd == 0] <- 1
  X_train_scaled <- scale(X_train_sel, center = X_train_mu, scale = X_train_sd)
  X_test_scaled <- scale(X_test_sel, center = X_train_mu, scale = X_train_sd)
  cv_ridge <- glmnet::cv.glmnet(x = as.matrix(X_train_scaled), y = y_outer_train, family = "binomial", alpha = 0)
  prob_test <- as.numeric(predict(cv_ridge, newx = as.matrix(X_test_scaled), s = cv_ridge$lambda.1se, type = "response"))
  outer_pred[test_idx] <- prob_test
  roc_tmp <- pROC::roc(response = y_outer_test, predictor = prob_test, quiet = TRUE)
  pred_class <- as.numeric(prob_test >= 0.5)
  tab <- table(factor(pred_class, levels = c(0, 1)), factor(y_outer_test, levels = c(0, 1)))
  outer_results <- rbind(outer_results, data.frame(fold = k, auc = as.numeric(pROC::auc(roc_tmp)), acc = sum(diag(tab)) / sum(tab), sens = tab[2, 2] / sum(tab[, 2]), spec = tab[1, 1] / sum(tab[, 1])))
}
outer_gene_counts <- table(unlist(outer_gene_sets))
outer_gene_freq <- outer_gene_counts / length(outer_gene_sets)
outer_selected_genes <- names(outer_gene_freq)[outer_gene_freq >= inner_freq_thresh]
tcga_roc <- pROC::roc(response = outer_true, predictor = outer_pred, quiet = TRUE)
tcga_roc_ci <- as.numeric(pROC::ci.auc(tcga_roc))
thr_youden <- as.numeric(pROC::coords(tcga_roc, x = "best", best.method = "youden", ret = "threshold", transpose = FALSE)[1, 1])
roc_dat_log <- data.frame(sample_id = rownames(datTraits), patient_id = substr(rownames(datTraits), 1, 12), status = outer_true, pred = outer_pred, HRD_score = as.numeric(datTraits$HRD_score), stringsAsFactors = FALSE)
roc_dat_log$score_oof <- roc_dat_log$pred
roc_dat_log$score_oof_z <- as.numeric(scale(qlogis(pmin(pmax(roc_dat_log$score_oof, 1e-6), 1 - 1e-6))))

##### final model #####
coldata1_full <- data.frame(row.names = colnames(counts1), condition = factor(y1))
dds1_full <- DESeqDataSetFromMatrix(countData = counts1, colData = coldata1_full, design = ~ 1)
geo_means_full <- exp(MatrixGenerics::rowMeans(log(counts(dds1_full))))
dds1_full <- estimateSizeFactors(dds1_full, geoMeans = geo_means_full)
dds1_full <- estimateDispersions(dds1_full)
vsd1_full <- vst(dds1_full, blind = TRUE)
X1_full_all <- t(assay(vsd1_full))
X1_full <- X1_full_all[, candidate_genes_model, drop = FALSE]
disp_fun_full <- dispersionFunction(dds1_full)
disp_full <- dispersions(dds1_full)
names(disp_full) <- rownames(dds1_full)
set.seed(456)
genes_cand <- colnames(X1_full)
gene_select_counts <- setNames(numeric(length(genes_cand)), genes_cand)
for (b in seq_len(B_boot)) {
  idx_b <- sample(seq_len(nrow(X1_full)), size = nrow(X1_full), replace = TRUE)
  X_b <- X1_full[idx_b, , drop = FALSE]
  y_b <- y1[idx_b]
  rf_data_b <- data.frame(y = factor(y_b), X_b, check.names = FALSE)
  rf_fit_b <- ranger::ranger(y ~ ., data = rf_data_b, importance = "permutation", num.trees = 10000, probability = FALSE)
  rf_imp_b <- rf_fit_b$variable.importance
  rf_genes_b <- names(sort(rf_imp_b, decreasing = TRUE))[seq_len(min(rf_top_n, length(rf_imp_b)))]
  X_b_scaled <- scale(X_b)
  cv_lasso_b <- glmnet::cv.glmnet(x = as.matrix(X_b_scaled), y = y_b, family = "binomial", alpha = 1)
  coef_lasso_b <- coef(cv_lasso_b, s = "lambda.1se")
  lasso_genes_b <- setdiff(rownames(coef_lasso_b)[which(as.vector(coef_lasso_b) != 0)], "(Intercept)")
  genes_b <- intersect(rf_genes_b, lasso_genes_b)
  if (length(genes_b) < min_overlap_genes) stop("Few genes were selected in one bootstrap resample.")
  gene_select_counts[genes_b] <- gene_select_counts[genes_b] + 1
}
gene_select_freq <- gene_select_counts / B_boot
stable_genes_df <- data.frame(gene = names(gene_select_freq), freq = as.numeric(gene_select_freq), stringsAsFactors = FALSE)
stable_genes_df <- stable_genes_df[order(-stable_genes_df$freq), ]
stable_genes <- intersect(stable_genes_df$gene[stable_genes_df$freq >= freq_thresh_boot], colnames(X1_full))
X1_final <- X1_full[, stable_genes, drop = FALSE]
mu_final <- colMeans(X1_final)
sd_final <- apply(X1_final, 2, sd)
sd_final[sd_final == 0] <- 1
X1_final_scaled <- scale(X1_final, center = mu_final, scale = sd_final)
set.seed(999)
cv_ridge_final <- glmnet::cv.glmnet(x = as.matrix(X1_final_scaled), y = y1, family = "binomial", alpha = 0)
lambda_final <- cv_ridge_final$lambda.1se
prob1 <- as.numeric(predict(cv_ridge_final, newx = as.matrix(X1_final_scaled), s = lambda_final, type = "response"))
coef_mat <- coef(cv_ridge_final, s = lambda_final)
model_coefficients <- data.frame(gene_id = rownames(coef_mat), beta = as.numeric(coef_mat[, 1]), stringsAsFactors = FALSE)
model_obj_boot <- list(stable_genes = stable_genes, mu_train = mu_final, sd_train = sd_final, glmnet_cv = cv_ridge_final, lambda_final = lambda_final, genes_norm = genes_norm_train, cand_genes = candidate_genes_model, geo_means_full = geo_means_full, disp_fun_full = disp_fun_full, disp_full = disp_full, threshold = thr_youden)
tcga_model_scores <- data.frame(sample_id = rownames(datTraits), patient_id = substr(rownames(datTraits), 1, 12), HRD_score = as.numeric(datTraits$HRD_score), HRD_class = y1, score_oof = outer_pred, score_final = prob1, stringsAsFactors = FALSE)
tcga_model_scores$score_oof_z <- as.numeric(scale(qlogis(pmin(pmax(tcga_model_scores$score_oof, 1e-6), 1 - 1e-6))))
tcga_model_scores$score_final_z <- as.numeric(scale(qlogis(pmin(pmax(tcga_model_scores$score_final, 1e-6), 1 - 1e-6))))
tcga_model_scores$predicted_HRD_oof <- as.numeric(tcga_model_scores$score_oof >= primary_threshold)
tcga_model_scores$predicted_HRD_final <- as.numeric(tcga_model_scores$score_final >= primary_threshold)

##### METABRIC validation #####
expr_mb <- read.delim(metabric_expression_file, check.names = FALSE, stringsAsFactors = FALSE)
expr_mb$ENSEMBL_ID <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = as.character(expr_mb$Entrez_Gene_Id), column = "ENSEMBL", keytype = "ENTREZID", multiVals = "first")
expr_mb_model <- expr_mb[!is.na(expr_mb$ENSEMBL_ID) & expr_mb$ENSEMBL_ID %in% stable_genes, ]
rownames(expr_mb_model) <- expr_mb_model$ENSEMBL_ID
expr_mb_model <- expr_mb_model[stable_genes, setdiff(colnames(expr_mb_model), c("Hugo_Symbol", "Entrez_Gene_Id", "ENSEMBL_ID")), drop = FALSE]
expr_mb_model <- t(expr_mb_model)
expr_mb_model_scaled <- scale(expr_mb_model, center = mu_final, scale = sd_final)
prob_mb <- as.numeric(predict(cv_ridge_final, newx = as.matrix(expr_mb_model_scaled), s = lambda_final, type = "response"))
metabric_hrd <- read.table(metabric_hrd_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
metabric_hrd$sample_id <- as.character(metabric_hrd$sample_id)
metabric_hrd$HRD_score <- as.numeric(metabric_hrd$HRD_score)
metabric_hrd$HRD_class <- as.numeric(metabric_hrd$HRD_score >= 42)
metabric_model_scores <- data.frame(sample_id = rownames(expr_mb_model), score_final = prob_mb, stringsAsFactors = FALSE)
metabric_model_scores <- merge(metabric_hrd, metabric_model_scores, by = "sample_id")
metabric_model_scores$score_final_z <- as.numeric(scale(qlogis(pmin(pmax(metabric_model_scores$score_final, 1e-6), 1 - 1e-6))))
metabric_model_scores$predicted_HRD_final <- as.numeric(metabric_model_scores$score_final >= primary_threshold)
roc_dat_mb <- data.frame(sample_id = metabric_model_scores$sample_id, status = metabric_model_scores$HRD_class, pred = metabric_model_scores$score_final, HRD_score = metabric_model_scores$HRD_score, stringsAsFactors = FALSE)
metabric_roc <- pROC::roc(response = roc_dat_mb$status, predictor = roc_dat_mb$pred, quiet = TRUE)
metabric_roc_ci <- as.numeric(pROC::ci.auc(metabric_roc))

##### figure 2 ROC objects #####
figure2_panel_A <- pROC::ggroc(tcga_roc, legacy.axes = TRUE, linewidth = 0.8, color = "steelblue") + geom_abline(intercept = 0, slope = 1, linetype = 2) + annotate("text", x = 0.72, y = 0.28, label = paste0("AUC = ", sprintf("%.2f", as.numeric(pROC::auc(tcga_roc))), "\nCI: ", sprintf("%.2f", tcga_roc_ci[1]), "-", sprintf("%.2f", tcga_roc_ci[3])), size = 3.2) + theme_classic(base_size = 9) + labs(title = "TCGA 5-fold cross-validation\n(stacked test predictions)", x = "1 - specificity", y = "True positive fraction")
figure2_panel_B <- pROC::ggroc(metabric_roc, legacy.axes = TRUE, linewidth = 0.8, color = "#A73030FF") + geom_abline(intercept = 0, slope = 1, linetype = 2) + annotate("text", x = 0.72, y = 0.28, label = paste0("AUC = ", sprintf("%.2f", as.numeric(pROC::auc(metabric_roc))), "\nCI: ", sprintf("%.2f", metabric_roc_ci[1]), "-", sprintf("%.2f", metabric_roc_ci[3])), size = 3.2) + theme_classic(base_size = 9) + labs(title = "METABRIC external validation", x = "1 - specificity", y = "True positive fraction")

##### survival analyses #####
tcga_cdr <- as.data.frame(readxl::read_excel(tcga_cdr_file))
colnames(tcga_cdr) <- janitor::make_clean_names(colnames(tcga_cdr))
tcga_cdr <- tcga_cdr[tcga_cdr$type == "BRCA", ]
tcga_clin <- data.frame(patient_id = tcga_cdr$bcr_patient_barcode, os_time = tcga_cdr$os_time / 30.44, os_event = tcga_cdr$os, dfi_time = tcga_cdr$dfi_time / 30.44, dfi_event = tcga_cdr$dfi, pfi_time = tcga_cdr$pfi_time / 30.44, pfi_event = tcga_cdr$pfi, dss_time = tcga_cdr$dss_time / 30.44, dss_event = tcga_cdr$dss, age = tcga_cdr$age_at_initial_pathologic_diagnosis, stage = tcga_cdr$ajcc_pathologic_tumor_stage, stringsAsFactors = FALSE)
tcga_clin$stage_clean <- toupper(trimws(tcga_clin$stage))
tcga_clin$stage_simple <- NA
tcga_clin$stage_simple[grepl("^STAGE I", tcga_clin$stage_clean)] <- "Stage I"
tcga_clin$stage_simple[grepl("^STAGE II", tcga_clin$stage_clean)] <- "Stage II"
tcga_clin$stage_simple[grepl("^STAGE III", tcga_clin$stage_clean)] <- "Stage III"
tcga_clin$stage_simple[grepl("^STAGE IV", tcga_clin$stage_clean)] <- "Stage IV"
tcga_clin$stage_simple <- factor(tcga_clin$stage_simple, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
tcga_metadata <- read.table(tcga_metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
tcga_surv_analysis <- merge(tcga_clin, unique(tcga_metadata[, c("patient_id", "PAM50")]), by = "patient_id", all.x = TRUE)
colnames(tcga_surv_analysis)[colnames(tcga_surv_analysis) == "PAM50"] <- "subtype"
tcga_surv_analysis <- merge(tcga_surv_analysis, tcga_model_scores, by = "patient_id")
tcga_surv_analysis$subtype <- factor(tcga_surv_analysis$subtype)
tcga_cox_models <- list()
tcga_cox_results <- data.frame()
tcga_survival_endpoints <- c("os", "dss", "pfi", "dfi")
for (endpoint in tcga_survival_endpoints) {
  dat <- tcga_surv_analysis[!is.na(tcga_surv_analysis[[paste0(endpoint, "_time")]]) & !is.na(tcga_surv_analysis[[paste0(endpoint, "_event")]]) & tcga_surv_analysis[[paste0(endpoint, "_time")]] > 0 & !is.na(tcga_surv_analysis$score_oof_z) & !is.na(tcga_surv_analysis$age) & !is.na(tcga_surv_analysis$stage_simple) & !is.na(tcga_surv_analysis$subtype), ]
  dat$stage_model <- droplevels(factor(dat$stage_simple))
  dat$subtype <- droplevels(factor(dat$subtype))
  tcga_cox_models[[endpoint]] <- survival::coxph(as.formula(paste0("survival::Surv(", endpoint, "_time, ", endpoint, "_event) ~ score_oof_z + age + stage_model + strata(subtype)")), data = dat, x = TRUE, y = TRUE)
  tcga_cox_results <- rbind(tcga_cox_results, data.frame(endpoint = endpoint, broom::tidy(tcga_cox_models[[endpoint]], exponentiate = TRUE, conf.int = TRUE), row.names = NULL))
}
metabric_clin <- read.table(metabric_clinical_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
metabric_metadata <- read.table(metabric_metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
metabric_surv_analysis <- merge(metabric_model_scores, metabric_clin, by = "sample_id")
metabric_surv_analysis <- merge(metabric_surv_analysis, unique(metabric_metadata[, c("sample_id", "PAM50")]), by = "sample_id", all.x = TRUE)
colnames(metabric_surv_analysis)[colnames(metabric_surv_analysis) == "PAM50"] <- "subtype"
metabric_surv_analysis$subtype <- factor(metabric_surv_analysis$subtype)
metabric_cox_models <- list()
metabric_cox_results <- data.frame()
metabric_survival_endpoints <- c("os", "rfs")
for (endpoint in metabric_survival_endpoints) {
  dat <- metabric_surv_analysis[!is.na(metabric_surv_analysis[[paste0(endpoint, "_time")]]) & !is.na(metabric_surv_analysis[[paste0(endpoint, "_event")]]) & metabric_surv_analysis[[paste0(endpoint, "_time")]] > 0 & !is.na(metabric_surv_analysis$score_final_z) & !is.na(metabric_surv_analysis$age) & !is.na(metabric_surv_analysis$npi) & !is.na(metabric_surv_analysis$subtype), ]
  dat$subtype <- droplevels(factor(dat$subtype))
  metabric_cox_models[[endpoint]] <- survival::coxph(as.formula(paste0("survival::Surv(", endpoint, "_time, ", endpoint, "_event) ~ score_final_z + age + npi + strata(subtype)")), data = dat, x = TRUE, y = TRUE)
  metabric_cox_results <- rbind(metabric_cox_results, data.frame(endpoint = endpoint, broom::tidy(metabric_cox_models[[endpoint]], exponentiate = TRUE, conf.int = TRUE), row.names = NULL))
}
figure2_panel_C <- list()
for (endpoint in names(tcga_cox_models)) {
  plot_df <- broom::tidy(tcga_cox_models[[endpoint]], exponentiate = TRUE, conf.int = TRUE)
  plot_df <- plot_df[!grepl("^strata", plot_df$term), ]
  plot_df$term <- gsub("score_oof_z", "Score", plot_df$term)
  plot_df$term <- gsub("stage_model", "", plot_df$term)
  plot_df$term <- factor(plot_df$term, levels = rev(plot_df$term))
  figure2_panel_C[[endpoint]] <- ggplot(plot_df, aes(x = estimate, y = term)) + geom_vline(xintercept = 1, linetype = 2) + geom_point(size = 1.5) + geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15) + scale_x_log10() + theme_classic(base_size = 7) + labs(title = paste("TCGA-BRCA", toupper(endpoint)), x = "Hazard ratio", y = NULL)
}
figure2_panel_D <- list()
for (endpoint in names(metabric_cox_models)) {
  plot_df <- broom::tidy(metabric_cox_models[[endpoint]], exponentiate = TRUE, conf.int = TRUE)
  plot_df <- plot_df[!grepl("^strata", plot_df$term), ]
  plot_df$term <- gsub("score_final_z", "Score", plot_df$term)
  plot_df$term <- factor(plot_df$term, levels = rev(plot_df$term))
  figure2_panel_D[[endpoint]] <- ggplot(plot_df, aes(x = estimate, y = term)) + geom_vline(xintercept = 1, linetype = 2) + geom_point(size = 1.5) + geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.15) + scale_x_log10() + theme_classic(base_size = 7) + labs(title = paste("METABRIC", toupper(endpoint)), x = "Hazard ratio", y = NULL)
}
figure2 <- cowplot::plot_grid(figure2_panel_A, figure2_panel_B, cowplot::plot_grid(plotlist = figure2_panel_C, ncol = 2), cowplot::plot_grid(plotlist = figure2_panel_D, ncol = 1), labels = c("A", "B", "C", "D"), ncol = 2, rel_heights = c(1, 1.2), label_size = 16)

##### model performance analyses #####
performance_datasets <- list(TCGA_outerCV = data.frame(dataset = "TCGA outer cross-validation", sample_id = tcga_model_scores$sample_id, HRD_score = tcga_model_scores$HRD_score, score = tcga_model_scores$score_oof), METABRIC = data.frame(dataset = "METABRIC external validation", sample_id = metabric_model_scores$sample_id, HRD_score = metabric_model_scores$HRD_score, score = metabric_model_scores$score_final))
classification_metrics_by_threshold <- data.frame()
confusion_matrices_primary_threshold <- data.frame()
for (nm in names(performance_datasets)) {
  dat <- performance_datasets[[nm]]
  for (cutoff in hrd_cutoffs) {
    status <- as.numeric(dat$HRD_score >= cutoff)
    roc_tmp <- pROC::roc(response = status, predictor = dat$score, quiet = TRUE)
    ci_tmp <- as.numeric(pROC::ci.auc(roc_tmp))
    pred_class <- as.numeric(dat$score >= primary_threshold)
    tab <- table(factor(pred_class, levels = c(0, 1)), factor(status, levels = c(0, 1)))
    TP <- tab[2, 2]; TN <- tab[1, 1]; FP <- tab[2, 1]; FN <- tab[1, 2]
    row_tmp <- data.frame(cohort = unique(dat$dataset), HRD_cutoff = cutoff, prediction_threshold = primary_threshold, n = length(status), n_HRD = sum(status == 1), n_HRP = sum(status == 0), AUC = as.numeric(pROC::auc(roc_tmp)), AUC_CI_low = ci_tmp[1], AUC_CI_high = ci_tmp[3], sensitivity = TP / (TP + FN), specificity = TN / (TN + FP), PPV = TP / (TP + FP), NPV = TN / (TN + FN), accuracy = (TP + TN) / sum(tab), balanced_accuracy = ((TP / (TP + FN)) + (TN / (TN + FP))) / 2, TP = TP, TN = TN, FP = FP, FN = FN)
    classification_metrics_by_threshold <- rbind(classification_metrics_by_threshold, row_tmp)
    if (cutoff == 42) confusion_matrices_primary_threshold <- rbind(confusion_matrices_primary_threshold, data.frame(cohort = unique(dat$dataset), TP = TP, TN = TN, FP = FP, FN = FN))
  }
}
calibration_plot_data <- data.frame()
for (nm in names(performance_datasets)) {
  dat <- performance_datasets[[nm]]
  dat$status <- as.numeric(dat$HRD_score >= 42)
  dat$bin <- cut(dat$score, breaks = unique(quantile(dat$score, probs = seq(0, 1, 0.1), na.rm = TRUE)), include.lowest = TRUE)
  cal <- aggregate(cbind(score, status) ~ bin, dat, mean, na.rm = TRUE)
  cal$n <- as.numeric(table(dat$bin)[as.character(cal$bin)])
  calibration_plot_data <- rbind(calibration_plot_data, data.frame(dataset = nm, bin = cal$bin, mean_predicted_probability = cal$score, observed_HRD_fraction = cal$status, n = cal$n))
}
decision_curve_analysis <- data.frame()
threshold_grid <- seq(0.01, 0.99, 0.01)
for (nm in names(performance_datasets)) {
  dat <- performance_datasets[[nm]]
  status <- as.numeric(dat$HRD_score >= 42)
  prevalence <- mean(status == 1)
  for (thr in threshold_grid) {
    pred <- as.numeric(dat$score >= thr)
    TP <- sum(pred == 1 & status == 1)
    FP <- sum(pred == 1 & status == 0)
    n <- length(status)
    decision_curve_analysis <- rbind(decision_curve_analysis, data.frame(dataset = nm, threshold_probability = thr, strategy = "Model", net_benefit = TP / n - FP / n * thr / (1 - thr)))
    decision_curve_analysis <- rbind(decision_curve_analysis, data.frame(dataset = nm, threshold_probability = thr, strategy = "Classify all as HRD", net_benefit = prevalence - (1 - prevalence) * thr / (1 - thr)))
    decision_curve_analysis <- rbind(decision_curve_analysis, data.frame(dataset = nm, threshold_probability = thr, strategy = "Classify none as HRD", net_benefit = 0))
  }
}
supplementary_figure1_panel_A <- ggplot(calibration_plot_data, aes(x = mean_predicted_probability, y = observed_HRD_fraction)) + geom_abline(intercept = 0, slope = 1, linetype = 2) + geom_point(aes(size = n)) + geom_line() + facet_wrap(~ dataset) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + theme_classic(base_size = 9) + labs(x = "Mean predicted probability", y = "Observed HRD fraction", size = "n")
supplementary_figure1_panel_B <- ggplot(decision_curve_analysis[decision_curve_analysis$threshold_probability <= 0.80, ], aes(x = threshold_probability, y = net_benefit, linetype = strategy)) + geom_line(linewidth = 0.8) + facet_wrap(~ dataset) + theme_classic(base_size = 9) + coord_cartesian(ylim = c(-0.05, 0.35)) + labs(x = "Threshold probability", y = "Net benefit", linetype = "")
supplementary_figure1 <- cowplot::plot_grid(supplementary_figure1_panel_A, supplementary_figure1_panel_B, labels = c("A", "B"), ncol = 2, label_size = 16)

##### subtype and continuous HRD analyses #####
tcga_eval <- merge(tcga_model_scores, tcga_metadata, by = c("sample_id", "patient_id"), all.x = TRUE)
tcga_eval$score_logit <- qlogis(pmin(pmax(tcga_eval$score_oof, 1e-6), 1 - 1e-6))
tcga_eval$score_z <- as.numeric(scale(tcga_eval$score_logit))
metabric_eval <- merge(metabric_model_scores, metabric_metadata, by = "sample_id", all.x = TRUE)
metabric_eval$score_logit <- qlogis(pmin(pmax(metabric_eval$score_final, 1e-6), 1 - 1e-6))
metabric_eval$score_z <- as.numeric(scale(metabric_eval$score_logit))
continuous_HRD_score_associations <- data.frame()
continuous_score_terms <- data.frame()
binary_score_terms <- data.frame()
subtype_adjustments <- list(PAM50 = "PAM50", ER_group = "ER_group", HER2_group = "HER2_group", TNBC_group = "TNBC_group")
for (cohort in c("TCGA_outerCV", "METABRIC")) {
  dat <- if (cohort == "TCGA_outerCV") tcga_eval else metabric_eval
  score_col <- if (cohort == "TCGA_outerCV") "score_oof" else "score_final"
  dat$score_logit <- qlogis(pmin(pmax(dat[[score_col]], 1e-6), 1 - 1e-6))
  dat$score_z <- as.numeric(scale(dat$score_logit))
  cor_tmp <- cor.test(dat$HRD_score, dat$score_logit, method = "spearman", exact = FALSE)
  fit_tmp <- lm(scale(HRD_score) ~ scale(score_logit), data = dat)
  continuous_HRD_score_associations <- rbind(continuous_HRD_score_associations, data.frame(cohort = cohort, n = sum(complete.cases(dat[, c("HRD_score", "score_logit")])), spearman_rho = unname(cor_tmp$estimate), spearman_p_value = cor_tmp$p.value, beta = coef(summary(fit_tmp))[2, 1], p_value = coef(summary(fit_tmp))[2, 4]))
  for (adj_name in names(subtype_adjustments)) {
    adj_col <- subtype_adjustments[[adj_name]]
    if (adj_col %in% colnames(dat)) {
      tmp <- dat[!is.na(dat$HRD_score) & !is.na(dat$score_z) & !is.na(dat[[adj_col]]), ]
      if (length(unique(tmp[[adj_col]])) > 1) {
        for (cutoff in hrd_cutoffs) {
          tmp$HRD_status_tmp <- as.numeric(tmp$HRD_score >= cutoff)
          fit_bin <- glm(HRD_status_tmp ~ score_z + factor(tmp[[adj_col]]), data = tmp, family = binomial())
          ci_bin <- suppressMessages(confint.default(fit_bin))
          binary_score_terms <- rbind(binary_score_terms, data.frame(cohort = cohort, HRD_cutoff = cutoff, subtype_adjustment = adj_name, model_type = "Subtype-adjusted", OR = exp(coef(fit_bin)["score_z"]), OR_CI_low = exp(ci_bin["score_z", 1]), OR_CI_high = exp(ci_bin["score_z", 2]), p_value = coef(summary(fit_bin))["score_z", 4]))
        }
        fit_cont <- lm(scale(HRD_score) ~ score_z + factor(tmp[[adj_col]]), data = tmp)
        ci_cont <- confint(fit_cont)
        continuous_score_terms <- rbind(continuous_score_terms, data.frame(cohort = cohort, subtype_adjustment = adj_name, model_type = "Subtype-adjusted", beta = coef(fit_cont)["score_z"], beta_CI_low = ci_cont["score_z", 1], beta_CI_high = ci_cont["score_z", 2], p_value = coef(summary(fit_cont))["score_z", 4]))
      }
    }
  }
}
supplementary_figure2_panel_A <- ggplot(tcga_eval, aes(x = score_logit, y = HRD_score, color = PAM50)) + geom_point(alpha = 0.65, size = 1.1) + geom_smooth(method = "lm", se = TRUE, color = "black") + geom_hline(yintercept = 42, linetype = 2) + geom_hline(yintercept = 33, linetype = 3) + scale_color_manual(values = pam50_cols, na.value = "grey70") + theme_classic(base_size = 9) + theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) + labs(title = "TCGA_outerCV: continuous HRD score vs model score", x = "Model score, logit-transformed predicted probability", y = "Continuous HRD score")
supplementary_figure2_panel_B <- ggplot(metabric_eval, aes(x = score_logit, y = HRD_score, color = PAM50)) + geom_point(alpha = 0.65, size = 1.1) + geom_smooth(method = "lm", se = TRUE, color = "black") + geom_hline(yintercept = 42, linetype = 2) + geom_hline(yintercept = 33, linetype = 3) + theme_classic(base_size = 9) + theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5)) + labs(title = "METABRIC: continuous HRD score vs model score", x = "Model score, logit-transformed predicted probability", y = "Continuous HRD score")
plot_df <- binary_score_terms
plot_df$label <- factor(paste0(plot_df$cohort, " | HRD≥", plot_df$HRD_cutoff, " | adj: ", plot_df$subtype_adjustment), levels = rev(paste0(plot_df$cohort, " | HRD≥", plot_df$HRD_cutoff, " | adj: ", plot_df$subtype_adjustment)))
supplementary_figure2_panel_C <- ggplot(plot_df, aes(x = OR, y = label)) + geom_vline(xintercept = 1, linetype = 2) + geom_point(size = 1.7) + geom_errorbarh(aes(xmin = OR_CI_low, xmax = OR_CI_high), height = 0.15) + scale_x_log10() + theme_classic(base_size = 8) + labs(x = "Odds ratio per 1 SD increase in model score", y = "")
plot_df <- continuous_score_terms
plot_df$label <- factor(paste0(plot_df$cohort, " | adj: ", plot_df$subtype_adjustment), levels = rev(paste0(plot_df$cohort, " | adj: ", plot_df$subtype_adjustment)))
supplementary_figure2_panel_D <- ggplot(plot_df, aes(x = beta, y = label)) + geom_vline(xintercept = 0, linetype = 2) + geom_point(size = 1.7) + geom_errorbarh(aes(xmin = beta_CI_low, xmax = beta_CI_high), height = 0.15) + theme_classic(base_size = 8) + labs(x = "Standardized beta for HRD score per 1 SD increase in model score", y = "")
supplementary_figure2 <- cowplot::plot_grid(supplementary_figure2_panel_A, supplementary_figure2_panel_B, supplementary_figure2_panel_C, supplementary_figure2_panel_D, labels = c("A", "B", "C", "D"), ncol = 2, label_size = 16)

##### BRCA subgroup analysis #####
brca_knijn <- read.csv(tcga_brca_status_file, stringsAsFactors = FALSE)
brca_patients <- aggregate(brca_knijn[, c("BRCA_any_mut", "BRCA_any_deepdel", "BRCA_any_silenced", "BRCA_any_Knijnenburg_alt", "BRCA_loss_conservative_DEL_or_SIL")], by = list(patient_id = brca_knijn$patient_id), FUN = max, na.rm = TRUE)
brca_patients$BRCA_analysis_group <- "BRCA-wildtype/no BRCA-loss evidence"
brca_patients$BRCA_analysis_group[brca_patients$BRCA_any_Knijnenburg_alt == 1] <- "BRCA1/2 mutation-only/ambiguous"
brca_patients$BRCA_analysis_group[brca_patients$BRCA_loss_conservative_DEL_or_SIL == 1] <- "BRCA1/2 loss/silencing"
roc_dat_log_knijn <- merge(roc_dat_log, brca_patients, by = "patient_id", all.x = TRUE)
roc_dat_log_knijn$BRCA_analysis_group[is.na(roc_dat_log_knijn$BRCA_analysis_group)] <- "Not annotated in Knijnenburg resource"
brca_group_summary <- aggregate(status ~ BRCA_analysis_group, roc_dat_log_knijn, length)
colnames(brca_group_summary)[2] <- "n"
brca_group_summary$n_HRD <- aggregate(status ~ BRCA_analysis_group, roc_dat_log_knijn, sum)$status
brca_group_summary$n_HRP <- brca_group_summary$n - brca_group_summary$n_HRD
brca_group_summary$percent_HRD <- 100 * brca_group_summary$n_HRD / brca_group_summary$n
brca_group_summary$median_prediction_score <- aggregate(pred ~ BRCA_analysis_group, roc_dat_log_knijn, median)$pred
roc_brca_wt_confident <- roc_dat_log_knijn[roc_dat_log_knijn$BRCA_analysis_group == "BRCA-wildtype/no BRCA-loss evidence", ]
roc_brca_wt <- pROC::roc(response = roc_brca_wt_confident$status, predictor = roc_brca_wt_confident$pred, quiet = TRUE)
auc_brca_wt_confident <- data.frame(group = "BRCA-wildtype/no BRCA-loss evidence", n = nrow(roc_brca_wt_confident), n_HRP = sum(roc_brca_wt_confident$status == 0), n_HRD = sum(roc_brca_wt_confident$status == 1), AUC = as.numeric(pROC::auc(roc_brca_wt)), AUC_CI_low = as.numeric(pROC::ci.auc(roc_brca_wt))[1], AUC_CI_high = as.numeric(pROC::ci.auc(roc_brca_wt))[3])

##### CRISPR screen analyses #####
square <- read.delim(crispr_square_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
crispr_mle <- read.delim(crispr_mle_file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
crispr_gene_col <- grep("^Gene$|gene|symbol", colnames(square), ignore.case = TRUE, value = TRUE)[1]
crispr_beta_col <- grep("beta|score", colnames(square), ignore.case = TRUE, value = TRUE)[1]
crispr_rank_col <- grep("rank", colnames(square), ignore.case = TRUE, value = TRUE)[1]
crispr_sl <- square[square$group == "bottomcenter", , drop = FALSE]
crispr_res <- square[square$group == "topcenter", , drop = FALSE]
crispr_sl <- crispr_sl[order(as.numeric(crispr_sl[[crispr_rank_col]])), , drop = FALSE]
crispr_res <- crispr_res[order(as.numeric(crispr_res[[crispr_rank_col]])), , drop = FALSE]
supplementary_table_s3 <- data.frame(Gene = crispr_sl[[crispr_gene_col]], Beta_Score_Treatment = as.numeric(crispr_sl[[crispr_beta_col]]), Rank = seq_len(nrow(crispr_sl)), stringsAsFactors = FALSE)
supplementary_table_s4 <- data.frame(Gene = crispr_res[[crispr_gene_col]], Beta_Score_Treatment = as.numeric(crispr_res[[crispr_beta_col]]), Rank = seq_len(nrow(crispr_res)), stringsAsFactors = FALSE)
crispr_68_genes <- unique(as.character(head(supplementary_table_s3$Gene, 68)))
crispr_universe_symbols <- unique(as.character(crispr_mle$Gene[!is.na(crispr_mle$Gene) & crispr_mle$Gene != ""]))
crispr_symbol_map <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(c(crispr_68_genes, crispr_universe_symbols)), keytype = "SYMBOL", columns = c("SYMBOL", "ENTREZID"))
crispr_symbol_map <- crispr_symbol_map[!is.na(crispr_symbol_map$ENTREZID), ]
crispr_68_entrez <- unique(as.character(crispr_symbol_map$ENTREZID[crispr_symbol_map$SYMBOL %in% crispr_68_genes]))
uni_crispr <- unique(as.character(crispr_symbol_map$ENTREZID[crispr_symbol_map$SYMBOL %in% crispr_universe_symbols]))
pooled_TERM2GENE_crispr <- pooled_TERM2GENE[pooled_TERM2GENE$gene %in% uni_crispr, ]
pooled_ora_crispr68 <- clusterProfiler::enricher(gene = intersect(crispr_68_entrez, uni_crispr), universe = uni_crispr, TERM2GENE = pooled_TERM2GENE_crispr, TERM2NAME = pooled_TERM2NAME, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500)
pooled_ora_crispr68_all <- as.data.frame(pooled_ora_crispr68)
pooled_ora_crispr68_all <- merge(pooled_ora_crispr68_all, term_metadata, by.x = "ID", by.y = "term", all.x = TRUE, sort = FALSE)
pooled_ora_crispr68_all$minus_log10_FDR <- -log10(pooled_ora_crispr68_all$p.adjust)
pooled_ora_crispr68_all <- pooled_ora_crispr68_all[order(pooled_ora_crispr68_all$p.adjust), ]
pooled_ora_crispr68_sig <- pooled_ora_crispr68_all[pooled_ora_crispr68_all$p.adjust <= 0.05, ]
crispr_plot_terms <- c("DNA Repair", "Recombinational Repair", "DNA Recombination", "DSB Repair via HR", "Homology Directed Repair")
ggpathways_crispr68 <- data.frame()
for (i in seq_along(crispr_plot_terms)) {
  hit <- grep(crispr_plot_terms[i], pooled_ora_crispr68_sig$Description, ignore.case = TRUE)
  if (length(hit) == 0) hit <- grep(crispr_plot_terms[i], pooled_ora_crispr68_sig$name, ignore.case = TRUE)
  if (length(hit) > 0) ggpathways_crispr68 <- rbind(ggpathways_crispr68, pooled_ora_crispr68_sig[hit[1], ])
}
ggpathways_crispr68$plot_label <- crispr_plot_terms[seq_len(nrow(ggpathways_crispr68))]
ggpathways_crispr68$Value <- -log10(ggpathways_crispr68$p.adjust)
figure3_panel_A <- ggplot(ggpathways_crispr68, aes(x = reorder(plot_label, Value), y = Value)) + geom_col(fill = "#0073C2FF", color = "black") + coord_flip() + theme_classic(base_size = 7) + theme(panel.border = element_rect(fill = NA, colour = "black"), axis.text = element_text(colour = "black")) + labs(title = "Synthetically lethal genes enrichment results", x = NULL, y = expression(-log[10] * " " * italic(q)))
figure3_panel_B_table <- supplementary_table_s3[supplementary_table_s3$Gene %in% c("USP1", "FANCA", "XRCC2", "RNASEH2A", "FANCD2", "RAD51", "FANCI", "FANCE", "MCM8", "AUNIP", "TRAIP", "UBE2T", "BLM"), c("Gene", "Beta_Score_Treatment", "Rank")]
figure3 <- cowplot::plot_grid(figure3_panel_A, ggplotify::as.ggplot(gridExtra::tableGrob(figure3_panel_B_table, rows = NULL)), labels = c("A", "B"), ncol = 2, label_size = 16)

##### published signature benchmarking #####
expr_mb_symbol <- expr_mb[!is.na(expr_mb$Entrez_Gene_Id), ]
expr_mb_symbol$gene_symbol <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = as.character(expr_mb_symbol$Entrez_Gene_Id), column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
expr_mb_symbol <- expr_mb_symbol[!is.na(expr_mb_symbol$gene_symbol), ]
expr_mb_symbol <- expr_mb_symbol[!duplicated(expr_mb_symbol$gene_symbol), ]
rownames(expr_mb_symbol) <- expr_mb_symbol$gene_symbol
expr_mb_symbol <- expr_mb_symbol[, intersect(colnames(expr_mb_symbol), metabric_eval$sample_id), drop = FALSE]
jacobson_templates <- as.data.frame(readxl::read_excel(jacobson_template_file))
jacobson_genes <- intersect(jacobson_templates$gene, rownames(expr_mb_symbol))
jacobson_scores <- data.frame(sample_id = colnames(expr_mb_symbol), Jacobson_TS228_score = NA_real_, stringsAsFactors = FALSE)
if (length(jacobson_genes) > 1) {
  jac_expr <- t(scale(t(as.matrix(expr_mb_symbol[jacobson_genes, , drop = FALSE]))))
  jac_hrd <- as.numeric(jacobson_templates$HRD_centroid[match(jacobson_genes, jacobson_templates$gene)])
  jac_hrp <- as.numeric(jacobson_templates$HRP_centroid[match(jacobson_genes, jacobson_templates$gene)])
  jac_delta <- jac_hrd - jac_hrp
  jacobson_scores$Jacobson_TS228_score <- as.numeric(t(jac_expr) %*% jac_delta) / length(jac_delta)
}
pan_scores_metabric <- read.table(pan_hrd200_scores_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
signature_compare_metabric <- merge(metabric_eval, jacobson_scores, by = "sample_id")
signature_compare_metabric <- merge(signature_compare_metabric, pan_scores_metabric[, c("sample_id", "Pan_HRD200_score")], by = "sample_id")
signature_compare_metabric$Our_5gene_score <- signature_compare_metabric$score_final
signature_compare_metabric$HRD_status <- as.numeric(signature_compare_metabric$HRD_score >= 42)
signature_compare_complete <- signature_compare_metabric[complete.cases(signature_compare_metabric[, c("HRD_status", "Our_5gene_score", "Jacobson_TS228_score", "Pan_HRD200_score")]), ]
signature_compare_tnbc <- signature_compare_complete[signature_compare_complete$TNBC_group == "TNBC", ]
roc_all_ours <- pROC::roc(signature_compare_complete$HRD_status, signature_compare_complete$Our_5gene_score, quiet = TRUE)
roc_all_jacobson <- pROC::roc(signature_compare_complete$HRD_status, signature_compare_complete$Jacobson_TS228_score, quiet = TRUE)
roc_all_pan <- pROC::roc(signature_compare_complete$HRD_status, signature_compare_complete$Pan_HRD200_score, quiet = TRUE)
roc_tnbc_ours <- pROC::roc(signature_compare_tnbc$HRD_status, signature_compare_tnbc$Our_5gene_score, quiet = TRUE)
roc_tnbc_jacobson <- pROC::roc(signature_compare_tnbc$HRD_status, signature_compare_tnbc$Jacobson_TS228_score, quiet = TRUE)
roc_tnbc_pan <- pROC::roc(signature_compare_tnbc$HRD_status, signature_compare_tnbc$Pan_HRD200_score, quiet = TRUE)
signature_auc_summary <- data.frame(analysis_set = rep(c("METABRIC: all analyzable tumors", "METABRIC: strict receptor-defined TNBC"), each = 3), signature = rep(c("Our 5-gene model", "Jacobson TS228", "Pan HRD200"), 2), AUC = c(as.numeric(pROC::auc(roc_all_ours)), as.numeric(pROC::auc(roc_all_jacobson)), as.numeric(pROC::auc(roc_all_pan)), as.numeric(pROC::auc(roc_tnbc_ours)), as.numeric(pROC::auc(roc_tnbc_jacobson)), as.numeric(pROC::auc(roc_tnbc_pan))))
roc_all_list <- stats::setNames(list(roc_all_ours, roc_all_jacobson, roc_all_pan), c(paste0("5-gene model (AUC = ", round(as.numeric(pROC::auc(roc_all_ours)), 2), ")"), paste0("Jacobson TS228 (AUC = ", round(as.numeric(pROC::auc(roc_all_jacobson)), 2), ")"), paste0("Pan HRD200 (AUC = ", round(as.numeric(pROC::auc(roc_all_pan)), 2), ")")))
roc_tnbc_list <- stats::setNames(list(roc_tnbc_ours, roc_tnbc_jacobson, roc_tnbc_pan), c(paste0("Our 5-gene model (AUC = ", round(as.numeric(pROC::auc(roc_tnbc_ours)), 2), ")"), paste0("Jacobson TS228 (AUC = ", round(as.numeric(pROC::auc(roc_tnbc_jacobson)), 2), ")"), paste0("Pan HRD200 (AUC = ", round(as.numeric(pROC::auc(roc_tnbc_pan)), 2), ")")))
supplementary_figure3_panel_A <- pROC::ggroc(roc_all_list, legacy.axes = TRUE, linewidth = 1) + geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey50") + scale_color_manual(values = c("#A73030FF", "#20854EFF", "#EFC000FF")) + theme_classic(base_size = 9) + theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom", legend.title = element_blank(), panel.border = element_rect(fill = NA, colour = "black")) + labs(title = "METABRIC: all included tumors", x = "1 - Specificity", y = "Sensitivity")
supplementary_figure3_panel_B <- pROC::ggroc(roc_tnbc_list, legacy.axes = TRUE, linewidth = 1) + geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey50") + scale_color_manual(values = c("#A73030FF", "#20854EFF", "#EFC000FF")) + theme_classic(base_size = 9) + theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom", legend.title = element_blank(), panel.border = element_rect(fill = NA, colour = "black")) + labs(title = "METABRIC: strict receptor-defined TNBC", x = "1 - Specificity", y = "Sensitivity")
supplementary_figure3 <- cowplot::plot_grid(supplementary_figure3_panel_A, supplementary_figure3_panel_B, labels = c("A", "B"), ncol = 2, label_size = 16)

##### supplementary table objects ##### 
supplementary_table_s2 <- classification_metrics_by_threshold
supplementary_tables <- list(S1 = supplementary_table_s1, S2 = supplementary_table_s2, S3 = supplementary_table_s3, S4 = supplementary_table_s4)

##### material objects #####
material_objects <- list(Figure1 = figure1, Figure2 = figure2, Figure3 = figure3, SupplementaryFigure1 = supplementary_figure1, SupplementaryFigure2 = supplementary_figure2, SupplementaryFigure3 = supplementary_figure3, SupplementaryTables = supplementary_tables, Model = model_obj_boot, Coefficients = model_coefficients, WGCNA = list(module_members = module.members, module_trait_correlation = moduleTraitCor, module_trait_pvalue = moduleTraitPvalue, module1_enrichment = pooled_ora_module1_all), Performance = list(classification_metrics = classification_metrics_by_threshold, confusion_matrices = confusion_matrices_primary_threshold, calibration = calibration_plot_data, decision_curve = decision_curve_analysis), Survival = list(TCGA = tcga_cox_results, METABRIC = metabric_cox_results), HRDSensitivity = list(binary = binary_score_terms, continuous = continuous_score_terms, continuous_associations = continuous_HRD_score_associations), PublishedSignatures = list(master = signature_compare_metabric, auc_summary = signature_auc_summary), CRISPR = list(synthetically_lethal = supplementary_table_s3, resistance = supplementary_table_s4, enrichment = pooled_ora_crispr68_all))
