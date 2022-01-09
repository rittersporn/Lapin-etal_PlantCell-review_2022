##### Script to draw a heatmap of gene expression
# input files (log2FC) and lists of genes of interest are in the same folder as the script heatmap-for-Figure2A/input

library("pheatmap")

# load the data, sets of genes, annotations, bring them to necessary formats
log2fc <- read.delim("Bjornson-etal_log2FC_STable1.txt", header = TRUE, sep = "\t")

log2fc_Ngou <- read.delim("Ngou-etal-2021_log2FC_V1_20220109.txt", header = TRUE, sep = "\t")
colnames(log2fc_Ngou) <- c("GeneID", "log2FC_e2ETI_4h", "padj_e2ETI_4h")

log2fc_Saile <- read.delim("Saile-etal-2020_log2FC_V1_20220109.txt", header = TRUE, sep = "\t")

log2fc <- merge(log2fc, log2fc_Ngou,
                by.x = "AGI", by.y = "GeneID")

log2fc <- merge(log2fc, log2fc_Saile,
                by.x = "AGI", by.y = "GeneID")

GOI <- as.character(read.table("unique-geneIDs_At-TIR-containing.txt", header = FALSE)[,1])


# subset and reorder columns
row.names(log2fc) <- log2fc$AGI
col_order <- c("Col_mock_090", "Col_flg22_090", "Col_elf18_090", "Col_nlp20_090", "Col_OGs_090", "Col_Pep1_090",
               "log2FC_e2ETI_4h",
               "log2FC_WT_EV_4h_vs_WT_untreated_0h", "log2FC_WT_avrRps4_4h_vs_WT_untreated_0h", "log2FC_WT_avrRpt2_4h_vs_WT_untreated_0h", "log2FC_WT_avrRpm1_4h_vs_WT_untreated_0h")
log2fc_GOI <- log2fc[log2fc$AGI %in% GOI, col_order]
log2fc_GOI <- as.matrix(log2fc_GOI)

# limit maximum absolute values to 5
log2fc_GOI[log2fc_GOI > 5] <- 5
log2fc_GOI[log2fc_GOI < -5] <- -5

# build the heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("#4dac26", "#252525", "#d01c8b"))(paletteLength)
myBreaks <- c(seq(-5, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(5/paletteLength, 5, length.out=floor(paletteLength/2)))

pheatmap(log2fc_GOI, cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, show_colnames = TRUE,
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2",
         color = myColor, breaks = myBreaks,
         border_color = NA,
         fontsize_col = 14,
         angle_col = 45, width = 9, height = 15,
         filename = "heatmap_GOI_V1_20220901.pdf")


