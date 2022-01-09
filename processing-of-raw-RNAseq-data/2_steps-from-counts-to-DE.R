###### Script to perform differential gene expression
# for RNAseq data Saile et al 2020 (PTI+ETI)
# doi.org/10.1371/journal.pbio.3000783

# for the ETI dataset Ngou et al 2021, adjust input files and names of output files

# libraries
library("tximport")
library("DESeq2")
library("gridExtra")


# load input files with transcript abundance information and summarize expression per gene using tximport
RunTable <- read.delim("SraRunTable_merged.txt", header = TRUE, sep = " ")
AccList <- RunTable$SRR
files <- file.path("quants", "mapped_cDNA", AccList, "quant.sf")
names(files) <- as.character(AccList)

tx2gene <- read.delim(file.path("reference", "tx2gene.txt"), stringsAsFactors = FALSE)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# merge data.frames to maintain order ("sort=FALSE" keeps order of names)
txiNames <- data.frame(SRR_txi = colnames(txi$counts))
colData <- merge(txiNames, RunTable, by.x = "SRR_txi", by.y = "SRR", sort = FALSE)
row.names(colData) <- colData$SRR_txi
colData$Sample_Name <- as.factor(paste0(colData$Genotype, "_", colData$Treatment, "_", colData$Timepoint))

# show summary of labels
table(colData$Sample_Name)

# import raw count estimates as DESeqDataSet (dds) object
dds <- DESeqDataSetFromTximport(txi = txi, colData = colData,
                                design = ~ Sample_Name)


# check the raw expression data
dim(dds) # dimensions
head(rownames(dds)) # gene names
head(colnames(dds)) # sample names

pdf("Saile-etal-2020_before-filter_V1_20220109.pdf", width = 7, height = 5)
hist((log2(counts(dds))), 50,
     xlab = "log2 raw counts",
     main = "Count Distribution before Filtering", col = "#bdd7e7")

hist((colSums(assay(dds)) / 1e6), 30,
     xlab = "M reads per sample",
     main = "Sample Size before Filtering", col = "#bae4b3")
dev.off()

# discard all genes that are not expressed in any of the samples
print(paste("Total number of genes:", nrow(dds)))
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
print(paste("Number of expressed genes:", nrow(dds)))


# additional filtering to discard genes with low expression
keep <- rowSums(counts(dds) >= 10) >= 30
dds <- dds[keep,]
print(paste("Number of genes after additional filtering:", nrow(dds)))

pdf("Saile-etal-2020_after-filter_V1_20220109.pdf", width = 7, height = 5)
hist((log2(counts(dds))), 50,
     xlab = "log2 raw counts",
     main = "Count Distribution after Filtering", col = "#6baed6")

hist((colSums(counts(dds)) / 1e6), 30,
     xlab = "M reads per sample",
     main = "Sample Size after Filtering", col = "#74c476")
dev.off()


# data transformation for analysis and visualization
# data for DE analysis stays untransformed

# both vst and rlog functions provide normalization:
# vst recommended for larger datasets (n > 30)
# rlog recommended for small datasets (n < 30)
trans_counts <- vst(dds, blind = TRUE)

pdf("Saile-etal-2020_after-normalization_V1_20220109.pdf", width = 7, height = 5)
hist(assay(trans_counts), 50,
     xlab = "log2 normalized counts",
     main = "normalized counts", col = "salmon")
dev.off()

# export transformed counts if needed
# write.table(assay(trans_counts), "Saile-etal-2020_normalized-counts_V1_20220109.txt", quote = FALSE)

# color samples based on different parameters using argument "intgroup"
pca <- plotPCA(trans_counts, intgroup = c("Sample_Name"), ntop = nrow(trans_counts))

# PCA plot for samples
pdf("Saile-etal-2020_PCA-for-normalized-counts_V1_20220109_.pdf", width = 12, height = 10, useDingbats = FALSE)
grid.arrange(pca)
dev.off()


### the actual DE analysis starts here ###
# calculate differential expression using Wald test
dds <- DESeq(dds, test = "Wald")

# extract results with contrast for two groups using WT untreated samples at 0 h as a baseline
comb_res_df <- data.frame()
ref_group <- "WT_untreated_0h"
for (treatment_group in c("WT_EV_4h", "WT_avrRps4_4h", "WT_avrRpt2_4h", "WT_avrRpm1_4h", "WT_EV_8h", "WT_avrRps4_8h", "WT_avrRpt2_8h", "WT_avrRpm1_8h")) {
        res <- results(dds, contrast = c("Sample_Name", treatment_group, ref_group))
        res_df <- data.frame(GeneID = rownames(res),
                             log2FC = res$log2FoldChange,
                             padj = res$padj)
        colnames(res_df) <- c("GeneID",
                              paste0("log2FC", "_", treatment_group, "_vs_", ref_group),
                              paste0("padj", "_", treatment_group, "_vs_", ref_group))
        if (nrow(comb_res_df) == 0){
                comb_res_df <- res_df   
        } else {
                comb_res_df <- merge(comb_res_df, res_df,
                                     by.x = "GeneID", by.y = "GeneID")
        }
        rm(res, res_df)
}



# write results of the differential expression analysis into a file
write.table(comb_res_df, "Saile-etal-2020_log2FC_V1_20220109.txt", sep = "\t",
            quote = FALSE, col.names = TRUE, row.names = FALSE)
