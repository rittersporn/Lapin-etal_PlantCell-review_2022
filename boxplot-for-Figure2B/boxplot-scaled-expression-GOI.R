##### script to viz. distribution of scaled expression levels of genes of interest


library("dplyr")
library("tidyr")
library("ggplot2")

#-------------------------------
# custom functions
#-------------------------------
  
#### function to load expression (tpm) values and standardize the gene expression levels ####
standardise_tpm <- function(tpm_df){
  tpm_df[tpm_df == 0] <- NA # remove 0 to prevent their influence on mean and sd
  tpm_df <- cbind(tpm_df$GeneID, as.data.frame(scale(log2(as.matrix(tpm_df[, 2:ncol(tpm_df)])))))
  colnames(tpm_df)[1] <- "GeneID"
  print(paste0("mean value per column (should be close to 0) - ", colMeans(as.matrix(tpm_df[, 2:ncol(tpm_df)]), na.rm = TRUE)))
  print(paste0("sd per column (should be 1) - ", apply(as.matrix(tpm_df[, 2:ncol(tpm_df)]), 2, sd, na.rm = TRUE)))
  return(tpm_df)
}

#### function to tidy up tpm dataframe
tidyup_df <- function(tpm_df, metadata_df, GOI){
  # tidy up the tpm dataframe
  tpm_tidy <- tpm_df %>% pivot_longer(cols = colnames(tpm_df)[2]:colnames(tpm_df)[length(colnames(tpm_df))],
                                      names_to = "SRR",
                                      values_to = "tpm")
  
  # remove rows without expression info (0 tpm initially, see standardise_tpm()) 
  tpm_tidy <- tpm_tidy[!is.na(tpm_tidy$tpm), ]
  # produce control messages
  print(paste0(length(unique(tpm_tidy$GeneID[tpm_tidy$GeneID %in% GOI])), " genes of interest have tpm>0"))
  print(paste0("only tpm>0 are used for plotting"))
  # annotate experiments as triggered or untriggered tissues
  tpm_tidy <- merge.data.frame(x = tpm_tidy,
                               y = metadata_df,
                               by.x = "SRR",
                               by.y = "SRA_accession",
                               all.x = TRUE,
                               all.y = FALSE)
  
  # return tidy df
  return(tpm_tidy)
}

#---------------------------------------------------
# perform visualization of scaled expression data #
#---------------------------------------------------
# load metadata about expression datasets
metadata_tpm <- read.delim("./TPM_values/metainfo.txt",
                             header = TRUE, stringsAsFactors = FALSE)

# load and standardize gene expression levels
tpm <- read.delim("./TPM_values/maize_quant_gene_level.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
#tpm <- read.delim("./TPM_values/rice_quant_gene_level.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
#tpm <- read.delim("./TPM_values/barley_quant_gene_level.txt", sep = " ", header = TRUE, stringsAsFactors = FALSE)
tpm <- standardise_tpm(tpm)

# define genes of interest
GOI <- c("Zm00001d021820", "Zm00001d006617", "Zm00001d007761")
#GOI <- c("Os07t0566800")
#GOI <- c("HORVU2Hr1G039670")

# tidy up the dataframe
tpm_tidy <- tidyup_df(tpm_df = tpm, metadata_df = metadata_tpm,
                      GOI = GOI)

# make boxplots
pdf("boxplot_GOI_barley_20220109.pdf")
df_for_boxplot <- tpm_tidy[tpm_tidy$GeneID %in% c(GOI),]

df_for_boxplot$Induction_status <- factor(df_for_boxplot$Induction_status, levels = c("untriggered", "triggered"))
p <- ggplot(data = df_for_boxplot, aes(x = Induction_status, y = tpm))
p + geom_boxplot(aes(fill = Induction_status), outlier.color = NA) +
  coord_cartesian(ylim = c(-3, 3))+
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(face = "bold", size = 8, angle = 45, hjust = 1))+
  labs(title = "barley",
       y = "expression Z-score")
dev.off()
