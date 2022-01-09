#--------------------------------------------------------------------------------------------------------
# Generate the bed files for gene sets of interest
#--------------------------------------------------------------------------------------------------------
rm(list = ls(all=TRUE))

# Working and output directories
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)
in_gff = file.path(wd, "TAIR10_GFF","TAIR10_GFF3_genes_transposons.gff")
out_data = file.path(wd, "bed_files")

# Load the TAIR10 gff file
gff = read.delim(in_gff, header=F, comment.char="#")
dim(gff) # 656309      9

# Retain necessary information
gff=gff[,c(1,3,4,5,7,9)]

# Get the gene regions
gff_gene = gff[which(gff$V3 == "gene"),]
dim(gff_gene) # 28775      6

# Extract the gene names
gff_gene$V9<-substr(gff_gene$V9, start = regexpr("AT", gff_gene$V9), stop = regexpr("AT", gff_gene$V9) + 8)
head(gff_gene)
dim(gff_gene) # 28775     

# Rename columns
colnames(gff_gene) = c("chrom","region","chromStart","chromEnd","strand","geneID")
head(gff_gene)

gff_gene=unique(gff_gene)
dim(gff_gene) # 28775     6

head(gff_gene)

# order df by geneID
gff_gene <- gff_gene[order(gff_gene$geneID),]


#Make new column named score, as this is required to attain the .bed file format (5th column)
#To let other systems know that column score is not relevant we fill it with "." as required for .bed files
#This column is necessary because BED6 files contain strands which we will include
gff_gene$score<- NA
gff_gene[is.na(gff_gene)]<-"."

head(gff_gene)

#Rename and reorder columns so that it can be read as a .bed file
gff_gene <- gff_gene[,c(1,3,4,6,7,5)]
head(gff_gene)

# Write table to file
filename =  file.path(out_data, "BED_gene_TAIR10.bed")
write.table(gff_gene, filename, sep="\t", row.names=FALSE, col.names = FALSE, quote=FALSE)

# write BED files for specific gene sets
filename = file.path(out_data, "BED_gene_ADR1_strict.bed")
df <- read.table(file.path(wd, "gene_sets", "DEGs_ADR1.txt"), stringsAsFactors = FALSE)
write.table(gff_gene[gff_gene$geneID %in% df$V1,], filename, sep="\t", row.names=FALSE, col.names = FALSE, quote=FALSE)

filename = file.path(out_data, "BED_gene_NRG1.bed")
df <- read.table(file.path(wd, "gene_sets", "DEGs_NRG1_helperless_DL.txt"), stringsAsFactors = FALSE)
write.table(gff_gene[gff_gene$geneID %in% df$V1,], filename, sep="\t", row.names=FALSE, col.names = FALSE, quote=FALSE)

