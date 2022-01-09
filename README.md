# Notes and code to generate figures in:
# Molecular innovations in plant TIR-based immunity signaling
by Dmitry Lapin, Oliver Johanndrees, Zhongshou Wu, Xin Li, Jane E. Parker

## Figure 2A
### Source data:

1. PTI - Bjornson et al. 2021 The transcriptional landscape of Arabidopsis thaliana pattern-triggered immunity; (https://www.nature.com/articles/s41477-021-00874-5); Supplementary Table 1 was used
2. ETI - Ngou et al. 2021 Mutual potentiation of plant immunity by cell-surface and intracellular receptors; (https://www.nature.com/articles/s41586-021-03315-7); raw RNAseq data were used, ENA/SRA ID - PRJEB34955
3. PTI+ETI - Saile et al. 2020 Two unequally redundant "helper" immune receptor families mediate Arabidopsis thaliana intracellular "sensor" immune receptor functions; (https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000783); raw data were used, ENA/SRA ID - PRJNA637584

### Related content of repository:

1. Code to process raw RNAseq data - ./processing-of-raw-RNAseq-data/1_steps-from-download-to-counts.txt
2. SRA accession codes for ETI and PTI+ETI datasets - ./processing-of-raw-RNAseq-data/PTI-and-ETI and ./processing-of-raw-RNAseq-data/ETI
3. Code to perform differential gene expression analysis - ./processing-of-raw-RNAseq-data/2_steps-from-counts-to-DE.R
4. Code to generate heatmap in Figure 2A from log2FC data - ./heatmap-for-Figure2A/heatmap-for-GOI.R with input log2FC data and geneIDs of interest in ./heatmap-for-Figure2A/input

## Figure 2B
### Source data:

Processed RNAseq data are from Johanndrees and Baggs et al. 2021 (preprint) Differential EDS1 requirement for cell death activities of plant TIR-domain proteins (https://www.biorxiv.org/content/10.1101/2021.11.29.470438v2.full)

Gene expression data and corresponding metadata (Supplementary-Table-5_V1_20211111.xlsx) can be found in Edmond collection
Supporting information for Johanndrees, Baggs et al 2021
https://edmond.mpdl.mpg.de/imeji/collection/Pv5t4gM8Sv0TOrCU

### Related content of repository:
1. Code to draw boxplots in Figure 2B - ./boxplot-for-Figure2B/boxplot-scaled-expression-GOI.R with input data from the above Edmond collection


## Figure 4C
### Source data:

Processed ChIP-seq data for SARD1 (Sun et al 2015, ChIP-seq reveals broad roles of SARD1 and CBP60g in regulating plant immunity https://www.nature.com/articles/ncomms10159) come from the preprint Griebel, Lapin et al. 2021 Arabidopsis Topless-related 1 mitigates physiological damage and growth penalties of induced immunity (https://www.biorxiv.org/content/10.1101/2021.07.07.451397v1.full)

SARD1 chromatin enrichment profiles (bigwig format) for Psm-infected plants can be found here:
https://edmond.mpdl.mpg.de/imeji/collection/U6N5zIOIWgjjMZCu

### Related content of repository:
1. Code to draw metaplot in Figure 4C - ./metaplot-for-Figure4C with input data from the above Edmond collection and geneset and bed files in respective folders in ./metaplot-for-Figure4C
