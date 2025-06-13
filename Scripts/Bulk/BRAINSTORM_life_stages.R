# packages loading ----
library("vsn")
library("hexbin")
library("IHW")
library(genefilter)
library(pasilla)
library(PasillaTranscriptExpr)
library("RColorBrewer")
library(DESeq2)
library(devtools)
library(ashr)
library(tidyverse)
library("dplyr")
library(ggplot2)
library(ggrepel)
library(limma)
library(pheatmap)

set.seed(1234)

#set wd in /home/irenerossi/RAW_DATA ----
# load reads and pheno objects ----
countData <- read.csv("/home/irenerossi/RAW_DATA/mating_bulk2021_trim2024/Ag_bulk_mating_raw_matrix2024.csv", header = TRUE, sep = "\t", row.names = 1) # count matrix
metaData <- read.csv("/home/irenerossi/RAW_DATA/phenofile_bulk_background.csv", header = TRUE, sep = "\t") # table of sample information

# check counts
# countData means
#mean_values <- colMeans(countData[sapply(countData, is.numeric)])
# countData sums
#sum_values <- colSums(countData[sapply(countData, is.numeric)])
#mean(countData$X3E.virgin.Brain.Ag_count) # 0.6544049
# Mean of all other columns (all values combined)
#mean(mean_values[names(mean_values) != "X3E.virgin.Brain.Ag_count"]) # 1339.56
#mean(sum_values[names(sum_values) != "X3E.virgin.Brain.Ag_count"]) # 18535491
# remove sample 3E because outlier
countData <- subset(countData, select = -X3E.virgin.Brain.Ag_count)
rownames(metaData) <- metaData$samplename
# ensure factors
metaData$mating.status <- as.factor(metaData$mating.status)
metaData$blood.feeding.status <- as.factor(metaData$blood.feeding.status)

metaData <- metaData[-13,-1]
# rename stuff
rownames(metaData) <- gsub("-Brain-Ag", "", as.character(rownames(metaData)))
rownames(metaData) <- gsub("-", "_", as.character(rownames(metaData))) # remove dashes (can interfere with R functions)
metaData$blood.feeding.status <- gsub("-", "_", as.character(metaData$blood.feeding.status)) # remove dashes (can interfere with R functions)

colnames(countData) <- rownames(metaData)
#is.data.frame(countData) # TRUE
# remove mated not blood-fed (blood meal offered) condition
countData <- countData %>% select(-contains('NBF'))
metaData <- metaData[- grep("NBF", rownames(metaData)),]

#ncol(countData) == nrow(metaData) # TRUE
colnames(countData)

# build DEseq2 objects ----
# create test vs ct objects ----
#colnames(metaData)
# data to create a mating object
countData_m <- countData %>% dplyr:: select(grep("virgin", names(countData)), grep("mated", names(countData)))
metaData_m <- subset(metaData, blood.feeding.status %in% "Non_blood.fed")
# order metadata to get it to match countData
metaData_m <- metaData_m[order(metaData_m$mating.status, decreasing = T),]
#tail(metaData_m) # correct
#length(colnames(countData_m)) # 19 samples
# check, it must be true!:
#colnames(countData_m) == rownames(metaData_m) # it is

# data to create a blood-feeding object
countData_bf <- countData %>% dplyr:: select(grep("mated", names(countData)), grep("BF", names(countData)))
metaData_bf <- subset(metaData, mating.status %in% "mated")
# order metadata to get it to match countData
metaData_bf <- metaData_bf[order(metaData_bf$blood.feeding.status, decreasing = T),]
#tail(metaData_bf) # correct
#length(colnames(countData_bf)) # 20 samples
# check, it must be true!:
#colnames(countData_bf) == rownames(metaData_bf) # it is

dds.m <- DESeqDataSetFromMatrix(countData = countData_m, 
                                colData = metaData_m, 
                                design = ~ mating.status)

dds.bf <- DESeqDataSetFromMatrix(countData = countData_bf, 
                                 colData = metaData_bf, 
                                 design = ~ blood.feeding.status)

# filter out by at least 3 samples with a count of 10 or higher
nrow(dds.m)
keep <- rowSums(counts(dds.m) >= 10) >= 3
dds.m <- dds.m[keep,]
nrow(dds.m) # 10,212; 3,625 genes removed

nrow(dds.bf)
keep <- rowSums(counts(dds.bf) >= 10) >= 3
dds.bf <- dds.bf[keep,]
nrow(dds.bf) # 10,296; 3,541 genes removed

####
# Set the wanted condition as reference level (otherwise dds chooses alphabetic
# order for the conditions) phenofile cathegories interactions
dds.m$group <- factor(paste0(dds.m$mating.status, dds.m$blood.feeding.status))
design(dds.m) <- ~ group
dds.bf$group <- factor(paste0(dds.bf$mating.status, dds.bf$blood.feeding.status))
design(dds.bf) <- ~ group
# set reference level
dds.m$group <- relevel(dds.m$group, "virginNon_blood.fed")
dds.m$group
dds.bf$group <- relevel(dds.bf$group, "matedNon_blood.fed")
dds.bf$group
####

# apply DESEQ function ----
dds.m <- DESeq(dds.m)
results(dds.m)

dds.bf <- DESeq(dds.bf)
results(dds.bf)


# EXPLORE EXPERIMENTAL CONDITIONS DEA ----
# MATING ----
# AAA: explore the two datasets (dds.m and dds.bf) separately and put them in the dds object not to re-write the script
dds <- dds.m
# shrinked results conditionwise ----
shres_MvsV <- lfcShrink(dds, type="apeglm", lfcThreshold=0, 
                        coef=c("group_matedNon_blood.fed_vs_virginNon_blood.fed"))
# padj < 0.05: subset results to return genes ----
sig_shres_MvsV <- shres_MvsV %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
#sig_shres_VvsM # look at FC


# PCA rld ----
# PC1,2 rld
rld <- rlog(dds, blind = FALSE)
# * plot PCA MATING for paper ----
rld_m <- rld
pcaData_m <- plotPCA(rld, intgroup=c("mating.status", "blood.feeding.status"), returnData=TRUE)
percentVar_m <- round(100 * attr(pcaData_m, "percentVar"))
ggplot(pcaData_m, aes(PC1, PC2, col = mating.status)) +
  geom_point(size=2, shape=16) +
  xlab(paste0("PC1: ",percentVar_m[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_m[2],"% variance")) +
  theme(
    axis.text = element_text(size = 10),  # Smaller axis labels
    axis.title = element_text(size = 10),  # Smaller PC_1 and PC_2 labels
    legend.title = element_text(size = 10),  # Smaller legend title
    legend.text = element_text(size = 10),  # Smaller legend text
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(colour = "black")
  ) +
  labs(colour="mating status") +
  xlim(-30, 30) +
  ylim(-30, 30) +
  scale_color_manual(labels = c("mated", "virgin"), values = c("#FF9966", "aquamarine4")) +
  coord_fixed()

# (LFC = 0.25 i.e. ~20% difference in gene expr; FC ~ 1.19)
LFChmMV <- (sig_shres_MvsV$log2FoldChange)
LFChmMVhigh <- which(LFChmMV>0.25)
LFChmMVlow <- which(LFChmMV<(-0.25))
heatmapMVhigh <- sig_shres_MvsV[LFChmMVhigh, ]
heatmapMVhigh <- heatmapMVhigh[order(heatmapMVhigh$log2FoldChange),]
heatmapMVlow <- sig_shres_MvsV[LFChmMVlow, ]
heatmapMVlow <- heatmapMVlow[order(heatmapMVlow$log2FoldChange),]

# GO analysis files ----
# Export gene list for GO TERMS analysis
# virgin vs mated
DE_MVall1 <- which(abs(sig_shres_MvsV$log2FoldChange)>0.25)
DE_MVall <- sig_shres_MvsV[DE_MVall1, ]
DE_MVall <- DE_MVall[order(DE_MVall$log2FoldChange),]
#write.table(DE_MVall, file = "/home/irenerossi/GOOGLE/Brainstorm_PhD/bulk-brain_results/Anopheles_bulk/R_results_bulk_mating/BMV_only/DE_lists_BMVonly_LFCplusminus025/DE_MVall", row.names = F, col.names = T)


# BLOOD-FEEDING ----
# AAA: explore the two datasets (dds.m and dds.bf) separately and put them in the dds object not to re-write the script
dds <- dds.bf
# shrinked results conditionwise ----
shres_BFvsM <- lfcShrink(dds, type="apeglm", lfcThreshold=0, 
                         coef=c("group_matedBlood.fed_vs_matedNon_blood.fed"))
# padj < 0.05: subset results to return genes ----
sig_shres_BFvsM <- shres_BFvsM %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
#sig_shres_BFvsM # look at FC

# PCA rld ----
# PC1,2 rld
rld <- rlog(dds, blind = FALSE)
plotPCA(rld, intgroup=c("mating.status", "blood.feeding.status"))# + stat_ellipse()

# * plot PCA BLOOD FEEDING for paper ----
rld_bf <- rld
pcaData_bf <- plotPCA(rld, intgroup=c("mating.status", "blood.feeding.status"), returnData=TRUE)
percentVar_bf <- round(100 * attr(pcaData_bf, "percentVar"))
ggplot(pcaData_bf, aes(PC1, PC2, col = blood.feeding.status, shape = blood.feeding.status)) +
  geom_point(size=2, shape=16) +
  xlab(paste0("PC1: ",percentVar_bf[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_bf[2],"% variance")) +
  theme(
    axis.text = element_text(size = 10),  # Smaller axis labels
    axis.title = element_text(size = 10),  # Smaller PC_1 and PC_2 labels
    legend.title = element_text(size = 10),  # Smaller legend title
    legend.text = element_text(size = 10),  # Smaller legend text
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(colour = "black")
  ) +
  labs(colour="condition") +
  xlim(-30, 30) +
  ylim(-30, 30) +
  scale_color_manual(labels = c("mated blood-fed", "mated sugar-fed"), values = c("#CC0066", "#FF9966")) +
  coord_fixed()


# (LFC = 0.25 i.e. ~20% difference in gene expr; FC ~ 1.19)
LFChbBFM <- (sig_shres_BFvsM$log2FoldChange)
LFChbBFMhigh <- which(LFChbBFM>0.25)
LFChbBFMlow <- which(LFChbBFM<(-0.25))
heatmapBFMhigh <- sig_shres_BFvsM[LFChbBFMhigh, ]
heatmapBFMhigh <- heatmapBFMhigh[order(heatmapBFMhigh$log2FoldChange),]
heatmapBFMlow <- sig_shres_BFvsM[LFChbBFMlow, ]
heatmapBFMlow <- heatmapBFMlow[order(heatmapBFMlow$log2FoldChange),]

# GO analysis files ----
# Export gene list for GO TERMS analysis
# virgin vs mated
DE_BFMall1 <- which(abs(sig_shres_BFvsM$log2FoldChange)>0.25)
DE_BFMall <- sig_shres_BFvsM[DE_BFMall1, ]
DE_BFMall <- DE_BFMall[order(DE_BFMall$log2FoldChange),]
#write.table(DE_BFMall, file = "/home/irenerossi/GOOGLE/Brainstorm_PhD/bulk-brain_results/Anopheles_bulk/R_results_bulk_mating/BMV_only/DE_lists_BMVonly_LFCplusminus025/DE_BFMall", row.names = F, col.names = T)



# * plots ----
# enrichment vaw curated based on VB GO, cytoscape ClueGO (including KEGG info) and manually (since some families are not well classified)
# !!! from dds<-dds.m rlog transformed
# heatmap mv ----
# plot trend of chosen genes of interest in the downloaded DEGs of mv
#library(reshape2)
#library(readODS)
# Define the list of genes
# the file has been ordered for each sub-cathegory to have the genes ordered by their FC (input order in VB)
mv_degs <- read_ods("/home/irenerossi/RAW_DATA/mating_bulk2021_trim2024/plotDEGs_bulk_mv.ods", sheet = 2, as_tibble = F, row_names = F)
mv_degs <- mv_degs %>%
  select(1,3,8:10)
gene_list <- mv_degs$`Gene ID`  # your actual gene list in the order you want
# in this case I put DEGs (upreg) in mv

# extract the assay matrix from the transformed object
mat.m <- assay(rld_m)
# Compute z-scores for each gene (row-wise scaling)
z_scores.m <- t(scale(t(mat.m)))

# extract data matrix
#expr_data <- assay(rld_m)
# or
expr_data <- z_scores.m


# Filter for genes and columns of interest
checked_data <- expr_data[rownames(expr_data) %in% gene_list, grep("_mated$|_virgin$", colnames(expr_data))]


# Check if checked_data has valid dimensions
if (nrow(checked_data) == 0 | ncol(checked_data) == 0) {
  stop("No matching genes or samples found.")
}

# Perform hierarchical clustering on genes
#distance_matrix <- dist(checked_data)            # Compute the distance matrix
#gene_clustering <- hclust(distance_matrix)        # Perform hierarchical clustering
#ordered_genes <- rownames(checked_data)[gene_clustering$order]  # Get gene order

# Reorder the data matrix based on clustering
#reordered_data <- checked_data[ordered_genes, ]

# Reorder the data matrix to match the order in gene_list
ordered_data <- checked_data[match(gene_list, rownames(checked_data)), ]

# get groups of virgin and mated samples
mated_data <- checked_data[, grep("_mated$", colnames(checked_data))]
virgin_data <- checked_data[, grep("_virgin$", colnames(checked_data))]

# to order samples based on their average expression ofo the different genes
# in the two groups, perform hierarchical clustering on each group separately
mated_clustering <- hclust(dist(t(mated_data)))
virgin_clustering <- hclust(dist(t(virgin_data)))

# combine in the clustering order
ordered_cols <- c(
  colnames(virgin_data)[virgin_clustering$order],
  colnames(mated_data)[mated_clustering$order]
)
# reorder the data matrix
ordered_data <- checked_data[, ordered_cols]
# reorder rows based on their log2FC (i.e. their order in gene_list)
ordered_data <- ordered_data[match(gene_list, rownames(ordered_data)), ]

# transform data for ggplot2
long_data <- as.data.frame(ordered_data) %>%
  rownames_to_column(var = "Gene") %>%
  melt(id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
# add sample group (Virgin/Mated) for annotation
long_data <- long_data %>%
  mutate(Group = ifelse(grepl("_virgin$", Sample), "Virgin", "Mated"))

# plot heatmap
library(viridis)
ggplot(long_data, aes(x = factor(Sample, levels = ordered_cols), y = factor(Gene, levels = rev(gene_list)), fill = Expression)) +
  geom_tile() +
  scale_fill_viridis(option = "B") +
  # alternative: 
  #scale_fill_gradient2(low = "white", high = "violetred4", mid = "yellow", 
  #                      midpoint = median(long_data$Expression), 
  #                      name = "Expression") +  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank(),
    legend.position = "right"
  ) +
  facet_grid(~Group, scales = "free_x", space = "free_x") +
  labs(title = "Gene Expression Heatmap (Virgin and Mated Groups)", fill = "Expression")


# heatmap bfm ----
# plot trend of chosen genes of interest in the downloaded DEGs of mv
#library(reshape2)
#library(readODS)
# Define the list of genes
bfm_degs_all <- read_ods("/home/irenerossi/RAW_DATA/mating_bulk2021_trim2024/plotDEGs_bulk_bfm.ods", sheet = 2, as_tibble = F, row_names = F)
# these genes have been obtained based on GO terms filtering (sheet1), and additional genes have been manually added


gene_list <- bfm_degs_all$gene_id  # your actual gene list in the order you want
# in this case I put DEGs (upreg) in bf


# extract the assay matrix from the transformed object
mat.bf <- assay(rld_bf)
# Compute z-scores for each gene (row-wise scaling)
z_scores.bf <- t(scale(t(mat.bf)))

# extract data matrix
#expr_data <- assay(rld_bf)
# or
expr_data <- z_scores.bf


# Filter for genes and columns of interest
checked_data <- expr_data[rownames(expr_data) %in% gene_list, grep("_mated$|_BF$", colnames(expr_data))]

# Check if checked_data has valid dimensions
if (nrow(checked_data) == 0 | ncol(checked_data) == 0) {
  stop("No matching genes or samples found.")
}

# Reorder the data matrix to match the order in gene_list
ordered_data <- checked_data[match(gene_list, rownames(checked_data)), ]

# get groups of mated and blood-fed samples
mated_data <- checked_data[, grep("_mated$", colnames(checked_data))]
bloodfed_data <- checked_data[, grep("_BF$", colnames(checked_data))]

# to order samples based on their average expression of the different genes
# in the two groups, perform hierarchical clustering on each group separately
mated_clustering <- hclust(dist(t(mated_data)))
bloodfed_clustering <- hclust(dist(t(bloodfed_data)))

# combine in the clustering order
ordered_cols <- c(
  colnames(mated_data)[mated_clustering$order],
  colnames(bloodfed_data)[bloodfed_clustering$order]
)
# reorder the data matrix
ordered_data <- checked_data[, ordered_cols]
# reorder rows based on their log2FC (i.e. their order in gene_list)
ordered_data <- ordered_data[match(gene_list, rownames(ordered_data)), ]

# transform data for ggplot2
long_data <- as.data.frame(ordered_data) %>%
  rownames_to_column(var = "Gene") %>%
  melt(id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
# add sample group (bloodfed/Mated) for annotation
long_data <- long_data %>%
  mutate(Group = ifelse(grepl("_BF$", Sample), "Bloodfed", "Mated"))

# plot heatmap
#library(viridis)
ggplot(long_data, aes(x = factor(Sample, levels = ordered_cols), y = factor(Gene, levels = rev(gene_list)), fill = Expression)) +
  geom_tile() +
  scale_fill_viridis(option = "B") +
  # alternative: 
  #scale_fill_gradient2(low = "white", high = "violetred4", mid = "yellow", 
  #                      midpoint = median(long_data$Expression), 
  #                      name = "Expression") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank(),
    legend.position = "right"
  ) +
  facet_grid(~Group, scales = "free_x", space = "free_x") +
  labs(title = "Gene Expression Heatmap (Bloodfed and Mated Groups)", fill = "Expression")
