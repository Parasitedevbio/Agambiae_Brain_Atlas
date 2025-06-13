library(readODS)
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
library(apeglm)
library(tidyverse)
library("dplyr")
library(ggplot2)
library(ggrepel)
library(limma)
library(pheatmap)

set.seed(1234)

# load reads and pheno objects ----
countData1 <- read.csv("/home/irenerossi/RAW_DATA/INFECTION_Ag_bulk_ALB0024_50-981620356/ALB0024_alignment/Ag_bulk_infections_raw_matrix.csv", header = T, sep = "\t")
metaData1 <- read_ods("/home/irenerossi/RAW_DATA/INFECTION_Ag_bulk_ALB0024_50-981620356/Ag_infections_phenofile.ods", sheet = 2)
# Convert wanted columns to factor
metaData1$batch <- as.factor(metaData1$batch)
metaData1$dissection_day <- as.factor(metaData1$dissection_day)


# AAAAAAAAAA::: CHECK BULK MATING 2024 SCRIPT FOR count and meta data objects creation and shape!
# rename stuff
colnames(countData1) <- gsub("_count", "", as.character(colnames(countData1)))
colnames(countData1) <- gsub("\\.", "_", as.character(colnames(countData1)))
# remove columns with variables that cannot be processed
metaData1 <- metaData1[-8] # i.e. intensity_sporozoites_+
# REMOVE OUTLIERS and sporo 2 and 4 ----
# GET AN OBJECT WITHOUT SPOROZOITE REPS 3 AND 4
countData <- countData1[-81:-96]
metaData <- metaData1[-80:-95,]
# check, it must be true!
colnames(countData)[2:80] == metaData$sample_name # it is

# check counts
# countData means
#mean_values <- colMeans(countData[sapply(countData, is.numeric)])
# countData sums
#sum_values <- colSums(countData[sapply(countData, is.numeric)])
#mean(countData$oocy2_6CT) # 0.6544049
# Mean of all other columns (all values combined)
#mean(mean_values[names(mean_values) != "oocy2_6CT"]) # 432.1687
#mean(sum_values[names(sum_values) != "oocy2_6CT"]) # 5970229
# you can check more, but for sample oocy2-6CT I already knew conc was off from Qubit measurements...
# oocyst2-6CT is an outlier, remove conts and metadata
countData <- countData %>% select(-c(oocy2_6CT))
metaData <- metaData[-which(metaData$sample_name == "oocy2_6CT"), ]
# check, it must be true!
colnames(countData)[2:79] == metaData$sample_name # it is

# create test vs ct objects ----
colnames(metaData)
# create an oocyst object
countData_oo <- countData %>% select(1, 2:39)
metaData_oo <- metaData[1:38,]
#tail(metaData_oo) # correct
# check, it must be true!:
#length(colnames(countData_oo)) # 39 samples
#colnames(countData_oo)[2:39] == metaData_oo$sample_name # it is
# create replicates objects
countData_oo1 <- countData_oo %>% select(1:21)
metaData_oo1 <- metaData_oo[1:20,]
countData_oo2 <- countData_oo %>% select(1, 22:39)
metaData_oo2 <- metaData_oo[21:38,]

# create a sporozoite object
countData_spz <- countData %>% select(1, 40:79)
metaData_spz <- metaData[39:78,]
#tail(metaData_spz) # correct
# check, it must be true!:
#length(colnames(countData_spz)) # 41 samples
#colnames(countData_spz)[2:41] == metaData_spz$sample_name # it is
# create replicates objects
countData_spz1 <- countData_spz %>% select(1:21)
metaData_spz1 <- metaData_spz[1:20,]
countData_spz2 <- countData_spz %>% select(1, 22:41)
metaData_spz2 <- metaData_spz[21:40,]


# * single experimental replicates ----
# batch 1 oo ----
dds_oo1 <- DESeqDataSetFromMatrix(countData = countData_oo1, 
                                  colData = metaData_oo1,
                                  design = ~ sample_type, tidy = TRUE) # warning: some variables in design formula are characters, converting to factors (it's ok)

# filter out by at least 3 samples with a count of 10 or higher
nrow(dds_oo1) # 13837
keep_oo1 <- rowSums(counts(dds_oo1) >= 10) >= 3
dds_oo1 <- dds_oo1[keep_oo1,]
nrow(dds_oo1) # 9223

# apply DESEQ function
dds_oo1 <- DESeq(dds_oo1)

# shrinked results conditionwise
shres_oo1 <- lfcShrink(dds_oo1, type="apeglm", lfcThreshold=0,
                       coef="sample_type_test_vs_control")

shres_oo1_ashr <- lfcShrink(dds_oo1, type="ashr", lfcThreshold=0,
                            coef="sample_type_test_vs_control")

# padj < 0.05: subset results to return genes ----
# oo2 ashr
sig_shres_oo1_ashr <- shres_oo1_ashr %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_oo1_ashr # look at FC # jut one gene, not present in res_oo_solo2

# padj < 0.05: subset results to return genes ----
# oo1 apeglm
sig_shres_oo1 <- shres_oo1 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_oo1 # look at FC # same gene as for not shrinked results
# ashr shrinkage gives the same significant gene

# batch 2 oo ----
dds_oo2 <- DESeqDataSetFromMatrix(countData = countData_oo2, 
                                  colData = metaData_oo2,
                                  design = ~ sample_type, tidy = TRUE) # warning: some variables in design formula are characters, converting to factors (it's ok)

# filter out by at least 3 samples with a count of 10 or higher
nrow(dds_oo2) # 13837
keep_oo2 <- rowSums(counts(dds_oo2) >= 10) >= 3
dds_oo2 <- dds_oo2[keep_oo2,]
nrow(dds_oo2) # 9083

# apply DESEQ function
dds_oo2 <- DESeq(dds_oo2)

# shrinked results conditionwise
shres_oo2 <- lfcShrink(dds_oo2, type="apeglm", lfcThreshold=0, # warning message suggesting possible overparametrization, so check with ashr
                       coef="sample_type_test_vs_control")

shres_oo2_ashr <- lfcShrink(dds_oo2, type="ashr", lfcThreshold=0,
                       coef="sample_type_test_vs_control")

# padj < 0.05: subset results to return genes ----
# oo2 ashr
sig_shres_oo2_ashr <- shres_oo2_ashr %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_oo2_ashr # look at FC # jut one gene, not present in res_oo_solo2

# oo2 apeglm
sig_shres_oo2 <- shres_oo2 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_oo2 # look at FC # jut one gene, not present in res_oo_solo2
# ashr gives the same significant gene


# batch 1 spz ----
dds_spz1 <- DESeqDataSetFromMatrix(countData = countData_spz1, 
                                   colData = metaData_spz1,
                                   design = ~ sample_type, tidy = TRUE) # warning: some variables in design formula are characters, converting to factors (it's ok)

# filter out by at least 3 samples with a count of 10 or higher
nrow(dds_spz1) # 13837
keep_spz1 <- rowSums(counts(dds_spz1) >= 10) >= 3
dds_spz1 <- dds_spz1[keep_spz1,]
nrow(dds_spz1) # 8889

# apply DESEQ function
dds_spz1 <- DESeq(dds_spz1)

# shrinked results conditionwise
shres_spz1 <- lfcShrink(dds_spz1, type="apeglm", lfcThreshold=0,
                        coef="sample_type_test_vs_control")

shres_spz1_ashr <- lfcShrink(dds_spz1, type="ashr", lfcThreshold=0,
                            coef="sample_type_test_vs_control")

# padj < 0.05: subset results to return genes ----
# spz1 ashr
sig_shres_spz1_ashr <- shres_spz1_ashr %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_spz1_ashr # look at FC # one gene, different form the two of res_spz_solo1

# spz1 apeglm
sig_shres_spz1 <- shres_spz1 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_spz1 # look at FC # one gene, different form the two of res_spz_solo1
# ashr gives the same significant gene

# batch 2 spz ----
dds_spz2 <- DESeqDataSetFromMatrix(countData = countData_spz2, 
                                   colData = metaData_spz2,
                                   design = ~ sample_type, tidy = TRUE) # warning: some variables in design formula are characters, converting to factors (it's ok)
str(metaData_spz2)
# filter out by at least 3 samples with a count of 10 or higher
nrow(dds_spz2) # 13837
keep_spz2 <- rowSums(counts(dds_spz2) >= 10) >= 3
dds_spz2 <- dds_spz2[keep_spz2,]
nrow(dds_spz2) # 9223

# apply DESEQ function
dds_spz2 <- DESeq(dds_spz2)

# shrinked results conditionwise
shres_spz2 <- lfcShrink(dds_spz2, type="apeglm", lfcThreshold=0,
                        coef="sample_type_test_vs_control")

shres_spz2_ashr <- lfcShrink(dds_spz2, type="ashr", lfcThreshold=0,
                             coef="sample_type_test_vs_control")

# padj < 0.05: subset results to return genes ----
# spz2 ashr
sig_shres_spz2_ashr <- shres_spz2_ashr %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_spz2_ashr # look at FC # no gene, while in res_spz_solo2 there was one

# spz2 apeglm
sig_shres_spz2 <- shres_spz2 %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_spz2 # look at FC # no gene, while in res_spz_solo2 there was one
# same results for ashr


# * infection stage replicates ----
# oocyst replicates ----
dds_oo <- DESeqDataSetFromMatrix(countData = countData_oo, 
                                 colData = metaData_oo,
                                 design = ~ sample_type, tidy = TRUE)

# filter out by at least 3 samples with a count of 10 or higher
nrow(dds_oo) # 13837
keep_oo <- rowSums(counts(dds_oo) >= 10) >= 3
dds_oo <- dds_oo[keep_oo,]
nrow(dds_oo) # 9492

# apply DESEQ function
dds_oo <- DESeq(dds_oo)
# shrinked results conditionwise
shres_oo <- lfcShrink(dds_oo, type="apeglm", lfcThreshold=0,
                      coef="sample_type_test_vs_control")

shres_oo_ashr <- lfcShrink(dds_oo, type="ashr", lfcThreshold=0,
                      coef="sample_type_test_vs_control")

# padj < 0.05: subset results to return genes ----
# oo-solo ashr
sig_shres_oo_ashr <- shres_oo_ashr %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_oo_ashr # look at FC # 3 genes 

# oo-solo apeglm
sig_shres_oo <- shres_oo %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_oo # look at FC # 3 genes
# ashr gives the same significant genes

# sporozoite replicates ----
dds_spz <- DESeqDataSetFromMatrix(countData = countData_spz, 
                                  colData = metaData_spz,
                                  design = ~ sample_type, tidy = TRUE)

# filter out by at least 3 samples with a count of 10 or higher
nrow(dds_spz) # 13837
keep_spz <- rowSums(counts(dds_spz) >= 10) >= 3
dds_spz <- dds_spz[keep_spz,]
nrow(dds_spz) # 9218

# apply DESEQ function
dds_spz <- DESeq(dds_spz)
# check whether the operations were successful
resultsNames(dds_spz)
dds_spz@design

# shrinked results conditionwise
shres_spz_solo <- lfcShrink(dds_spz, type="apeglm", lfcThreshold=0,
                            coef="sample_type_test_vs_control")

shres_spz_solo_ashr <- lfcShrink(dds_spz, type="ashr", lfcThreshold=0,
                            coef="sample_type_test_vs_control")

# padj < 0.05: subset results to return genes ----
# spz-solo ashr
sig_shres_spz_ashr <- shres_spz_solo_ashr %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_spz_ashr # look at FC # zero genes

# spz-solo apeglm
sig_shres_spz <- shres_spz_solo %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_spz # look at FC # zero genes
# ashr also gives zero sig gene

# no DE gene


# * Preliminary inspection of aging effect: build DEseq2 object for all the conditions ----
# make DESeq matrix
dds <- DESeqDataSetFromMatrix(countData = countData, 
                              colData = metaData,
                              design = ~ sample_type, tidy = TRUE)

# set reference level
dds$sample_type <- relevel(dds$sample_type, ref = "control")
# add phenofile cathegories interactions ("design" argument in the above function "DESeqDataSetFromMatrix")
dds$group <- factor(paste0(dds$sample_type, dds$infection_stage))
dds$group <- relevel(dds$group, ref = "controloocyst")
design(dds) <- ~ group
str(dds)
# this way later the result levels would be all contrasted to oocyst controls.

# filter out by at least 3 samples with a count of 10 or higher ----
nrow(dds) # 13837
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
nrow(dds) # 9725

# inf not inf
# oo inf not inf (subs)
# spz inf not inf (subs)

# apply DESEQ function ----
dds <- DESeq(dds)
# check DSeq2 object structure
#colData(dds)
#rowRanges(dds)
#dds$sample_name
#head(countData$feature)

# check whether the operations were successful
resultsNames(dds)
#dds$group

# Use Case:
# - apeglm: ideal for small datasets or when you want highly reliable LFC estimates. Tends to shrink LFCs toward zero more aggressively than ashr, especially for genes with low confidence.
# - ashr: better customization contrasts. Suitable for larger datasets or when custom contrasts are needed. Moderate shrinkage, less aggressive than apeglm.
# So I went on analyzing just the oocyst samples versus their controls in this setup, that give two significant genes.
shres_ooall <- lfcShrink(dds, type="apeglm", lfcThreshold=0, 
                      coef=c("group_testoocyst_vs_controloocyst"))
# apeglm methos was not used because it was rising this error:
# In nbinomGLM(x = x, Y = YNZ, size = size, weights = weightsNZ, offset = offsetNZ, : the line search routine failed, unable to sufficiently decrease the function value
# usually this rises for extremely low counts - already filtered out - or problematic genes with zero logFC - that were not found in my data.
# Basically it seems that this warning is raised for no reason, so apeglm results are kept for omogeneity with the previous analyses.
# In addition, ashr method is also applied later without raising any error, giving the same DEGs (2 antimicrobial peptides),
# so I will consider these results for homogeneity.
# OO
sig_shres_ooall <- shres_ooall %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_ooall # look at FC

# When setting: dds$group <- relevel(dds$group, ref = "controlsporozoite")
# after the creation of dds$group this would still give 0 genes in the sporozoite vs their controls lfcShrink (apeglm) results.




# ashr shrinked results conditionwise (ashr is used because apeglm cannot handle complex setups) ----
#shres_oo <- lfcShrink(dds, type="ashr", lfcThreshold=0, 
#                      contrast=c("group", "testoocyst", "controloocyst"))
#shres_spz <- lfcShrink(dds, type="ashr", lfcThreshold=0, 
#                       contrast=c("group", "testsporozoite", "controlsporozoite"))
shres_spzVSoo <- lfcShrink(dds, type="ashr", lfcThreshold=0, 
                           contrast=c("group", "testsporozoite", "testoocyst"))
shres_cts <- lfcShrink(dds, type="ashr", lfcThreshold=0, 
                       contrast=c("group", "controlsporozoite", "controloocyst"))


# padj < 0.05: subset results to return genes ----
# OO
#sig_shres_oo <- shres_oo %>%
#  data.frame() %>%
#  rownames_to_column(var="gene") %>% 
#  as_tibble() %>% 
#  filter(padj < 0.05)
#sig_shres_oo # look at FC # 2 genes

# SPZ
#sig_shres_spz <- shres_spz %>%
#  data.frame() %>%
#  rownames_to_column(var="gene") %>% 
#  as_tibble() %>% 
#  filter(padj < 0.05)
#sig_shres_spz # look at FC # no DEG

# SPZ vs OO
sig_shres_spzVSoo <- shres_spzVSoo %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_spzVSoo # look at FC

# CTs (spz_ct vs oo_ct)
sig_shres_cts <- shres_cts %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_cts # look at FC

# drop lines that don't contain gene information
tail(sig_shres_spzVSoo)
tail(sig_shres_cts)
sig_shres_spzVSoo <- head(sig_shres_spzVSoo,-4)
sig_shres_cts <- head(sig_shres_cts,-2)


# PCA rld ----
# PC1,2 rld
rld <- rlog(dds, blind = FALSE)
#plotPCA(rld, intgroup=c("sample_type", "infection_stage")) + stat_ellipse()
#plotPCA(rld, intgroup = "sample_condition_subsetting") + stat_ellipse()

#plotPCA(rld, intgroup = "sample_type") + stat_ellipse()# + geom_label(aes(label = colnames(rld)))

# get PC1,3
alternative.pca1_3 <- function (object, ...) 
{
  .local <- function (object, intgroup = "condition", ntop = 500, 
                      returnData = FALSE) 
  {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                       length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
      stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                                 drop = FALSE])
    group <- if (length(intgroup) > 1) {
      factor(apply(intgroup.df, 1, paste, collapse = ":"))
    }
    else {
      colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 3], group = group, 
                    intgroup.df, name = colnames(object))
    if (returnData) {
      attr(d, "percentVar") <- percentVar[1:3]
      return(d)
    }
    ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
      geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
                                                          100), "% variance")) + ylab(paste0("PC3: ", round(percentVar[3] * 
                                                                                                              100), "% variance")) + coord_fixed() +
      stat_ellipse()
  }
  .local(object, ...)
}


# * for the paper: ggplot PCAs ----
#pcaData <- plotPCA(rld, intgroup= "sample_type", returnData=TRUE)
#percentVar <- round(100 * attr(pcaData, "percentVar"))
# PC1,2
#ggplot(pcaData, aes(PC1, PC2, col = sample_type, shape = sample_type)) +
#  scale_color_manual(values = c("magenta4", "#DC3220")) +
#  geom_point(size=3) +
#  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#  labs(colour="infection state") + 
#  scale_color_manual(labels = c("not infected", "infected"), values = c("magenta4", "#DC3220")) +
#  coord_fixed()
# PC1,3
pcaData_inf <- alternative.pca1_3(rld, intgroup= "sample_type", returnData=TRUE)
percentVar_inf <- round(100 * attr(pcaData_inf, "percentVar"))
ggplot(pcaData_inf, aes(PC1, PC2, col = sample_type, shape = sample_type)) +
  geom_point(size=2, shape=16) +
  xlab(paste0("PC1: ",percentVar_inf[1],"% variance")) +
  ylab(paste0("PC3: ",percentVar_inf[2],"% variance")) +
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
  labs(colour="infection state") +
  xlim(-30, 30) +
  ylim(-30, 30) +
  scale_color_manual(labels = c("not infected", "infected"), values = c("magenta4", "#DC3220")) +
  coord_fixed()


#plotPCA(rld, intgroup="infection_stage") + stat_ellipse()
# ggplot it
# PC1,2
#pcaData <- plotPCA(rld, intgroup= "infection_stage", returnData=TRUE)
#percentVar <- round(100 * attr(pcaData, "percentVar"))
#ggplot(pcaData, aes(PC1, PC2, col = infection_stage, shape = infection_stage)) +
#  scale_color_manual(values = c("orange2", "royalblue")) +
#  geom_point(size=3) +
#  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
#  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
#  labs(colour="age") + 
#  scale_color_manual(labels = c("14 days old", "22 days old"), values = c("orange2", "royalblue")) +
#  coord_fixed()
# PC1,3
pcaData_age <- alternative.pca1_3(rld, intgroup= "infection_stage", returnData=TRUE)
percentVar_age <- round(100 * attr(pcaData_age, "percentVar"))
ggplot(pcaData_age, aes(PC1, PC2, col = infection_stage, shape = infection_stage)) +
  geom_point(size=2, shape=16) +
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
  xlab(paste0("PC1: ",percentVar_age[1],"% variance")) +
  ylab(paste0("PC3: ",percentVar_age[2],"% variance")) +
  labs(colour="age") +
  xlim(-30, 30) +
  ylim(-30, 30) +
  scale_color_manual(labels = c("14 days old", "22 days old"), values = c("orange2", "royalblue")) +
  coord_fixed()

#z <- plotPCA(rld, intgroup = "sample_name")
#z + geom_label(aes(label = name)) # oocy2_3 seems an outlier but I don't think it is


# * DEA AGING ALL SAMPLES: ----
# build DEseq2 object for all the conditions
# make DESeq matrix
dds.aging <- DESeqDataSetFromMatrix(countData = countData, 
                                    colData = metaData,
                                    design = ~ infection_stage, tidy = TRUE)

# filter out by at least 3 samples with a count of 10 or higher ----
nrow(dds.aging) # 13837
keep <- rowSums(counts(dds.aging) >= 10) >= 3
dds.aging <- dds.aging[keep,]
nrow(dds.aging) # 9725

# apply DESEQ function ----
dds.aging <- DESeq(dds.aging)
#dds.aging1 <- DESeq(dds.aging, test="LRT", reduced=~1)
# check DSeq2 object structure
#colData(dds.aging)
#rowRanges(dds.aging)
#dds.aging$sample_name
#head(countData$feature)

# check whether the operations were successful (same for dds.aging and dds.aging1)
resultsNames(dds.aging)

# shrinked results conditionwise
shres_aging <- lfcShrink(dds.aging, type="apeglm", lfcThreshold=0,
                         coef=c("infection_stage_sporozoite_vs_oocyst"))

# old VS young
sig_shres_aging <- shres_aging %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble() %>% 
  filter(padj < 0.05)
sig_shres_aging # look at FC

# up
up_sig_shres_aging <- sig_shres_aging %>%
  data.frame() %>%
  as_tibble() %>% 
  filter(log2FoldChange > 0)
up_sig_shres_aging # look at FC

# down
down_sig_shres_aging <- sig_shres_aging %>%
  data.frame() %>%
  as_tibble() %>% 
  filter(log2FoldChange < 0)
down_sig_shres_aging # look at FC

# percentages
length(sig_shres_aging$gene)*100/length(rownames(dds.aging)) # 45.2545
length(up_sig_shres_aging$gene)*100/length(rownames(dds.aging)) # 22.29306
length(down_sig_shres_aging$gene)*100/length(rownames(dds.aging)) # 22.96144

# download files
#write.table(sig_shres_aging, file = "/home/irenerossi/drive_ird/Brainstorm_PhD/INFECTION_EXPERIMENTS/INFECTION_dataanalysis_RESULTS/AGING_ALL_results/aging_degs_lists_and_annotation/degs_aging", row.names = F, col.names = T)
#write.table(up_sig_shres_aging, file = "/home/irenerossi/drive_ird/Brainstorm_PhD/INFECTION_EXPERIMENTS/INFECTION_dataanalysis_RESULTS/AGING_ALL_results/aging_degs_lists_and_annotation,degs_up_aging", row.names = F, col.names = T)
#write.table(down_sig_shres_aging, file = "/home/irenerossi/drive_ird/Brainstorm_PhD/INFECTION_EXPERIMENTS/INFECTION_dataanalysis_RESULTS/AGING_ALL_results/aging_degs_lists_and_annotation/degs_down_aging", row.names = F, col.names = T)


# * FOR THE PAPER check overlapping aging and infection lists ----
# merge gene expression and annotation info and sort by FC
# VectorBase gene annotation added to degs_aging genes
#all_degs_aging_annotat <- read.csv("/home/irenerossi/drive_ird/Brainstorm_PhD/INFECTION_EXPERIMENTS/INFECTION_dataanalysis_RESULTS/AGING_ALL_results/aging_degs_lists_and_annotation/all_degs_aging_info_VB.csv", header = T, sep = ",")
# keep wanted columns
#all_degs_aging_annotat <- all_degs_aging_annotat[, c(1, 3, 5)]
# merge by gene
#all_aging_info <- sig_shres_aging %>%
#  left_join(all_degs_aging_annotat, by = c("gene" = "Gene.ID"))
# sort FCs
#all_aging_info <- all_aging_info %>%
#  arrange(log2FoldChange)
#download
#write.table(all_aging_info, file = "/home/irenerossi/drive_ird/Brainstorm_PhD/INFECTION_EXPERIMENTS/INFECTION_dataanalysis_RESULTS/AGING_ALL_results/FINAL_aging_genes", row.names = F, col.names = T)
# upload final table (not aligned and ambiguous reads terms (3 terms) values removed)
all_aging_info <- read_ods("/home/irenerossi/drive_ird/Brainstorm_PhD/INFECTION_EXPERIMENTS/INFECTION_dataanalysis_RESULTS/AGING_ALL_results/FINAL_aging_genes.ods")

# top FC and FDR genes to plot ----
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)
all_aging_info$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
all_aging_info$diffexpressed[all_aging_info$log2FoldChange > 1.5 & all_aging_info$padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
all_aging_info$diffexpressed[all_aging_info$log2FoldChange < -1.5 & all_aging_info$padj < 0.05] <- "DOWN"
# Explore a bit
head(all_aging_info[order(all_aging_info$padj) & all_aging_info$diffexpressed == 'DOWN', ])

# which genes are there
diffexpressed_top_genes <- all_aging_info %>% 
  filter(diffexpressed != "NO")

# get genes to plot
# Filter the rows in diffexpressed_top_genes where p-value is < 10^20
genes_to_plot_vp <- diffexpressed_top_genes[diffexpressed_top_genes$padj < 10^(-19), ]
# inspect
genes_to_plot_vp$Product.Description
# * check if they are all present in all_aging_info ----
#length(genes_to_plot_vp$gene)
plot_aging_strength <- Reduce(intersect, list(all_aging_info$gene, genes_to_plot_vp$gene))
identical(sort(genes_to_plot_vp$gene), sort(plot_aging_strength)) # TRUE


# heatmap top DEGs ----
# plot trend of chosen genes of interest in the downloaded DEGs of aging
# sometimes VB messes up with the gene IDs, recognition, that's why the next step is important:
library(reshape2)

# Define the list of genes
gene_list <- genes_to_plot_vp$gene
# in this case I put DEGs interesting in aging

# check if they are all present in all_aging_info ----
#length(gene_list)
plot_aging_strength2 <- Reduce(intersect, list(all_aging_info$gene, gene_list))
identical(sort(gene_list), sort(plot_aging_strength2)) # TRUE

# extract the assay matrix from the transformed object
mat.age <- assay(rlog(dds.aging, blind = FALSE))
# Compute z-scores for each gene (row-wise scaling)
z_scores.age <- t(scale(t(mat.age)))

# extract data matrix
expr_data <- z_scores.age

# Filter for genes and columns of interest
checked_data <- expr_data[rownames(expr_data) %in% gene_list, grep("^oocy|^sporo", colnames(expr_data))]

# Check if checked_data has valid dimensions
if (nrow(checked_data) == 0 | ncol(checked_data) == 0) {
  stop("No matching genes or samples found.")
}

# Reorder the data matrix to match the order in gene_list
ordered_data <- checked_data[match(gene_list, rownames(checked_data)), ]

# get groups of virgin and mated samples
oo_data <- checked_data[, grep("^oocy", colnames(checked_data))]
spz_data <- checked_data[, grep("^sporo", colnames(checked_data))]

# to order samples based on their average expression ofo the different genes
# in the two groups, perform hierarchical clustering on each group separately
oo_clustering <- hclust(dist(t(oo_data)))
spz_clustering <- hclust(dist(t(spz_data)))

# combine in the clustering order
ordered_cols <- c(
  colnames(oo_data)[oo_clustering$order],
  colnames(spz_data)[spz_clustering$order]
)
# reorder the data matrix
ordered_data <- checked_data[, ordered_cols]
# reorder rows based on their log2FC (i.e. their order in gene_list)
ordered_data <- ordered_data[match(gene_list, rownames(ordered_data)), ]

# transform data for ggplot2
long_data <- as.data.frame(ordered_data) %>%
  rownames_to_column(var = "Gene") %>%
  melt(id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
# Create a dataframe for sample type
sample_type_df <- data.frame(
  Sample = colnames(dds.aging),
  sample_type = dds.aging$sample_type
)
# Create a dataframe for sample type
sample_type_df <- data.frame(
  Sample = colnames(dds.aging),
  sample_type = dds.aging$sample_type,
  Group = ifelse(grepl("^oocy", colnames(dds.aging)), "14 d.o.", "22 d.o.")
)

# Add sample_type to long_data
long_data <- merge(long_data, sample_type_df, by = "Sample")
# add sample group (young/old) for annotation
long_data <- long_data %>%
  mutate(Group = ifelse(grepl("^oocy", Sample), "14 d.o.", "22 d.o."))


# plot heatmap
library(viridis)
# Required to use two different fill scales
library(ggnewscale)
# ggplot with gene names and sample bar
ggplot() +
  # Sample type color bar
  geom_tile(data = distinct(sample_type_df),
            aes(x = factor(Sample, levels = ordered_cols), 
                y = "Sample Type", 
                fill = sample_type)) +
  scale_fill_manual(name = "Sample Type",
                    values = c("control" = "magenta4", "test" = "#DC3220")) +
  ggnewscale::new_scale_fill() +
  
  # Expression heatmap
  geom_tile(data = long_data,
            aes(x = factor(Sample, levels = ordered_cols), 
                y = factor(Gene, levels = rev(gene_list)), 
                fill = Expression)) +
  scale_fill_viridis(option = "B", name = "Expression") +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 9),
    axis.title = element_blank(),
    legend.position = "right"
  ) +
  facet_grid(~Group, scales = "free_x", space = "free_x") +
  labs(title = "Gene Expression Heatmap (14 d.o. and 22 d.o. Groups)") +
  scale_y_discrete(labels = c("Sample Type" = "", gene_list))  # Hide label on the annotation bar
# put gene_name_map instead of gene_list to have gene names


# * distance plot (OBJECT CREATION SAME TO SOME PREVIOUS CODE) ----
# build DEseq2 object for all the conditions
# make DESeq matrix
#dds.distance <- DESeqDataSetFromMatrix(countData = countData, 
#                                    colData = metaData,
#                                    design = ~ 1, tidy = T)
# add phenofile cathegories interactions ("design" argument in the above function "DESeqDataSetFromMatrix")
#dds.distance$group <- factor(paste0(dds.distance$sample_type, dds.distance$infection_stage, dds.distance$batch))
#dds.distance$group <- relevel(dds.distance$group, ref = "controloocyst1")
#design(dds.distance) <- ~ group
#str(dds.distance)
#unique(dds.distance$group) # see the values
#dds.distance$group <- gsub("3", "1", dds.distance$group)
#dds.distance$group <- gsub("4", "2", dds.distance$group)
# Convert group to factor AFTER string substitution
#dds.distance$group <- factor(dds.distance$group)

# filter out by at least 3 samples with a count of 10 or higher ----
#nrow(dds.distance) # 13837
#keep <- rowSums(counts(dds.distance) >= 10) >= 3
#dds.distance <- dds.distance[keep,]
#nrow(dds.distance) # 9725

# apply DESEQ function
#dds.distance <- DESeq(dds.distance)


#unique(dds.distance$group) # see the final values
# Define your conditions
#conditions <- rev(c("controloocyst", "controlsporozoite", "testoocyst", "testsporozoite"))
#replicates <- 1:2  # assuming you have replicates 1 and 2

# Initialize a list to store all results
#results_list <- list()

# Generate all unique non-self pairs
#for (cond1 in conditions) {
#  for (rep1 in replicates) {
#    for (cond2 in conditions) {
#      for (rep2 in replicates) {
# Skip self-comparisons
#        if (!(cond1 == cond2 && rep1 == rep2)) {
# Create condition names
#          name1 <- paste0(cond1, rep1)
#          name2 <- paste0(cond2, rep2)

# Create a meaningful name for the result
#          result_name <- paste0(gsub("control", "c", gsub("test", "", cond2)),
#                                rep2,
#                                gsub("control", "c", gsub("test", "", cond1)),
#                                rep1)

# Perform the shrinkage
#          cat("Processing contrast:", name2, "vs", name1, "...\n")
#          res <- lfcShrink(dds.distance, type="ashr", lfcThreshold=0,
#                           contrast=c("group", name2, name1))

# Store the result
#          results_list[[result_name]] <- res
#        }
#      }
#    }
#  }
#}
#length(results_list)
#names(results_list)

# Apply filtering and formatting to all results
#sig_results <- lapply(results_list, function(res) {
#  res %>% 
#    data.frame() %>%
#    rownames_to_column(var = "gene") %>% 
#    as_tibble() %>% 
#    filter(padj < 0.05)
#})

# Access individual results using:
#sig_results$csporozoite1sporozoite1
# ... and so on for other comparisons

# To check number of significant genes per comparison:
#sapply(sig_results, nrow)

# Define all conditions (same order as your matrix)
#conditions <- c("coocyst1", "coocyst2", 
#                "csporozoite1", "csporozoite2",
#                "oocyst1", "oocyst2",
#                "sporozoite1", "sporozoite2")

# initialize distance matrix
#distance_matrix <- matrix(0, nrow=length(conditions), ncol=length(conditions),
#                          dimnames=list(conditions, conditions))

# Fill the matrix in a loop
#for (i in 1:length(conditions)) {
#  for (j in 1:length(conditions)) {
#    if (i == j) {
#      distance_matrix[i,j] <- 0  # Diagonal is zero
#    } else {
# Construct the result object name
#      result_name <- paste0(conditions[j], conditions[i])

# Check if this comparison exists in sigresults
#      if (exists(result_name, where = sig_results)) {
#        res <- sig_results[[result_name]]
#        distance_matrix[i,j] <- length(res$gene)  # Or your actual DEG count metric
#      } else {
# Try reverse comparison if this one doesn't exist
#        reverse_name <- paste0(conditions[i], conditions[j])
#        if (exists(reverse_name, where = sig_results)) {
#          res <- sig_results[[reverse_name]]
#          distance_matrix[i,j] <- length(res$gene)  # Or your actual DEG count metric
#        } else {
#          distance_matrix[i,j] <- NA  # No comparison available
#        }
#      }
#    }
#  }
#}
# View the completed distance matrix
#print(distance_matrix)

# Remove bottom part of the matrix (mirrored)
#distance_matrix[lower.tri(distance_matrix)] <- 0

# save matrix 
#sep.reps.distance <- distance_matrix
#write.table(sep.reps.distance, file = "/home/irenerossi/drive_ird/Brainstorm_PhD/INFECTION_EXPERIMENTS/INFECTION_dataanalysis_RESULTS/distance_matrix_sep_reps_inf.csv", row.names = T, col.names = T)
# upload matrix
sep.reps.distance <- read.csv("/home/irenerossi/drive_ird/Brainstorm_PhD/INFECTION_EXPERIMENTS/INFECTION_dataanalysis_RESULTS/distance_matrix_sep_reps_inf.csv", row.names = 1, header = T, sep = " ")


# Melt the matrix
library(reshape2)
melted_matrix <- as.data.frame(as.table(as.matrix(sep.reps.distance))) # for a distance matrix with empty diagonal

# 5. Force column names to be exact
names(melted_matrix) <- c("Var1", "Var2", "Value")

# 6. Plot with guaranteed correct column names
library(ggplot2)

ggplot(melted_matrix, aes(x=Var1, y=Var2, fill=Value)) +
  geom_tile(color="white") +
  scale_fill_viridis_c(name="# DEGs", option="plasma", na.value="white") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=10),
    axis.text.y = element_text(size=10)
  ) +
  labs(x="", y="") +
  coord_fixed()


# DO IT FOR MERGED SAMPLE REPLICATES (see code - much - above, data obtained the same way, just re-run here to have all the lines one after the other)
# build DEseq2 object for all the conditions
# make DESeq matrix
#dds.distance <- DESeqDataSetFromMatrix(countData = countData, 
#                              colData = metaData,
#                              design = ~ sample_type, tidy = TRUE)
# set reference level
#dds.distance$sample_type <- relevel(dds.distance$sample_type, ref = "control")
# add phenofile cathegories interactions ("design" argument in the above function "DESeqDataSetFromMatrix")
#dds.distance$group <- factor(paste0(dds.distance$sample_type, dds.distance$infection_stage))
#dds.distance$group <- relevel(dds.distance$group, ref = "controloocyst")
#design(dds.distance) <- ~ group
#str(dds.distance)
# this way later the result levels would be all contrasted to oocyst controls.

# filter out by at least 3 samples with a count of 10 or higher ----
#nrow(dds.distance) # 13837
#keep <- rowSums(counts(dds.distance) >= 10) >= 3
#dds.distance <- dds.distance[keep,]
#nrow(dds.distance) # 9725

# apply DESEQ function ----
#dds.distance <- DESeq(dds.distance)

# ashr shrinked results conditionwise (ashr is used because apeglm cannot handle complex setups) ----
#shres_oo <- lfcShrink(dds.distance, type="ashr", lfcThreshold=0, 
#                      contrast=c("group", "testoocyst", "controloocyst"))
#shres_spz <- lfcShrink(dds.distance, type="ashr", lfcThreshold=0, 
#                       contrast=c("group", "testsporozoite", "controlsporozoite"))
#shres_spzVSoo <- lfcShrink(dds.distance, type="ashr", lfcThreshold=0, 
#                           contrast=c("group", "testsporozoite", "testoocyst"))
#shres_cts <- lfcShrink(dds.distance, type="ashr", lfcThreshold=0, 
#                       contrast=c("group", "controlsporozoite", "controloocyst"))

# padj < 0.05: subset results to return genes ----
# OO
#sig_shres_oo <- shres_oo %>%
#  data.frame() %>%
#  rownames_to_column(var="gene") %>% 
#  as_tibble() %>% 
#  filter(padj < 0.05)
#sig_shres_oo # look at FC # 2 genes

# SPZ
#sig_shres_spz <- shres_spz %>%
#  data.frame() %>%
#  rownames_to_column(var="gene") %>% 
#  as_tibble() %>% 
#  filter(padj < 0.05)
#sig_shres_spz # look at FC # no DEG

# SPZ vs OO
#sig_shres_spzVSoo <- shres_spzVSoo %>%
#  data.frame() %>%
#  rownames_to_column(var="gene") %>% 
#  as_tibble() %>% 
#  filter(padj < 0.05)
#sig_shres_spzVSoo
# check and fix
#tail(sig_shres_spzVSoo)
#sig_shres_spzVSoo <- head(sig_shres_spzVSoo,-4)
#tail(sig_shres_spzVSoo)

# CTs (spz_ct vs oo_ct)
#sig_shres_cts <- shres_cts %>%
#  data.frame() %>%
#  rownames_to_column(var="gene") %>% 
#  as_tibble() %>% 
#  filter(padj < 0.05)
#sig_shres_cts
# check and fix
#tail(sig_shres_cts)
#sig_shres_cts <- head(sig_shres_cts,-2)
#tail(sig_shres_cts)

# Create a distance matrix based on the absolute difference in list lengths
#distance_matrix <- matrix(0, nrow = 4, ncol = 4)
#rownames(distance_matrix) <- c("ct_oo", "oo", "ct_spz", "spz")
#colnames(distance_matrix) <- c("ct_oo", "oo", "ct_spz", "spz")
# Fill the matrix with the absolute difference in lengths
#distance_matrix[1, 2] <- length(sig_shres_oo$gene)
#distance_matrix[1, 3] <- length(sig_shres_cts$gene)
#distance_matrix[1, 4] <- 0

#distance_matrix[2, 3] <- 0
#distance_matrix[2, 4] <- length(sig_shres_spzVSoo$gene)

#distance_matrix[3, 4] <- length(sig_shres_spz$gene)

# save matrix 
#aggr.reps.distance <- distance_matrix
#write.table(aggr.reps.distance, file = "/home/irenerossi/drive_ird/Brainstorm_PhD/INFECTION_EXPERIMENTS/INFECTION_dataanalysis_RESULTS/distance_matrix_aggr_reps_inf.csv", row.names = T, col.names = T)
# upload matrix
aggr.reps.distance <- read.csv("/home/irenerossi/drive_ird/Brainstorm_PhD/INFECTION_EXPERIMENTS/INFECTION_dataanalysis_RESULTS/distance_matrix_aggr_reps_inf.csv", row.names = 1, header = T, sep = " ")

# Now let's plot the distance matrix
#library(ggplot2)
#library(reshape2)

# Convert the matrix to a long format for ggplot
melted_matrix <- as.data.frame(as.table(as.matrix(aggr.reps.distance))) # for a distance matrix with empty diagonal
names(melted_matrix) <- c("Row", "Column", "Value")

# Create the heatmap
# Reverse axes in the heatmap plot
ggplot(melted_matrix, aes(x = Row, y = Column, fill = Value)) +  # Swap x and y
  geom_tile(color = "white") +
  scale_fill_viridis_c(limits = c(0, 3000), 
                       breaks = seq(0, 3000, by = 500),
                       option = "viridis") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  labs(
    fill = "",
    x = "",
    y = ""
  ) +
  coord_fixed()  # Makes the tiles square


# * aging immune DEGs family enrichment ----
# p-values are calculated using the upper-tail (right-tail) probability 
# P(X ≥ a). The FDR is adjusted across all families using the Benjamini-Hochberg method.

# * all DEGs
# Parameters
N <- 13845   # Total genes in Ag PEST genome
n <- 4398    # Number of DEGs in aging analysis
# Gene families: m = number in gene family, a = observed overlap
# families taken form Christophides 2002 (doi: 10.1126/science.1077136), some families mentioned by Christophides
# are absent because their gene name could not be found in the DEGs or in both degs and total gene list
# the counts have been done manually due to impossibility to always assign univocal strings to the gene names of a given family
families <- c("PGRP", "TEP", "GNBP", "SCR",
              "CTL", "GALE", "CLIP", # "FBN",
              "SRPN", "IAP", "TOLL + cactus", # "MyD88", "tube", "pelle",
              "REL", "Imd",
              "STAT", "PPO", "DEF", "CEC",
              "CASP")
m <- c(7, 14, 6, 18,
       21, 9, 57,
       14, 7, 8,
       2, 1,
       2, 9, 4, 3,
       12) # total number of genes in the family in Ag PEST genome (accessed 2025-05-02, GCA_000005575.2, AgamP4.14)

# NB: considered DEGs with logFC > |0.25|
a <- c(3, 4, 0, 8,
       8, 2, 29,
       8, 1, 1,
       0, 0,
       1, 4, 0, 0,
       0) # number of DEGs belonging to the family

a_up <- c(2, 3, 0, 3,
          4, 1, 27,
          8, 0, 0,
          0, 0,
          1, 3, 0, 0,
          0) # number of DEGs belonging to the family
a_down <- c(1, 1, 0, 5,
            4, 1, 2,
            0, 1, 1,
            0, 0,
            0, 1, 0, 0,
            0) # number of DEGs belonging to the family


# ALL GENES
# Calculate hypergeometric p-values
p_values <- phyper(a-1, m, N-m, n, lower.tail = FALSE)
# Calculate FDR (Benjamini-Hochberg)
fdr <- p.adjust(p_values, method = "fdr")
# Create a results table
results <- data.frame(
  Family = families,
  m = m,
  a = a,
  p_value = p_values,
  fdr = fdr
)
# Sort by FDR (ascending)
results <- results[order(results$fdr), ]
# Print the results
print(results)
# Optional: Save to CSV
# write.csv(results, "gene_family_enrichment_results.csv", row.names = FALSE)


# * just upregulated genes
# Parameters
N <- 13845   # Total genes in Ag PEST genome
n <- 2165    # Number of urpeg DEGs in aging analysis
# ALL GENES
# Calculate hypergeometric p-values
p_values <- phyper(a_up-1, m, N-m, n, lower.tail = FALSE)
# Calculate FDR (Benjamini-Hochberg)
fdr <- p.adjust(p_values, method = "fdr")
# Create a results table
results <- data.frame(
  Family = families,
  m = m,
  a = a_up,
  p_value = p_values,
  fdr = fdr
)
# Sort by FDR (ascending)
results_up <- results[order(results$fdr), ]
# Print the results
print(results_up)
# Optional: Save to CSV
# write.csv(results_up, "gene_family_enrichment_results_up.csv", row.names = FALSE)


# * just downregulated genes
# Parameters
N <- 13845   # Total genes in Ag PEST genome
n <- 2233    # Number of dowreg DEGs in aging analysis
# ALL GENES
# Calculate hypergeometric p-values
p_values <- phyper(a_down-1, m, N-m, n, lower.tail = FALSE)
# Calculate FDR (Benjamini-Hochberg)
fdr <- p.adjust(p_values, method = "fdr")
# Create a results table
results <- data.frame(
  Family = families,
  m = m,
  a = a_down,
  p_value = p_values,
  fdr = fdr
)
# Sort by FDR (ascending)
results_down <- results[order(results$fdr), ]
# Print the results
print(results_down)
# Optional: Save to CSV
# write.csv(results_down, "gene_family_enrichment_results_down.csv", row.names = FALSE)
