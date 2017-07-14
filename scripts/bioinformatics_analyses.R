# Eskew et al. wood frog/American bullfrog transcriptomics study

# Bioinformatics Analyses

# This script generates: Figures 6, 7, 8a, 8b, 8c, 9, 10, 11, 12

#==============================================================================


# Load packages


# Load the tidyverse
library(tidyverse)

# Data visualization packages
library(cowplot)
library(gplots)
library(VennDiagram)
library(gridExtra)

# Differential expression package
library(edgeR)

# Annotation packages
library(biomaRt)
library(GO.db)
library(GOstats)
library(GSEABase)


# Source functions
source("../R/bioinformatics_functions.R")

#==============================================================================


# Prepare a Lithobates clamitans annotation dataframe using biomaRt


# The primary goal of this section is to supplement the previously existing
# L. clamitans transcriptome annotations with GO terms to provide more 
# functional information

# Import L. clamitans annotation file supplied by Robertson and Cornman 2014, 
# which includes RefSeq, Entrez, and UniProt IDs
clamitans.data <- 
  read_csv("../data/annotation/Lithobates_clamitans_annotation.csv") %>%
  as.data.frame()

# Set up the BioMart database, using Ensembl data
ensembl.mart <- 
  useMart("ENSEMBL_MART_ENSEMBL", 
          dataset = "xtropicalis_gene_ensembl",
          host = "jul2015.archive.ensembl.org")

# Get annotation information (particularly GO terms)
# for the L. clamitans contigs, given the supplied Entrez gene IDs
anno.table <- 
  getBM(attributes = c("description", "refseq_peptide", 
                       "refseq_peptide_predicted", "entrezgene", "go_id",
                       "name_1006", "definition_1006", "go_linkage_type",
                       "namespace_1003"),
        filters = "entrezgene", 
        values = clamitans.data$'Entrez gene ID', 
        mart = ensembl.mart)

# Populate a contig_identifier column so we can identify the L. clamitans
# contigs easily. Note that some Entrez Gene IDs mapped to multiple L. 
# clamitans contigs (i.e., multiple contigs had the same Entrez ID), so in 
# these cases I'm just taking the first matching contig
anno.table$contig_identifier <- 
  sapply(1:nrow(anno.table), function(x) as.character(clamitans.data[
    match(anno.table$entrezgene[x], clamitans.data$'Entrez gene ID')[1], 1])) 

# Reorder and clean up the dataframe
anno.table$go_linkage_type <- as.factor(anno.table$go_linkage_type)
anno.table$namespace_1003 <- as.factor(anno.table$namespace_1003)
anno.table$contig_identifier <- as.factor(anno.table$contig_identifier)
anno.table <- anno.table[c(10, 1:9)]

str(anno.table)

# How many unique genes are present in the biomaRt annotation table?
# Should match the number of levels for anno.table$contig_identifier
unique(anno.table$entrezgene) %>% length()
levels(anno.table$contig_identifier) %>% length()


# Import Trinotate annotations
trinotate.table <-
  read_csv("../data/annotation/trinotate_annotation_report.csv") %>%
  as.data.frame()

#==============================================================================


# Data packaging for differential expression analyses


# Import metadata on RNA-seq samples
d.samples <- read_csv("../data/RNAseq_metadata.csv") %>%
  as.data.frame() %>%
  filter(Well != "")
d.samples$Species <- 
  as.factor(ifelse(d.samples$Species == "Wood Frog", "WoodFrog", "Bullfrog"))
d.samples$Treatment <- relevel(as.factor(d.samples$Treatment), ref = "Control")
d.samples$DayOfSacrifice <- as.factor(d.samples$DayOfSacrifice)
d.samples <- droplevels(d.samples)

# Create species-specific dataframes
d.bull <- droplevels(filter(d.samples, Species == "Bullfrog"))
d.wood <- droplevels(filter(d.samples, Species == "WoodFrog"))

# Create sample well labels
labels <- paste("Sample", d.samples$Well, sep = "")
labels.bull <- paste("Sample", d.bull$Well, sep = "")
labels.wood <- paste("Sample", d.wood$Well, sep = "")

# Create a Group factor representing unique Species, Treatment, and 
# day of tissue harvesting combinations
d.samples$Group <- 
  factor(paste(d.samples$Species, d.samples$Treatment, 
               d.samples$DayOfSacrifice, sep = "."))
d.bull$Group <- 
  factor(paste(d.bull$Treatment, d.bull$DayOfSacrifice, sep = "."))
d.wood$Group <- 
  factor(paste(d.wood$Treatment, d.wood$DayOfSacrifice, sep = "."))


# Create a vector of file names to import Salmon outputs

# Full pipeline:
# Preprocessing - Trimmomatic, ConDeTri
# Reference - Lithobates clamitans + Bd + AMPs transcripts concatenated
# Alignment - Salmon
files <- 
  sapply(d.samples$Well, 
         function(z) paste("../data/salmon/to_concatenated_transcriptome", 
                           "_w_amps/Sample", z, 
                           "_noadapt_hq15trim_clamitans_bd_amps_concatenated",
                           "_salmon_quant/quant.sf", sep = ""))
files.bull <- 
  sapply(d.bull$Well, 
         function(z) paste("../data/salmon/to_concatenated_transcriptome", 
                           "_w_amps/Sample", z, 
                           "_noadapt_hq15trim_clamitans_bd_amps_concatenated",
                           "_salmon_quant/quant.sf", sep = ""))
files.wood <- 
  sapply(d.wood$Well, 
         function(z) paste("../data/salmon/to_concatenated_transcriptome", 
                           "_w_amps/Sample", z, 
                           "_noadapt_hq15trim_clamitans_bd_amps_concatenated",
                           "_salmon_quant/quant.sf", sep = ""))


# Or create a vector of file names to import Salmon outputs following 
# expHTS preprocessing

# Full pipeline:
# Preprocessing - expHTS
# Reference - Lithobates clamitans + Bd + AMPs transcripts concatenated
# Alignment - Salmon
files2 <- 
  sapply(d.samples$Well, 
         function(z) paste("../data/salmon/to_concatenated_transcriptome", 
                           "_w_amps_expHTS/Sample", z, 
                           "_expHTS_clamitans_bd_amps_concatenated",
                           "_salmon_quant/quant.sf", sep = ""))
files2.bull <- 
  sapply(d.bull$Well, 
         function(z) paste("../data/salmon/to_concatenated_transcriptome", 
                           "_w_amps_expHTS/Sample", z, 
                           "_expHTS_clamitans_bd_amps_concatenated",
                           "_salmon_quant/quant.sf", sep = ""))
files2.wood <- 
  sapply(d.wood$Well, 
         function(z) paste("../data/salmon/to_concatenated_transcriptome", 
                           "_w_amps_expHTS/Sample", z, 
                           "_expHTS_clamitans_bd_amps_concatenated",
                           "_salmon_quant/quant.sf", sep = ""))


# Read in data to a DGEList object. Use sample names as labels and only take 
# the 1st and 5th columns in the original data files (representing gene IDs 
# and counts from eXpress)

# data <- readDGE(files, columns = c(1, 5), labels = labels)
data <- readDGE(files2, columns = c(1, 5), labels = labels)

# data.bull <- readDGE(files.bull, columns = c(1, 5), labels = labels.bull)
data.bull <- readDGE(files2.bull, columns = c(1, 5), labels = labels.bull)
                              
# data.wood <- readDGE(files.wood, columns = c(1, 5), labels = labels.wood)
data.wood <- readDGE(files2.wood, columns = c(1, 5), labels = labels.wood)


# Round the counts and view the data
data <- DGEList(counts = round(data$counts))

nrow(data)
head(data$counts)
data$samples

data.bull <- DGEList(counts = round(data.bull$counts))

nrow(data.bull)
head(data.bull$counts)
data.bull$samples

data.wood <- DGEList(counts = round(data.wood$counts))

nrow(data.wood)
head(data.wood$counts)
data.wood$samples


# Add sample info to the dataframes
data$samples <- cbind(data$samples, d.samples)
data.bull$samples <- cbind(data.bull$samples, d.bull)
data.wood$samples <- cbind(data.wood$samples, d.wood)

# Relevel treatment factors
data$samples$Treatment <- 
  relevel(as.factor(data$samples$Treatment), ref = "Control")
data.bull$samples$Treatment <- 
  relevel(as.factor(data.bull$samples$Treatment), ref = "Control")
data.wood$samples$Treatment <- 
  relevel(as.factor(data.wood$samples$Treatment), ref = "Control")

#==============================================================================


# Data pre-processing for differential expression analyses


# Add basic Lithobates clamitans annotation data to the dataframe
data$genes <- 
  data.frame(clamitans.data$Contig, clamitans.data$description,
             clamitans.data$'RefSeq ID', clamitans.data$'Entrez gene ID', 
             clamitans.data$'UniProt ID')
colnames(data$genes) <- 
  c("Contig_Name", "Description", "RefSeq_ID", "Entrez_ID", "UniProt_ID")


# Create datasets representing host- and pathogen-associated contigs
keep <- grepl("BDET", rownames(data$counts))
data.bd <- data[keep, , keep.lib.sizes=FALSE]

keep <- grepl("L", rownames(data$counts))
data <- data[keep, , keep.lib.sizes=FALSE]

nrow(data) # Should be 50,249
nrow(data.bd) # Should be 8,819


# Filter out low expression contigs, keeping any transcript that has 
# a count > 5 in seven or more samples
keep <- rowSums(data$counts > 5) >= 7
data <- data[keep, , keep.lib.sizes=FALSE]

keep <- rowSums(data.bd$counts > 5) >= 7
data.bd <- data.bd[keep, , keep.lib.sizes=FALSE]


# Report the total number of transcripts remaining in the datasets
nrow(data)
nrow(data.bd)

# Report number of transcripts mapping to L. clamitans in the host dataset
sum(grepl("Lithobates", rownames(data$counts)))

# Report number of transcripts mapping to AMPs in the host dataset
sum(grepl("Lcla", rownames(data$counts)))


# Further filter Bd dataset to include only contigs with counts of < 10
# in all control samples
keep <-
  rowSums(data.bd$counts[, data.bd$samples$Treatment == "Control"] < 10) == 30
data.bd <- data.bd[keep, , keep.lib.sizes=FALSE]

nrow(data.bd)


# Create DGEList objects and calculate normalization factors.
dge <- DGEList(counts = data)
dge <- calcNormFactors(dge)

dge.bd <- DGEList(counts = data.bd)
dge.bd <- calcNormFactors(dge.bd)

#==============================================================================


# MDS Plots


# Assign vectors for plotting
treatment.vec <- d.samples$Treatment
species.vec <- d.samples$Species
time.vec <- d.samples$DayOfSacrifice

# Make sure colors are assigned to treatment groups appropriately
levels(treatment.vec)
colors <- c("green", "blue", "orange", "red")


# Plot of all RNA-seq data with labels
plotMDS(dge, top = 100, xlim = c(-3, 3), ylim = c(-3, 3),
        cex.axis = 1.2, cex.lab = 1.4, las = 1, col = colors[treatment.vec],
        cex = 1.5, xlab = "Dimension 1", ylab = "Dimension 2",
        labels = ifelse(species.vec == "WoodFrog", "W", "B"))
legend("bottom", cex = 1.2,
       c("Control", "Carter Meadow", "Section Line",
         "PE Section Line"), 
       fill = c("green", "blue", "red", "orange"), bty = "n")


# Create figure

tiff("../figures/Figure_6.tiff", width = 1000, height = 1100, res = 96)

mat <- matrix(c(1, 1, 1, 2, 3, 4), 2, 3, byrow = TRUE)
layout(mat, heights = c(2, 1))

# Plot of all RNA-seq data with shapes

plotMDS(dge, top = 100, xlim = c(-3, 3), ylim = c(-3, 3),
        cex.axis = 1.2, cex.lab = 1.4, las = 1, col = colors[treatment.vec],
        cex = 2.3, xlab = "Dimension 1", ylab = "Dimension 2",
        pch = c(19, 17)[species.vec])
legend("bottom", cex = 2,
       c("Control", "Carter Meadow", "Section Line",
         "PE Section Line"), 
       fill = c("green", "blue", "red", "orange"), bty = "n")

# Tryptch of MDS plots

transparency <- c(1, 0, 0)[time.vec]

plotMDS(dge, top = 100, xlim = c(-3, 3), ylim = c(-3, 3),
        cex.axis = 1.2, cex.lab = 1.4, las = 1, 
        col = alpha(colors[treatment.vec], alpha = transparency),
        cex = 1.5, xlab = "Dimension 1", ylab = "Dimension 2",
        pch = c(19, 17)[species.vec],
        main = "Day 3 Only", cex.main = 1.5)

transparency <- c(0, 1, 0)[time.vec]

plotMDS(dge, top = 100, xlim = c(-3, 3), ylim = c(-3, 3),
        cex.axis = 1.2, cex.lab = 1.4, las = 1, 
        col = alpha(colors[treatment.vec], alpha = transparency),
        cex = 1.5, xlab = "Dimension 1", ylab = "Dimension 2",
        pch = c(19, 17)[species.vec],
        main = "Day 7 Only", cex.main = 1.5)

transparency <- c(0, 0, 1)[time.vec]

plotMDS(dge, top = 100, xlim = c(-3, 3), ylim = c(-3, 3),
        cex.axis = 1.2, cex.lab = 1.4, las = 1, 
        col = alpha(colors[treatment.vec], alpha = transparency),
        cex = 1.5, xlab = "Dimension 1", ylab = "Dimension 2",
        pch = c(19, 17)[species.vec],
        main = "Day 10 Only", cex.main = 1.5)

dev.off()

#==============================================================================


# Differential expression analyses in edgeR


# Define design matrices
design <- model.matrix(~0 + Group, data = data$samples)
colnames(design) <- levels(data$samples$Group)
design


# Define contrasts
contrasts <- makeContrasts(
  Bull_Con_7V3 = Bullfrog.Control.7 - Bullfrog.Control.3,
  Bull_Con_10V3 = Bullfrog.Control.10 - Bullfrog.Control.3,
  Bull_CarVCon_3 = Bullfrog.Carter.3 - Bullfrog.Control.3,
  Bull_SecVCon_3 = Bullfrog.SectionLine.3 - Bullfrog.Control.3,
  Bull_SecVCon_7 = Bullfrog.SectionLine.7 - Bullfrog.Control.7,
  Bull_SecVCon_10 = Bullfrog.SectionLine.10 - Bullfrog.Control.10,
  Wood_Con_7V3 = WoodFrog.Control.7 - WoodFrog.Control.3,
  Wood_Con_10V3 = WoodFrog.Control.10 - WoodFrog.Control.3,
  Wood_CarVCon_3 = WoodFrog.Carter.3 - WoodFrog.Control.3,
  Wood_CarVCon_7 = WoodFrog.Carter.7 - WoodFrog.Control.7,
  Wood_CarVCon_10 = WoodFrog.Carter.10 - WoodFrog.Control.10,
  Wood_SecVCon_3 = WoodFrog.SectionLine.3 - WoodFrog.Control.3,
  Wood_SecVCon_7 = WoodFrog.SectionLine.7 - WoodFrog.Control.7,
  Wood_SecVCon_10 = WoodFrog.SectionLine.10 - WoodFrog.Control.10,
  Wood_PESecVCon_3 = 
    WoodFrog.PreviouslyExposedSectionLine.3 - WoodFrog.Control.3,
  Wood_PESecVCon_7 = 
    WoodFrog.PreviouslyExposedSectionLine.7 - WoodFrog.Control.7,
  Wood_PESecVSec_3 =
    WoodFrog.PreviouslyExposedSectionLine.3 - WoodFrog.SectionLine.3,
  Wood_PESecVSec_7 =
    WoodFrog.PreviouslyExposedSectionLine.7 - WoodFrog.SectionLine.7,
  WoodVBull_CarVCon_3 = 
    (WoodFrog.Carter.3 - WoodFrog.Control.3) -
    (Bullfrog.Carter.3 - Bullfrog.Control.3),
  WoodVBull_SecVCon_3 = 
    (WoodFrog.SectionLine.3 - WoodFrog.Control.3) -
    (Bullfrog.SectionLine.3 - Bullfrog.Control.3),
  WoodVBull_SecVCon_7 = 
    (WoodFrog.SectionLine.7 - WoodFrog.Control.7) -
    (Bullfrog.SectionLine.7 - Bullfrog.Control.7),
  WoodVBull_SecVCon_10 = 
    (WoodFrog.SectionLine.10 - WoodFrog.Control.10) -
    (Bullfrog.SectionLine.10 - Bullfrog.Control.10),
  WoodVBull = 
    ((WoodFrog.Control.3 + WoodFrog.Control.7 + WoodFrog.Control.10)/3) -
    ((Bullfrog.Control.3 + Bullfrog.Control.7 + Bullfrog.Control.10)/3),
  levels = design)


# Estimate dispersions
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
plotBCV(dge)


# voom fitting
v <- voom(dge, design, plot = T)

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts = contrasts)
efit <- eBayes(vfit)
plotSA(efit)

dt <- decideTests(efit)
summary(dt)

tfit <- treat(vfit, lfc = log(2, 2))

dt2 <- decideTests(tfit)
summary(dt2)


# Compute quasi-likelihood F-tests
qlf_fit <- glmQLFit(dge, design)

qlf_Bull_Con_7V3 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Bull_Con_7V3"])
qlf_Bull_Con_10V3 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Bull_Con_10V3"])
qlf_Bull_CarVCon_3 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Bull_CarVCon_3"])
qlf_Bull_SecVCon_3 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Bull_SecVCon_3"])
qlf_Bull_SecVCon_7 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Bull_SecVCon_7"])
qlf_Bull_SecVCon_10 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Bull_SecVCon_10"])

qlf_Wood_Con_7V3 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Wood_Con_7V3"])
qlf_Wood_Con_10V3 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Wood_Con_10V3"])
qlf_Wood_CarVCon_3 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Wood_CarVCon_3"])
qlf_Wood_CarVCon_7 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Wood_CarVCon_7"])
qlf_Wood_CarVCon_10 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Wood_CarVCon_10"])
qlf_Wood_SecVCon_3 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Wood_SecVCon_3"])
qlf_Wood_SecVCon_7 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Wood_SecVCon_7"])
qlf_Wood_SecVCon_10 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Wood_SecVCon_10"])
qlf_Wood_PESecVCon_3 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Wood_PESecVCon_3"])
qlf_Wood_PESecVCon_7 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Wood_PESecVCon_7"])

qlf_Wood_PESecVSec_3 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Wood_PESecVSec_3"])
qlf_Wood_PESecVSec_7 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "Wood_PESecVSec_7"])

qlf_WoodVBull_CarVCon_3 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "WoodVBull_CarVCon_3"])
qlf_WoodVBull_SecVCon_3 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "WoodVBull_SecVCon_3"])
qlf_WoodVBull_SecVCon_7 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "WoodVBull_SecVCon_7"])
qlf_WoodVBull_SecVCon_10 <- 
  glmQLFTest(qlf_fit, contrast = contrasts[, "WoodVBull_SecVCon_10"])

summary(decideTestsDGE(qlf_Bull_Con_7V3))
summary(decideTestsDGE(qlf_Bull_Con_10V3))
summary(decideTestsDGE(qlf_Bull_CarVCon_3))
summary(decideTestsDGE(qlf_Bull_SecVCon_3))
summary(decideTestsDGE(qlf_Bull_SecVCon_7))
summary(decideTestsDGE(qlf_Bull_SecVCon_10))

summary(decideTestsDGE(qlf_Wood_Con_7V3))
summary(decideTestsDGE(qlf_Wood_Con_10V3))
summary(decideTestsDGE(qlf_Wood_CarVCon_3))
summary(decideTestsDGE(qlf_Wood_CarVCon_7))
summary(decideTestsDGE(qlf_Wood_CarVCon_10))
summary(decideTestsDGE(qlf_Wood_SecVCon_3))
summary(decideTestsDGE(qlf_Wood_SecVCon_7))
summary(decideTestsDGE(qlf_Wood_SecVCon_10))
summary(decideTestsDGE(qlf_Wood_PESecVCon_3))
summary(decideTestsDGE(qlf_Wood_PESecVCon_7))

summary(decideTestsDGE(qlf_WoodVBull_CarVCon_3))
summary(decideTestsDGE(qlf_WoodVBull_SecVCon_3))
summary(decideTestsDGE(qlf_WoodVBull_SecVCon_7))
summary(decideTestsDGE(qlf_WoodVBull_SecVCon_10))

topTags(qlf_Bull_Con_7V3)
topTags(qlf_Bull_Con_10V3)
topTags(qlf_Bull_CarVCon_3)
topTags(qlf_Bull_SecVCon_3)
topTags(qlf_Bull_SecVCon_7)
topTags(qlf_Bull_SecVCon_10)

topTags(qlf_Wood_Con_7V3)
topTags(qlf_Wood_Con_10V3)
topTags(qlf_Wood_CarVCon_3)
topTags(qlf_Wood_CarVCon_7)
topTags(qlf_Wood_CarVCon_10)
topTags(qlf_Wood_SecVCon_3)
topTags(qlf_Wood_SecVCon_7)
topTags(qlf_Wood_SecVCon_10)
topTags(qlf_Wood_PESecVCon_3)
topTags(qlf_Wood_PESecVCon_7)

topTags(qlf_WoodVBull_CarVCon_3)
topTags(qlf_WoodVBull_SecVCon_3)
topTags(qlf_WoodVBull_SecVCon_7)
topTags(qlf_WoodVBull_SecVCon_10)

#==============================================================================


# Define gene sets


# Pairwise treatment vs. control contrasts
genes_Bull_Con_7V3 <- get_de(qlf_Bull_Con_7V3)
genes_Bull_Con_10V3 <- get_de(qlf_Bull_Con_10V3)
genes_Bull_CarVCon_3 <- get_de(qlf_Bull_CarVCon_3)
genes_Bull_SecVCon_3 <- get_de(qlf_Bull_SecVCon_3)
genes_Bull_SecVCon_7 <- get_de(qlf_Bull_SecVCon_7)
genes_Bull_SecVCon_10 <- get_de(qlf_Bull_SecVCon_10)

genes_Wood_Con_7V3 <- get_de(qlf_Wood_Con_7V3)
genes_Wood_Con_10V3 <- get_de(qlf_Wood_Con_10V3)
genes_Wood_CarVCon_3 <- get_de(qlf_Wood_CarVCon_3)
genes_Wood_CarVCon_7 <- get_de(qlf_Wood_CarVCon_7)
genes_Wood_CarVCon_10 <- get_de(qlf_Wood_CarVCon_10)
genes_Wood_SecVCon_3 <- get_de(qlf_Wood_SecVCon_3)
genes_Wood_SecVCon_7 <- get_de(qlf_Wood_SecVCon_7)
genes_Wood_SecVCon_10 <- get_de(qlf_Wood_SecVCon_10)
genes_Wood_PESecVCon_3 <- get_de(qlf_Wood_PESecVCon_3)
genes_Wood_PESecVCon_7 <- get_de(qlf_Wood_PESecVCon_7)

genes_Wood_PESecVSec_3 <- get_de(qlf_Wood_PESecVSec_3)
genes_Wood_PESecVSec_7 <- get_de(qlf_Wood_PESecVSec_7)


genes_Bull_Con_7V3_up <- get_de_up(qlf_Bull_Con_7V3)
genes_Bull_Con_10V3_up <- get_de_up(qlf_Bull_Con_10V3)
genes_Bull_CarVCon_3_up <- get_de_up(qlf_Bull_CarVCon_3)
genes_Bull_SecVCon_3_up <- get_de_up(qlf_Bull_SecVCon_3)
genes_Bull_SecVCon_7_up <- get_de_up(qlf_Bull_SecVCon_7)
genes_Bull_SecVCon_10_up <- get_de_up(qlf_Bull_SecVCon_10)

genes_Wood_Con_7V3_up <- get_de_up(qlf_Wood_Con_7V3)
genes_Wood_Con_10V3_up <- get_de_up(qlf_Wood_Con_10V3)
genes_Wood_CarVCon_3_up <- get_de_up(qlf_Wood_CarVCon_3)
genes_Wood_CarVCon_7_up <- get_de_up(qlf_Wood_CarVCon_7)
genes_Wood_CarVCon_10_up <- get_de_up(qlf_Wood_CarVCon_10)
genes_Wood_SecVCon_3_up <- get_de_up(qlf_Wood_SecVCon_3)
genes_Wood_SecVCon_7_up <- get_de_up(qlf_Wood_SecVCon_7)
genes_Wood_SecVCon_10_up <- get_de_up(qlf_Wood_SecVCon_10)
genes_Wood_PESecVCon_3_up <- get_de_up(qlf_Wood_PESecVCon_3)
genes_Wood_PESecVCon_7_up <- get_de_up(qlf_Wood_PESecVCon_7)


genes_Bull_Con_7V3_down <- get_de_down(qlf_Bull_Con_7V3)
genes_Bull_Con_10V3_down <- get_de_down(qlf_Bull_Con_10V3)
genes_Bull_CarVCon_3_down <- get_de_down(qlf_Bull_CarVCon_3)
genes_Bull_SecVCon_3_down <- get_de_down(qlf_Bull_SecVCon_3)
genes_Bull_SecVCon_7_down <- get_de_down(qlf_Bull_SecVCon_7)
genes_Bull_SecVCon_10_down <- get_de_down(qlf_Bull_SecVCon_10)

genes_Wood_Con_7V3_down <- get_de_down(qlf_Wood_Con_7V3)
genes_Wood_Con_10V3_down <- get_de_down(qlf_Wood_Con_10V3)
genes_Wood_CarVCon_3_down <- get_de_down(qlf_Wood_CarVCon_3)
genes_Wood_CarVCon_7_down <- get_de_down(qlf_Wood_CarVCon_7)
genes_Wood_CarVCon_10_down <- get_de_down(qlf_Wood_CarVCon_10)
genes_Wood_SecVCon_3_down <- get_de_down(qlf_Wood_SecVCon_3)
genes_Wood_SecVCon_7_down <- get_de_down(qlf_Wood_SecVCon_7)
genes_Wood_SecVCon_10_down <- get_de_down(qlf_Wood_SecVCon_10)
genes_Wood_PESecVCon_3_down <- get_de_down(qlf_Wood_PESecVCon_3)
genes_Wood_PESecVCon_7_down <- get_de_down(qlf_Wood_PESecVCon_7)


# Explicit species-level response comparison
genes_WoodVBull_CarVCon_3 <- get_de(qlf_WoodVBull_CarVCon_3)
genes_WoodVBull_SecVCon_3 <- get_de(qlf_WoodVBull_SecVCon_3)
genes_WoodVBull_SecVCon_7 <- get_de(qlf_WoodVBull_SecVCon_7)
genes_WoodVBull_SecVCon_10 <- get_de(qlf_WoodVBull_SecVCon_10)

genes_WoodVBull_CarVCon_3_up <- get_de_up(qlf_WoodVBull_CarVCon_3)
genes_WoodVBull_SecVCon_3_up <- get_de_up(qlf_WoodVBull_SecVCon_3)
genes_WoodVBull_SecVCon_7_up <- get_de_up(qlf_WoodVBull_SecVCon_7)
genes_WoodVBull_SecVCon_10_up <- get_de_up(qlf_WoodVBull_SecVCon_10)

genes_WoodVBull_CarVCon_3_down <- get_de_down(qlf_WoodVBull_CarVCon_3)
genes_WoodVBull_SecVCon_3_down <- get_de_down(qlf_WoodVBull_SecVCon_3)
genes_WoodVBull_SecVCon_7_down <- get_de_down(qlf_WoodVBull_SecVCon_7)
genes_WoodVBull_SecVCon_10_down <- get_de_down(qlf_WoodVBull_SecVCon_10)


# Wood frog genes that were upregulated at any SL timepoint
genes_Wood_SecVCon_up <- Reduce(union, list(
  genes_Wood_SecVCon_3_up, genes_Wood_SecVCon_7_up, genes_Wood_SecVCon_10_up))

# Wood frog genes that were downregulated at any SL timepoint
genes_Wood_SecVCon_down <- Reduce(union, list(
  genes_Wood_SecVCon_3_down, genes_Wood_SecVCon_7_down, 
  genes_Wood_SecVCon_10_down))

# Wood frog genes that were upregulated at any PE SL timepoint
genes_Wood_PESecVCon_up <- Reduce(union, list(
  genes_Wood_PESecVCon_3_up, genes_Wood_PESecVCon_7_up))

# Wood frog genes that were downregulated at any PE SL timepoint
genes_Wood_PESecVCon_down <- Reduce(union, list(
  genes_Wood_PESecVCon_3_down, genes_Wood_PESecVCon_7_down))

# Wood frog genes that were upregulated at any CM timepoint
genes_Wood_CarVCon_up <- Reduce(union, list(
  genes_Wood_CarVCon_3_up, genes_Wood_CarVCon_7_up, genes_Wood_CarVCon_10_up))

# Wood frog genes that were downregulated at any CM timepoint
genes_Wood_CarVCon_down <- Reduce(union, list(
  genes_Wood_CarVCon_3_down, genes_Wood_CarVCon_7_down, 
  genes_Wood_CarVCon_10_down))


# Bullfrog genes that were upregulated at any SL timepoint
genes_Bull_SecVCon_up <- Reduce(union, list(
  genes_Bull_SecVCon_3_up, genes_Bull_SecVCon_7_up, genes_Bull_SecVCon_10_up))

# Bullfrog genes that were downregulated at any SL timepoint
genes_Bull_SecVCon_down <- Reduce(union, list(
  genes_Bull_SecVCon_3_down, genes_Bull_SecVCon_7_down, 
  genes_Bull_SecVCon_10_down))


# Wood frog genes that are upregulated in all Bd exposure treatment groups
genes_Wood_disease_up <- Reduce(intersect, list(
  genes_Wood_CarVCon_3_up, genes_Wood_CarVCon_7_up, genes_Wood_CarVCon_10_up, 
  genes_Wood_SecVCon_3_up, genes_Wood_SecVCon_7_up, genes_Wood_SecVCon_10_up,
  genes_Wood_PESecVCon_3_up, genes_Wood_PESecVCon_7_up))

# Wood frog genes that are downregulated in all Bd exposure treatment groups
genes_Wood_disease_down <- Reduce(intersect, list(
  genes_Wood_CarVCon_3_down, genes_Wood_CarVCon_7_down, 
  genes_Wood_CarVCon_10_down, 
  genes_Wood_SecVCon_3_down, genes_Wood_SecVCon_7_down, 
  genes_Wood_SecVCon_10_down,
  genes_Wood_PESecVCon_3_down, genes_Wood_PESecVCon_7_down))


# Gene sets making species comparisons


# Genes that are differentially regulated in either species at any time point
# for Carter Meadow
genes_WoodOrBull_CarVCon_DE <- Reduce(union, list(
  genes_Wood_CarVCon_up, genes_Wood_CarVCon_down,
  genes_Bull_CarVCon_3_up, genes_Bull_CarVCon_3_down))

# Genes regulated in the same direction for both species at any time point 
# for Carter Meadow
genes_WoodAndBull_CarVCon_DE <- Reduce(union, list(
  intersect(genes_Wood_CarVCon_3_up, genes_Bull_CarVCon_3_up),
  intersect(genes_Wood_CarVCon_3_down, genes_Bull_CarVCon_3_down)))

genes_WoodAndBull_CarVCon_DE_up <- Reduce(union, list(
  intersect(genes_Wood_CarVCon_3_up, genes_Bull_CarVCon_3_up)))

genes_WoodAndBull_CarVCon_DE_down <- Reduce(union, list(
  intersect(genes_Wood_CarVCon_3_down, genes_Bull_CarVCon_3_down)))

# Genes regulated differently in the two species at any time point 
# for Carter Meadow
genes_WoodVBull_CarVCon_DE <- genes_WoodVBull_CarVCon_3


# Genes that are differentially regulated in either species at any time point
# for Section Line
genes_WoodOrBull_SecVCon_DE <- Reduce(union, list(
  genes_Wood_SecVCon_up, genes_Wood_SecVCon_down,
  genes_Bull_SecVCon_up, genes_Bull_SecVCon_down))

# Genes regulated in the same direction for both species at any time point 
# for Section Line
genes_WoodAndBull_SecVCon_DE <- Reduce(union, list(
  intersect(genes_Wood_SecVCon_3_up, genes_Bull_SecVCon_3_up),
  intersect(genes_Wood_SecVCon_3_down, genes_Bull_SecVCon_3_down),
  intersect(genes_Wood_SecVCon_7_up, genes_Bull_SecVCon_7_up),
  intersect(genes_Wood_SecVCon_7_down, genes_Bull_SecVCon_7_down),
  intersect(genes_Wood_SecVCon_10_up, genes_Bull_SecVCon_10_up),
  intersect(genes_Wood_SecVCon_10_down, genes_Bull_SecVCon_10_down)))

genes_WoodAndBull_SecVCon_DE_up <- Reduce(union, list(
  intersect(genes_Wood_SecVCon_3_up, genes_Bull_SecVCon_3_up),
  intersect(genes_Wood_SecVCon_7_up, genes_Bull_SecVCon_7_up),
  intersect(genes_Wood_SecVCon_10_up, genes_Bull_SecVCon_10_up)))

genes_WoodAndBull_SecVCon_DE_down <- Reduce(union, list(
  intersect(genes_Wood_SecVCon_3_down, genes_Bull_SecVCon_3_down),
  intersect(genes_Wood_SecVCon_7_down, genes_Bull_SecVCon_7_down),
  intersect(genes_Wood_SecVCon_10_down, genes_Bull_SecVCon_10_down)))

# Genes regulated differently in the two species at any time point 
# for Section Line
genes_WoodVBull_SecVCon_DE <- Reduce(union, list(
  genes_WoodVBull_SecVCon_3, genes_WoodVBull_SecVCon_7, 
  genes_WoodVBull_SecVCon_10))

genes_WoodVBull_SecVCon_DE_up <- Reduce(union, list(
  genes_WoodVBull_SecVCon_3_up, genes_WoodVBull_SecVCon_7_up, 
  genes_WoodVBull_SecVCon_10_up))

genes_WoodVBull_SecVCon_DE_down <- Reduce(union, list(
  genes_WoodVBull_SecVCon_3_down, genes_WoodVBull_SecVCon_7_down, 
  genes_WoodVBull_SecVCon_10_down))


# Genes that are differentially regulated in wood frogs in either SL or 
# PE SL at any time point
genes_Wood_SecOrPESec_DE <- Reduce(union, list(
  genes_Wood_SecVCon_up, genes_Wood_SecVCon_down, 
  genes_Wood_PESecVCon_up, genes_Wood_PESecVCon_down))

#==============================================================================


# Differential expression analyses in edgeR for Bd contigs only


# Define design matrices
design.bd <- model.matrix(~0 + Group, data = data.bd$samples)
colnames(design.bd) <- levels(data.bd$samples$Group)
design.bd


# Define contrasts
contrasts.bd <- makeContrasts(
  Bd_Wood_SecVCar_3 = (WoodFrog.SectionLine.3) - (WoodFrog.Carter.3),
  Bd_Wood_SecVCar_7 = (WoodFrog.SectionLine.7) - (WoodFrog.Carter.7),
  Bd_Wood_SecVCar_10 = (WoodFrog.SectionLine.10) - (WoodFrog.Carter.10),
  Bd_Wood_PESecVSec_3 = 
    (WoodFrog.PreviouslyExposedSectionLine.3) - (WoodFrog.SectionLine.3),
  Bd_Wood_PESecVSec_7 = 
    (WoodFrog.PreviouslyExposedSectionLine.7) - (WoodFrog.SectionLine.7),
  Bd_Bull_SecVCar_3 = (Bullfrog.SectionLine.3) - (Bullfrog.Carter.3),
  levels = design.bd)


# Estimate dispersions
dge.bd <- estimateGLMCommonDisp(dge.bd, design.bd) 
dge.bd <- estimateGLMTrendedDisp(dge.bd, design.bd)
dge.bd <- estimateGLMTagwiseDisp(dge.bd, design.bd)
plotBCV(dge.bd)


# Compute quasi-likelihood F-tests
qlf_fit <- glmQLFit(dge.bd, design.bd)

qlf_Bd_Wood_SecVCar_3 <-
  glmQLFTest(qlf_fit, contrast = contrasts.bd[, "Bd_Wood_SecVCar_3"])
qlf_Bd_Wood_SecVCar_7 <-
  glmQLFTest(qlf_fit, contrast = contrasts.bd[, "Bd_Wood_SecVCar_7"])
qlf_Bd_Wood_SecVCar_10 <-
  glmQLFTest(qlf_fit, contrast = contrasts.bd[, "Bd_Wood_SecVCar_10"])
qlf_Bd_Wood_PESecVSec_3 <-
  glmQLFTest(qlf_fit, contrast = contrasts.bd[, "Bd_Wood_PESecVSec_3"])
qlf_Bd_Wood_PESecVSec_7 <-
  glmQLFTest(qlf_fit, contrast = contrasts.bd[, "Bd_Wood_PESecVSec_7"])
qlf_Bd_Bull_SecVCar_3 <-
  glmQLFTest(qlf_fit, contrast = contrasts.bd[, "Bd_Bull_SecVCar_3"])

summary(decideTestsDGE(qlf_Bd_Wood_SecVCar_3))
summary(decideTestsDGE(qlf_Bd_Wood_SecVCar_7))
summary(decideTestsDGE(qlf_Bd_Wood_SecVCar_10))
summary(decideTestsDGE(qlf_Bd_Wood_PESecVSec_3))
summary(decideTestsDGE(qlf_Bd_Wood_PESecVSec_7))
summary(decideTestsDGE(qlf_Bd_Bull_SecVCar_3))

topTags(qlf_Bd_Wood_SecVCar_3)
topTags(qlf_Bd_Wood_SecVCar_7)
topTags(qlf_Bd_Wood_SecVCar_10)
topTags(qlf_Bd_Wood_PESecVSec_3)
topTags(qlf_Bd_Wood_PESecVSec_7)
topTags(qlf_Bd_Bull_SecVCar_3)


# Define gene sets
genes_Bd_Wood_SecVCar_3_up <- 
  rownames(data.bd$counts)[decideTestsDGE(qlf_Bd_Wood_SecVCar_3) == 1]
genes_Bd_Wood_SecVCar_3_down <- 
  rownames(data.bd$counts)[decideTestsDGE(qlf_Bd_Wood_SecVCar_3) == -1]
genes_Bd_Wood_SecVCar_7_up <- 
  rownames(data.bd$counts)[decideTestsDGE(qlf_Bd_Wood_SecVCar_7) == 1]
genes_Bd_Wood_SecVCar_7_down <- 
  rownames(data.bd$counts)[decideTestsDGE(qlf_Bd_Wood_SecVCar_7) == -1]
genes_Bd_Wood_SecVCar_10_up <- 
  rownames(data.bd$counts)[decideTestsDGE(qlf_Bd_Wood_SecVCar_10) == 1]
genes_Bd_Wood_SecVCar_10_down <- 
  rownames(data.bd$counts)[decideTestsDGE(qlf_Bd_Wood_SecVCar_10) == -1]
genes_Bd_Wood_PESecVSec_3_up <- 
  rownames(data.bd$counts)[decideTestsDGE(qlf_Bd_Wood_PESecVSec_3) == 1]
genes_Bd_Wood_PESecVSec_3_down <- 
  rownames(data.bd$counts)[decideTestsDGE(qlf_Bd_Wood_PESecVSec_3) == -1]
genes_Bd_Wood_PESecVSec_7_up <- 
  rownames(data.bd$counts)[decideTestsDGE(qlf_Bd_Wood_PESecVSec_7) == 1]
genes_Bd_Wood_PESecVSec_7_down <- 
  rownames(data.bd$counts)[decideTestsDGE(qlf_Bd_Wood_PESecVSec_7) == -1]

#==============================================================================


# Other design matrix structures


# Full interaction formula (with intercept)

# Create DGEList objects and calculate normalization factors
# dge2 <- DGEList(counts = data)
# dge2 <- calcNormFactors(dge2)


# Define design matrix
# design2 <- model.matrix(~Species * Treatment * DayOfSacrifice, 
#                        data = data$samples)
# dim(design2)
# apply(design2, 2, sum)

# Remove empty columns and those that are not estimable
# Column 9: "SpeciesWoodFrog:TreatmentPreviouslyExposedSectionLine" 
# unestimable because it is synonymous with 
# "TreatmentPreviouslyExposedSectionLine"
# Column 17: Empty
# Column 19: "SpeciesWoodFrog:TreatmentCarter:DayOfSacrifice7" unestimable 
# because it is synonymous with "TreatmentCarter:DayOfSacrifice7"
# Column 20: 
# "SpeciesWoodFrog:TreatmentPreviouslyExposedSectionLine:DayOfSacrifice7"
# unestimable because it is synonymous with 
# "TreatmentPreviouslyExposedSectionLine:DayOfSacrifice10"
# Column 22: "SpeciesWoodFrog:TreatmentCarter:DayOfSacrifice10" unestimable
# because it is synonymous with "TreatmentCarter:DayOfSacrifice10"
# Column 23: Empty
# design2 <- design2[ ,-c(9, 17, 19, 20, 22, 23)]
# dim(design2)
# apply(design2, 2, sum)


# Estimate dispersions
# dge2 <- estimateGLMCommonDisp(dge2, design2)
# dge2 <- estimateGLMTrendedDisp(dge2, design2)
# dge2 <- estimateGLMTagwiseDisp(dge2, design2)
# plotBCV(dge2)


# Fit model
# qlf_fit2 <- glmQLFit(dge2, design2)


# Comparing differential expression calls between design matrix formulations

# summary(decideTestsDGE(qlf_Bull_Con_7V3))
# summary(decideTestsDGE(glmQLFTest(qlf_fit2, 
#                               contrast = 
#                                 c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0))))
# 
# summary(decideTestsDGE(qlf_Bull_Con_10V3))
# summary(decideTestsDGE(glmQLFTest(qlf_fit2,
#                               contrast = 
#                                 c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0))))

# Genes that respond to Carter at day 3 (i.e., in bullfrogs)
# summary(decideTestsDGE(glmQLFTest(qlf_fit2,
#                                   contrast = 
#                                     c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))))
# Equivalent to:
# summary(decideTestsDGE(qlf_Bull_CarVCon_3))

# Genes that respond to Section Line at day 3 (i.e., in bullfrogs)
# summary(decideTestsDGE(glmQLFTest(qlf_fit2,
#                                   contrast = 
#                                     c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0))))
# Equivalent to:
# summary(decideTestsDGE(qlf_Bull_SecVCon_3))

# Genes that respond to Section Line at day 7 (i.e., in bullfrogs)
# summary(decideTestsDGE(glmQLFTest(qlf_fit2,
#                                   contrast = 
#                                     c(0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0))))
# Equivalent to:
# summary(decideTestsDGE(qlf_Bull_SecVCon_7))

# Genes that respond to Section Line at day 10 (i.e., in bullfrogs)
# summary(decideTestsDGE(glmQLFTest(qlf_fit2,
#                                   contrast = 
#                                     c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,0,0))))
# Equivalent to:
# summary(decideTestsDGE(qlf_Bull_SecVCon_10))


# Genes that respond to Carter only in wood frogs at day 3
# summary(decideTestsDGE(glmQLFTest(qlf_fit2,
#                                   contrast = 
#                                     c(0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0))))
# Equivalent to:
# summary(decideTestsDGE(qlf_Wood_CarVCon_3))

# Genes that respond to Section Line only in wood frogs at day 3
# summary(decideTestsDGE(glmQLFTest(qlf_fit2,
#                                   contrast = 
#                                     c(0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0))))
# Equivalent to:
# summary(decideTestsDGE(qlf_Wood_SecVCon_3))

# Genes that respond to Section Line only in wood frogs at day 7
# summary(decideTestsDGE(glmQLFTest(qlf_fit2,
#                                   contrast = 
#                                     c(0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,1,0))))
# Equivalent to:
# summary(decideTestsDGE(qlf_Wood_SecVCon_7))

# Genes that respond DIFFERENTLY to Section Line in wood frogs and
# bullfrogs at day 3
# summary(decideTestsDGE(glmQLFTest(qlf_fit2,
#                                   contrast = 
#                                     c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0))))
# Equivalent to:
# summary(decideTestsDGE(qlf_WoodVBull_SecVCon_3))

#==============================================================================


# GOstats annotation analysis


# Compile necessary GO data for custom GOstats analysis
# Need GO terms, evidence codes, and gene identifiers
goframeData <- 
  data.frame(anno.table$go_id, anno.table$go_linkage_type, 
             anno.table$contig_identifier)
colnames(goframeData) <- c("GO_ID", "Evidence_Code", "Contig")
nrow(goframeData)


# Or get GO data from Trinotate annotation
go.table <- 
  read_tsv("../data/annotation/go_annotations.txt") %>%
  as.data.frame()
colnames(go.table) <- c("Contig", "GO_ID")

go.table <- go.table %>% 
  mutate(V3 = strsplit(GO_ID, ",")) %>% 
  unnest(V3) %>% 
  dplyr::select(Contig, V3) %>%
  dplyr::rename(GO_ID = V3)

goframeData <- 
  data.frame(go.table$GO_ID, rep("IEA", nrow(go.table)), go.table$Contig)
colnames(goframeData) <- c("GO_ID", "Evidence_Code", "Contig")
nrow(goframeData)


# Subset the dataframe so that it doesn't have GO terms with missing
# evidence codes
goframeData <- goframeData[goframeData$Evidence_Code != "", ]
nrow(goframeData)


# Construct GO dataframes
goFrame <- GOFrame(goframeData, organism = "Xenopus tropicalis")
goAllFrame <- GOAllFrame(goFrame)


# Convert to geneSetCollection object
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())


# Set up parameter object
universe <- rownames(data$counts)


# Perform GO testing

# Wood frog contigs
params_Wood_SecVCon_3_up <- get.go.params(genes_Wood_SecVCon_3_up)
enriched_Wood_SecVCon_3_up <- hyperGTest(params_Wood_SecVCon_3_up)

params_Wood_SecVCon_7_up <- get.go.params(genes_Wood_SecVCon_7_up)
enriched_Wood_SecVCon_7_up <- hyperGTest(params_Wood_SecVCon_7_up)

params_Wood_SecVCon_10_up <- get.go.params(genes_Wood_SecVCon_10_up)
enriched_Wood_SecVCon_10_up <- hyperGTest(params_Wood_SecVCon_10_up)

params_Wood_SecVCon_up <- get.go.params(genes_Wood_SecVCon_up) 
enriched_Wood_SecVCon_up <- hyperGTest(params_Wood_SecVCon_up)

params_Wood_SecVCon_3_down <- get.go.params(genes_Wood_SecVCon_3_down)
enriched_Wood_SecVCon_3_down <- hyperGTest(params_Wood_SecVCon_3_down)

params_Wood_SecVCon_7_down <- get.go.params(genes_Wood_SecVCon_7_down)
enriched_Wood_SecVCon_7_down <- hyperGTest(params_Wood_SecVCon_7_down)

params_Wood_SecVCon_10_down <- get.go.params(genes_Wood_SecVCon_10_down)
enriched_Wood_SecVCon_10_down <- hyperGTest(params_Wood_SecVCon_10_down)

params_Wood_SecVCon_down <- get.go.params(genes_Wood_SecVCon_down)
enriched_Wood_SecVCon_down <- hyperGTest(params_Wood_SecVCon_down)

enriched_Wood_SecVCon_3_up
enriched_Wood_SecVCon_7_up
enriched_Wood_SecVCon_10_up
enriched_Wood_SecVCon_up

enriched_Wood_SecVCon_3_down
enriched_Wood_SecVCon_7_down
enriched_Wood_SecVCon_10_down
enriched_Wood_SecVCon_down

summary(enriched_Wood_SecVCon_3_up)$Term
summary(enriched_Wood_SecVCon_7_up)$Term
summary(enriched_Wood_SecVCon_10_up)$Term
summary(enriched_Wood_SecVCon_up)$Term

summary(enriched_Wood_SecVCon_3_down)$Term
summary(enriched_Wood_SecVCon_7_down)$Term
summary(enriched_Wood_SecVCon_10_down)$Term
summary(enriched_Wood_SecVCon_down)$Term


params_Wood_PESecVCon_3_up <- get.go.params(genes_Wood_PESecVCon_3_up)
enriched_Wood_PESecVCon_3_up <- hyperGTest(params_Wood_PESecVCon_3_up)

params_Wood_PESecVCon_7_up <- get.go.params(genes_Wood_PESecVCon_7_up)
enriched_Wood_PESecVCon_7_up <- hyperGTest(params_Wood_PESecVCon_7_up)

params_Wood_PESecVCon_up <- get.go.params(genes_Wood_PESecVCon_up)
enriched_Wood_PESecVCon_up <- hyperGTest(params_Wood_PESecVCon_up)

params_Wood_PESecVCon_3_down <- get.go.params(genes_Wood_PESecVCon_3_down)
enriched_Wood_PESecVCon_3_down <- hyperGTest(params_Wood_PESecVCon_3_down)

params_Wood_PESecVCon_7_down <- get.go.params(genes_Wood_PESecVCon_7_down)
enriched_Wood_PESecVCon_7_down <- hyperGTest(params_Wood_PESecVCon_7_down)

params_Wood_PESecVCon_down <- get.go.params(genes_Wood_PESecVCon_down)
enriched_Wood_PESecVCon_down <- hyperGTest(params_Wood_PESecVCon_down)

enriched_Wood_PESecVCon_3_up
enriched_Wood_PESecVCon_7_up
enriched_Wood_PESecVCon_up

enriched_Wood_PESecVCon_3_down
enriched_Wood_PESecVCon_7_down
enriched_Wood_PESecVCon_down

summary(enriched_Wood_PESecVCon_3_up)$Term
summary(enriched_Wood_PESecVCon_7_up)$Term
summary(enriched_Wood_PESecVCon_up)$Term

summary(enriched_Wood_PESecVCon_3_down)$Term
summary(enriched_Wood_PESecVCon_7_down)$Term
summary(enriched_Wood_PESecVCon_down)$Term


params_Wood_CarVCon_3_up <- get.go.params(genes_Wood_CarVCon_3_up)
enriched_Wood_CarVCon_3_up <- hyperGTest(params_Wood_CarVCon_3_up)

params_Wood_CarVCon_7_up <- get.go.params(genes_Wood_CarVCon_7_up)
enriched_Wood_CarVCon_7_up <- hyperGTest(params_Wood_CarVCon_7_up)

params_Wood_CarVCon_10_up <- get.go.params(genes_Wood_CarVCon_10_up)
enriched_Wood_CarVCon_10_up <- hyperGTest(params_Wood_CarVCon_10_up)

params_Wood_CarVCon_up <- get.go.params(genes_Wood_CarVCon_up)
enriched_Wood_CarVCon_up <- hyperGTest(params_Wood_CarVCon_up)

params_Wood_CarVCon_3_down <- get.go.params(genes_Wood_CarVCon_3_down)
enriched_Wood_CarVCon_3_down <- hyperGTest(params_Wood_CarVCon_3_down)

params_Wood_CarVCon_7_down <- get.go.params(genes_Wood_CarVCon_7_down)
enriched_Wood_CarVCon_7_down <- hyperGTest(params_Wood_CarVCon_7_down)

params_Wood_CarVCon_10_down <- get.go.params(genes_Wood_CarVCon_10_down)
enriched_Wood_CarVCon_10_down <- hyperGTest(params_Wood_CarVCon_10_down)

params_Wood_CarVCon_down <- get.go.params(genes_Wood_CarVCon_down)
enriched_Wood_CarVCon_down <- hyperGTest(params_Wood_CarVCon_down)

enriched_Wood_CarVCon_3_up
enriched_Wood_CarVCon_7_up
enriched_Wood_CarVCon_10_up
enriched_Wood_CarVCon_up

enriched_Wood_CarVCon_3_down
enriched_Wood_CarVCon_7_down
enriched_Wood_CarVCon_10_down
enriched_Wood_CarVCon_down

summary(enriched_Wood_CarVCon_3_up)$Term
summary(enriched_Wood_CarVCon_7_up)$Term
summary(enriched_Wood_CarVCon_10_up)$Term
summary(enriched_Wood_CarVCon_up)$Term

summary(enriched_Wood_CarVCon_3_down)$Term
summary(enriched_Wood_CarVCon_7_down)$Term
summary(enriched_Wood_CarVCon_10_down)$Term
summary(enriched_Wood_CarVCon_down)$Term


# American bullfrog contigs
params_Bull_SecVCon_3_up <- get.go.params(genes_Bull_SecVCon_3_up)
enriched_Bull_SecVCon_3_up <- hyperGTest(params_Bull_SecVCon_3_up)

params_Bull_SecVCon_7_up <- get.go.params(genes_Bull_SecVCon_7_up)
enriched_Bull_SecVCon_7_up <- hyperGTest(params_Bull_SecVCon_7_up)

params_Bull_SecVCon_10_up <- get.go.params(genes_Bull_SecVCon_10_up)
enriched_Bull_SecVCon_10_up <- hyperGTest(params_Bull_SecVCon_10_up)

params_Bull_SecVCon_up <- get.go.params(genes_Bull_SecVCon_up)
enriched_Bull_SecVCon_up <- hyperGTest(params_Bull_SecVCon_up)

params_Bull_SecVCon_3_down <- get.go.params(genes_Bull_SecVCon_3_down)
enriched_Bull_SecVCon_3_down <- hyperGTest(params_Bull_SecVCon_3_down)

params_Bull_SecVCon_7_down <- get.go.params(genes_Bull_SecVCon_7_down)
enriched_Bull_SecVCon_7_down <- hyperGTest(params_Bull_SecVCon_7_down)

params_Bull_SecVCon_10_down <- get.go.params(genes_Bull_SecVCon_10_down)
enriched_Bull_SecVCon_10_down <- hyperGTest(params_Bull_SecVCon_10_down)

params_Bull_SecVCon_down <- get.go.params(genes_Bull_SecVCon_down)
enriched_Bull_SecVCon_down <- hyperGTest(params_Bull_SecVCon_down)

enriched_Bull_SecVCon_3_up
enriched_Bull_SecVCon_7_up
enriched_Bull_SecVCon_10_up
enriched_Bull_SecVCon_up

enriched_Bull_SecVCon_3_down
enriched_Bull_SecVCon_7_down
enriched_Bull_SecVCon_10_down
enriched_Bull_SecVCon_down

summary(enriched_Bull_SecVCon_3_up)$Term
summary(enriched_Bull_SecVCon_7_up)$Term
summary(enriched_Bull_SecVCon_10_up)$Term
summary(enriched_Bull_SecVCon_up)$Term

summary(enriched_Bull_SecVCon_3_down)$Term
summary(enriched_Bull_SecVCon_7_down)$Term
summary(enriched_Bull_SecVCon_10_down)$Term
summary(enriched_Bull_SecVCon_down)$Term


params_Bull_CarVCon_3_up <- get.go.params(genes_Bull_CarVCon_3_up)
enriched_Bull_CarVCon_3_up <- hyperGTest(params_Bull_CarVCon_3_up)

params_Bull_CarVCon_3_down <- get.go.params(genes_Bull_CarVCon_3_down)
enriched_Bull_CarVCon_3_down <- hyperGTest(params_Bull_CarVCon_3_down)

enriched_Bull_CarVCon_3_up
enriched_Bull_CarVCon_3_down

summary(enriched_Bull_CarVCon_3_up)$Term
summary(enriched_Bull_CarVCon_3_down)$Term


# Species shared responses to Bd
params_WoodAndBull_SecVCon_DE_up <- 
  get.go.params(genes_WoodAndBull_SecVCon_DE_up)
enriched_WoodAndBull_SecVCon_DE_up <- 
  hyperGTest(params_WoodAndBull_SecVCon_DE_up)

params_WoodAndBull_SecVCon_DE_down <- 
  get.go.params(genes_WoodAndBull_SecVCon_DE_down)
enriched_WoodAndBull_SecVCon_DE_down <- 
  hyperGTest(params_WoodAndBull_SecVCon_DE_down)

enriched_WoodAndBull_SecVCon_DE_up
enriched_WoodAndBull_SecVCon_DE_down

summary(enriched_WoodAndBull_SecVCon_DE_up)$Term
summary(enriched_WoodAndBull_SecVCon_DE_down)$Term


# Species DIVERGENT responses to Bd
params_WoodVBull_SecVCon_DE_up <- 
  get.go.params(genes_WoodVBull_SecVCon_DE_up)
enriched_WoodVBull_SecVCon_DE_up <- 
  hyperGTest(params_WoodVBull_SecVCon_DE_up)

params_WoodVBull_SecVCon_DE_down <- 
  get.go.params(genes_WoodVBull_SecVCon_DE_down)
enriched_WoodVBull_SecVCon_DE_down <- 
  hyperGTest(params_WoodVBull_SecVCon_DE_down)

enriched_WoodVBull_SecVCon_DE_up
enriched_WoodVBull_SecVCon_DE_down

summary(enriched_WoodVBull_SecVCon_DE_up)$Term
summary(enriched_WoodVBull_SecVCon_DE_down)$Term


# Contigs broadly induced during disease
params_Wood_disease_up <- get.go.params(genes_Wood_disease_up)
enriched_Wood_disease_up <- hyperGTest(params_Wood_disease_up)

params_Wood_disease_down <- get.go.params(genes_Wood_disease_down)
enriched_Wood_disease_down <- hyperGTest(params_Wood_disease_down)

enriched_Wood_disease_up
enriched_Wood_disease_down

summary(enriched_Wood_disease_up)$Term
summary(enriched_Wood_disease_down)$Term


# Info used in tables for manuscript
enriched_Wood_CarVCon_up
length(genes_Wood_CarVCon_up)
summary(enriched_Wood_CarVCon_up)$Term[1:10]

enriched_Wood_CarVCon_down
length(genes_Wood_CarVCon_down)
summary(enriched_Wood_CarVCon_down)$Term[1:10]

enriched_Wood_SecVCon_up
length(genes_Wood_SecVCon_up)
summary(enriched_Wood_SecVCon_up)$Term[1:10]

enriched_Wood_SecVCon_down
length(genes_Wood_SecVCon_down)
summary(enriched_Wood_SecVCon_down)$Term[1:10]

enriched_Wood_PESecVCon_up
length(genes_Wood_PESecVCon_up)
summary(enriched_Wood_PESecVCon_up)$Term[1:10]

enriched_Wood_PESecVCon_down
length(genes_Wood_PESecVCon_down)
summary(enriched_Wood_PESecVCon_down)$Term[1:10]

enriched_Bull_CarVCon_3_up
length(genes_Bull_CarVCon_3_up)
summary(enriched_Bull_CarVCon_3_up)$Term[1:10]

enriched_Bull_CarVCon_3_down
length(genes_Bull_CarVCon_3_down)
summary(enriched_Bull_CarVCon_3_down)$Term[1:10]

enriched_Bull_SecVCon_up
length(genes_Bull_SecVCon_up)
summary(enriched_Bull_SecVCon_up)$Term[1:10]

enriched_Bull_SecVCon_down
length(genes_Bull_SecVCon_down)
summary(enriched_Bull_SecVCon_down)$Term[1:10]

# Contrasts comparing species responses to Bd
enriched_WoodAndBull_SecVCon_DE_up
length(genes_WoodAndBull_SecVCon_DE_up)
summary(enriched_WoodAndBull_SecVCon_DE_up)$Term[1:10]

enriched_WoodAndBull_SecVCon_DE_down
length(genes_WoodAndBull_SecVCon_DE_down)
summary(enriched_WoodAndBull_SecVCon_DE_down)$Term[1:10]

enriched_WoodVBull_SecVCon_DE_up
length(genes_WoodVBull_SecVCon_DE_up)
summary(enriched_WoodVBull_SecVCon_DE_up)$Term[1:10]

enriched_WoodVBull_SecVCon_DE_down
length(genes_WoodVBull_SecVCon_DE_down)
summary(enriched_WoodVBull_SecVCon_DE_down)$Term[1:10]

#==============================================================================


# GOstats annotation analysis

# Bd contigs only


universe <- rownames(data.bd$counts)


params_Bd_Wood_SecVCar_3_up <- get.go.params(genes_Bd_Wood_SecVCar_3_up)
enriched_Bd_Wood_SecVCar_3_up <- hyperGTest(params_Bd_Wood_SecVCar_3_up)

params_Bd_Wood_SecVCar_3_down <- get.go.params(genes_Bd_Wood_SecVCar_3_down)
enriched_Bd_Wood_SecVCar_3_down <- hyperGTest(params_Bd_Wood_SecVCar_3_down)

params_Bd_Wood_PESecVSec_3_up <- get.go.params(genes_Bd_Wood_PESecVSec_3_up)
enriched_Bd_Wood_PESecVSec_3_up <- hyperGTest(params_Bd_Wood_PESecVSec_3_up)

params_Bd_Wood_PESecVSec_3_down <- 
  get.go.params(genes_Bd_Wood_PESecVSec_3_down)
enriched_Bd_Wood_PESecVSec_3_down <- 
  hyperGTest(params_Bd_Wood_PESecVSec_3_down)

params_Bd_Wood_PESecVSec_7_up <- get.go.params(genes_Bd_Wood_PESecVSec_7_up)
enriched_Bd_Wood_PESecVSec_7_up <- hyperGTest(params_Bd_Wood_PESecVSec_7_up)

params_Bd_Wood_PESecVSec_7_down <- 
  get.go.params(genes_Bd_Wood_PESecVSec_7_down)
enriched_Bd_Wood_PESecVSec_7_down <- 
  hyperGTest(params_Bd_Wood_PESecVSec_7_down)

enriched_Bd_Wood_SecVCar_3_up
enriched_Bd_Wood_SecVCar_3_down
enriched_Bd_Wood_PESecVSec_3_up
enriched_Bd_Wood_PESecVSec_3_down
enriched_Bd_Wood_PESecVSec_7_up
enriched_Bd_Wood_PESecVSec_7_down

summary(enriched_Bd_Wood_SecVCar_3_up)$Term
summary(enriched_Bd_Wood_SecVCar_3_down)$Term
summary(enriched_Bd_Wood_PESecVSec_3_up)$Term
summary(enriched_Bd_Wood_PESecVSec_3_down)$Term
summary(enriched_Bd_Wood_PESecVSec_7_up)$Term
summary(enriched_Bd_Wood_PESecVSec_7_down)$Term


# Info used in tables for manuscript
enriched_Bd_Wood_SecVCar_3_up
length(genes_Bd_Wood_SecVCar_3_up)
summary(enriched_Bd_Wood_SecVCar_3_up)$Term[1:10]

enriched_Bd_Wood_SecVCar_3_down
length(genes_Bd_Wood_SecVCar_3_down)
summary(enriched_Bd_Wood_SecVCar_3_down)$Term[1:10]

enriched_Bd_Wood_PESecVSec_3_up
length(genes_Bd_Wood_PESecVSec_3_up)
summary(enriched_Bd_Wood_PESecVSec_3_up)$Term[1:10]

enriched_Bd_Wood_PESecVSec_3_down
length(genes_Bd_Wood_PESecVSec_3_down)
summary(enriched_Bd_Wood_PESecVSec_3_down)$Term[1:10]

enriched_Bd_Wood_PESecVSec_7_up
length(genes_Bd_Wood_PESecVSec_7_up)
summary(enriched_Bd_Wood_PESecVSec_7_up)$Term[1:10]

enriched_Bd_Wood_PESecVSec_7_down
length(genes_Bd_Wood_PESecVSec_7_down)
summary(enriched_Bd_Wood_PESecVSec_7_down)$Term[1:10]

#==============================================================================


# Differential expression plots


dayspost <- c(3, 7, 10)

# Wood frog differential expression plot
w.car.de <- 
  c(get_num_de(qlf_Wood_CarVCon_3), get_num_de(qlf_Wood_CarVCon_7),
    get_num_de(qlf_Wood_CarVCon_10))
w.sec.de <- 
  c(get_num_de(qlf_Wood_SecVCon_3), get_num_de(qlf_Wood_SecVCon_7),
    get_num_de(qlf_Wood_SecVCon_10))
w.pe.sec.de <-
  c(get_num_de(qlf_Wood_PESecVCon_3), get_num_de(qlf_Wood_PESecVCon_7), NA)

plot(w.car.de ~ dayspost, xlim = c(0, 12), ylim = c(0, 4000), las = 1,
     bty = "l", cex.axis = 1.2, cex.lab = 1.4, cex.main = 1.5, type = "n",
     xlab = "Days Post-Exposure", 
     ylab= "No. Differentially Expressed Genes")
legend(x = 5, y = 4000, c("Carter Meadow", "Section Line", "PE Section Line"), 
       fill = c("blue", "red", "orange"), bty = "n", cex = 1.4)

points(w.car.de ~ dayspost, col = "blue", pch = 16, cex = 2)
lines(w.car.de ~ dayspost, col = "blue", lty = 1)

points(w.sec.de ~ dayspost, col= "red", pch = 16, cex = 2)
lines(w.sec.de ~ dayspost, col= "red", lty = 1)

points(w.pe.sec.de ~ dayspost, col= "orange", pch = 16, cex = 2)
lines(w.pe.sec.de ~ dayspost, col= "orange", lty = 1)


# Bullfrog differential expression plot
b.car.de <- c(get_num_de(qlf_Bull_CarVCon_3), NA, NA)
b.sec.de <- 
  c(get_num_de(qlf_Bull_SecVCon_3), get_num_de(qlf_Bull_SecVCon_7),
    get_num_de(qlf_Bull_SecVCon_10))

plot(b.car.de ~ dayspost, xlim = c(0, 12), ylim = c(0,200), las = 1,
     bty = "l", cex.axis = 1.2, cex.lab = 1.4, cex.main = 1.5, type = "n",
     xlab = "Days Post-Exposure", 
     ylab= "No. Differentially Expressed Genes")
legend(x = 0, y = 200, c("Carter Meadow", "Section Line"), 
       fill = c("blue", "red"), bty = "n", cex = 1.4)

points(b.car.de ~ dayspost, col = "blue", pch = 16, cex = 2)
lines(b.car.de ~ dayspost, col = "blue", lty = 1)

points(b.sec.de ~ dayspost, col= "red", pch = 16, cex = 2)
lines(b.sec.de ~ dayspost, col= "red", lty = 1)

#==============================================================================


# Plot comparing species responses to Bd exposure


tiff("../figures/Figure_12.tiff", width = 1000, height = 700, res = 96)

par(mfrow = c(2, 3))

col.vec <- ifelse(decideTestsDGE(qlf_Bull_CarVCon_3) != 0 |
                    decideTestsDGE(qlf_Wood_CarVCon_3) != 0,
                  "blue", "darkgrey")

plot(qlf_Bull_CarVCon_3$table$logFC ~ qlf_Wood_CarVCon_3$table$logFC,
     xlim = c(-8, 8), ylim = c(-8, 8),
     main = "Day 3", cex.main = 2,
     xlab = "Wood Frog Response", ylab = "American Bullfrog Response",
     col = col.vec, pch = 19, cex.lab = 1.3)
points(qlf_Bull_CarVCon_3$table$logFC[col.vec == "blue"] ~
         qlf_Wood_CarVCon_3$table$logFC[col.vec == "blue"],
       col = "blue", pch = 19)
abline(lm(qlf_Wood_CarVCon_3$table$logFC ~ qlf_Bull_CarVCon_3$table$logFC),
       lty = 2)
text(0, -7, 
     paste("Correlation = ", 
           round(cor(qlf_Wood_CarVCon_3$table$logFC, 
                     qlf_Bull_CarVCon_3$table$logFC), digits = 3)), cex = 1.3)

plot(0, pch = '',
     xlim = c(-8, 8), ylim = c(-8, 8),
     main = "Day 7", cex.main = 2, cex.lab = 1.3,
     xlab = "Wood Frog Response", ylab = "American Bullfrog Response")
text(0, 0, "No Relevant Data", cex = 1.5)

plot(0, pch = '',
     xlim = c(-8, 8), ylim = c(-8, 8),
     main = "Day 10", cex.main = 2, cex.lab = 1.3,
     xlab = "Wood Frog Response", ylab = "American Bullfrog Response")
text(0, 0, "No Relevant Data", cex = 1.5)


col.vec <- ifelse(decideTestsDGE(qlf_Bull_SecVCon_3) != 0 |
         decideTestsDGE(qlf_Wood_SecVCon_3) != 0,
         "red", "darkgrey")

plot(qlf_Bull_SecVCon_3$table$logFC ~ qlf_Wood_SecVCon_3$table$logFC,
     xlim = c(-8, 8), ylim = c(-8, 8),
     #main = "Day 3", cex.main = 2,
     xlab = "Wood Frog Response", ylab = "American Bullfrog Response",
     col = col.vec, pch = 19, cex.lab = 1.3)
points(qlf_Bull_SecVCon_3$table$logFC[col.vec == "red"] ~
         qlf_Wood_SecVCon_3$table$logFC[col.vec == "red"],
       col = "red", pch = 19)
abline(lm(qlf_Wood_SecVCon_3$table$logFC ~ qlf_Bull_SecVCon_3$table$logFC),
       lty = 2)
text(0, -7, 
     paste("Correlation = ", 
           round(cor(qlf_Wood_SecVCon_3$table$logFC, 
                     qlf_Bull_SecVCon_3$table$logFC), digits = 3)), cex = 1.3)

col.vec <- ifelse(decideTestsDGE(qlf_Bull_SecVCon_7) != 0 |
                    decideTestsDGE(qlf_Wood_SecVCon_7) != 0,
                  "red", "darkgrey")

plot(qlf_Bull_SecVCon_7$table$logFC ~ qlf_Wood_SecVCon_7$table$logFC,
     xlim = c(-8, 8), ylim = c(-8, 8),
     #main = "Day 7", cex.main = 2,
     xlab = "Wood Frog Response", ylab = "American Bullfrog Response",
     col = col.vec, pch = 19, cex.lab = 1.3)
points(qlf_Bull_SecVCon_7$table$logFC[col.vec == "red"] ~
         qlf_Wood_SecVCon_7$table$logFC[col.vec == "red"],
       col = "red", pch = 19)
abline(lm(qlf_Wood_SecVCon_7$table$logFC ~ qlf_Bull_SecVCon_7$table$logFC),
       lty = 2)
text(0, -7, 
     paste("Correlation = ", 
           round(cor(qlf_Wood_SecVCon_7$table$logFC, 
                     qlf_Bull_SecVCon_7$table$logFC), digits = 3)), cex = 1.3)

col.vec <- ifelse(decideTestsDGE(qlf_Bull_SecVCon_10) != 0 |
                    decideTestsDGE(qlf_Wood_SecVCon_10) != 0,
                  "red", "darkgrey")

plot(qlf_Bull_SecVCon_10$table$logFC ~ qlf_Wood_SecVCon_10$table$logFC,
     xlim = c(-8, 8), ylim = c(-8, 8),
     #main = "Day 10", cex.main = 2,
     xlab = "Wood Frog Response", ylab = "American Bullfrog Response",
     col = col.vec, pch = 19, cex.lab = 1.3)
points(qlf_Bull_SecVCon_10$table$logFC[col.vec == "red"] ~
         qlf_Wood_SecVCon_10$table$logFC[col.vec == "red"],
       col = "red", pch = 19)
abline(lm(qlf_Wood_SecVCon_10$table$logFC ~ qlf_Bull_SecVCon_10$table$logFC),
       lty = 2)
text(0, -7, 
     paste("Correlation = ", 
           round(cor(qlf_Wood_SecVCon_10$table$logFC, 
                     qlf_Bull_SecVCon_10$table$logFC), digits = 3)), cex = 1.3)

dev.off()

#==============================================================================


# Venn diagrams

venn(list("Car 3" = genes_Wood_CarVCon_3, "Car 7" = genes_Wood_CarVCon_7,
          "Car 10" = genes_Wood_CarVCon_10))

venn(list("SL 3" = genes_Wood_SecVCon_3, "SL 7" = genes_Wood_SecVCon_7,
          "SL 10" = genes_Wood_SecVCon_10))

venn1 <- venn(list("CM" = genes_Wood_CarVCon_3,
                   "SL" = genes_Wood_SecVCon_3,
                   "PE SL" = genes_Wood_PESecVCon_3))

venn2 <- venn(list("CM" = genes_Wood_CarVCon_7,
                   "SL" = genes_Wood_SecVCon_7,
                   "PE SL" = genes_Wood_PESecVCon_7))

venn3 <- venn(list("CM" = genes_Wood_CarVCon_10,
                   "SL" = genes_Wood_SecVCon_10))

venn4 <- venn(list("CM" = genes_Bull_CarVCon_3,
                   "SL" = genes_Bull_SecVCon_3))

venn5 <- venn(list("SL" = genes_Bull_SecVCon_7))

venn6 <- venn(list("SL" = genes_Bull_SecVCon_10))


# Create Venn diagrams for every species by treatment group combination

cat.cex <- 0
cex <- 1.3
main.cex <- 1.5
alpha <- 0.7
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

venn1 <- venn.diagram(list("CM" = genes_Wood_CarVCon_3,
                           "SL" = genes_Wood_SecVCon_3,
                           "PE SL" = genes_Wood_PESecVCon_3), NULL,
                      fill = c("blue", "red", "orange"),
                      fontfamily = "Arial",
                      # main = "Wood Frog, Day 3", main.cex = 1.3,
                      cat.cex = cat.cex, cex = cex, alpha = alpha)

venn2 <- venn.diagram(list("CM" = genes_Wood_CarVCon_7,
                           "SL" = genes_Wood_SecVCon_7,
                           "PE SL" = genes_Wood_PESecVCon_7), NULL,
                      fill = c("blue", "red", "orange"),
                      fontfamily = "Arial",
                      # main = "Wood Frog, Day 7", main.cex = 1.3,
                      cat.cex = cat.cex, cex = cex, alpha = alpha)

venn3 <- venn.diagram(list("CM" = genes_Wood_CarVCon_10,
                           "SL" = genes_Wood_SecVCon_10), NULL,
                      fill = c("blue", "red"),
                      fontfamily = "Arial",
                      # main = "Wood Frog, Day 10", main.cex = 1.3,
                      cat.cex = cat.cex, cex = cex, alpha = alpha,
                      rotation.degree = 180, scaled = F)

venn4 <- venn.diagram(list("CM" = genes_Bull_CarVCon_3,
                           "SL" = genes_Bull_SecVCon_3), NULL,
                      fill = c("blue", "red"),
                      fontfamily = "Arial",
                      # main = "Day 3",
                      cat.cex = cat.cex, cex = cex, alpha = alpha,
                      rotation.degree = 180, scaled = F)

venn5 <- venn.diagram(list("SL" = genes_Bull_SecVCon_7), NULL,
                      fill = c("red"),
                      fontfamily = "Arial", margin = 0.105,
                      # main = "Day 7",
                      cat.cex = cat.cex, cex = cex, alpha = alpha)

venn6 <- venn.diagram(list("SL" = genes_Bull_SecVCon_10), NULL,
                      fill = c("red"),
                      fontfamily = "Arial",
                      # main = "Day 10",
                      cat.cex = cat.cex, cex = cex, alpha = alpha)


# Put Venn diagrams together on one plot

tiff("../figures/Figure_7.tiff", width = 800, height = 700, res = 96)

grid.arrange(textGrob("Wood Frog\nDay 3", gp = gpar(cex = cex)),
             textGrob("Wood Frog\nDay 7", gp = gpar(cex = cex)),
             textGrob("Wood Frog\nDay 10", gp = gpar(cex = cex)),
             gTree(children = venn1), 
             gTree(children = venn2),
             gTree(children = venn3),
             textGrob("American Bullfrog\nDay 3", gp = gpar(cex = cex)),
             textGrob("American Bullfrog\nDay 7", gp = gpar(cex = cex)),
             textGrob("American Bullfrog\nDay 10", gp = gpar(cex = cex)),
             gTree(children = venn4),
             gTree(children = venn5),
             textGrob("No\nDifferential Expression\nDetected"),
             ncol = 3, nrow = 4, heights = c(1, 3, 1, 3))

dev.off()

#==============================================================================


# Heatmaps, part 1


# Colors could indicate Treatment
colors <- c("blue", "green", "orange", "red")
color.vec <- colors[treatment.vec]

# Or colors could indicate Species
colors <- c("forestgreen", "brown")
color.vec <- colors[species.vec]


# Select genes to plot and put these gene rows in a matrix
genes.to.plot <- rownames(topTags(qlf_Wood_SecVCon_7, 500))
h <- as.matrix(data$counts[rownames(data$counts) %in% genes.to.plot, ])

var.genes <- 
  head(order(genefilter::rowVars(data$counts), decreasing = T), 500)
h <- data$counts[var.genes, ]


# Create a vector for labelling columns
ylabels <- treatment.vec
ylabels <- 
  ifelse(ylabels == "Carter", "Carter Meadow", as.character(ylabels))
ylabels <- 
  ifelse(ylabels == "SectionLine", "Section Line", as.character(ylabels))
ylabels <- 
  ifelse(ylabels == "PreviouslyExposedSectionLine", 
         "PE Section Line", as.character(ylabels))


# Plot regular heatmap
heatmap(h, labRow = NA, labCol = ylabels, Rowv = NA)


# Plot heatmap.2
heatmap.2(h, col = redgreen(100), labRow = NA, dendrogram = c("column"),
          labCol = ylabels, ColSideColors = color.vec,
          density.info = "none", trace = "none", scale = "row",
          margins = c(8, 4), lwid = c(1, 5))

# heatmap.2, only wood frogs
heatmap.2(h[ , 1:52], 
          col = redgreen(100), labRow = NA, dendrogram = c("column"),
          labCol = ylabels[1:52], ColSideColors = color.vec[1:52],
          density.info = "none", trace = "none", scale = "row",
          margins = c(8, 4), lwid = c(1, 5))

# heatmap.2, only bullfrogs
heatmap.2(h[ , 53:87], 
          col = redgreen(100), labRow = NA, dendrogram = c("column"),
          labCol = ylabels[53:87], ColSideColors = color.vec[53:87],
          density.info = "none", trace = "none", scale = "row",
          margins = c(8, 4), lwid = c(1, 5))

#==============================================================================


# Heatmaps, part 2


my.palette <- colorRampPalette(c("steelblue2", "black", "yellow2"))
width <- 900
height <- 580


# Heatmap of genes differentially expressed in Carter Meadow in either species
# at any time point

h <- as.data.frame(cbind(qlf_Wood_CarVCon_3$table$logFC,
                         qlf_Wood_CarVCon_7$table$logFC,
                         qlf_Wood_CarVCon_10$table$logFC,
                         qlf_Bull_CarVCon_3$table$logFC),
                   row.names = rownames(data$counts))
genes_to_plot <- genes_WoodOrBull_CarVCon_DE
h <- as.matrix(h[rownames(h) %in% genes_to_plot, ])

ylabels <- c("WF, Day 3", "WF, Day 7", "WF, Day 10", "AB, Day 3")

length(genes_to_plot)

# Plot heatmap.2

tiff("../figures/Figure_9.tiff", width = width, height = height)

main <- 
  paste("\n\n\nContigs differentially expressed at any time point (n = ",
        length(genes_to_plot), ")", sep = "")

heatmap.2(h, col = my.palette(200), breaks = seq(-3, 3, length.out = 201),
          labRow = NA, labCol = ylabels, Colv = NA, dendrogram = c("none"),
          srtCol = 0, adjCol = c(0.5, -36), cexCol = 1.5,
          density.info = "none", trace = "none",
          colsep = c(3), sepcolor = "white",
          lwid = c(1, 5), main = main) 

dev.off()


# Heatmap of genes differentially expressed in Section Line in either species
# at any time point

h <- as.data.frame(cbind(qlf_Wood_SecVCon_3$table$logFC,
                         qlf_Wood_SecVCon_7$table$logFC,
                         qlf_Wood_SecVCon_10$table$logFC,
                         qlf_Bull_SecVCon_3$table$logFC,
                         qlf_Bull_SecVCon_7$table$logFC,
                         qlf_Bull_SecVCon_10$table$logFC),
                   row.names = rownames(data$counts))
genes_to_plot <- genes_WoodOrBull_SecVCon_DE
h <- as.matrix(h[rownames(h) %in% genes_to_plot, ])

ylabels <- c("WF, Day 3", "WF, Day 7", "WF, Day 10", 
             "AB, Day 3", "AB, Day 7", "AB, Day 10")

length(genes_to_plot)

# Plot heatmap.2

tiff("../figures/Figure_8a.tiff", width = width, height = height)

main <- 
  paste("\nA\n\nContigs differentially expressed at any time point (n = ",
        length(genes_to_plot), ")", sep = "")

# Have to hard code number of genes
main <- expression(bold(atop(
    (italic(a)),
    "Contigs differentially expressed at any time point (n = 5767)"
)))

heatmap.2(h, col = my.palette(200), breaks = seq(-3, 3, length.out = 201),
          labRow = NA, labCol = ylabels, Colv = NA, dendrogram = c("none"),
          srtCol = 0, adjCol = c(0.5, -36), cexCol = 1.5,
          density.info = "none", trace = "none",
          colsep = c(3), sepcolor = "white",
          lwid = c(1, 5),
          main = main)

dev.off()


# Heatmap of genes COMMONLY expressed in wood frogs relative to 
# American bullfrogs following Section Line exposure

h <- as.data.frame(cbind(qlf_Wood_SecVCon_3$table$logFC,
                         qlf_Wood_SecVCon_7$table$logFC,
                         qlf_Wood_SecVCon_10$table$logFC,
                         qlf_Bull_SecVCon_3$table$logFC,
                         qlf_Bull_SecVCon_7$table$logFC,
                         qlf_Bull_SecVCon_10$table$logFC),
                   row.names = rownames(data$counts))
genes_to_plot <- genes_WoodAndBull_SecVCon_DE
h <- as.matrix(h[rownames(h) %in% genes_to_plot, ])

ylabels <- c("WF, Day 3", "WF, Day 7", "WF, Day 10", 
             "AB, Day 3", "AB, Day 7", "AB, Day 10")

length(genes_to_plot)

# Plot heatmap.2

tiff("../figures/Figure_8b.tiff", width = width, height = height)

main <- 
  paste("\nB\n\nContigs with a common response between species (n = ",
        length(genes_to_plot), ")", sep = "")

# Have to hard code number of genes
main <- expression(bold(atop(
  (italic(b)),
  "Contigs with a common response between species (n = 35)"
)))

heatmap.2(h, col = my.palette(200), breaks = seq(-3, 3, length.out = 201),
          labRow = NA, labCol = ylabels, Colv = NA, dendrogram = c("none"),
          srtCol = 0, adjCol = c(0.5, -36), cexCol = 1.5,
          density.info = "none", trace = "none",
          colsep = c(3), sepcolor = "white",
          lwid = c(1, 5),
          main = main)

dev.off()


# Heatmap of genes DIFFERENTLY expressed in wood frogs relative to 
# American bullfrogs following Section Line exposure

h <- as.data.frame(cbind(qlf_Wood_SecVCon_3$table$logFC,
                         qlf_Wood_SecVCon_7$table$logFC,
                         qlf_Wood_SecVCon_10$table$logFC,
                         qlf_Bull_SecVCon_3$table$logFC,
                         qlf_Bull_SecVCon_7$table$logFC,
                         qlf_Bull_SecVCon_10$table$logFC),
                   row.names = rownames(data$counts))
genes_to_plot <- genes_WoodVBull_SecVCon_DE
h <- as.matrix(h[rownames(h) %in% genes_to_plot, ])

ylabels <- c("WF, Day 3", "WF, Day 7", "WF, Day 10", 
             "AB, Day 3", "AB, Day 7", "AB, Day 10")

length(genes_to_plot)

# Plot heatmap.2

tiff("../figures/Figure_8c.tiff", width = width, height = height)

main <- 
  paste("\nC\n\nContigs with a different response between species (n = ",
        length(genes_to_plot), ")", sep = "")

# Have to hard code number of genes
main <- expression(bold(atop(
  (italic(c)),
  "Contigs with a different response between species (n = 813)"
)))

heatmap.2(h, col = my.palette(200), breaks = seq(-3, 3, length.out = 201),
          labRow = NA, labCol = ylabels, Colv = NA, dendrogram = c("none"),
          srtCol = 0, adjCol = c(0.5, -36), cexCol = 1.5,
          density.info = "none", trace = "none",
          colsep = c(3), sepcolor = "white",
          lwid = c(1, 5),
          main = main)

dev.off()


# Heatmap of genes differentially expressed in Section Line or PE Section
# Line in wood frogs at any time point

h <- as.data.frame(cbind(qlf_Wood_SecVCon_3$table$logFC,
                         qlf_Wood_SecVCon_7$table$logFC,
                         qlf_Wood_SecVCon_10$table$logFC,
                         qlf_Wood_PESecVCon_3$table$logFC,
                         qlf_Wood_PESecVCon_7$table$logFC),
                   row.names = rownames(data$counts))
genes_to_plot <- genes_Wood_SecOrPESec_DE
h <- as.matrix(h[rownames(h) %in% genes_to_plot, ])

ylabels <- c("SL, Day 3", "SL, Day 7", "SL, Day 10", 
             "PE SL, Day 3", "PE SL, Day 7")

length(genes_to_plot)

# Plot heatmap.2

heatmap.2(h, col = my.palette(200), breaks = seq(-3, 3, length.out = 201),
          labRow = NA, labCol = ylabels, Colv = NA, dendrogram = c("none"),
          srtCol = 0, adjCol = c(0.5, -36), cexCol = 1.5,
          density.info = "none", trace = "none",
          colsep = c(3), sepcolor = "white",
          lwid = c(1, 5),
          main = paste("\nSection Line or PE Section Line Bd-exposed Wood Frogs
          \nContigs differentially expressed at any time point (n = ",
          length(genes_to_plot), ")", sep = ""))

#==============================================================================


# Bd among sample comparison


# Re-import raw data, so that we're working with the full dataset
data <- readDGE(files2, columns = c(1, 5), labels = labels)
data <- DGEList(counts = round(data$counts))
nrow(data) 
# Should be 59,068 (representing all reference transcriptome contigs)

# Calculate a Bd CPM metric
bd_cpm_matrix <- cpm(data$counts)[grepl("BDET", rownames(data$counts)), ]
d.samples$bd_cpm <- apply(bd_cpm_matrix, 2, mean)

# Code to help make the plots
d.samples$Treatment <- relevel(d.samples$Treatment, ref = "Control")
d.samples$Species <- relevel(d.samples$Species, ref = "WoodFrog")
label_list <- c("Bullfrog" = "American Bullfrog", "WoodFrog" = "Wood Frog")

ggplot(aes(x = Treatment, y = bd_cpm, col = Treatment), 
       data = d.samples) +
  geom_jitter(size = 3, width = 0.3) +
  facet_wrap(~Species, labeller = as_labeller(label_list), scales = "free_x") +
  scale_x_discrete(labels = 
                     c("Control" = "Control", "Carter" = "CM", 
                       "SectionLine" = "SL", 
                       "PreviouslyExposedSectionLine" = "PE SL")) +
  ylab("Mean Counts per Million Mapping to Bd") +
  scale_color_manual(values = c("green", "blue", "orange", "red")) +
  guides(color = F) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_blank(),
        panel.border = element_rect(size = 3),
        text = element_text(size = 24))

#==============================================================================


# Gene-level among sample expression plots


# Re-import raw data, so that we're working with the full dataset
data <- readDGE(files2, columns = c(1, 5), labels = labels)
data <- DGEList(counts = round(data$counts))
nrow(data) 
# Should be 59,068 (representing all reference transcriptome contigs)

# Set up dataframe
gene.based.df <- 
  data.frame(
    d.samples$Species,
    d.samples$ID,
    d.samples$Treatment,
    d.samples$DayOfSacrifice,
    cpm(data$counts)[which(rownames(data$counts) == 
                             "Lcla-E"), ])
colnames(gene.based.df) <- c("Species", "ID", "Treatment", "Day", "Count")
gene.based.df$Treatment <- relevel(gene.based.df$Treatment, ref = "Control")
gene.based.df$Species <- relevel(gene.based.df$Species, ref = "WoodFrog")
label_list <- c("Bullfrog" = "American Bullfrog", "WoodFrog" = "Wood Frog")

ggplot(aes(x = Treatment, y = Count, color = Treatment, shape = Day), 
       data = gene.based.df) + #ggtitle("Contig 6394 - hemostasis") +
  geom_jitter(size = 3, width = 0.3) +
  ylab("Counts per Million") +
  scale_color_manual(values = c("green", "blue", "orange", "red")) +
  facet_wrap(~Species, labeller = as_labeller(label_list), scales = "free_x") +
  scale_x_discrete(labels =
    c("Control" = "Control", "Carter" = "CM", 
      "SectionLine" = "SL", "PreviouslyExposedSectionLine" = "PE SL")) +
  guides(color = F) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_blank(),
        panel.border = element_rect(size = 3),
        text = element_text(size = 24, family = "Helvetica"))

#==============================================================================


# Irg1 Figure


# Re-import raw data, so that we're working with the full dataset
data <- readDGE(files2, columns = c(1, 5), labels = labels)
data <- DGEList(counts = round(data$counts))
nrow(data) 
# Should be 59,068 (representing all reference transcriptome contigs)


gene.based.df <- 
  data.frame(
    d.samples$Species, 
    d.samples$Treatment,
    d.samples$DayOfSacrifice,
    cpm(data$counts)[which(rownames(data$counts) == 
                             "Lithobates.clamitans_contig_252"), ])
colnames(gene.based.df) <- c("Species", "Treatment", "Day", "Count")
gene.based.df$Treatment <- relevel(gene.based.df$Treatment, ref = "Control")
gene.based.df$Species <- relevel(gene.based.df$Species, ref = "WoodFrog")
label_list <- c("Bullfrog" = "American Bullfrog", "WoodFrog" = "Wood Frog")


tiff("../figures/Figure_11.tiff", width = 800, height = 600, res = 96)

set.seed(8)

ggplot(aes(x = Treatment, y = Count, color = Treatment, shape = Day), 
         data = gene.based.df) + ggtitle("Immune-responsive Gene 1") +
  geom_jitter(size = 3, width = 0.2) +
  ylab("Counts per Million") +
  scale_color_manual(values = c("green", "blue", "orange", "red")) +
  facet_wrap(~Species, labeller = as_labeller(label_list), scales = "free_x") +
  scale_x_discrete(labels =
                     c("Control" = "Control", "Carter" = "CM", 
                       "SectionLine" = "SL", 
                       "PreviouslyExposedSectionLine" = "PE SL")) +
  guides(color = F, shape = F) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_blank(),
        panel.border = element_rect(size = 3),
        text = element_text(size = 24, family = "Helvetica"),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))

dev.off()

#==============================================================================


# AMP figure

# Re-import raw data, so that we're working with the full dataset
data <- readDGE(files3, columns = c(1, 5), labels = labels)
data <- DGEList(counts = round(data$counts))
nrow(data) 
# Should be 59,068 (representing all reference transcriptome contigs)


gene.based.df <- 
  data.frame(
    d.samples$Species, 
    d.samples$Treatment,
    d.samples$DayOfSacrifice,
    cpm(data$counts)[which(rownames(data$counts) == 
                             "Lcla-B"), ])
colnames(gene.based.df) <- c("Species", "Treatment", "Day", "Count")
gene.based.df$Treatment <- relevel(gene.based.df$Treatment, ref = "Control")
gene.based.df$Species <- relevel(gene.based.df$Species, ref = "WoodFrog")
label_list <- c("Bullfrog" = "American Bullfrog", "WoodFrog" = "Wood Frog")

set.seed(8)

plot1 <-
  ggplot(aes(x = Treatment, y = Count, color = Treatment, shape = Day), 
       data = gene.based.df) + ggtitle("Temporin-1Cb-like Peptide") +
  geom_jitter(size = 3, width = 0.2) +
  ylab("Counts per Million") +
  scale_color_manual(values = c("green", "blue", "orange", "red")) +
  facet_wrap(~Species, labeller = as_labeller(label_list), scales = "free_x") +
  scale_x_discrete(labels =
                     c("Control" = "Control", "Carter" = "CM", 
                       "SectionLine" = "SL", 
                       "PreviouslyExposedSectionLine" = "PE SL")) +
  guides(color = F, shape = F) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_blank(),
        panel.border = element_rect(size = 3),
        text = element_text(size = 24, family = "Helvetica"),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))

gene.based.df <- 
  data.frame(
    d.samples$Species, 
    d.samples$Treatment,
    d.samples$DayOfSacrifice,
    cpm(data$counts)[which(rownames(data$counts) == 
                             "Lcla-H"), ])
colnames(gene.based.df) <- c("Species", "Treatment", "Day", "Count")
gene.based.df$Treatment <- relevel(gene.based.df$Treatment, ref = "Control")
gene.based.df$Species <- relevel(gene.based.df$Species, ref = "WoodFrog")
label_list <- c("Bullfrog" = "American Bullfrog", "WoodFrog" = "Wood Frog")

set.seed(8)

plot2 <-
ggplot(aes(x = Treatment, y = Count, color = Treatment, shape = Day), 
       data = gene.based.df) + ggtitle("Palustrin-like Peptide") +
  geom_jitter(size = 3, width = 0.2) +
  ylab("Counts per Million") +
  scale_color_manual(values = c("green", "blue", "orange", "red")) +
  facet_wrap(~Species, labeller = as_labeller(label_list), scales = "free_x") +
  scale_x_discrete(labels =
                     c("Control" = "Control", "Carter" = "CM", 
                       "SectionLine" = "SL", 
                       "PreviouslyExposedSectionLine" = "PE SL")) +
  guides(color = F, shape = F) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_blank(),
        panel.border = element_rect(size = 3),
        text = element_text(size = 24, family = "Helvetica"),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = 0.5))


tiff("../figures/Figure_10.tiff", width = 800, height = 1200, res = 96)

plot_grid(plot1, plot2, nrow = 2, scale = 0.95)

dev.off()
