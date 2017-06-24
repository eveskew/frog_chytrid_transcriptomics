# Eskew et al. wood frog/American bullfrog transcriptomics study

# expHTS Preprocessing Metadata Analyses

# This script uses metadata produced by expHTS software to check for any
# biases in read processing that may have been induced due to treatment or
# host species. In other words, you don't want samples to be systematically
# processed differently based on some experimental factor.

# "expHTS_preprocess_summary.log" is the original file output by expHTS (not
# used in this script).
# "expHTS_preprocess_summary_modified.log" is a manually modified version of
# that file that attempts to provide a more comprehensive set of data on the
# read processing pipeline.
# "expHTS_preprocess_summary_cleaned.log" is a pared down version of the
# "expHTS_preprocess_summary_modified.log" file that is used for the MDS
# plotting here. In contrast to "expHTS_preprocess_summary_modified.log", it
# does not include absolute data (such as absolute numbers of reads) because
# these metrics do vary considerably across samples and therefore might
# falsely distinguish samples. Rather, it focuses only on proportional data 
# related to the various expHTS filtering steps (i.e., percentage of reads
# filtered by Sickle).

#==============================================================================


# Load packages

library(dplyr)

#==============================================================================


input.file <- 
  "../data/expHTS_preprocess/expHTS_preprocess_summary_cleaned.log"

samples <- 
  read.table("../data/expHTS_preprocess/samples.txt", sep = "\t", header = T)

rna.seq.info <- 
  read.csv("../data/RNAseq_metadata.csv") %>%
  filter(Well != "") %>%
  arrange(Well)

samples <- 
  cbind(samples, rna.seq.info$Species, rna.seq.info$Treatment,
        paste(rna.seq.info$Species, rna.seq.info$Treatment),
        rna.seq.info$DayOfSacrifice)

colnames(samples) <- c("SEQUENCE_ID", "SAMPLE_ID", "Species", "Treatment",
                       "S_by_T", "DayOfSacrifice")

samples <- droplevels(samples)

#==============================================================================


# Plotting MDS of expHTS metadata


# Color by treatment to look for any treatment effects

col <- as.numeric(as.factor(samples$Treatment)) + 1
names(col) <- samples$SAMPLE_ID

# Read in the expHTS summary table

tb <- read.table(input.file, sep = "\t", header = T, row.names = 1, 
                 as.is = T, quote = "", comment.char = "", fill = T)

# Scale the data

scaled.tb <- apply(tb, 2, function(x) scale(as.numeric(sub("%", "",  x))))

# Perform multi-dimentional scaling

mds <- cmdscale(dist(scaled.tb))

# Plot the result

plot(mds[, 1], mds[, 2], type = "n", xlab = "MDS1", ylab = "MDS2", 
     asp = 1, axes = TRUE, main = "cmdscale (scaled results)",
     xlim = c(min(mds[, 1] - 5), max(mds[, 1] + 5)))
text(mds[, 1], mds[, 2], rownames(tb), cex = 0.6, 
     col = col[match(rownames(tb), names(col))])


# Repeat the same analysis, looking for species effects

col <- as.numeric(as.factor(samples$Species)) + 1
names(col) <- samples$SAMPLE_ID

plot(mds[, 1], mds[, 2], type = "n", xlab = "MDS1", ylab = "MDS2", 
     asp = 1, axes = TRUE, main = "cmdscale (scaled results)",
     xlim = c(min(mds[, 1] - 5), max(mds[, 1] + 5)))
text(mds[, 1], mds[, 2], rownames(tb), cex = 0.6, 
     col = col[match(rownames(tb), names(col))])

#==============================================================================


# Calculations of post-processing read counts


input.file2 <- 
  "../data/expHTS_preprocess/expHTS_preprocess_summary_modified.log"

tb2 <- read.table(input.file2, sep = "\t", header = T, row.names = 1, 
                  as.is = T, quote = "", comment.char = "", fill = T)

# How many post-processing reads total across samples?

sum(tb2$Short_PE_Kept)*2 + sum(tb2$Short_SE_Kept)

# What is the range of post-processing reads across samples?

range((tb2$Short_PE_Kept*2) + tb2$Short_SE_Kept)
