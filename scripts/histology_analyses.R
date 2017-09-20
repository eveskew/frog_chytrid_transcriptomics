# Eskew et al. wood frog/American bullfrog transcriptomics study

# Histology Analyses

# This script generates: Figure 5

#==============================================================================


# Load packages

library(ggplot2)
library(dplyr)

#==============================================================================


# Import histology dataframe

d <- read.csv("../data/histology_data.csv")


# Filter to only include records with histology scores

d <- filter(d, !is.na(HistoScore))
nrow(d)

# Round histology scores to reflect the precision of original field-by-field
# quantifications

d$HistoScore <- round(d$HistoScore)

# How many samples with histology come from sacrificed animals?

d %>%
  filter(!is.na(HistoScore) & !is.na(DayOfSacrifice)) %>%
  nrow()

# How many samples with histology come from monitored animals
# (i.e., experiment-induced euthanasia/mortality)?

d %>%
  filter(!is.na(HistoScore) & is.na(DayOfSacrifice)) %>%
  nrow()

# Summarize histology scores by Species and Treatment

group_by(d, Species, Treatment) %>%
  summarize(SampleSize = n(),
            Avg = mean(HistoScore, na.rm = T),
            Min = min(HistoScore, na.rm = T),
            Max = max(HistoScore, na.rm = T))

# What is the correlation coefficient between qPCR loads and histology scores?

cor(d$AdjustedBdLoad, d$HistoScore)

# Look only at records that had a histology score of 0
# How many qPCR records also report 0 vs. > 0?

filter(d, HistoScore == 0) %>%
  summarize(NumTotal = n(), NumMatching = sum(AdjustedBdLoad == 0),
            NumNotMatching = sum(AdjustedBdLoad > 0))

# Look only at records that had a qPCR load of 0
# How many Histology records also report 0 vs. > 0?

filter(d, AdjustedBdLoad == 0) %>%
  summarize(NumTotal = n(), NumMatching = sum(HistoScore == 0),
            NumNotMatching = sum(HistoScore > 0))

#==============================================================================


# Plotting histology vs. qPCR comparison

labels <- c("Bullfrog" = "American Bullfrog", "Wood Frog" = "Wood Frog")

d$Species <- factor(d$Species, levels = c("Wood Frog", "Bullfrog"))
d$Treatment <- factor(d$Treatment,
                      levels = c("Control", "Carter", "SectionLine", 
                                 "PreviouslyExposedSectionLine",
                                 "Untreated"))


tiff("../figures/Figure_5.tiff", width = 800, height = 700, res = 96)

filter(d) %>%
  ggplot(aes(x = AdjustedBdLoad, y = HistoScore)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.justification = c(0, 1),
        legend.position = c(0.01, 0.99),
        legend.key = element_rect(color = "grey"),
        text = element_text(size = 16, face = "bold"),
        axis.text = element_text(color = "black"),
        strip.background = element_blank()) +
  geom_point(aes(size = 1.5, color = Treatment, shape = Species)) +
  geom_smooth(method = lm, se = F, color = "black", fullrange = T) + 
  xlim(0, 15000) + ylim(0, 300) +
  xlab("Bd Load (Zoospore Equivalents)") + ylab("Histology Infection Score") +
  scale_shape_manual(values = c(4, 1)) +
  scale_color_manual(labels = c("Control", "Carter Meadow", 
                                "Section Line", "PE Section Line",
                                 "Untreated"),
                     values = c("green", "blue", "red", "orange", "purple")) +
  guides(size = F, shape = F)
  #facet_wrap(~ Species, labeller = labeller(Species = labels))

dev.off()
