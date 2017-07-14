# Eskew et al. wood frog/American bullfrog transcriptomics study

# Body Measurement Analyses

# This script generates: Figures 2a, 2b

#==============================================================================


# Load packages

library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(lme4)
library(rethinking)

#==============================================================================


# Import the raw data

d <- read.csv("../data/body_measurement_data.csv")

# Add a body condition column

d$BodyCondition <- (d$Mass/(d$SUL^3))*1000000

# Fill in a "Status" column using the infection dataset
# Status = 1; the animal was euthanized during the experiment
# Status = 2; the animal survived until the end of the experiment
# Status = 3: the animal was sacrificed during the experiment for tissues

d$Status <- rep(NA, nrow(d))
d.infec <- read.csv("../data/infection_data.csv")

for (i in 1:nrow(d)) {
  species_query <- as.character(d[i, ]$Species)
  ID_query <- d[i, ]$ID
  d[i, ]$Status <- 
    filter(d.infec, Species == species_query & ID == ID_query)[1, "Status"]
}

#==============================================================================


# Subset dataset to get rid of NAs in DaysPost as these represent data taken
# multiple weeks prior to experimental treatments

d <- droplevels(d[!is.na(d$DaysPost), ])

# Reorder the DaysPost and Treatment factors

d$DaysPost <- factor(d$DaysPost, 
                     levels = c("Pre", 3, 4, 7, 10, 11, 18, 25, 32, 39, 46))

d$Treatment <- factor(d$Treatment, 
                      levels = c("Control", "Carter", "SectionLine", 
                                 "PreviouslyExposedSectionLine"))

# Get rid of harvesting days to make the data strictly weekly

d <- droplevels(filter(d, DaysPost != 3 & DaysPost != 7 & DaysPost != 10))


# Subset data by species

d.wood <- filter(d, Species == "Wood Frog") %>% droplevels()
d.bull <- filter(d, Species == "Bullfrog") %>% droplevels()


# Subset data by treatments

w.control <- filter(d.wood, Treatment == "Control") %>% droplevels()
w.carter <- filter(d.wood, Treatment == "Carter") %>% droplevels()
w.section <- filter(d.wood, Treatment == "SectionLine") %>% droplevels()
w.pesection <- 
  filter(d.wood, Treatment == "PreviouslyExposedSectionLine") %>%
  droplevels()

b.control <- filter(d.bull, Treatment == "Control") %>% droplevels()
b.carter <- filter(d.bull, Treatment == "Carter") %>% droplevels()
b.section <- filter(d.bull, Treatment == "SectionLine") %>% droplevels()

#==============================================================================


# Were there any differences in frog body mass prior to treatments?

anova(lm(Mass ~ Treatment, data = filter(d.wood, DaysPost == "Pre")))
# Answer: No
plot(Mass ~ Treatment, data = filter(d.wood, DaysPost == "Pre"))

anova(lm(Mass ~ Treatment, data = filter(d.bull, DaysPost == "Pre")))
# Answer: No
plot(Mass ~ Treatment, data = filter(d.bull, DaysPost == "Pre"))

#==============================================================================


# Tests and plots on Mass using all the data for each species

anova(lm(Mass ~ DaysPost * Treatment, data = d.wood))
plot(Mass ~ DaysPost * Treatment, data = d.wood) 

anova(lm(Mass ~ DaysPost * Treatment, data = d.bull))
plot(Mass ~ DaysPost * Treatment, data = d.bull)


# Tests and plots on Body Condition using all the data for each species

anova(lm(BodyCondition ~ DaysPost * Treatment, data = d.wood))
plot(BodyCondition ~ DaysPost * Treatment, data = d.wood)

anova(lm(BodyCondition ~ DaysPost * Treatment, data = d.bull))
plot(BodyCondition ~ DaysPost * Treatment, data = d.bull)

#==============================================================================


# Tests and plots on data subsetted by species and treatment

anova(lm(Mass ~ DaysPost, data = w.control))
plot(Mass ~ DaysPost, data = w.control, col = "green", ylim = c(0.5,3),
     cex = 1.3, cex.axis = 1.5, cex.lab = 2, las = 1,
     xlab = "Days Post-Exposure", ylab = "Mass (g)")
TukeyHSD(aov(Mass ~ DaysPost, data = w.control))

anova(lm(BodyCondition ~ DaysPost, data = w.control))
plot(BodyCondition ~ DaysPost, data = w.control, col = "green", 
     ylim = c(50,200), cex = 1.3, cex.axis = 1.5, cex.lab = 2, las = 1, 
     xlab = "Days Post-Exposure", ylab = "Body Condition")
TukeyHSD(aov(BodyCondition ~ DaysPost, data = w.control))


anova(lm(Mass ~ DaysPost, data = w.carter))
plot(Mass~ DaysPost, data = w.carter, col = "blue", ylim = c(0.5,3),
     cex = 1.3, cex.axis = 1.5, cex.lab = 2, las = 1, 
     xlab = "Days Post-Exposure", ylab = "Mass (g)")
TukeyHSD(aov(Mass ~ DaysPost, data = w.carter))

anova(lm(BodyCondition ~ DaysPost, data = w.carter))
plot(BodyCondition ~ DaysPost, data = w.carter, col = "blue", 
     ylim = c(50,200), cex = 1.3, cex.axis = 1.5, cex.lab = 2, las = 1,
     xlab = "Days Post-Exposure", ylab = "Body Condition")
TukeyHSD(aov(BodyCondition ~ DaysPost, data = w.carter))


anova(lm(Mass ~ DaysPost, data = w.section))
plot(Mass~ DaysPost, data = w.section, col = "red", ylim = c(0.5,3),
     cex = 1.3, cex.axis = 1.5, cex.lab = 2, las = 1,
     xlab = "Days Post-Exposure", ylab = "Mass (g)")
TukeyHSD(aov(Mass ~ DaysPost, data = w.section))

anova(lm(BodyCondition ~ DaysPost, data = w.section))
plot(BodyCondition ~ DaysPost, data = w.section, col = "red", 
     ylim = c(50,200), cex = 1.3, cex.axis = 1.5, cex.lab = 2, las = 1,
     xlab = "Days Post-Exposure", ylab = "Body Condition")
TukeyHSD(aov(BodyCondition ~ DaysPost, data = w.section))


anova(lm(Mass ~ DaysPost, data = w.pesection))
plot(Mass~ DaysPost, data = w.pesection, col = "orange", ylim = c(0.5,3),
     cex = 1.3, cex.axis = 1.5, cex.lab = 2, las = 1,
     xlab = "Days Post-Exposure", ylab = "Mass (g)")
TukeyHSD(aov(Mass ~ DaysPost, data = w.pesection))

anova(lm(BodyCondition ~ DaysPost, data = w.pesection))
plot(BodyCondition ~ DaysPost, data = w.pesection, col = "orange", 
     ylim = c(50,200), cex = 1.3, cex.axis = 1.5, cex.lab = 2, las = 1,
     xlab = "Days Post-Exposure", ylab = "Body Condition")
TukeyHSD(aov(BodyCondition ~ DaysPost, data = w.pesection))


anova(lm(Mass ~ DaysPost, data = b.control))
plot(Mass ~ DaysPost, data = b.control, col = "green", ylim = c(0,20),
     cex = 1.3, cex.axis = 1.5, cex.lab = 2, las = 1,
     xlab = "Days Post-Exposure", ylab = "Mass (g)")
TukeyHSD(aov(Mass ~ DaysPost, data = b.control))

anova(lm(BodyCondition ~ DaysPost, data = b.control))
plot(BodyCondition ~ DaysPost, data = b.control, col = "green", 
     ylim = c(50,150), cex = 1.3, cex.axis = 1.5, cex.lab = 2, las = 1,
     xlab = "Days Post-Exposure", ylab = "Body Condition")
TukeyHSD(aov(BodyCondition ~ DaysPost, data = b.control))


anova(lm(Mass ~ DaysPost, data = b.carter))
plot(Mass~ DaysPost, data = b.carter, col = "blue", ylim = c(0,20),
     cex = 1.3, cex.axis = 1.5, cex.lab = 2, las = 1,
     xlab = "Days Post-Exposure", ylab = "Mass (g)")
TukeyHSD(aov(Mass ~ DaysPost, data = b.carter))

anova(lm(BodyCondition ~ DaysPost, data = b.carter))
plot(BodyCondition ~ DaysPost, data = b.carter, col = "blue", 
     ylim = c(50,150))
TukeyHSD(aov(BodyCondition ~ DaysPost, data = b.carter))


anova(lm(Mass ~ DaysPost, data = b.section))
plot(Mass ~ DaysPost, data = b.section, col = "red", ylim = c(0,20),
     cex = 1.3, cex.axis = 1.5, cex.lab = 2, las = 1,
     xlab = "Days Post-Exposure", ylab = "Mass (g)")
TukeyHSD(aov(Mass ~ DaysPost, data = b.section))

anova(lm(BodyCondition ~ DaysPost, data = b.section))
plot(BodyCondition ~ DaysPost, data = b.section, col = "red", 
     ylim = c(50,150), cex = 1.3, cex.axis = 1.5, cex.lab = 2, las = 1,
     xlab = "Days Post-Exposure", ylab = "Body Condition")
TukeyHSD(aov(BodyCondition ~ DaysPost, data = b.section))

#==============================================================================


# Fit linear mixed models within treatments to get at changes in Mass
# over time. Can't fit a model with Carter Meadow bullfrogs since they only
# have pre-treatment Mass data.

m1 <- lmer(Mass ~ DaysPost + (1|ID), data = w.control)
m2 <- lmer(Mass ~ DaysPost + (1|ID), data = w.carter)
m3 <- lmer(Mass ~ DaysPost + (1|ID), data = w.section)
m4 <- lmer(Mass ~ DaysPost + (1|ID), data = w.pesection)
m5 <- lmer(Mass ~ DaysPost + (1|ID), data = b.control)
m6 <- lmer(Mass ~ DaysPost + (1|ID), data = b.section)

m1
m2
m3
m4
m5
m6

precis(m1, prob = 0.95)
precis(m2, prob = 0.95)
precis(m3, prob = 0.95)
precis(m4, prob = 0.95)
precis(m5, prob = 0.95)
precis(m6, prob = 0.95)


# Create a new dataframe summarizing whether or not Masses at each time
# period within each treatment were different from the "Pre" mass within
# that treatment. This "Direction" information is manually input using 
# the "precis" output.

sig.df <- as.data.frame(rep(c("Pre", "4", "11", "18", "25", "32",
                                       "39", "46"), 5))
colnames(sig.df)[1] <- "DaysPost"

sig.df$Species <- c(rep("Wood Frog", 24), rep("Bullfrog", 16))

sig.df$Treatment <- c(rep("Control", 8), rep("Carter", 8), 
                      rep("SectionLine", 8), rep("Control", 8),
                      rep("SectionLine", 8))

sig.df$Direction <- 
  c("none", "none", "none", "none", "none", "up", "up", "up",
    "none", "down", "down", "down", "down", "down", "down", "down",
    "none", "down", "down", "down", "down", "none", "none", "none",
    "none", "none", "none", "up", "up", "up", "up", "up",
    "none", "down", "none", "none", "none", "up", "up", "up")

# Add a row to the dataframe representing Carter Meadow bullfrogs

sig.df <- rbind(sig.df, c("Pre", "Bullfrog", "Carter", "none")) 

# Add rows to the dataframe representing previously exposed
# Section Line wood frogs

sig.df <- 
  rbind(sig.df, c("Pre", "Wood Frog", "PreviouslyExposedSectionLine", "none"))
sig.df <- 
  rbind(sig.df, c("4", "Wood Frog", "PreviouslyExposedSectionLine", "down"))


# Merge the species-specific dataframes and the "sig.df" dataframes
# so that significance information is a column in those dataframes

d.wood <- merge(d.wood, sig.df, by = c("DaysPost", "Species", "Treatment"))
d.bull <- merge(d.bull, sig.df, by = c("DaysPost", "Species", "Treatment"))

#==============================================================================


# Creating violin plots for body mass over time


# Create a list of treatment labels

treatment_names <- c("Control" = "Control", "Carter" = "Carter Meadow",
                     "SectionLine" = "Section Line",
                     "PreviouslyExposedSectionLine" = "PE Section Line")


# Generate wood frog plot, using set.seed() to maintain jittered point locations

set.seed(8)

wood_plot <- ggplot(data = d.wood, aes(x = DaysPost, y = Mass)) +
  xlab("Days Post-Exposure") + ylab("Mass (g)") + ylim(0, 3) +
  geom_violin(aes(fill = Direction)) +
  geom_jitter(width = 0.1, size = 0.8) +
  facet_wrap(~Treatment, nrow = 1, labeller = as_labeller(treatment_names)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 20, color = "black"),
        axis.title = element_text(size = 28, face = "bold"),
        strip.text.x = element_text(size = 20, face = "bold")) +
  scale_fill_manual(values = c("red", "darkgrey", "green")) +
  guides(fill = F)

tiff("../figures/Figure_2a.tiff", width = 1200, height = 600, res = 96)

grid.arrange(wood_plot, 
             top = textGrob(expression((italic(a))), 
                            gp = gpar(fontsize = 30)))

dev.off()


# Generate bullfrog plot, using set.seed() to maintain jittered point locations

set.seed(8)

bull_plot <- ggplot(data = d.bull, aes(x = DaysPost, y = Mass)) +
  xlab("Days Post-Exposure") + ylab("Mass (g)") + ylim(0, 20) +
  geom_violin(aes(fill = Direction)) +
  geom_jitter(width = 0.1, size = 0.8) +
  facet_wrap(~Treatment, nrow = 1, labeller = as_labeller(treatment_names)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        axis.text = element_text(size = 20, color = "black"),
        axis.title = element_text(size = 28, face = "bold"),
        strip.text.x = element_text(size = 20, face = "bold")) +
  scale_fill_manual(values = c("red", "darkgrey", "green")) +
  guides(fill = F)

tiff("../figures/Figure_2b.tiff", width = 1200, height = 600, res = 96)

grid.arrange(bull_plot, 
             top = textGrob(expression((italic(b))), 
                            gp = gpar(fontsize = 30)))

dev.off()

#==============================================================================


# Facetted mass plots by individual frog


w.control <- filter(w.control, Status != 3)

ggplot(data = w.control, aes(x = DaysPost, y = Mass, group = ID)) +
  xlab("Days Post-Exposure") + ylab("Mass (g)") + ylim(0,3) +
  ggtitle("Control Wood Frogs") +
  geom_point(size = 1) +
  geom_smooth(se = F, size = 2, color = "green", span = 1) +
  facet_wrap(~ID)


w.carter <- filter(w.carter, Status != 3)

ggplot(data = w.carter, aes(x = DaysPost, y = Mass, group = ID)) +
  xlab("Days Post-Exposure") + ylab("Mass (g)") + ylim(0,3) +
  ggtitle("Carter Meadow Wood Frogs") +
  geom_point(size = 1) +
  geom_smooth(se = F, size = 2, color = "blue", span = 1) +
  facet_wrap(~ID)


w.section <- filter(w.section, Status != 3)

ggplot(data = w.section, aes(x = DaysPost, y = Mass, group = ID)) +
  xlab("Days Post-Exposure") + ylab("Mass (g)") + ylim(0,3) +
  ggtitle("Section Line Wood Frogs") +
  geom_point(size = 1) +
  geom_smooth(se = F, size = 2, color = "red", span = 1) +
  facet_wrap(~ID)


w.pesection <- filter(w.pesection, Status != 3)

ggplot(data = w.pesection, aes(x = DaysPost, y = Mass, group = ID)) +
  xlab("Days Post-Exposure") + ylab("Mass (g)") + ylim(0,3) +
  ggtitle("Previously Exposed Section Line Wood Frogs") +
  geom_point(size = 1) +
  geom_smooth(se = F, size = 2, color = "red", span = 1) +
  facet_wrap(~ID)


b.control <- filter(b.control, Status != 3)

ggplot(data = b.control, aes(x = DaysPost, y = Mass, group = ID)) +
  xlab("Days Post-Exposure") + ylab("Mass (g)") + ylim(0,20) +
  ggtitle("Control Bullfrogs") +
  geom_point(size = 1) +
  geom_smooth(se = F, size = 2, color = "green", span = 1) +
  facet_wrap(~ID)


b.carter <- filter(b.carter, Status != 3)

ggplot(data = b.carter, aes(x = DaysPost, y = Mass, group = ID)) +
  xlab("Days Post-Exposure") + ylab("Mass (g)") + ylim(0,20) +
  ggtitle("Carter Meadow Bullfrogs") +
  geom_point(size = 1) +
  geom_smooth(se = F, size = 2, color = "blue", span = 1) +
  facet_wrap(~ID)


b.section <- filter(b.section, Status != 3)

ggplot(data = b.section, aes(x = DaysPost, y = Mass, group = ID)) +
  xlab("Days Post-Exposure") + ylab("Mass (g)") + ylim(0,20) +
  ggtitle("Section Line Bullfrogs") +
  geom_point(size = 1) +
  geom_smooth(se = F, size = 2, color = "red", span = 1) +
  facet_wrap(~ID)
