# Eskew et al. wood frog/American bullfrog transcriptomics study

# Survival Analyses

# This script generates: Figures 1a, 1b

#==============================================================================


# Load packages

library(survival)
library(dplyr)

#==============================================================================


# Load survival data

# Note: 
# Status == "1" indicates an animal that was euthanized because of health 
# status during the course of the study;
# Status == "2" indicates an animal that survived until the end of the study;
# Status == "3" indicates an animal that was sacrificed for tissue collection
# (and was therefore censored in survival analyses)

d <- read.csv("../data/survival_data.csv")
d$DaysSurviving <- as.numeric(d$DaysSurviving)

# Select relevant treatment groups

d <- d %>% 
  filter(Treatment != "Untreated") %>%
  droplevels()

# Subset dataframe by species for plotting purposes

dwood <- filter(d, Species == "Wood Frog")
dbull <- filter(d, Species == "Bullfrog")

#==============================================================================


# Fit survival objects with full dataframe. Status is set to "1" because 
# those represent animals that were euthanized during the course of the 
# experiment (most presumably suffering bad health from Bd infection). 
# The model fit is just an intercept model

mSurv <- Surv(d$DaysSurviving, d$Status == 1)
summary(mSurv)

mfit <- survfit(Surv(d$DaysSurviving, d$Status == 1) ~ 1)
mfit
summary(mfit)


# Fit survival objects with differences by treatment and species

mfitTreat <- survfit(
  Surv(d$DaysSurviving, d$Status == 1) ~ d$Treatment + d$Species)
mfitTreat
summary(mfitTreat)

# Calculate Chi square on these results

survdiff(Surv(d$DaysSurviving, d$Status == 1) ~ d$Treatment + d$Species)


# Fit survival objects with differences by treatment for wood frogs only

mfitTreat_wood <- survfit(
  Surv(dwood$DaysSurviving, dwood$Status == 1) ~ dwood$Treatment)
mfitTreat_wood
summary(mfitTreat_wood)

# Calculate Chi square on these results

survdiff(Surv(dwood$DaysSurviving, dwood$Status == 1) ~ dwood$Treatment)


# Verify that control wood frog survival differs significantly from 
# Section Line wood frog survival

dwood.sub1 <- filter(dwood, Treatment %in% c("Control", "SectionLine"))

survdiff(Surv(dwood.sub1$DaysSurviving, dwood.sub1$Status == 1) ~ 
           dwood.sub1$Treatment)


# Verify that control wood frog survival differs significantly from 
# previously exposed Section Line wood frog survival

dwood.sub2 <- 
  filter(dwood, Treatment %in% c("Control", "PreviouslyExposedSectionLine"))

survdiff(Surv(dwood.sub2$DaysSurviving, dwood.sub2$Status == 1) ~ 
           dwood.sub2$Treatment)


# Fit survival objects with differences by treatment for bullfrogs only

mfitTreat_bull <- survfit(
  Surv(dbull$DaysSurviving, dbull$Status == 1) ~ dbull$Treatment)
mfitTreat_bull
summary(mfitTreat_bull)

# Calculate Chi square on these results

survdiff(Surv(dbull$DaysSurviving, dbull$Status == 1) ~ dbull$Treatment)

#==============================================================================


# Plotting


# Plot all wood frog survival data at once

plot(mfitTreat_wood, col = c("blue", "green", "orange", "red"), 
     xlab = "Days Post-Exposure", ylab = "Survival")
legend(x = "bottomleft", 
       c("Control", "Carter Meadow", "Section Line", "PE Section Line"), 
       fill = c("green", "blue", "red", "orange"), bty = "n")

# Plot coded to show treatments sequentially

tiff("../figures/Figure_1a.tiff", width = 800, height = 600, res = 96)

plot(mfitTreat_wood[2, ], col = "green", lwd = 6, 
     cex.axis = 1.2, cex.lab = 1.4, las = 1, bty = "l", 
     xlab = "Days Post-Exposure", ylab = "Survival", conf.int = T)

legend(x = "bottomleft", 
       c("Control", "Carter Meadow", "Section Line", "PE Section Line"), 
       fill = c("green", "blue", "red", "orange"), bty = "n", cex = 1.2)

lines(mfitTreat_wood[1, ], col = "blue", lwd = 6, conf.int = F)

lines(mfitTreat_wood[4, ], col = "red", lwd = 6, conf.int = F)

lines(mfitTreat_wood[3, ], col = "orange", lwd = 6, conf.int = F)

mtext(expression((italic(a))), cex = 2)

dev.off()


# Plot all bullfrog survival data at once

plot(mfitTreat_bull, col = c("blue", "green", "red"), 
     xlab = "Days Post-Exposure", ylab = "Survival") 
legend(x = "bottomleft", 
       c("Control", "Carter Meadow", "Section Line"), 
       fill = c("green", "blue", "red"), bty = "n")

# Plot coded to show treatments sequentially

tiff("../figures/Figure_1b.tiff", width = 800, height = 600, res = 96)

plot(mfitTreat_bull[2, ], col = "green", lwd = 6, 
     cex.axis = 1.2, cex.lab = 1.4, las = 1, bty = "l", 
     xlab = "Days Post-Exposure", ylab = "Survival", conf.int = F) 

legend(x = "bottomleft", c("Control", "Carter Meadow", "Section Line"), 
       fill = c("green", "blue", "red"), bty = "n", cex = 1.2)

lines(mfitTreat_bull[1, ], col = "blue", lwd = 6, conf.int = F)

lines(mfitTreat_bull[3, ], col = "red", lwd = 6, conf.int = F)

mtext(expression((italic(b))), cex = 2)

dev.off()
