# Eskew et al. wood frog/American bullfrog transcriptomics study

# Infection Prevalence and Load Analyses

# This script generates: Figures 3a, 3b, 4a, 4b

#==============================================================================


# Load packages

library(dplyr)
library(ggplot2)
library(reshape)

#==============================================================================


# Read in the qPCR zoospore load data

d <- read.csv("../data/infection_data.csv", check.names = F)

# Convert qPCR data to raw zoospore counts per swab 
# To do so, multiply the relevant data by 160

d[, 9:ncol(d)] <- d[, 9:ncol(d)]*160

# Convert zoospore counts. Add 1 and take the log base 10 of that measurement.

d[, 9:ncol(d)] <- log10(d[, 9:ncol(d)] + 1)


# Subset the data based on treatments

w.control <- filter(d, Treatment == "Control" & Species == "Wood Frog")
w.section <- filter(d, Treatment == "SectionLine" & Species == "Wood Frog")
w.presection <- 
  filter(d, Treatment == "PreviouslyExposedSectionLine" & 
           Species == "Wood Frog")
w.carter <- filter(d, Treatment == "Carter" & Species == "Wood Frog")
b.control <- filter(d, Treatment == "Control" & Species == "Bullfrog")
b.section <- filter(d, Treatment == "SectionLine" & Species == "Bullfrog")
b.carter <- filter(d, Treatment == "Carter" & Species == "Bullfrog")


# Choose the columns of data to be included in the analysis. 
# Since Bd exposures occurred on different days for different treatments, 
# this became more complicated. My solution was to create these "dates" 
# vectors that include the names of the columns for the dates on which 
# data was collected for the various treatment groups. Then I use 
# sapply() statements to create the "times" vectors which list the columns in 
# the original dataframes where this data can be found. These vectors then 
# get referenced extensively in the following code. Finally, I record the 
# days post-exposure at which data was collected in the "dayspost" vector. 
# This vector is used later to plot the data at an appropriate scale.

dates <- c("10.16.14", "10.17.14", "10.20.14", "10.23.14", "10.24.14", 
           "10.31.14", "11.07.14", "11.14.14", "11.21.14", "11.28.14")
dates.wpresection <- c("10.16.14", "10.17.14", "10.20.14")
dates.wcarter <- c("10.23.14", "10.24.14", "10.27.14", "10.30.14", 
                  "10.31.14", "11.07.14", "11.14.14", "11.21.14", 
                  "11.28.14", "12.05.14")
dates.bcarter <- c("11.06.14")

# Or if you want dates without the harvesting days...

dates <- c("10.17.14", "10.24.14", "10.31.14", "11.07.14", 
           "11.14.14", "11.21.14", "11.28.14")
dates.wpresection <- c("10.17.14")
dates.wcarter <- c("10.24.14", "10.31.14", "11.07.14", "11.14.14", 
                   "11.21.14", "11.28.14", "12.05.14")
dates.bcarter <- NA

times <- sapply(dates, function(z) which(colnames(d) == z))
times.wpresection <- 
  sapply(dates.wpresection, function(z) which(colnames(d) == z))
times.wcarter <- sapply(dates.wcarter, function(z) which(colnames(d) == z))
times.bcarter <- sapply(dates.bcarter, function(z) which(colnames(d) == z))


# Set up a dayspost vector, and if you're using the dates vector WITHOUT 
# harvesting days then you need the second vector listed here

dayspost <- c(3, 4, 7, 10, 11, 18, 25, 32, 39, 46)
dayspost <- c(4, 11, 18, 25, 32, 39, 46)

#==============================================================================


# Calculate the number of swab samples collected at each time point

w.control.sample <- 
  sapply(times, function(z) sum(w.control[, z] >= 0, na.rm = T))
w.section.sample <- 
  sapply(times, function(z) sum(w.section[, z] >= 0, na.rm = T))
w.presection.sample <- 
  sapply(times.wpresection, function(z) sum(w.presection[, z] >= 0, na.rm = T))
w.carter.sample <- 
  sapply(times.wcarter, function(z) sum(w.carter[, z] >= 0, na.rm = T))

b.control.sample <- 
  sapply(times, function(z) sum(b.control[, z] >= 0, na.rm = T))
b.section.sample <- 
  sapply(times, function(z) sum(b.section[, z] >= 0, na.rm = T))
b.carter.sample <- 
  sapply(times.bcarter, function(z) sum(b.carter[, z] >= 0, na.rm = T))

w.control.sample
w.section.sample
w.presection.sample
w.carter.sample
b.control.sample
b.section.sample
b.carter.sample


# Calculate average zoospore load over time by treatment

w.control.avg <- 
  sapply(times, function(z) mean(w.control[,z], na.rm = T))
w.section.avg <- 
  sapply(times, function(z) mean(w.section[, z], na.rm = T))
w.presection.avg <- 
  sapply(times.wpresection, function(z) mean(w.presection[, z], na.rm = T))
w.carter.avg <- 
  sapply(times.wcarter, function(z) mean(w.carter[, z], na.rm = T))

b.control.avg <- 
  sapply(times, function(z) mean(b.control[, z], na.rm = T))
b.section.avg <- 
  sapply(times, function(z) mean(b.section[, z], na.rm = T))
b.carter.avg <- 
  sapply(times.bcarter, function(z) mean(b.carter[, z], na.rm = T))

w.control.avg
w.section.avg
w.presection.avg
w.carter.avg
b.control.avg
b.section.avg
b.carter.avg


# Or average zoospore load over time by treatment with only positive swabs

w.control.avg.pos <- sapply(times, function(z) 
  sum(w.control[, z], na.rm = T)/sum(w.control[, z] > 0, na.rm = T))
w.section.avg.pos <- sapply(times, function(z) 
  sum(w.section[, z], na.rm = T)/sum(w.section[, z] > 0, na.rm = T))
w.presection.avg.pos <- sapply(times.wpresection, function(z) 
  sum(w.presection[, z], na.rm = T)/sum(w.presection[, z] > 0, na.rm = T))
w.carter.avg.pos <- sapply(times.wcarter, function(z) 
  sum(w.carter[, z], na.rm = T)/sum(w.carter[, z] > 0, na.rm = T))

b.control.avg.pos <- sapply(times, function(z) 
  sum(b.control[, z], na.rm = T)/sum(b.control[, z] > 0, na.rm = T))
b.section.avg.pos <- sapply(times, function(z) 
  sum(b.section[, z], na.rm = T)/sum(b.section[, z] > 0, na.rm = T))
b.carter.avg.pos <- sapply(times.bcarter, function(z) 
  sum(b.carter[, z], na.rm = T)/sum(b.carter[, z] > 0, na.rm = T))

w.control.avg.pos <- 
  sapply(w.control.avg.pos, function(z) ifelse(is.nan(z), 0, z))
w.section.avg.pos <- 
  sapply(w.section.avg.pos, function(z) ifelse(is.nan(z), 0, z))
w.presection.avg.pos <- 
  sapply(w.presection.avg.pos, function(z) ifelse(is.nan(z), 0, z))
w.carter.avg.pos <- 
  sapply(w.carter.avg.pos, function(z) ifelse(is.nan(z), 0, z))

b.control.avg.pos <- 
  sapply(b.control.avg.pos, function(z) ifelse(is.nan(z), 0, z))
b.section.avg.pos <- 
  sapply(b.section.avg.pos, function(z) ifelse(is.nan(z), 0, z))
b.carter.avg.pos <- 
  sapply(b.carter.avg.pos, function(z) ifelse(is.nan(z), 0, z))

w.control.avg.pos
w.section.avg.pos
w.presection.avg.pos
w.carter.avg.pos
b.control.avg.pos
b.section.avg.pos
b.carter.avg.pos


# Calculate infection prevalence over time by treatment

w.control.prev <- sapply(times, function(z) 
  sum(w.control[, z] > 0, na.rm = T)/sum(w.control[, z] >= 0, na.rm = T))
w.section.prev <- sapply(times, function(z) 
  sum(w.section[, z] > 0, na.rm = T)/sum(w.section[, z] >= 0, na.rm = T))
w.presection.prev <- sapply(times.wpresection, function(z) 
  sum(w.presection[, z] > 0, na.rm = T)/sum(w.presection[, z] >= 0, na.rm = T))
w.carter.prev <- sapply(times.wcarter, function(z) 
  sum(w.carter[, z] > 0, na.rm = T)/sum(w.carter[, z] >= 0, na.rm = T))

b.control.prev <- sapply(times, function(z) 
  sum(b.control[, z] > 0, na.rm = T)/sum(b.control[, z] >= 0, na.rm = T))
b.section.prev <- sapply(times, function(z) 
  sum(b.section[, z] > 0, na.rm = T)/sum(b.section[, z] >= 0, na.rm = T))
b.carter.prev <- sapply(times.bcarter, function(z) 
  sum(b.carter[, z] > 0, na.rm = T)/sum(b.carter[, z] >= 0, na.rm = T))

w.control.prev
w.section.prev
w.presection.prev
w.carter.prev
b.control.prev
b.section.prev
b.carter.prev


# Calculate lower confidence interval on the prevalence value

w.control.prev.lower <- sapply(times, function(z) 
  binom.test(sum(w.control[, z] > 0, na.rm = T), 
             sum(w.control[, z] >= 0, na.rm = T))$conf.int[1])
w.section.prev.lower <- sapply(times, function(z) 
  binom.test(sum(w.section[, z] > 0, na.rm = T), 
             sum(w.section[, z] >= 0, na.rm = T))$conf.int[1])
w.presection.prev.lower <- sapply(times.wpresection, function(z) 
  binom.test(sum(w.presection[, z] > 0, na.rm = T), 
             sum(w.presection[, z] >= 0, na.rm = T))$conf.int[1])
w.carter.prev.lower <- sapply(times.wcarter, function(z) 
  binom.test(sum(w.carter[, z] > 0, na.rm = T), 
             sum(w.carter[, z] >= 0, na.rm = T))$conf.int[1])

b.control.prev.lower <- sapply(times, function(z) 
  binom.test(sum(b.control[, z] > 0, na.rm = T), 
             sum(b.control[, z] >= 0, na.rm = T))$conf.int[1])
b.section.prev.lower <- sapply(times, function(z) 
  binom.test(sum(b.section[, z] > 0, na.rm = T), 
             sum(b.section[, z] >= 0, na.rm = T))$conf.int[1])
b.carter.prev.lower <- sapply(times.bcarter, function(z) 
  binom.test(sum(b.carter[, z] > 0, na.rm = T), 
             sum(b.carter[, z] >= 0, na.rm = T))$conf.int[1])


# Calculate upper confidence interval on the prevalence value

w.control.prev.upper <- sapply(times, function(z) 
  binom.test(sum(w.control[, z] > 0, na.rm = T), 
             sum(w.control[, z] >= 0, na.rm = T))$conf.int[2])
w.section.prev.upper <- sapply(times, function(z) 
  binom.test(sum(w.section[, z] > 0, na.rm = T), 
             sum(w.section[, z] >= 0, na.rm = T))$conf.int[2])
w.presection.prev.upper <- sapply(times.wpresection, function(z) 
  binom.test(sum(w.presection[, z] > 0, na.rm = T), 
             sum(w.presection[, z] >= 0, na.rm = T))$conf.int[2])
w.carter.prev.upper <- sapply(times.wcarter, function(z) 
  binom.test(sum(w.carter[, z] > 0, na.rm = T), 
             sum(w.carter[, z] >= 0, na.rm = T))$conf.int[2])

b.control.prev.upper <- sapply(times, function(z) 
  binom.test(sum(b.control[, z] > 0, na.rm = T), 
             sum(b.control[, z] >= 0, na.rm = T))$conf.int[2])
b.section.prev.upper <- sapply(times, function(z) 
  binom.test(sum(b.section[, z] > 0, na.rm = T), 
             sum(b.section[, z] >= 0, na.rm = T))$conf.int[2])
b.carter.prev.upper <- sapply(times.bcarter, function(z) 
  binom.test(sum(b.carter[, z] > 0, na.rm = T), 
             sum(b.carter[, z] >= 0, na.rm = T))$conf.int[2])

#==============================================================================


# Infection prevalence plots


# Infection prevalence plot for wood frogs

tiff("../figures/Figure_3a.tiff", width = 800, height = 600, res = 96)

plot(w.control.prev ~ dayspost, xlim = c(0, 50), ylim = c(0,1), 
     bty = "l", xlab = "Days Post-Exposure", ylab = "Infection Prevalence", 
     cex.axis = 1.2, cex.lab = 1.4, las = 1, type = "n")

points(w.control.prev ~ dayspost, col = "green", pch = 16, cex = 2)
polygon.x <- c(dayspost, rev(dayspost))
polygon.y <- c(w.control.prev.lower, rev(w.control.prev.upper))
lines(w.control.prev ~ dayspost, col = "green", lty = 1)

points(w.carter.prev ~ dayspost, col = "blue", pch = 16, cex = 2)
polygon.x <- c(dayspost, rev(dayspost))
polygon.y <- c(w.carter.prev.lower, rev(w.carter.prev.upper))
lines(w.carter.prev ~ dayspost, col = "blue", lty = 1)

points(w.section.prev ~ dayspost, col= "red", pch = 16, cex = 2)
polygon.x <- c(dayspost, rev(dayspost))
polygon.y <- c(w.section.prev.lower, rev(w.section.prev.upper))
lines(w.section.prev ~ dayspost, col = "red", lty = 1)

points(w.presection.prev ~ dayspost[1:length(w.presection.prev)], 
       col= "orange", pch = 16, cex = 2)
lines(w.presection.prev ~ dayspost[1:length(w.presection.prev)],
      col = "orange", lty = 1)

mtext(expression((italic(a))), cex = 2)
legend(x = 0, y = 0.35, 
       c("Control", "Carter Meadow", "Section Line", "PE Section Line"), 
       fill = c("green", "blue", "red", "orange"), bty = "n", cex = 1.2)

dev.off()


# Infection prevalence plot for bullfrogs

tiff("../figures/Figure_3b.tiff", width = 800, height = 600, res = 96)

plot(b.control.prev ~ dayspost, xlim = c(0, 50), ylim = c(0,1), 
     bty = "l", xlab = "Days Post-Exposure", ylab = "Infection Prevalence", 
     cex.axis = 1.2, cex.lab = 1.4, las = 1, type = "n")

points(b.control.prev ~ dayspost, col = "green", pch = 16, cex = 2)
polygon.x <- c(dayspost, rev(dayspost))
polygon.y <- c(b.control.prev.lower, rev(b.control.prev.upper))
lines(b.control.prev ~ dayspost, col = "green", lty = 1)

points(b.carter.prev ~ dayspost[1:length(b.carter.prev)], 
       col = "blue", pch = 16, cex = 2)
lines(b.carter.prev ~ dayspost[1:length(b.carter.prev)], 
       col = "blue", lty = 1)

points(b.section.prev ~ dayspost, col= "red", pch = 16, cex = 2)
polygon.x <- c(dayspost, rev(dayspost))
polygon.y <- c(b.section.prev.lower, rev(b.section.prev.upper))
lines(b.section.prev ~ dayspost, col= "red", lty = 1)

mtext(expression((italic(b))), cex = 2)
legend(x = 35, y = 1, c("Control", "Section Line"), 
       fill = c("green", "red"), bty = "n", cex = 1.2)

dev.off()

#==============================================================================


# Infection load plots


# Infection load plot for wood frogs

plot(w.control.avg ~ dayspost, xlim = c(0, 50), ylim = c(0,4), las = 1,
     bty = "l", cex.axis = 1.2, cex.lab = 1.4, type = "n",
     xlab = "Days Post-Exposure", ylab= "Log(Zoospore Equivalents + 1)")

points(w.control.avg ~ dayspost, col = "green", pch = 16, cex = 2)
lines(w.control.avg ~ dayspost, col = "green", lty = 1)

points(w.carter.avg ~ dayspost, col = "blue", pch = 16, cex = 2)
lines(w.carter.avg ~ dayspost, col = "blue", lty = 1)

points(w.section.avg ~ dayspost, col = "red", pch = 16, cex = 2)
lines(w.section.avg ~ dayspost, col = "red", lty = 1)

points(w.presection.avg ~ dayspost[1:length(w.presection.avg)],
      col = "orange", pch = 16, cex = 2)
lines(w.presection.avg ~ dayspost[1:length(w.presection.avg)], 
      col = "orange", lty = 1)

mtext(expression((italic(a))), cex = 2)
legend(x = 35, y = 4, 
       c("Control", "Carter Meadow", "Section Line", "PE Section Line"), 
       fill = c("green", "blue", "red", "orange"), bty = "n", cex = 1.4)


# Infection load plot for wood frogs (only positives)

tiff("../figures/Figure_4a.tiff", width = 800, height = 600, res = 96)

plot(w.control.avg.pos ~ dayspost, xlim = c(0, 50), ylim = c(0,4), las = 1,
     bty = "l", cex.axis = 1.2, cex.lab = 1.4, type = "n",
     xlab = "Days Post-Exposure", ylab= "Log(Zoospore Equivalents + 1)")

points(w.control.avg.pos ~ dayspost, col = "green", pch = 16, cex = 2)
lines(w.control.avg.pos ~ dayspost, col = "green", lty = 1)

points(w.carter.avg.pos ~ dayspost, col = "blue", pch = 16, cex = 2)
lines(w.carter.avg.pos ~ dayspost, col = "blue", lty = 1)

points(w.section.avg.pos ~ dayspost, col= "red", pch = 16, cex = 2)
lines(w.section.avg.pos ~ dayspost, col= "red", lty = 1)

points(w.presection.avg.pos ~ dayspost[1:length(w.presection.avg.pos)],
       col = "orange", pch = 16, cex = 2)
lines(w.presection.avg.pos ~ dayspost[1:length(w.presection.avg.pos)], 
      col = "orange", lty = 1)

mtext(expression((italic(a))), cex = 2)
legend(x = 35, y = 4, 
       c("Control", "Carter Meadow", "Section Line", "PE Section Line"), 
       fill = c("green", "blue", "red", "orange"), bty = "n", cex = 1.4)

dev.off()


# Infection load plot for bullfrogs

plot(b.control.avg ~ dayspost, xlim = c(0, 50), ylim = c(0,4), las = 1,
     bty = "l", cex.axis = 1.2, cex.lab = 1.4, type = "n",
     xlab = "Days Post-Exposure", ylab= "Log(Zoospore Equivalents + 1)")

points(b.control.avg ~ dayspost, col = "green", pch = 16, cex = 2)
lines(b.control.avg ~ dayspost, col = "green", lty = 1)

points(b.carter.avg ~ dayspost[1:length(b.carter.avg)], 
       col = "blue", pch = 16, cex = 2)
lines(b.carter.avg ~ dayspost[1:length(b.carter.avg)], 
       col = "blue", lty = 1)

points(b.section.avg ~ dayspost, col= "red", pch = 16, cex = 2)
lines(b.section.avg ~ dayspost, col = "red", lty = 1)

mtext(expression((italic(b))), cex = 2)
legend(x = 35, y = 4, c("Control", "Section Line"), 
       fill = c("green", "red"), bty = "n", cex = 1.4)


# Infection load plot for bullfrogs (only positives)

tiff("../figures/Figure_4b.tiff", width = 800, height = 600, res = 96)

plot(b.control.avg.pos ~ dayspost, xlim = c(0, 50), ylim = c(0,4), las = 1,
     bty = "l", cex.axis = 1.2, cex.lab = 1.4, type = "n",
     xlab = "Days Post-Exposure", ylab= "Log(Zoospore Equivalents + 1)")

points(b.control.avg.pos ~ dayspost, col = "green", pch = 16, cex = 2)
lines(b.control.avg.pos ~ dayspost, col = "green", lty = 1)

points(b.carter.avg.pos ~ dayspost[1:length(b.carter.avg.pos)], 
       col = "blue", pch = 16, cex = 2)
lines(b.carter.avg.pos ~ dayspost[1:length(b.carter.avg.pos)], 
      col = "blue", lty = 1)

points(b.section.avg.pos ~ dayspost, col= "red", pch = 16, cex = 2)
lines(b.section.avg.pos ~ dayspost, col= "red", lty = 1)

mtext(expression((italic(b))), cex = 2)
legend(x = 35, y = 4, c("Control", "Section Line"), 
       fill = c("green", "red"), bty = "n", cex = 1.4)

dev.off()

#==============================================================================


# Facetted infection load plots by individual frog


# Reload date, time, and dayspost vectors

dates <- c("10.16.14", "10.17.14", "10.20.14", "10.23.14", "10.24.14", 
           "10.31.14", "11.07.14", "11.14.14", "11.21.14", "11.28.14")
dates.wpresection <- c("10.16.14", "10.17.14", "10.20.14")
dates.wcarter <- c("10.23.14", "10.24.14", "10.27.14", "10.30.14", 
                   "10.31.14", "11.07.14", "11.14.14", "11.21.14", 
                   "11.28.14", "12.05.14")
dates.bcarter <- c("11.06.14")

times <- sapply(dates, function(z) which(colnames(d) == z))
times.wpresection <- 
  sapply(dates.wpresection, function(z) which(colnames(d) == z))
times.wcarter <- sapply(dates.wcarter, function(z) which(colnames(d) == z))
times.bcarter <- sapply(dates.bcarter, function(z) which(colnames(d) == z))

dayspost <- c(3, 4, 7, 10, 11, 18, 25, 32, 39, 46)


# Control wood frogs

new.w.control <- cbind(w.control[, 1:8], w.control[, times])
colnames(new.w.control)[9:ncol(new.w.control)] <- as.character(dayspost)
new.w.control <- melt(new.w.control, id = colnames(new.w.control)[1:8])
colnames(new.w.control)[9:10] <- c("Dayspost", "ZoosporeLoad")
new.w.control$Dayspost <- 
  as.numeric(levels(new.w.control$Dayspost))[new.w.control$Dayspost]

# To exclude sacrificed animals
new.w.control <- filter(new.w.control, Status != 3)

ggplot(data = new.w.control, aes(x = Dayspost, y = ZoosporeLoad)) +
  xlab("Days Post-Exposure") + ylab("Zoospore Load") + 
  xlim(0, 50) + ylim(0, 5) +
  ggtitle("Control Wood Frogs") +
  geom_point(size = 2) +
  geom_smooth(se = F, size = 2, color = "green", span = 1) +
  facet_wrap(~ID)


# Carter Meadow wood frogs

new.w.carter <- cbind(w.carter[, 1:8], w.carter[, times.wcarter])
colnames(new.w.carter)[9:ncol(new.w.carter)] <- as.character(dayspost)
new.w.carter <- melt(new.w.carter, id = colnames(new.w.carter)[1:8])
colnames(new.w.carter)[9:10] <- c("Dayspost", "ZoosporeLoad")
new.w.carter$Dayspost <- 
  as.numeric(levels(new.w.carter$Dayspost))[new.w.carter$Dayspost]

# To exclude sacrificed animals
new.w.carter <- filter(new.w.carter, Status != 3)

ggplot(data = new.w.carter, aes(x = Dayspost, y = ZoosporeLoad)) +
  xlab("Days Post-Exposure") + ylab("Zoospore Load") + 
  xlim(0, 50) + ylim(0, 5) +
  ggtitle("Carter Meadow Wood Frogs") +
  geom_point(size = 2) +
  geom_smooth(se = F, size = 2, color = "blue", span = 1) +
  facet_wrap(~ID)


# Section Line wood frogs

new.w.section <- cbind(w.section[, 1:8], w.section[, times])
colnames(new.w.section)[9:ncol(new.w.section)] <- as.character(dayspost)
new.w.section <- melt(new.w.section, id = colnames(new.w.section)[1:8])
colnames(new.w.section)[9:10] <- c("Dayspost", "ZoosporeLoad")
new.w.section$Dayspost <- 
  as.numeric(levels(new.w.section$Dayspost))[new.w.section$Dayspost]

# To exclude sacrificed animals
new.w.section <- filter(new.w.section, Status != 3)

ggplot(data = new.w.section, aes(x = Dayspost, y = ZoosporeLoad)) +
  xlab("Days Post-Exposure") + ylab("Zoospore Load") + 
  xlim(0, 50) + ylim(0, 5) +
  ggtitle("Section Line Wood Frogs") +
  geom_point(size = 2) +
  geom_smooth(se = F, size = 2, color = "red", span = 1) +
  facet_wrap(~ID)


# PE Section Line wood frogs

new.w.presection <- 
  cbind(w.presection[, 1:8], w.presection[, times.wpresection])
colnames(new.w.presection)[9:ncol(new.w.presection)] <- 
  as.character(dayspost[1:length(times.wpresection)])
new.w.presection <- 
  melt(new.w.presection, id = colnames(new.w.presection)[1:8])
colnames(new.w.presection)[9:10] <- c("Dayspost", "ZoosporeLoad")
new.w.presection$Dayspost <- 
  as.numeric(levels(new.w.presection$Dayspost))[new.w.presection$Dayspost]

# To exclude sacrificed animals
new.w.presection <- filter(new.w.presection, Status != 3)

ggplot(data = new.w.presection, aes(x = Dayspost, y = ZoosporeLoad)) +
  xlab("Days Post-Exposure") + ylab("Zoospore Load") + 
  xlim(0, 50) + ylim(0, 5) +
  ggtitle("Previously Exposed Section Line Wood Frogs") +
  geom_point(size = 2) +
  geom_smooth(se = F, size = 2, color = "red", span = 1) +
  facet_wrap(~ID)


# Control bullfrogs

new.b.control <- cbind(b.control[, 1:8], b.control[, times])
colnames(new.b.control)[9:ncol(new.b.control)] <- as.character(dayspost)
new.b.control <- melt(new.b.control, id = colnames(new.b.control)[1:8])
colnames(new.b.control)[9:10] <- c("Dayspost", "ZoosporeLoad")
new.b.control$Dayspost <- 
  as.numeric(levels(new.b.control$Dayspost))[new.b.control$Dayspost]

# To exclude sacrificed animals
new.b.control <- filter(new.b.control, Status != 3)

ggplot(data = new.b.control, aes(x = Dayspost, y = ZoosporeLoad)) +
  xlab("Days Post-Exposure") + ylab("Zoospore Load") + 
  xlim(0, 50) + ylim(0, 5) +
  ggtitle("Control Bullfrogs") +
  geom_point(size = 2) +
  geom_smooth(se = F, size = 2, color = "green", span = 1) +
  facet_wrap(~ID)


# Carter Meadow bullfrogs

new.b.carter <- cbind(b.carter[, 1:8], b.carter[, times.bcarter])
colnames(new.b.carter)[9:ncol(new.b.carter)] <- 
  as.character(dayspost[1:length(times.bcarter)])
new.b.carter <- melt(new.b.carter, id = colnames(new.b.carter)[1:8])
colnames(new.b.carter)[9:10] <- c("Dayspost", "ZoosporeLoad")
new.b.carter$Dayspost <- 
  as.numeric(levels(new.b.carter$Dayspost))[new.b.carter$Dayspost]

# To exclude sacrificed animals
new.b.carter <- filter(new.b.carter, Status != 3)

ggplot(data = new.b.carter, aes(x = Dayspost, y = ZoosporeLoad)) +
  xlab("Days Post-Exposure") + ylab("Zoospore Load") + 
  xlim(0, 50) + ylim(0, 5) +
  ggtitle("Carter Meadow Bullfrogs") +
  geom_point(size = 2) +
  geom_smooth(se = F, size = 2, color = "blue", span = 1) +
  facet_wrap(~ID)


# Section Line bullfrogs

new.b.section <- cbind(b.section[, 1:8], b.section[, times])
colnames(new.b.section)[9:ncol(new.b.section)] <- as.character(dayspost)
new.b.section <- melt(new.b.section, id = colnames(new.b.section)[1:8])
colnames(new.b.section)[9:10] <- c("Dayspost", "ZoosporeLoad")
new.b.section$Dayspost <- 
  as.numeric(levels(new.b.section$Dayspost))[new.b.section$Dayspost]

# To exclude sacrificed animals
new.b.section <- filter(new.b.section, Status != 3)

ggplot(data = new.b.section, aes(x = Dayspost, y = ZoosporeLoad)) +
  xlab("Days Post-Exposure") + ylab("Zoospore Load") + 
  xlim(0, 50) + ylim(0, 5) +
  ggtitle("Section Line Bullfrogs") +
  geom_point(size = 2) +
  geom_smooth(se = F, size = 2, color = "red", span = 1) +
  facet_wrap(~ID)

#==============================================================================


# Create plot for Grad Slam talk


new.b.section.2 <- filter(new.b.section, ID == 29 | ID == 108)

ggplot(data = new.b.section.2, aes(x = Dayspost, y = ZoosporeLoad)) +
  xlab("\nDays After Exposure") + ylab("Infection Level\n") + 
  xlim(0,50) + ylim(0,5) +
  geom_point(size = 5) +
  geom_smooth(se = F, size = 5, color = "red", span = 1) +
  #geom_hline(yintercept = -0.2) +
  facet_grid(ID~.) +
  theme_bw() +
  theme(strip.background = element_blank(), strip.text.y = element_blank(),
        panel.border = element_rect(size = 3),
        axis.text = element_text(size = 32, family = "Franklin Gothic Medium"),
        axis.title = element_text(size = 40, family = "Franklin Gothic Medium"))
