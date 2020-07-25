## Scatterplot
# Analyzing NYC Marathon Data
library(UsingR)
data(nym.2002, package="UsingR")
library(dplyr)

males <- filter(nym.2002, gender == "Male")
females <- filter(nym.2002, gender == "Female")
maleage <- select(males, age) %>% unlist()
maletime <- select(males, time) %>% unlist()
femaleage <- select(females, age) %>% unlist()
femaletime <- select(females, time) %>% unlist()

cor(maleage, maletime)
cor(femaleage, femaletime)
# Correlations suggest time increases with age

# Boxplots reveal they are approx. constant until 50-60
groups <- split(maletime, round(maleage, digits = -1))
boxplot(groups)

## Log Ratios
library(UsingR)
data(nym.2002, package="UsingR")
library(dplyr)

time <- sort(nym.2002$time)
time[1]/median(time)
time[length(time)]/median(time)

# Plot on linear scale at 0.5X and 2X median
plot(time/median(time), ylim=c(1/4,4))
abline(h=c(1/2,1,2))

# Plot the same on a log scale
plot(log2(time/median(time)),ylim=c(-2,2))
abline(h=-1:1)
