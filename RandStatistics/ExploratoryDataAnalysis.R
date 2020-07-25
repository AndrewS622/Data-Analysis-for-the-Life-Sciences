## QQ (Quantile-Quantile) Plots

load("skew.RData")
dim(dat)

par(mfrow = c(3,3)) # multifigure grid filled row-by-row
for (i in 1:9) {
  qqnorm(dat[,i])
}

## Columns 4 and 9 are skewed right and left, respectively

par(mfrow=c(1,2))
hist(dat[,4])
hist(dat[,9])

## Boxplots
## EXamine counts of insects for different sprays

dat <- InsectSprays
library(dplyr)
values <- select(dat, count) %>% unlist()
factor <- select(dat, spray) %>% unlist()

boxplot(values ~ factor)
boxplot(split(values, factor))

## Spray C is most effective (smallest median)

## Explore NYC Marathon 2002 data set

library(UsingR)
data(nym.2002, package="UsingR")
times <- dplyr::select(nym.2002, time) %>% unlist()
genders <- dplyr::select(nym.2002, gender) %>% unlist()

males <- filter(nym.2002, gender == "Male")
females <- filter(nym.2002, gender == "Female")

par(mfrow=c(1,1))
boxplot(times ~ genders)

par(mfrow=c(1,2))
hist(dplyr::select(males, time)$time)
hist(dplyr::select(females,time)$time)

mean(dplyr::select(males, time) %>% unlist())
mean(dplyr::select(females, time) %>% unlist())
