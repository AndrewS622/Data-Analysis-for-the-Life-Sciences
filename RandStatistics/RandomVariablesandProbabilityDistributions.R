## Examine weight dataset and random sampling

library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- basename(url)
download(url, destfile=filename)
x <- unlist(read.csv(filename))

mean(x)
set.seed(1)
xs <- sample(x, 5)
mean(xs)
abs(mean(xs) - mean(x))

set.seed(5)
xs5 <- sample(x, 5)
mean(xs5)
abs(mean(xs5) - mean(x))

## Null Distribution Generation

xm <- mean(x)

set.seed(1)
count <- 0
N <- 1000
n <- 5
for (i in 1:N) {
  xs <- sample(x, n)
  xsm <- mean(xs)
  if (abs(xsm - xm) > 1) {
    count <- count + 1
  }
}
count/N

N <- 10000
set.seed(1)
count <- 0
for (i in 1:N) {
  xs <- sample(x, n)
  xsm <- mean(xs)
  if (abs(xsm - xm) > 1) {
    count <- count + 1
  }
}
count/N

set.seed(1)
count <- 0
N <- 1000
n <- 50
for (i in 1:N) {
  xs <- sample(x, n)
  xsm <- mean(xs)
  if (abs(xsm - xm) > 1) {
    count <- count + 1
  }
}
count/N

## Probability Distributions
## Exploring country and life expectancy in the gapminder dataset

library(gapminder)
data(gapminder)
head(gapminder)
library(dplyr)

x <- filter(gapminder, year == "1952") %>% dplyr::select(lifeExp) %>% unlist()
(hist(x))
mean(x <= 40)
mean(x <= 60) - mean(x <= 40)

# empirical cumulative distribution function
plot(ecdf(x))

# building ecdf manually
prop = function(q) {
  mean(x <= q)
}

prop(40)
qs = seq(from=min(x), to=max(x), length=20)
props = sapply(qs, prop)
plot(qs, props)

# or in a condensed version
props = sapply(qs, function(q) mean(x <= q))
