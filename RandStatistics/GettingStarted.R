## Go through swirl Basic Building Blocks
## Basics - version, vectors, loops

R.version
mean(c(2.23, 3.45, 1.87, 2.11, 7.33, 18.34, 19.23))

sumval <- 0
for (i in 1:25) {
  sumval = sumval + i^2
}
sumval

## Look at cars dataset

class(cars)
nrow(cars)
names(cars)
mean(cars[,2])
which(cars[,2] == 85)

## Basics of data manipulation, random sampling

library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename <- "femaleMiceWeights.csv"
download(url, destfile=filename)

data <- read.csv("femaleMiceWeights.csv")
head(data)
data[12,2]
weights <- data$Bodyweight
weights[11]
length(weights)

idx <- which(data$Diet == "hf")
mean(weights[idx])
set.seed(1)
sample(weights[idx], 1)


## dyplr basics

library(downloader)
url="https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
filename <- basename(url)
download(url,filename)

data <- read.csv(filename)
class(data)

dim(data)
names(data)

library(dplyr)

primates <- filter(data, order=="Primates")
nrow(primates)
class(primates)

primateSleep <- select(primates, sleep_total) %>% unlist()
class(primateSleep)
mean(primateSleep)

primateSleepAvg <- filter(data, order=="Primates") %>% summarize(mean(sleep_total))
primateSleepAvg
