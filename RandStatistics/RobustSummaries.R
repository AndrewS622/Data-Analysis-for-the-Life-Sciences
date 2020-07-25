## Median, MAD, and Spearman Correlation
# Exploring weight of chicks over time
data("ChickWeight")
head(ChickWeight)
plot(ChickWeight$Time, ChickWeight$weight, col=ChickWeight$Diet)

# reshape to wide so each chick is a row
chick <- reshape(ChickWeight, idvar=c("Chick","Diet"), timevar="Time", 
                 direction="wide")
head(chick)
chick <- na.omit(chick)

d4 <- chick$weight.4

# Effect of outliers on mean, median, SD, and MAD
mean(c(d4,3000))/mean(d4)
median(c(d4,3000))/mean(d4)

sd(c(d4,3000))/sd(d4)
mad(c(d4,3000))/mad(d4)

d21 <- chick$weight.21

plot(d4, d21)

# Effect of outliers on correlation
cor(c(d4, 3000), c(d21, 3000))/cor(d4, d21)
cor(c(d4, 3000), c(d21, 3000), method="spearman")/cor(d4, d21, method="spearman")

## Mann-Whitney-Wilcoxon Test
library(dplyr)
x <- chick %>% filter(Diet == 1) %>%
  select(weight.4) %>% 
  unlist()
y <- chick %>% filter(Diet == 4) %>%
  select(weight.4) %>% 
  unlist()
t.test(x,y)$p.value
wilcox.test(x,y)$p.value

# Effect of outliers
t.test(c(x,200), y)$p.value
wilcox.test(c(x,200), y)$p.value

# But how effective is it when samples have no overlap?
library(rafalib)
mypar(1,3)
boxplot(x,y)
boxplot(x,y+10)
boxplot(x,y+100)

t.test(x, y+10)$statistic - t.test(x, y+100)$statistic
wilcox.test(x, y+10)$statistic - wilcox.test(x, y+100)$statistic
# Entirely unaffected when samples are non-overlapping

wilcox.test(c(1,2,3), c(4,5,6))
wilcox.test(c(1,2,3), c(400,500,600))
