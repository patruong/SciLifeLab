"MANOVA EXAMPLE"
"http://www.sthda.com/english/wiki/manova-test-in-r-multivariate-analysis-of-variance"
"IRIS Data"
"https://www.r-statistics.com/tag/iris-data-set/"


#IRIS EXAMPLE CODE

# source("https://www.r-statistics.com/wp-content/uploads/2012/01/source_https.r.txt") # Making sure we can source code from github
# source_https("https://raw.github.com/talgalili/R-code-snippets/master/clustergram.r")
# data(iris)
# library("dplyr")
# #random sampling of iris-data set
# set.seed(1234)
# dplyr::sample_n(iris, 10)
# 
# sepl <- iris$Sepal.Length
# petl <- iris$Petal.Length
# 
# # MANOVA test
# res.man <- manova(cbind(Sepal.Length, Petal.Length) ~ Species, data = iris)
# summary(res.man)
# summary.aov(res.man)
# 



rm(list = ls())
par(mfrow=c(1,1))
library(zoo) # time series object
# READIN DATA

setwd("C:/Users/Patrick/Desktop/Thesis In/LukasK")
df <- read.csv("P1902_V.csv", header = TRUE, row.names = 1)

# MANOVA test

res.man <- manova(cbind(V2, V3, V4, V5, V6, V7, V8, V9, V10, V11) ~ V1, data = df)
summary(res.man)
summary(res.man, test = "Pillai")
summary(res.man, test = "Wilks")
summary(res.man, test = "Hotelling-Lawley")
summary(res.man, test = "Roy")

summary.aov(res.man)

