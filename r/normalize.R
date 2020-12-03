getwd()
setwd("C:/1Rob_Documents/EUFL bonneted bats APAFR/EUFL_occupancy/data/")
library(dplyr)

minT<- read.csv("EUFL_Tmin.csv", header = TRUE)
minTnorm <- minT %>% 
  mutate_if (is.numeric, scale)
summary(minTnorm)

maxT<- read.csv("EUFL_Tmax.csv", header = TRUE)
maxTnorm <- maxT %>% 
  mutate_if (is.numeric, scale)
summary(maxTnorm)

meanT<- read.csv("EUFL_Tmean.csv", header = TRUE)
meanTnorm <- meanT %>% 
  mutate_if (is.numeric, scale)
summary(meanTnorm)

noise<- read.csv("EUFL_noise.csv", header = TRUE)
noisenorm <- noise %>% 
  mutate_if (is.numeric, scale)
summary(noisenorm)

TABR<- read.csv("EUFL_TABRnoise.csv", header = TRUE)
TABRnorm <- TABR %>% 
  mutate_if (is.numeric, scale)
summary(TABRnorm)

scoto<- read.csv("EUFL_scoto.csv", header = TRUE)
scotonorm <- scoto %>% 
  mutate_if (is.numeric, scale)
summary(scotonorm)



