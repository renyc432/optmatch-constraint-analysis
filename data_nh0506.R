library(tidyverse)
library(bigmatch)
library(optmatch)
setwd("C:/Users/rs/Desktop/Academics/Research/Optimal Matching/optmatch_code/optmatch_constraint")
source("bigmatch_cost.R")



#################### bigmatch caliper+constant constraint #######################

covariates = c("female", "black", "hispanic", "education", "povertyr","bmi")
caliper = optcal(nh0506$z,nh0506$propens, rank=FALSE)
constant = optconstant(nh0506$z,nh0506$propens,caliper=caliper$caliper,rank=FALSE)

dist.nc.nh0506 = bigmatch_cost("nc",nh0506,covariates)
dist.c.nh0506 = bigmatch_cost("c",nh0506,covariates,caliper$caliper)
dist.cc.nh0506 = bigmatch_cost("cc",nh0506,covariates,caliper$caliper,constant$constant)

# data=nh0506
# data$female = as.factor(female)
# data$black= as.factor(black)
# data$hispanic= as.factor(hispanic)
# 
# 
# 
# ggplot(data)+
#   geom_bar(aes(x=education))
# 
# ggplot(data)+
#   geom_histogram(aes(x=bmi),binwidth=2)
