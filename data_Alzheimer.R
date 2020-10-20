library(AppliedPredictiveModeling)
library(tidyverse)
library(optmatch)
library(bigmatch)
setwd("C:/Users/rs/Desktop/Academics/Research/Optimal Matching/optmatch_code/optmatch_constraint")
source("bigmatch_cost.R")

# Source dataset
data(AlzheimerDisease)
subjects = data.frame(predictors)
set.seed(1234)
z = sample(c(0,1), replace = TRUE, size=333, prob = c(.75,.25))
subjects$z = z

mod.propens = glm(z~.,family=binomial("logit"),data=subjects)
subjects = mutate(subjects, propens = mod.propens$fitted.values)


# subjects: optmatch with no constraint & optmatch with caliper
covariates = c("AXL", "Adiponectin","CD40","BMP_6","Eotaxin_3","male","age","Genotype")
caliper = optcal(subjects$z,subjects$propens, rank=FALSE)
constant = optconstant(subjects$z,subjects$propens,caliper=caliper$caliper,rank=FALSE)

dist.nc.AD = bigmatch_cost("nc",subjects,covariates)
dist.c.AD = bigmatch_cost("c",subjects,covariates,caliper$caliper)
dist.cc.AD = bigmatch_cost("cc",subjects,covariates,caliper$caliper,constant$constant)



# data=subjects
# data$male = as.factor(data$male)
# 
# ggplot(data)+
#   geom_bar(aes(x=Genotype))
# 
# ggplot(data)+
#   geom_histogram(aes(x=Genotype))


