library(tidyverse)
library(optmatch)
library(bigmatch)
setwd("C:/Users/rs/Desktop/Academics/Research/Optimal Matching/optmatch_code/optmatch_constraint")
source("bigmatch_cost.R")

##################################################################
################### Dataset Generator ############################
##################################################################

# Same features and distributions as subjects.sim.L2

sim.generator = function(n,seed) {
  set.seed(seed)
  Age = sample(20:80, replace = TRUE, size = n)
  Height.cm = rnorm(n, mean=175, sd = 10)
  Weight.kg = rnorm(n, mean = 166, sd = 15)
  Gender = sample(c(0,1), replace = TRUE, size=n)
  measure.1 = rpois(n, 5.43)
  measure.2 = rpois(n, 8)
  measure.3 = rnorm(n, 123, 5)
  measure.4 = rnorm(n, 56, 2)
  measure.5 = rnorm(n, 98, 8)
  z = sample(c(0,1), replace = TRUE, size=n, prob = c(.80,.20))
  
  subjects.sim.L2 = data.frame(z, Age, Height.cm, Weight.kg, Gender, measure.1, measure.2, measure.3,measure.4,measure.5)
  mod.propens.sim.L2 = glm(z~., family = binomial("logit"), data=subjects.sim.L2)
  subjects.sim.L2$propens = mod.propens.sim.L2$fitted.values
  
  return(subjects.sim.L2)
}


covariates = c("Age","Height.cm","Weight.kg","Gender","measure.1","measure.2","measure.3")
results = data.frame(spl.size=NA,cost.nc=NA,cost.c=NA,cost.caliper=NA)[numeric(0),]
spl.sizes = c(1:20)*100


for (size in spl.sizes) {
  data.sim = sim.generator(size,182123)
  caliper = optcal(data.sim$z,data.sim$propens, rank = FALSE)

  dist.nc = bigmatch_cost("nc",data.sim,covariates)
  dist.c = bigmatch_cost("c",data.sim,covariates,caliper$caliper)
  diff = (dist.c-dist.nc)/dist.c
  results = results %>% add_row(spl.size=size,cost.nc=dist.nc,cost.c=dist.c,cost.caliper=diff)
}
ggplot(results) +
  geom_point(mapping=aes(x=spl.size,y=cost.caliper))
