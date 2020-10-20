library(tidyverse)
library(optmatch)
library(bigmatch)
setwd("C:/Users/rs/Desktop/Academics/Research/Optimal Matching/optmatch_code/optmatch_constraint")
source("bigmatch_cost.R")

##################################################################
################ Simulated Dataset Small #########################
##################################################################

# Create dataset
set.seed(890)
Age = sample(30:80, replace = TRUE, size = 1500)
Height.cm = rnorm(1500, mean=170, sd = 20)
Weight.kg = rnorm(1500, mean = 150, sd = 20)
Gender = sample(c(0,1), replace = TRUE, size=1500)
measure.1 = runif(1500, 1, 10)
measure.2 = rpois(1500, 3)
measure.3 = runif(1500, 100, 1000)
z = sample(c(0,1), replace = TRUE, size=1500, prob = c(.75,.25))

subjects.sim = data.frame(z, Age, Height.cm, Weight.kg, Gender, measure.1, measure.2, measure.3)
mod.propens.sim = glm(z~., family = binomial("logit"), data=subjects.sim)
subjects.sim$propens = mod.propens.sim$fitted.values

# subjects.sim: optmatch with no constraint & optmatch with caliper
covariates = c("Age","Height.cm","Weight.kg","measure.1","measure.2","measure.3")
caliper = optcal(subjects.sim$z,subjects.sim$propens, rank = FALSE)
constant = optconstant(subjects.sim$z,subjects.sim$propens,caliper=caliper$caliper,rank=FALSE)

dist.nc = bigmatch_cost("nc",subjects.sim,covariates)
dist.c = bigmatch_cost("c",subjects.sim,covariates,caliper$caliper)
dist.cc = bigmatch_cost("cc",subjects.sim,covariates,caliper$caliper,constant$constant)


##################################################################
################ Simulated Dataset Large #########################
##################################################################

# Create dataset
set.seed(643)
n = 15000
Age = sample(30:80, replace = TRUE, size = n)
Height.cm = rnorm(n, mean=170, sd = 20)
Weight.kg = rnorm(n, mean = 150, sd = 20)
Gender = sample(c(0,1), replace = TRUE, size=n)
measure.1 = runif(n, 1, 10)
measure.2 = rpois(n, 3)
measure.3 = runif(n, 100, 1000)
z = sample(c(0,1), replace = TRUE, size=n, prob = c(.75,.25))

subjects.sim.L = data.frame(z, Age, Height.cm, Weight.kg, Gender, measure.1, measure.2, measure.3)
mod.propens.sim.L = glm(z~., family = binomial("logit"), data=subjects.sim.L)
subjects.sim.L$propens = mod.propens.sim.L$fitted.values


# subjects.sim.L: optmatch with no constraint & optmatch with caliper
covariates = c("Age","Height.cm","Gender","measure.1","measure.2","measure.3")
caliper = optcal(subjects.sim.L$z,subjects.sim.L$propens, rank = FALSE)
constant = optconstant(subjects.sim.L$z,subjects.sim.L$propens,caliper=caliper$caliper,rank=FALSE)

options("optmatch_max_problem_size" = Inf)

dist.nc.L.1 = bigmatch_cost("nc",subjects.sim.L,covariates)
dist.c.L.1 = bigmatch_cost("c",subjects.sim.L,covariates,caliper$caliper)




# data=subjects.sim.L
# data$Gender = as.factor(data$Gender)
# 
# ggplot(data)+
#   geom_bar(aes(x=Gender))
# 
# ggplot(data)+
#   geom_histogram(aes(x=measure.3))
# 
# ,binwidth = 5

###################################
## optmatch with covariate set 2 ##
###################################
covariates = c("Age","Gender","Height.cm","Weight.kg","measure.1","measure.2")
dist.nc.L.2 = bigmatch_cost("nc",subjects.sim.L,covariates)
dist.c.L.2 = bigmatch_cost("c",subjects.sim.L,covariates,caliper$caliper)

# dist.nc = 1741.792
# dist.c  = 1994.001

# This takes too long to finish (longer than 8 hours on my laptop and still not finished)
# dist.cc = bigmatch_cost("cc",subjects.sim.L,covariates,caliper$caliper,constant$constant)


##################################################################
################ Simulated Dataset Large 2 #######################
##################################################################

# Create dataset
set.seed(98236)
n = 20000
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

# subjects.sim.L: optmatch with no constraint & optmatch with caliper
## ## optmatch with covariate set 1
covariates = c("Age","Height.cm","Weight.kg","Gender","measure.1","measure.2","measure.3","measure.5")
caliper = optcal(subjects.sim.L2$z,subjects.sim.L2$propens, rank = FALSE)
#constant = optconstant(subjects.sim.L2$z,subjects.sim.L2$propens,caliper=caliper$caliper,rank=FALSE)

options("optmatch_max_problem_size" = Inf)

dist.nc.L2.1 = bigmatch_cost("nc",subjects.sim.L2,covariates)
dist.c.L2.1 = bigmatch_cost("c",subjects.sim.L2,covariates,caliper$caliper)

# dist.nc = 3566.933
# dist.c  = 3597.968

# This takes too long to finish (longer than 8 hours on my laptop and still not finished)
# dist.cc = bigmatch_cost("cc",subjects.sim.L,covariates,caliper$caliper,constant$constant)


## optmatch with covariate set 2
covariates = c("Age","Height.cm","Weight.kg","measure.1","measure.2","measure.4")
dist.nc.L2.2 = bigmatch_cost("nc",subjects.sim.L2,covariates)
dist.c.L2.2 = bigmatch_cost("c",subjects.sim.L2,covariates,caliper$caliper)

# dist.nc = 2410.067
# dist.c  = 2455.539


# "Age","Height.cm","Weight.kg","Gender","measure.1","measure.2","measure.3","measure.4","measure.5"
# 
# data=subjects.sim.L2
# data$Gender = as.factor(data$Gender)
# 
# ggplot(data)+
#   geom_bar(aes(x=Gender))
# 
# ggplot(data)+
#   geom_histogram(aes(x=Weight.kg))
# 
# ,binwidth = 5





##################################################################
################### Bootstrap Dataset ############################
##################################################################

bootreps = 100
n = 2475
covariates = c("female", "black", "hispanic", "education", "povertyr","bmi")
summary.boot = data.frame(boot=NA,cost.nc=NA,cost.c=NA)[numeric(0),]
counter = 1
while (counter <= bootreps) {
  sample.boot = sample(nh0506, size = n, replace=TRUE)
  sample.boot.caliper = optcal(sample.boot$z, sample.boot$propens, rank=FALSE)
  cost.nc = bigmatch_cost("nc", sample.boot, covariates)
  cost.c = bigmatch_cost("c", sample.boot, covariates,sample.boot.caliper$caliper)
  summary.boot = summary.boot %>% add_row(boot=counter,cost.nc=cost.nc,cost.c=cost.c)
  counter = counter + 1
}
mean(summary.boot$cost.nc)
mean(summary.boot$cost.c)

# bootreps = 100
# -> mean(cost.nc) = 174.4363
# -> mean(cost.c) = 205.2499
