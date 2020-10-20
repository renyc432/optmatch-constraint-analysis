library(tidyverse)
library(bigmatch)
library(optmatch)
library(effectsize)

# constraint: nc, c, cc
# caliper: optional, default to 0
# constant: optioanl, default to 0
# data: must contain columns z and propens
# covariates: a vector of covariate (string)


bigmatch_cost = function(constraint, data, covariates, caliper=-1, constant=-1) {
  
  
  # concatenate covariates
  cov = paste(covariates,collapse="+")
  cov = str_c("z~",cov)
  
  # If no constraint
  if (constraint == "nc") {
    dist.nc = match_on(as.formula(cov), data=data)
    prm.nc = pairmatch(dist.nc, data=data)
    prm.nc.dist = prm.nc %>% 
      matched.distances(dist.nc) %>% 
      sum()
    return (prm.nc.dist)
  }

  # If caliper
  if (constraint == "c") {
    if (caliper < 0) {return ("No caliper supplied")}
    pooledsd = sd_pooled(data$propens,as.factor(data$z))
    caliper = caliper/pooledsd
    
    mhd = match_on(as.formula(cov), data = data) + caliper(match_on(z~propens, data=data), width=caliper)
    
    prm.c = pairmatch(mhd, data = data)
    prm.c.dist = prm.c %>% 
      matched.distances(mhd) %>% 
      sum()
    return (prm.c.dist)
  }
  
  # If caliper+constant
  if (constraint == "cc") {
    if (caliper < 0) {return ("No caliper supplied")}
    if (constant < 0) {return ("No constant supplied")}
    index = 1
    X= matrix()
    while (index <= length(covariates)) {
      if (index == 1) {X = cbind(data[[covariates[index]]])}
      X = cbind(X,data[[covariates[index]]])
      index = index + 1
    }
    
    mDist = smahal(data$z, data$propens, caliper = caliper, constant = constant, X=X)

    bigmatch.t = mDist$start
    bigmatch.c = mDist$end
  
    dist.cc = match_on(as.formula(cov), caliper=400, data=data)
    dist.optmatch = data.frame(cols = dist.cc@cols,rows = dist.cc@rows,dist = dist.cc@.Data)
  
    dist.temp = c()
    index = 1
    numT = nrow(filter(data, z==1))
    while (index <= length(bigmatch.c)) {
      value = filter(dist.optmatch, cols == bigmatch.c[index]-numT, rows == bigmatch.t[index])$dist
      dist.temp = c(dist.temp, value)
      index = index+1
    }
    
    dist.cc@.Data = dist.temp
    dist.cc@cols = as.integer(bigmatch.c - numT)
    dist.cc@rows = as.integer(bigmatch.t)
  
    prm.cc = pairmatch(dist.cc, data=data)
    prm.cc.dist = prm.cc %>% 
      matched.distances(dist.cc) %>% 
      sum()
  return (prm.cc.dist)
  }
}

# ##### TESTING + DEBUGGING
# # Simulated dataset
# covariates = c("Age","Height.cm","Weight.kg","Gender","measure.1","measure.2","measure.3")
# 
# caliper = optcal(subjects.sim$z,subjects.sim$propens,rank=F)
# constant = optconstant(subjects.sim$z,subjects.sim$propens)
# 
# near(bigmatch_cost("nc",subjects.sim,covariates),379.18,tol=0.1)
# near(bigmatch_cost("c",subjects.sim,covariates, caliper$caliper),416.74,tol=0.1)
# near(bigmatch_cost("cc",subjects.sim,covariates,caliper$caliper,constant$constant),996.39,tol=0.1)
# 
# # nh0506 dataset
# covariates = c("female", "black", "hispanic", "education", "povertyr","bmi")
# 
# caliper = optcal(nh0506$z,nh0506$propens, rank=F)
# constant = optconstant(nh0506$z,nh0506$propens)
# 
# near(bigmatch_cost("nc",nh0506,covariates),174.44,tol=0.1)
# near(bigmatch_cost("c",nh0506,covariates,caliper$caliper),205.25,tol=0.1)
# near(bigmatch_cost("cc",nh0506,covariates,caliper$caliper,constant$constant),381.68,tol=0.1)
# 
# 
# # Alzheimer
# covariates = c("AXL", "Adiponectin","CD40","BMP_6","Eotaxin_3","male","age","Genotype")
# 
# caliper = optcal(subjects$z,subjects$propens, rank=F)
# constant = optconstant(subjects$z,subjects$propens)
# 
# near(bigmatch_cost("nc",subjects,covariates),149.78,tol=0.1)
# near(bigmatch_cost("c",subjects,covariates,caliper$caliper),228.2,tol=0.1)
# near(bigmatch_cost("cc",subjects,covariates,caliper$caliper,constant$constant),256.96,tol=0.1)

##### caliper+constant
# Algorithm
# 1. Calculate caliper and constant
# 2. Use smahal to eliminate possible pairings with caliper and constant
# 3. Use match_on to create an InfinitySparseMatrix to feed to pairmatch
# 4. Update the InfinitySparseMatrix so that it only includes pairings from smahal
# 5. Pairmatch using updated InfinitySparseMatrix
# 6. Caclulate cost with matched.distances

