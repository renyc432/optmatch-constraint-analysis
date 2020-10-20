#library(devtools)
#install_github("markmfredrickson/optmatch", ref="issue54-hinting")

library(RItools)
library(optmatch)
library(tidyverse)
library(bigmatch)
library(effectsize)
source("sim_generator.R")

# trace(utils:::unpackPkgZip, edit=TRUE)
setwd("C:/Users/rs/Desktop/Academics/Research/Optimal Matching/optmatch_code/optmatch_constraint")
source("gurm_match.R")




## gurm_match.R

regrets  <- function(matchres, ISM)
{
  stopifnot(isTRUE(inherits(matchres, "optmatch")) || is(matchres, "MCFSolutions"))
  if (inherits(matchres, "optmatch")) matchres  <- attr(matchres, "MCFSolutions")
  if (is.null(matchres)) stop("No MCF info")
  
  
  
  regret_origdist  <- matchres@subproblems$lagrangian_value -
    matchres@subproblems$dual_value
  
  
  # regret_origdist  <-
  #   with(matchres@subproblems,
  #        sum(lagrangian_value) - sum(dual_value)
  #   )
  
  
  #  lagrangeval  <- optmatch:::evaluate_lagrangian(ISM, matchres)
  #  dualval  <- optmatch:::evaluate_dual(ISM, matchres)
  objective  <- sum(matchres@subproblems[,'lagrangian_value'])
  dualval  <- evaluate_dual(ISM, matchres)
  regret_newdist  <-  (objective - dualval)
  regret_diff = regret_newdist - regret_origdist
  
  if (regret_diff > objective) {regret_diff=objective} 
  #   print("#########regrets() debug#############")
  #   print(paste0("dist.c: ", 205.25))
  #   print(paste0("objective: ", objective))
  #   print(paste0("regret_diff: ", regret_diff))
     
  # #  print(paste0("lagrangeval: ", lagrangeval))
  #   print(paste0("dualval: ", dualval))
  #   print("#####################################")
  
  return (c(discretization=regret_origdist, newdist=regret_newdist,
    diff=regret_diff)/objective)
}


diff = function(data, covariates) {
  
  cov = paste(covariates,collapse="+")
  cov = str_c("z~",cov)
  
  caliper = optcal(data$z, data$propens,rank=FALSE)
  pooledsd = sd_pooled(data$propens,as.factor(data$z))
  caliper.adj = caliper$caliper/pooledsd
  
  dist.dense = match_on(as.formula(cov), data = data)
  dist.sparse = match_on(as.formula(cov), data = data) + 
    caliper(match_on(z~propens, data=data), width=caliper.adj)
  
  prm.c = pairmatch(dist.sparse, data = data)  
  diff = regrets(prm.c, dist.dense)['diff']
  
  return (diff)
}

diff.dct = function(data,covariates) {
  cov = paste(covariates,collapse="+")
  cov = str_c("z~",cov)
  
  dist.mt = match_on(as.formula(cov),data = data)
  prm.c = pairmatch(dist.mt,data=data)
  diff = regrets(prm.c,dist.mt)['discretization']
  return (diff)
}




# nh0506
data=nh0506
covariates = c("female", "black", "hispanic", "education", "povertyr","bmi")
regrets.nh0506 = diff(nh0506,covariates)
regrets.dct.nh0506 = diff.dct(nh0506,covariates)

# AD
data=subjects
covariates = c("AXL", "Adiponectin","CD40","BMP_6","Eotaxin_3","male","age","Genotype")
regrets.AD = diff(subjects,covariates)
regrets.dct.AD = diff.dct(subjects,covariates)


# subjects.sim
data = subjects.sim
covariates = c("Age","Height.cm","Weight.kg","measure.1","measure.2","measure.3")
regrets.sim = diff(subjects.sim,covariates)
regrets.dct.sim = diff.dct(subjects.sim,covariates)


















# subjects.sim
covariates = c("Age","Height.cm","Weight.kg","measure.1","measure.2","measure.3")
cov = paste(covariates,collapse="+")
cov = str_c("z~",cov)

results = data.frame(spl.size=NA,cost.caliper=NA)[numeric(0),]
spl.sizes = c(1:20)*100
for (size in spl.sizes) {
  data = sim.generator(size,182123)
  caliper = optcal(data$z, data$propens,rank=FALSE)
  pooledsd = sd_pooled(data$propens,as.factor(data$z))
  caliper.adj = caliper$caliper/pooledsd
  
  dist.dense = match_on(as.formula(cov), data = data)
  dist.sparse = match_on(as.formula(cov), data = data) + 
    caliper(match_on(z~propens, data=data), width=caliper.adj)
  
  prm.c = pairmatch(dist.sparse, data = data)  
  diff = regrets(prm.c, dist.dense)["diff"]
  results = results %>% add_row(spl.size=size,cost.caliper=diff)
}

ggplot(results) +
  geom_point(mapping=aes(x=spl.size,y=cost.caliper))


data= subjects.sim








## Compute regrets for different sample seeds

covariates = c("Age","Height.cm","Weight.kg","Gender","measure.1","measure.2","measure.3")
diff.regrets = c()
seeds = sample(1:10000, replace = FALSE,size=20)
for (i in seeds) {
  data=sim.generator(1500,i)
  diff.regrets = c(diff.regrets, diff(data,covariates))
}
ggplot()+
  geom_histogram(aes(x=diff.regrets),binwidth = 0.01)

diff.regrets2 = c()
spl.sizes = c(2:40)*50

## Compute regrets for different sizes
## NOTE: Try a different seed and see if the result deplicates
seed = 234241

for (size in spl.sizes) {
  data=sim.generator(size,129423)
  diff.regrets2 = c(diff.regrets2, diff(data,covariates))
}
ggplot()+
  geom_point(aes(y=diff.regrets2,x=spl.sizes))


## Compute the difference between regrets and actual cost
spl.sizes = c(1:40)*100
diff.summary = data.frame(spl.size=NA,diff.regrets3=NA,diff.regfromactual=NA)[numeric(0),]
diff.regfromactual = c()
diff.regrets3 = c()
for (size in spl.sizes) {
  data=sim.generator(size,129423)
  diff = diff2(data,covariates)
#  print(diff)
#  diff.regrets3 = c(diff.regrets3,diff[1])
#  diff.regfromactual = c(diff.regfromactual, diff[2])
  diff.summary = diff.summary %>% add_row(spl.size = size, diff.regrets3 = diff[1], diff.regfromactual = diff[2])
}
ggplot(diff.summary)+
  geom_point(aes(x=spl.sizes,y=diff.regrets3))

ggplot(diff.summary)+
  geom_point(aes(x=spl.sizes,y=diff.regfromactual))

filter(diff.summary, spl.sizes > 1000) %>% filter(diff.regfromactual<0,diff.regfromactual>-.2) %>% nrow()

diff2 = function(data, covariates) {
  
  cov = paste(covariates,collapse="+")
  cov = str_c("z~",cov)
  
  caliper = optcal(data$z, data$propens,rank=FALSE)
  pooledsd = sd_pooled(data$propens,as.factor(data$z))
  caliper.adj = caliper$caliper/pooledsd
  
  dist.dense = match_on(as.formula(cov), data = data)
  dist.sparse = match_on(as.formula(cov), data = data) + 
    caliper(match_on(z~propens, data=data), width=caliper.adj)
  
  prm.c = pairmatch(dist.sparse, data = data)
  prm.nc = pairmatch(dist.dense, data = data)
  
  dist.nc = prm.nc %>%
    matched.distances(dist.dense) %>% 
    sum()
  dist.c = prm.c %>% 
    matched.distances(dist.sparse) %>% 
    sum()
  
  diff.actual = (dist.c - dist.nc)/dist.c
  
  diff.regret = regrets2(prm.c, dist.dense)
  
#  print("***********************DEBUG***********************")
#  print(paste0("actual cost: ",diff.actual))
#  print(paste0("regret cost: ",diff.regret))
  
  
  return (c(diff.regret, (diff.regret - diff.actual)))
}


#MCFS = attr(prm.c,"MCFSolutions")
#attr(prm.c,"exceedances")

#regrets(prm.c, dist.sparse)



regrets2 = function(matchres, ISM){
  stopifnot(isTRUE(inherits(matchres, "optmatch")) || is(matchres, "MCFSolutions"))
  if (inherits(matchres, "optmatch")) matchres  <- attr(matchres, "MCFSolutions")
  if (is.null(matchres)) stop("No MCF info")
  objective  <- sum(matchres@subproblems[,'lagrangian_value'])
  dualval  <- evaluate_dual(ISM, matchres)
  
  lagrangeval  <- evaluate_lagrangian(ISM, matchres)
  
#  print(paste0("objective: ", objective))
#  print(paste0("dual: ", dualval))
#  print(paste0("lagrangianeval: ", lagrangeval))
  return ((objective-dualval)/objective)
}

regrets2(prm.c,dist.dense)

