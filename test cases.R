library(tidyverse)
library(bigmatch)
library(optmatch)

# TEST 1 with caliper

f1 = c(1,6,2,7,8)
f2 = c(12,29,10,40,120)
z = c(1,1,0,0,0)
test1 = data.frame(f1,f2,z)
mod.propens.sim = glm(z~., family = binomial("logit"), data=test1)
test1$propens = mod.propens.sim$fitted.values

dist.test1 = match_on(z ~ f1+f2, caliper=400, data=test1)
temp = data.frame(c = dist.test1@cols, t = dist.test1@rows, cost = dist.test1@.Data)
test1.cost = (filter(temp,c==1&t==1) + filter(temp,c==2&t==2))$cost

covariates = c("f1","f2")
near(bigmatch_cost("nc",test1,covariates),test1.cost,tol=0.1)

# TEST 2 with caliper
## nc/c return the same cost
caliper = 0.15

f1 = c(1,2,3,1,3,1,2,2,1,1,3,3,3,2,3)
f2 = c(0.1,0.9,0.5,0.3,0.4,0.05,0.69,0.7,0.19,0.3,0.5,0.61,0.41,0.99,0.45)
f3 = c(155,190,178,188,172,158,185,180,162,187,179,177,169,192,178)
z=c(rep(1,5),rep(0,10))
set.seed(12432)
propens= runif(15,0,1)
test2 = data.frame(f1,f2,f3,propens,z)

dist.test2 = match_on(z ~ f1+f2+f3, caliper=400, data=test2)
temp = data.frame(c = dist.test2@cols, t = dist.test2@rows, cost = dist.test2@.Data)
test2.cost = (filter(temp,t==1&c==1) + 
                filter(temp,t==2&c==9) + 
                filter(temp,t==3&c==6) + 
                filter(temp,t==4&c==5) +
                filter(temp,t==5&c==8))$cost

covariates = c("f1","f2","f3")
bigmatch_cost("nc",test2,covariates)
near(bigmatch_cost("c",test2,covariates,caliper=caliper),test2.cost,tol=0.1)

# TEST 3 with caliper and constant
