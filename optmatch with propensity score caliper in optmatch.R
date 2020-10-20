# ## Approach 1
# aGlm <- glm(z~.-propens, family=binomial(), data=subjects.sim)
# ps1 <- match_on(aGlm)
# 
# dist.c.2 = match_on(z ~ .-propens, 
#                   caliper=caliper(ps1, .38), data=subjects.sim)
# 
# 
# # this is an optmatch caliper, not bigmatch caliper
# pscal <- ps1 + caliper(ps1, .38)
# prm.c.2 = pairmatch(pscal, data = subjects.sim)
# prm.c.2.dist = matched.distances(prm.c.2,dist.c.2) %>% 
#   sum()

## Approach 2
## how to set width? how to connect width with bigmatch caliper


ISM = match_on(z~propens,data=subjects.sim)
caliper = caliper(x=ISM, width=.38)
ISM = match_on(z~.-propens, within = caliper, data=subjects.sim)
prm.ISM = pairmatch(ISM,data=subjects.sim)
prm.c.3.dist = matched.distances(prm.ISM,ISM) %>% 
  sum()
# temp = caliper


## Approach 3
library(effectsize)

#data=subjects.sim
#caliper = optcal(subjects.sim$z,subjects.sim$propens,rank=F)
#covariates = c("Age","Height.cm","Weight.kg","Gender","measure.1","measure.2","measure.3")

data=nh0506
caliper = optcal(nh0506$z,nh0506$propens, rank=F)
covariates = c("female", "black", "hispanic", "education", "povertyr","bmi")

data = subjects
caliper = optcal(subjects$z,subjects$propens, rank=F)
covariates = c("AXL", "Adiponectin","CD40","BMP_6","Eotaxin_3","male","age","Genotype")


cov = paste(covariates,collapse="+")
cov = str_c("z~",cov)

### Pair matching within a propensity score caliper.
pooledsd = sd_pooled(data$propens,as.factor(data$z))
caliper = caliper$caliper/pooledsd

# ppty <- glm(z ~ .-(propens+z), family = binomial(), data = data)

mhd <- match_on(as.formula(cov), data = data) + caliper(match_on(z~propens, data=data), width=caliper)
#mhd <- match_on(as.formula(cov), data = subjects.sim) + a

pm2 <- pairmatch(mhd, data = data)
pm2.dist = pm2 %>% 
  matched.distances(mhd) %>% 
  sum()

### The warning messages may be due to linear separation, considering the number of variables in
### nh0506 and Alzheimer, this is possible.
### Try glmnet (adds penalty) instead of glm to resolve this and avoid warning:
### 1. glm.fit: algorithm did not converge
### 2. glm.fit: fitted probabilities numerically 0 or 1 occurred
