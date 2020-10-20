# Learning to use bigmathc, sample codes

library(tidyverse)
library(bigmatch)
data(nh0506)
attach(nh0506)

# start = proc.time()

caliper = optcal(z, propens)
constant = optconstant(z, propens)

X=cbind(female, age, black, hispanic, education, povertyr, bmi)

# Matching with caliper, constant, and fine balance
match_sparse = nfmatch(z, propens, caliper=caliper$caliper, constant = constant$constant, dat=nh0506, rank=FALSE, X=X, fine=education)

# sparse_time = proc.time() - start

# start = proc.time()

# Matching with no constraint
match_dense = nfmatch(z, rep(1,length(z)), caliper = 1, dat=nh0506, rank=FALSE,X=X)

# dense_time = proc.time() - start


library(optmatch)

###########################################################################################################
######################### Distance comparison between sparse and dense matching ###########################
###########################################################################################################

######################### use optmatch (match_on) distance matrix ###########################


# Set up a dataframe for the distances among all nodes
distance.caliper = match_on(z ~ female + black + hispanic + education + povertyr + bmi, caliper=caliper$caliper, data=nh0506)

c = distance.caliper@cols
t = distance.caliper@rows
distance = distance.caliper@.Data
optmatch.distance = data.frame(c, t, distance)

## Divide the original table into treated and control so that treated rownum is 1:512 and control rownum is 1:1963
## If the InfinitySparseMatrix is not organized in this fashion, then this will not work
control.index = distance.caliper@colnames
treated.index = distance.caliper@rownames


## Calculate cost of sparse matching
match.sparse = match_sparse$data
# Get the row numbers of each row
match.sparse.index = as.numeric(rownames(match.sparse))

## Calculate cost of sparse matching
match.dense = match_dense$data
# Get the row numbers of each row
match.dense.index = as.numeric(rownames(match.dense))


row = 1
cost.sparse.optmatch = 0
while (row < nrow(match.sparse)){
  treated = match.sparse.index[row]
  treated.rownum = which(treated.index==treated)
  control = match.sparse.index[row+1]
  control.rownum = which(control.index==control)
  
  cost.sparse.optmatch = cost.sparse.optmatch + filter(optmatch.distance, c == control.rownum&t == treated.rownum)$distance
  row = row + 2
}


## Calculate cost of dense matching
row = 1
cost.dense.optmatch = 0
while (row < nrow(match.dense)){
  treated = match.dense.index[row]
  treated.rownum = which(treated.index==treated)
  control = match.dense.index[row+1]
  control.rownum = which(control.index==control)
  
  cost.dense.optmatch = cost.dense.optmatch + filter(optmatch.distance, c == control.rownum&t == treated.rownum)$distance
  row = row + 2
}



############################### use bigmatch (smhal) distance matrix ###########################33

## bigmatch's function for calculating distance
mDist = smahal(z,propens,caliper = caliper$caliper, X=X, constant = constant$constant)
distance = mDist$d
bigmatch.t = mDist$start
bigmatch.c = mDist$end
bigmatch.distance = data.frame(bigmatch.t,bigmatch.c,distance)


# cost of sparse matching by bigmatch distance
row = 1
cost.sparse.bigmatch = 0
while (row < nrow(match.sparse)){
  treated = match.sparse.index[row]
  treated.rownum = which(treated.index==treated)
  control = match.sparse.index[row+1]
  control.rownum = which(control.index==control) + 512
  
  cost.sparse.bigmatch = cost.sparse.bigmatch + filter(bigmatch.distance, bigmatch.c == control.rownum&bigmatch.t == treated.rownum)$distance
  row = row + 2
}


########################

# Recalculate distance matrix, with no caliper (set to large so no effect)
mDist = smahal(z,propens, X=X, caliper = 400)
distance = mDist$d
bigmatch.t = mDist$start
bigmatch.c = mDist$end
bigmatch.distance = data.frame(bigmatch.t,bigmatch.c,distance)


# cost of dense matching by bigmatch distance
row = 1
cost.dense.bigmatch = 0
while (row < nrow(match.dense)){
  treated = match.dense.index[row]
  treated.rownum = which(treated.index==treated)
  control = match.dense.index[row+1]
  control.rownum = which(control.index==control) + 512
  
  cost.dense.bigmatch = cost.dense.bigmatch + filter(bigmatch.distance, bigmatch.c == control.rownum&bigmatch.t == treated.rownum)$distance
  row = row + 2
}









###########################################################################################################
########## Determine the caliper such that sparse and dense matrix has same results #######################
###########################################################################################################

