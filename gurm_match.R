#' 
#+ setup, include=FALSE
knitr::opts_chunk$set(echo = TRUE, warning=FALSE)
#' 
#+ packages-and-data, message = FALSE
##setwd("~/Dropbox")
library(tidyverse)
tmp_libs  <- '~/tmp/Rlibs'
#with a nonstandard RItools, have to load before optmatch
stopifnot( packageVersion("RItools", lib.loc = tmp_libs)
          >= "0.2.0.9003" ) # i.e. i90-broom branch 
withr::with_libpaths(tmp_libs, library('RItools'), "prefix")
stopifnot( packageVersion("optmatch", lib.loc = tmp_libs)
          >= "0.9.11.9007" ) 
withr::with_libpaths(tmp_libs, library('optmatch'), "prefix")
withr::with_libpaths(tmp_libs, requireNamespace('bigmatch'), "prefix")
requireNamespace("sessioninfo")
options("optmatch_max_problem_size" = Inf)
load("~/Documents/data/Gurm-pci2007/propensity_scoring-scores.rda")
names(pscores)[names(pscores)=='em.subclasses']  <- 'em_subclasses'
#' 
#' Gurm et al 2013 discussed a bunch of matches.  Their ultimate, "chosen" match
#' used `em_subclasses` for exact matching, within these classes full matching
#' on a Euclidean distance combination of `PS.incl.hospvolume` and `PS.glmer` 
#' favoring the former, i.e.
#' 
#+ PS1
head(pscores)
pscores$ps <- pscores$PS.incl.hospvolume

pdf(file="figure/ps1.pdf", width=4, height=6)
boxplot(ps ~ factor(pr_vcd, levels=0:1, labels=c("control", "treatment")),
        data=pscores,
        xlab=NULL, main=NULL, ylab=expression(paste(X, symbol("\242"), hat(beta))),
        varwidth=TRUE, horizontal=FALSE)
dev.off()
#' Right off the bat there's a strong hint that we should think about
#' calipering out treatment group members as well as controls. 
#' Outlier information:
## ----PS2-----------------------------------------------------------------
pscores %>% filter(ps>5)
hi_outlier <- row.names(pscores)[pscores$ps>5]
stopifnot(as.logical(pscores[hi_outlier, "pr_vcd"]),
          pscores[hi_outlier, "em_subclasses"]=='med.0.0.1.1.0.0')
hi_outlier_closest_C  <-
    pscores %>% filter(em_subclasses=='med.0.0.1.1.0.0', pr_vcd==0) %>%
    with(min(abs(ps - pscores[hi_outlier, "ps"])))
#'
lo_outlier  <- row.names(pscores)[pscores$ps< -15]
stopifnot(!(pscores[lo_outlier, "pr_vcd"]),
          pscores[lo_outlier, "em_subclasses"]=='med.0.1.1.1.0.1')
pscores %>% filter(em_subclasses=='med.0.0.1.1.0.0') %>% with(table(pr_vcd))
#' 
#' For reference, the canonical caliper of one fifth of a pooled SD, 
#' calculated as in optmatch using `mad()` rather than `sd()`, for outlier 
#' resistance, is:
## ------------------------------------------------------------------------
(canon_cal <- with(pscores, optmatch:::match_on_szn_scale(ps, pr_vcd))/5)

#' From separate calculations, after excluding the `hi_outlier` the
#' narrowest full matching-compatible caliper is:
source("algorithm_1.R")
minimum_cals <- algorithm_1(as.matrix(pscores[,"ps", drop=F])[,"ps"],
                            pscores$pr_vcd)
## drop the outlier
(our_cal <- sort(minimum_cals, decreasing =  TRUE)[2])
#' 

#' Most exact matching categories have more control subjects than treatments, but i 
#' a handful that isn't the case.  
#+ emcategories0-------------------------------------------------------
pscores$pr_vcd %>% tapply(pscores$em_subclasses, mean)  ->
    zbar_by_em_subclasses
stem(zbar_by_em_subclasses)

#' Takeaways for matching and related calculations below:
tx_exceed_ctls  <- levels(pscores$em_subclasses)[zbar_by_em_subclasses>.5]
pairmatch_OK  <- !(pscores$em_subclasses %in% tx_exceed_ctls)
mean_c  <- zbar_by_em_subclasses %>%
    (function(p) {(1-p)/p -.Machine$double.eps^.5 }) %>%   pmin(1)
#' 
#' # Calipers
#' 
#' ## Optimal caliper
#'
#' The optimal match in the sense of Yu, Silber & Rosenbaum can be no less than
hi_outlier_closest_C
#' 
#' The `bigmatch` package implementation of their algoritm gives
#+ optcal, eval=TRUE, cache=2-----------------------------------------------------
system.time({
  pscores[-which(row.names(pscores)==hi_outlier),] %>% 
  with( bigmatch::optcal(z=pr_vcd, p=ps, exact=em_subclasses)
       ) -> o
})
o

#' ## Cutting down problem size for demo
#'
#' <!-- Let's restrict attn to the 2 subgroups containing extreme outliers.-->
#+ echo=FALSE, eval=FALSE
pscores  <- subset(pscores,em_subclasses=='med.0.0.1.1.0.0' |
                           em_subclasses=='med.0.1.1.1.0.1' )
#' 
#' ## Ideal distance
#'
#' Here's a matching distance combining the caliper of Algorithm 1, after
#' removing the `hi_outlier`, with `em_subclasses`, which were
#' *not* folded into the application of Algorithm 1.
#' (Without the subclasses, every treatment group member except
#' the `hi_outlier` would be permitted a match according to this distance;
#' but due to exact matching within subclasses it would omit larger subset of the
#' treatment group, to be termed the `excluded`.
#+ firstISM, eval=TRUE, error=FALSE, cache=2
system.time({ 
 smaller_ism_full <- match_on(pr_vcd ~ ps + strata(em_subclasses),
                    caliper=our_cal,
                    method="euclidean", data =pscores)
})

#' sparsity:
length(smaller_ism_full)/prod(dim(smaller_ism_full))
excluded  <- setdiff(1L:nrow(smaller_ism_full),smaller_ism_full@rows)
length(excluded)
excluded  <- smaller_ism_full@rownames[excluded]

smaller_ism_full <- match_on(pr_vcd ~ ps + strata(em_subclasses),
                    caliper=our_cal, exclude=excluded,
                    method="euclidean", data =pscores)
#' Because of the `exclude=`, full matching with this second 
#' version includes the whole treatment
#' group.  The following full match is calibrated
#' (via its `mean.controls=` argument) to exclude about as much of
#' the control group as pair matching would).
#'

full_unr  <- smaller_ism_full %>% fullmatch(mean.c=mean_c, data=pscores)
summary(full_unr)
full_unr_primal <- sum(attr(full_unr, "MCFSolutions")@subproblems[,"lagrangian_value"])

#'
#' Full matching w/ restrictions.  We use a loop to determine the smallest
#' `m` making it feasible to full match under the restriction that matched
#' set sizes range from m:1 to 1:m.  (Exceptions being made for strata in
#' which treatments outnumber controls enough that the overall ratio of
#' treatments to controls exceeds m:1.)
#+ fullmatchloop, error=FALSE, cache=2
full_r  <- NULL
m  <- 1L 
matchtimes1  <- list()
while (m <10 & (is.null(full_r) || any(matchfailed(full_r), na.rm=TRUE)))
    {
        m  <- m+1L
        cat(paste("Trying m of", m, "...\n"))
        matchtime  <- system.time({
            full_r  <- fullmatch(smaller_ism_full, min.c=pmin(mean_c,1/m),
                                 mean.c=mean_c,
                                 max.c=m,
                                 data=pscores, tol=1,
                                 hint=full_r) #unclear whether this hint should help
        })
        matchtime  <- list(matchtime)
        names(matchtime)  <- m
        matchtimes1  <- c(matchtimes1, matchtime)
        }
sapply(matchtimes1, summary)
m
system.time({
full_r  <- fullmatch(smaller_ism_full, min.c=pmin(mean_c,1/m),
                     mean.c=mean_c,
                     max.c=m,
                     data=pscores, tol=0.001, #optmatch 0.9x default
                                 hint=full_r) #this hint should help
})
summary(full_r)
#' Summary info about the match
mds  <- unlist(matched.distances(full_r, smaller_ism_full))
cut(mds/our_cal, breaks=c(0, 2^(-5:1), Inf), include.lowest=T) %>%
    (function(x) table(x)/length(x))
rm(smaller_ism_full)
#' ## regret
#'
#+ echo=FALSE
regrets  <- function(matchres, ISM)
    {
        stopifnot(isTRUE(inherits(matchres, "optmatch")) || is(matchres, "MCFSolutions"))
        if (inherits(matchres, "optmatch")) matchres  <- attr(matchres, "MCFSolutions")
        if (is.null(matchres)) stop("No MCF info")
        regret_origdist  <-
                with(matchres@subproblems,
                     sum(lagrangian_value) - sum(dual_value)
                     )
        lagrangeval  <- optmatch:::evaluate_lagrangian(ISM, matchres)
        dualval  <- optmatch:::evaluate_dual(ISM, matchres)
        regret_newdist  <-  (lagrangeval - dualval)
        objective  <- sum(matchres@subproblems[,'lagrangian_value'])
        c(discretization=regret_origdist, newdist=regret_newdist,
          diff=regret_newdist - regret_origdist)/objective
        }
#' (Definition of `regrets()` suppressed from output.)
#+ echo=TRUE
args(regrets)
#' Regret against caliper constraints embedded in ISM used for matching, and
#' against an ISM w/ same exact matching categories but leaving out the caliper.
#+ bigISM, echo=TRUE, eval=TRUE, cache=2
big_ism_full <- match_on(pr_vcd ~ ps + strata(em_subclasses),
                    method="euclidean", data =pscores)
(full_r_regret <- regrets(full_r, big_ism_full))
#' (As percentages of the optimization objective value.)
#'
system.time({
big_fm_r <- fullmatch(big_ism_full, min.c=pmin(mean_c,1/m),
                      mean.c=mean_c, max.c=m,
                      data = pscores, hint=full_r)
})
rm(big_ism_full)

full_r_primal <- sum(attr(full_r, "MCFSolutions")@subproblems[,"lagrangian_value"])
big_fm_primal <- sum(attr(big_fm_r, "MCFSolutions")@subproblems[,"lagrangian_value"])
full_r_sum <- summary(full_r)
(big_fm_r_sum  <- summary(big_fm_r))
#' Improvement in restricted full matching objective from including
#' additional arcs, as compared to corresp regret bounds.  (But observe
#' that `big_fm_r` and `full_r` may have been optimized under slightly
#' different constraints if internal rounding led to differences in how
#' `mean_c` was translated into target sizes of the control group by
#' `em_subclasses` category.)
1- big_fm_primal/full_r_primal
full_r_regret

#' ## Addendum: results for a pairmatch with modified "optimal" caliper
#'
#' I.e., the Yu-Silber-Rosenbaum optimal caliper concept, but w/ informal implementation
#' and with an exception to include the `hi_outlier`. 
#'
#+ matchloop, error=FALSE, cache=2
pairs  <- NULL
multiplier  <- 1
matchtimes  <- list()
medium_ism_for_pairs  <-
    match_on(pr_vcd ~ ps + strata(em_subclasses),
                    method="euclidean", data =pscores[pairmatch_OK,],
             caliper=5*our_cal, exclude=excluded)
while (multiplier <5 & (is.null(pairs) || any(matchfailed(pairs), na.rm=TRUE)))
    {
        multiplier  <- multiplier+.1
        cat(paste("Trying multiplier of", multiplier, "...\n"))
        matchtime  <- system.time({
        pairs  <- medium_ism_for_pairs %>%
            caliper(width=multiplier*our_cal,
                    exclude=excluded, values=TRUE) %>%
            pairmatch(data=pscores[pairmatch_OK,],
                      tol=1,
                      hint=pairs
                      )
        })
        matchtime  <- list(matchtime)
        names(matchtime)  <- multiplier
        matchtimes  <- c(matchtimes, matchtime)
        }
sapply(matchtimes, summary) 
multiplier
system.time({
        pairs  <- medium_ism_for_pairs %>%
            caliper(width=multiplier*our_cal,
                    exclude=excluded, values=TRUE) %>%
            pairmatch(data=pscores[pairmatch_OK,],
                      tol=0.001,# optmatch v0.9x default
                      hint=pairs
                      )
})
summary(pairs)
rm(medium_ism_for_pairs)
#'
#' Regret, against distance used for matching and also against
#' distance w/o caliper restrictions. As fraction of optimized
#' objective value.
#'
#+ bigISMpairs, cache=2
big_ism_pairs  <-  match_on(pr_vcd ~ ps + strata(em_subclasses),
                    method="euclidean", data =pscores[pairmatch_OK,])
(pairs_regret <- regrets(pairs, big_ism_pairs))
#' 
#' # Reference information 
#'
#+ bigpm, cache=2
big_pm <-
    pairmatch(big_ism_pairs, data = pscores[pairmatch_OK,], tol=0.001, hint=pairs)
save(big_ism_pairs, file="gurm_match_big_ism_pairs.RData")
rm(big_ism_pairs)

save(full_unr, full_r, pairs, big_fm_r, big_pm, file="gurm_match_matches.RData")

pairs_primal  <- sum(attr(pairs, "MCFSolutions")@subproblems[,"lagrangian_value"])
big_pm_primal <- sum(attr(big_pm, "MCFSolutions")@subproblems[,"lagrangian_value"])


pairs_sum <- summary(pairs)
#' 
#'  Improvement in pair matching objective from including
#' additional arcs, as compared to corresp regret bounds.
#'
1- big_pm_primal/pairs_primal
pairs_regret
#'
#' Primal objective value, full vs pairs.
full_r_primal
pairs_primal
#' Average matched distances, pairs vs restricted full matching.
#'
full_r_primal/(sum(matched(full_r))-nlevels(full_r))
pairs_primal/(sum(matched(pairs))-nlevels(pairs))

save(full_r_sum, big_fm_r_sum, pairs_sum, pairs_regret, full_r_regret,
     full_r_primal, pairs_primal, big_fm_primal, big_pm_primal,
     file = "gurm_summary_regret.rda")

sessioninfo::session_info()
