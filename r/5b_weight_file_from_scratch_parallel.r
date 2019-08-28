


#****************************************************************************************************
#                Libraries ####
#****************************************************************************************************
library("magrittr")
library("plyr") # needed for ldply; must be loaded BEFORE dplyr
library("tidyverse")
options(tibble.print_max = 60, tibble.print_min = 60) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats

library("scales")
library("hms") # hms, for times.
library("lubridate") # lubridate, for date/times.
library("readxl") # readxl, for .xls and .xlsx files.
library("haven") # haven, for SPSS, SAS and Stata files.
library("vctrs")
library("precis")

library("tibbletime") # https://business-science.github.io/tibbletime/

library("grDevices")
library("knitr")

library("zoo") # for rollapply

# devtools::install_github("donboyd5/btools")
library("btools") # library that I created (install from github)

library("ipoptr")
library("nloptr")

# library("synthpop") # note: masks select in dplyr


#****************************************************************************************************
#                Globals ####
#****************************************************************************************************


#****************************************************************************************************
#                Includes ####
#****************************************************************************************************
source("./r/includes/globals_system_specific_boyd.r") # use a different version of this file if changing systems
source("./r/includes/globals_other.r")

source("./r/includes/functions_general.r")
source("./r/includes/functions_target_setup_and_analysis.r")
source("./r/includes/functions_ipopt.r")

# functions specific to the weighting from scratch approach:
source("./r/includes/functions_weight_from_scratch.r")


#****************************************************************************************************
#                Get puf variable names ####
#****************************************************************************************************
puf.vnames <- get_puf_vnames()


#******************************************************************************************************************
#  TEST: Get previously-prepared synfile-PUF and tax output, merge, and separate PUF and synfile ####
#******************************************************************************************************************
# CAUTION: temporary fix below to make sure the file has a variable called wt !!!! ----

#.. My files start with synthpop ----
# sfname <- "synthpop3"
# sfname <- "synthpop7"
# sfname <- "synthpop8"
sfname <- "synthpop10"
synprep <- readRDS(paste0(globals$tc.dir, sfname, "_rwprep.rds"))


#.. Max's files start with synpuf ----
sfname <- "synpuf17"
# get weight variable as it is needed below
synprep <- readRDS(paste0(globals$tc.dir, sfname, "_rwprep.rds"))
#.... here is the fix ----
synprep$tc.base <- synprep$tc.base %>% mutate(wt=wt.puf) 
#.... end of temporary fix ----

#.. Get the data and test ----
# synprep <- readRDS(paste0(globals$tc.dir, sfname, "_rwprep.rds"))
# synprep <- readRDS(paste0(globals$synd, sfname, "_all_rwprep.rds")) # this file is now in synpuf in Google Drive
names(synprep)

count(synprep$tc.base, ftype)


# now get the reforms and merge in taxbc

# merge and then split
tcvars <- c("c00100", "taxbc") # taxcalc vars
mrgdf <- left_join(synprep$tc.base, synprep$tc.output %>% dplyr::select(RECID, tcvars))
  # backwards compatibility for synthpop3
  #  %>%mutate(ftype=case_when(ftype=="puf.full" ~ "puf",
  #                       ftype=="synthpop3" ~ "syn"))
glimpse(mrgdf)
count(mrgdf, ftype)


#******************************************************************************************************************
#  pick a subset ####
#******************************************************************************************************************
mrgdf2 <- mrgdf %>%
  filter(MARS==2, c00100>=0, c00100<=25e3)

puffile <- mrgdf2 %>% filter(ftype=="puf")
synfile <- mrgdf2 %>% filter(ftype=="syn", m==1)

# synfile <- puffile %>% mutate(ftype="syn")


#******************************************************************************************************************
#  set up the "constraint" components of the objective function ####
#******************************************************************************************************************
# I refer to "constraint" components of the objective function as the components for which we want to minimize the squared difference of sum vs target
# for each component, we need:
#   the variable involved
#   the "constraint" type, i.e., one of:  n.all, sum.all, n.pos, n.neg, sum.pos, sum.neg
#   the "constraint" priority -- a multiplier of the squared diff -- the larger the mult, the more important this becomes
#   the target
#   the coefficient for each record determining how it enters into the sum that will be compared to the target
#     for value sums:
#       for the weight variable it will be 1 x the weight
#       for other variables it will be the variable times the weight
#     for numbers of returns it is simply the weight variable

# maybe make a list of objective function elements, or a data frame -- with an associated list of the coefficients
# the df would have:
#   elname, elvar, eltype, priority; link to a vector of coefficients based on elname, link to target based on elname

# we may want to see how far off we are on each constraint to help us determine priorities
# or even automate the priority setting process

# wt, c00100.val, taxbc.val, e00200.val
recipe <- read_csv("vname, fn
                   wt, n.sum
                   c00100, val.sum
                   taxbc, val.sum
                   e00200, val.sum,
                   p23250, val.neg")

# or...
recipe <- get_recipe_long(get_weighting_recipe("recipe5")) %>%
  filter(vname %in% names(puffile)) %>%
  dplyr::select(vname, vname, fn)


# start here to adjust a previously created recipe ----
tscale <- 1
recipe <- recipe %>%
  rowwise() %>%
  mutate(target=do.call(fn, list(puffile, vname, puffile$wt))) %>%
  ungroup %>%
  mutate(scale=ifelse(target!=0, abs(target / tscale), 1/ tscale),
         obj.element=paste0(vname, "_", fn)) %>%
  dplyr::select(obj.element, vname, fn, scale, target) %>%
  arrange(vname, fn)
recipe

#..weed out unnecessary elements of the recipe ----
# if the target is 0 for negative values we can drop the neg versions AS LONG AS WE HAVE NO SYNTH NEG VERSIONS
# if the val.pos and val.sum versions are identical then we can drop the val.sum version
# can drop the "neg" and "pos" versions
recipe.flagged <- recipe %>%
  rowwise() %>%
  mutate(syn.unwtd=do.call(fn, list(synfile, vname, rep(1, nrow(synfile))))) %>% # so we can check if negs!
  group_by(vname) %>%
  mutate(flag.dropneg=ifelse(str_detect(fn, "neg") & target==0 & syn.unwtd==0, 1, 0),
         flag.dropdupsum=ifelse(target==target[match("val.sum", fn)] & (fn=="val.pos"), 1, 0),
         flag.dropdupn=ifelse(target==target[match("n.sum", fn)] & (fn=="n.pos"), 1, 0)) %>%
  mutate_at(vars(starts_with("flag")), ~naz(.)) %>%
  ungroup %>%
  arrange(vname, fn)
recipe.flagged

# remove recipe elements where the target is zero
recipe.use <- recipe.flagged %>%
  filter(!(flag.dropneg | flag.dropdupsum | flag.dropdupn)) %>%
  filter(target!=0) %>%
  dplyr::select(obj.element, vname, fn, scale, target)
recipe.use

# finally, add priority weights
recipe.use <- recipe.use %>%
  mutate(priority.weight=case_when(vname %in% c("wt", "c00100", "e00200", "taxbc") ~ 100,
                                   fn %in% c("n.sum", "val.sum") ~ 100,
                                   TRUE ~ 1))  %>% 
  left_join(puf.vnames %>% dplyr::select(vname, vdesc))
recipe.use

recipe.use %>% arrange(-priority.weight)

# What would our objective function be if each targeted variable was off by a given % (as a decimal)?
# pct <- .01
# sum(recipe.use$priority.weight * (pct^2))
# # what if they were off by that same pct on average but with a random variation?
# pctv <- rnorm(nrow(recipe.use), pct, sd=.05)
# sum(recipe.use$priority.weight * (pctv^2))
# (pctv *100) %>% round(., 1)


#******************************************************************************************************************
#  prepare the input list ####
#******************************************************************************************************************
inputs <- list()
inputs$recipe <- recipe.use
inputs$synsub <- synfile[, unique(inputs$recipe$vname)] %>% mutate(wt=1)
synlong <- inputs$synsub %>%
  mutate(wtnum=row_number()) %>%
  gather(vname, value, -wtnum)

# create a data frame with one row for each weight and obj.element combination
coeffs <- expand.grid(wtnum=1:nrow(inputs$synsub), obj.element=inputs$recipe$obj.element, stringsAsFactors = FALSE) %>%
  ungroup %>%
  left_join(inputs$recipe %>% dplyr::select(obj.element, vname, fn, scale, priority.weight, target)) %>%
  left_join(synlong) %>%
  mutate(coeff=case_when(fn=="val.sum" ~ value,
                         fn=="val.pos" ~ value*(value>0),
                         fn=="val.neg" ~ value*(value<0),
                         fn=="n.sum" ~ 1,
                         fn=="n.pos" ~ 1*(value>0),
                         fn=="n.neg" ~ 1*(value<0),
                         TRUE  ~ 0)) %>%
  dplyr::select(obj.element, vname, fn, wtnum, scale, priority.weight, value, coeff, target)
# glimpse(coeffs)
# ht(coeffs)
inputs$coeffs <- coeffs


#******************************************************************************************************************
#  run optimizer(s) ####
#******************************************************************************************************************
# inputs$recipe

# bounds on the weights
xlb <- rep(1, nrow(synfile))
xub <- rep(1.5*max(puffile$wt), nrow(synfile))

# starting point:
# x0 <- (xlb + xub) / 2
x0 <- rep(mean(puffile$wt), nrow(synfile))

# alternatively:
# set.seed(1234)
# x0 <- rnorm(nrow(synfile), mean(puffile$wt), sd(puffile$wt)); x0 <- pmax(x0, xlb); x0 <- pmin(x0, xub); quantile(x0)
# set.seed(1234)
# x0 <- runif(nrow(synfile), min=1, max=1.5*max(puffile$wt))


# PRE-CHECK: Take a look at the values at the starting point
start <- inputs$recipe %>%
  rowwise() %>%
  mutate(calc=do.call(fn, list(synfile, vname, x0)),
         diff=calc - target,
         pdiff=diff / target * 100,
         apdiff=abs(pdiff),
         sdiffsq=(diff / scale)^2,
         objfn=sdiffsq * priority.weight) %>%
  ungroup
start %>% arrange(-apdiff)
start %>% arrange(-sdiffsq)
start %>% arrange(-objfn)
# END PRE-CHECK

opts <- list("print_level" = 5,
             "file_print_level" = 5, # integer
             "linear_solver" = "ma57", # mumps pardiso ma27 ma57 ma77 ma86 ma97
             "max_iter"=1000,
             # "derivative_test"="first-order",
             # "derivative_test_print_all"="yes",
             "output_file" = "scratch3.out")

a <- proc.time()
result <- ipoptr(x0 = x0,
                 lb = xlb,
                 ub = xub,
                 eval_f = eval_f_wtfs, 
                 eval_grad_f = eval_grad_f_wtfs,
                 opts = opts,
                 inputs = inputs)
b <- proc.time()
b - a


a <- proc.time()
t2 <- mma(x0, 
          fn=eval_f_wtfs, 
          gr = eval_grad_f_wtfs,
          lower=xlb, upper=xub,
          nl.info = FALSE, 
          control=list(maxeval=500),
          inputs=inputs)
b <- proc.time()
b - a # 16 mins

names(t2)
# t2$par
t2$value
t2$iter
t2$convergence
t2$message

w.sol <- t2$par




# set.seed(1234)
# x0 <- runif(nrow(synfile), min=1, max=1.5*max(puffile$wt))
# x0 <- seq(1, 1.5*max(puffile$wt), length.out=nrow(synfile))
# x0 <- rep(mean(puffile$wt), nrow(synfile))


# nloptr.print.options()
# control = list(xtol_rel = 1e-8)
nloptr.print.options(opts.show = c("algorithm"))

# iter: obj, secs -- with x0=mean
# NLOPT_LD_MMA 20: 4.7, 8 secs; 50: 0.5, 18 secs; 100: 0.07, 36 secs; 200: 0.006, 67 secs; 500: 0.00004, 172 secs
# MMA 50 iter with runif and seed: 1 0.39, 2 0.38
# MMA 100 iter with runif and seed: 1 0.06, 2 0.08, 3 .067
# MMA 200 iter with runif and seed: 1 0.0036, 2 , 3 0.005
# obj 0.0036 looks like it gets apdiff close enough for priority variables

# NLOPT_LD_LBFGS 20: 121, 8; 50: 47, 18 secs; 100, 18.9, 36 secs
# NLOPT_LD_VAR1 20: 136, 8 secs; 100 21, 37 secs
# NLOPT_LD_VAR2 20: 136, 7 secs
# NLOPT_LD_TNEWTON 20: 284, 8 secs
# NLOPT_LD_TNEWTON_PRECOND_RESTART 20: 256, 7; 100: 96, 35
# NLOPT_LD_AUGLAG:
#   NLOPT_LD_MMA 20, 4.7, 8
#   NLOPT_LN_BOBYQA 20, 737, 11 secs
#   NLOPT_LN_COBYLA
#   DO NOT DO THIS NLOPT_LD_SLSQP

# no derivative
# NLOPT_LN_BOBYQA # 20: 737, 11 secs
# NLOPT_LN_COBYLA 20: 737 9 secs

# bad
# DO NOT DO THIS NLOPT_LD_SLSQP # RStudio BLOWS UP
# NLOPT_LD_LBFGS_NOCEDAL DOESN'T WORK

# uncomment this for auglag, but don't think it's needed
# local_opts <- list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel"  = 1.0e-8) # 
# opts <- list("algorithm" = "NLOPT_LD_AUGLAG",
#              "xtol_rel"  = 1.0e-8,
#              "maxeval"   = 20,
#              "local_opts" = local_opts)


# set.seed(1)
# x0 <- runif(nrow(synfile), min=1, max=1.5*max(puffile$wt))
opts <- list("algorithm"="NLOPT_LD_MMA",
             "xtol_rel"=1.0e-8,
             "maxeval"=200)

a <- proc.time()
t3 <- nloptr(x0, 
             eval_f=eval_f_wtfs,
             eval_grad_f = eval_grad_f_wtfs,
             lb = xlb, ub = xub,
             opts = opts, inputs=inputs)
b <- proc.time()
b - a
t3$objective

# format(3.864633e-05, digits=9, scientific=FALSE)

# names(t3)
t3$objective
t3$iterations
t3$message

w.sol <- t3$solution

# a <- proc.time()
# t4 <- slsqp(x0, 
#              fn=eval_f_wtfs,
#              gr = eval_grad_f_wtfs,
#              lower = xlb, upper = xub,
#              control = list(maxeval=20), inputs=inputs)
# b <- proc.time()
# b - a
# t4$objective





#******************************************************************************************************************
#  Examine results ####
#******************************************************************************************************************
# retrieve a previous run or else use the results from above


# ------------------
w.sol <- result$solution
# w.sol <- val$solution

inputs$recipe %>% arrange(-priority.weight) %>% head(10)

comp <- inputs$recipe %>%
  rowwise() %>%
  mutate(calc=do.call(fn, list(synfile, vname, w.sol)),
         diff=calc - target,
         pdiff=diff / target * 100,
         apdiff=abs(pdiff),
         sdiffsq=(diff / scale)^2, # scaled diff sq
         objfn=sdiffsq * priority.weight) %>% # weighted sdiffsq -- the element in the objective function
  select(obj.element, vname, fn, scale, priority.weight, target, calc, diff, pdiff, apdiff, sdiffsq, objfn, vdesc)

sum(comp$objfn)
result$objective 

comp %>%
  arrange(-sdiffsq)

comp %>%
  arrange(-objfn)

comp %>%
  arrange(-apdiff)

comp %>%
  arrange(apdiff)

comp %>% filter(var %in% c("wt", 'c00100', "e00200", "taxbc"))

comp %>%
  select(obj.element, target, calc, diff, pdiff, apdiff, vdesc) %>%
  arrange(-apdiff) %>%
  kable(digits=c(0, 0, 0, 0, 3, 3, 0), format.args=list(big.mark = ','))

comp %>%
  select(obj.element, target, calc, diff, pdiff, apdiff, vdesc, priority.weight) %>%
  filter(obj.element %in% inputs$recipe$obj.element[inputs$recipe$priority.weight > 1]) %>%
  arrange(-apdiff) %>%
  kable(digits=c(0, 0, 0, 0, 3, 3, 0, 0), format.args=list(big.mark = ','))
  
quantile(w.sol, probs=0:10/10)
quantile(synfile$wt, probs=0:10/10)
quantile(puffile$wt, probs=0:10/10)

p <- bind_rows(tibble(w=puffile$wt, type="1_puf"),
          tibble(w=w.sol, type="2_weights_from_scratch"),
          tibble(w=synfile$wt, type="3_synthesized")) %>%
  ggplot(aes(w)) +
  geom_histogram(binwidth=25, fill="blue") +
  geom_vline(aes(xintercept = median(w))) +
  scale_x_continuous(breaks=seq(0, 5000, 250)) +
  theme(axis.text.x=element_text(size=8, angle=30)) +
  facet_wrap(~type, nrow=3) +
  ggtitle("Distribution of weights")
p

# ggsave("./results/optim_example_hist.png", plot=p)


#******************************************************************************************************************
# FULL file - solve problem in pieces ####
#******************************************************************************************************************
#.. My files start with synthpop ----
# sfname <- "synthpop5"
# sfname <- "synthpop6"
# sfname <- "synthpop7"
# sfname <- "synthpop8"
sfname <- "synthpop10"
synprep <- readRDS(paste0(globals$tc.dir, sfname, "_rwprep.rds"))

#.. Max's files start with synpuf ----
# sfname <- "synpuf17"
sfname <- "synpuf20"
sfname <- "synpuf20_lowmatch"
# get weight variable as it is needed below
synprep <- readRDS(paste0(globals$tc.dir, sfname, "_rwprep.rds"))
#.... here is the fix ----
synprep$tc.base <- synprep$tc.base %>%
  mutate(wt=ifelse(ftype=="puf", wt.puf, wt.syn),
         wt=ifelse(ftype=="syn", wt / 5, wt))

synprep$tc.base %>%
  group_by(ftype) %>%
  do(qtiledf(.$wt))
#.... end of temporary fix ----


names(synprep)

# now get the reforms and merge in taxbc, then merge and then split
tcvars <- c("c00100", "taxbc") # taxcalc vars
mrgdf <- left_join(synprep$tc.base, synprep$tc.output %>% dplyr::select(RECID, tcvars)) %>%
  mutate(mgroup=ifelse(MARS %in% 1:2, MARS, 3))
glimpse(mrgdf)
count(mrgdf, ftype)
count(mrgdf, ftype, mgroup, msname)
names(mrgdf) %>% sort

mrgdf %>%
  group_by(ftype, m) %>%
  summarise(n=n(), wtsum=sum(wt))

# tmp <- mrgdf %>% filter(is.na(c00100))
# count(tmp, ftype, mgroup, m)
# names(tmp)
# summary(tmp)
# 
# mrgdf %>%
#   group_by(ftype) %>%
#   summarise(vmax=max(e01200))
# anyDuplicated(mrgdf$RECID)  


# get full versions of the two files
m.num <- max(mrgdf$m)
puf <- mrgdf %>% filter(ftype=="puf")
syn <- mrgdf %>% filter(ftype=="syn") %>% mutate(wt=wt / m.num)
sum(puf$wt)
sum(syn$wt)
# note that RECID is sequential from 1 on puf to the highest value on syn and is unique


# get groups that are split by marital status and by agi in thousand dollar increments,
# then collapse so that they have at least 1,000 in both puf and syn
low.g1 <- seq(0, 80e3, 1e3)
low.g2 <-  seq(80e3, 100e3, 5e3)
low.g <- c(low.g1, low.g2)
mid.g <- seq(low.g[length(low.g)], 1e6, 10e3)
high.g <- seq(mid.g[length(mid.g)], 10e6, 100e3)
agibreaks <- c(-Inf, -1e5, -5e3, low.g, mid.g, high.g, Inf) %>% unique %>% sort
agibreaks

summary(mrgdf$c00100)

groups <- mrgdf %>%
  mutate(agibreak=cut(c00100, agibreaks),
         ibreak=as.integer(agibreak)) %>%
  group_by(ftype, mgroup, agibreak, ibreak) %>%
  summarise(n=n()) %>%
  mutate(imin=agibreaks[ibreak])
groups
groups %>% filter(ftype=="puf")
tmp <- groups %>% filter(ftype=="syn")

groupit <- function(df){
  # if group puf or syn is less than min.gcount, put the rec into the prior group
  min.gcount <- 500
  
  puf.gcount <- 0
  df$puf.gcount <- 0
  
  df$imin.new <- df$imin
  
  for(i in 1:nrow(df)){
    if(puf.gcount < min.gcount){
      puf.gcount <- puf.gcount + df$puf[i] # add this rec to group count
      df$puf.gcount[i] <- puf.gcount
      if(i > 1) df$imin.new[i] <- df$imin.new[i - 1]
    } else {
      df$puf.gcount[i] <- df$puf[i]
      puf.gcount <- df$puf[i] # set the group counter to start with this record
    }
  }
  return(df)
}

g2 <- groups %>% 
  spread(ftype, n) %>%
  group_by(mgroup) %>%
  mutate_at(vars(puf, syn), funs(naz)) %>%
  do(groupit(.))
glimpse(g2)

g3 <- g2 %>%
  group_by(mgroup, imin.new) %>%
  summarise(puf=sum(puf), syn=sum(syn))
g3

g3 %>% filter(mgroup==1)
g3 %>% filter(mgroup==2) %>% as.data.frame
g3 %>% filter(mgroup==3)

g3 %>%
  group_by(mgroup) %>%
  summarise(puf.min=min(puf), puf.max=max(puf), syn.min=min(syn), syn.max=max(syn))
  

#.. create split rules ----
split.rules <- g3 %>% 
  rename(imin=imin.new) %>%
  group_by(mgroup) %>%
  mutate(imax=lead(imin),
         imax=ifelse(is.na(imax), Inf, imax)) %>%
  ungroup %>%
  mutate(group=row_number()) %>%
  dplyr::select(group, mgroup, imin, imax, pufcount=puf, syncount=syn)
sum(split.rules$pufcount)
sum(split.rules$syncount)


# prepare the files for splitting
getgroup <- function(mgroup.in, c00100){
  split <- split.rules %>% filter(mgroup==mgroup.in[1])
  igroup.element <- function(c00100) min(which(c00100 < split$imax))
  group <- split$group[sapply(c00100, igroup.element)]
  # split$group[min(which(c00100 < split$imax))]
  return(group)
}
getgroup(1, c(-100, -1, 0, 1))

# decide on the sample ----

idfile <- mrgdf %>%
  mutate(mgroup=ifelse(MARS %in% 1:2, MARS, 3)) %>%
  group_by(ftype, mgroup) %>%
  mutate(group=getgroup(mgroup, c00100)) %>%
  ungroup %>%
  dplyr::select(ftype, RECID, mgroup, group) %>%
  arrange(RECID)
ht(idfile)
count(idfile, mgroup, ftype, group) %>% spread(ftype, n) %>% mutate(diff=syn - puf, sum=puf + syn)

count(idfile, group) %>% ht(20)

# now we are ready to run the file in pieces


rungroup <- function(group.ind){
  a <- proc.time()
  
  getrec <- function(puf, syn, recipe, puf.vnames){
    # adjust a previously created recipe
    tscale <- 1
    
    # use get rather than do.call to find a function as do.call does not seem to work in parallel
    recipe$target <- NA_real_
    for(i in 1:nrow(recipe)){
      recipe$target[i] <- get(recipe$fn[i])(df=puf, var=recipe$vname[i], puf$wt)
    }
    
    recipe <- recipe %>%
      mutate(scale=ifelse(target!=0, abs(target / tscale), 1/ tscale),
             obj.element=paste0(vname, "_", fn)) %>%
      dplyr::select(obj.element, vname, fn, scale, target) %>%
      arrange(vname, fn)

    #..weed out unnecessary elements of the recipe ----
    # if the target is 0 for negative values we can drop the neg versions AS LONG AS WE HAVE NO SYNTH NEG VERSIONS
    # if the val.pos and val.sum versions are identical then we can drop the val.sum version
    # can drop the "neg" and "pos" versions
    recipe.flagged <- recipe 
    recipe.flagged$syn.unwtd <- NA_real_
    for(i in 1:nrow(recipe.flagged)){
      recipe.flagged$syn.unwtd[i] <- get(recipe.flagged$fn[i])(df=syn, var=recipe.flagged$vname[i], weight=rep(1, nrow(syn)))
    }
    recipe.flagged <- recipe.flagged %>%
      group_by(vname) %>%
      mutate(flag.dropneg=ifelse(str_detect(fn, "neg") & target==0 & syn.unwtd==0, 1, 0),
             flag.dropdupsum=ifelse(target==target[match("val.sum", fn)] & (fn=="val.pos"), 1, 0),
             flag.dropdupn=ifelse(target==target[match("n.sum", fn)] & (fn=="n.pos"), 1, 0)) %>%
      mutate_at(vars(starts_with("flag")), funs(naz)) %>%
      ungroup %>%
      arrange(vname, fn)

    # remove recipe elements where the target is zero
    recipe.use <- recipe.flagged %>%
      filter(!(flag.dropneg | flag.dropdupsum | flag.dropdupn)) %>%
      filter(target!=0) %>%
      dplyr::select(obj.element, vname, fn, scale, target)

    # finally, add priority weights
    recipe.use <- recipe.use %>%
      mutate(priority.weight=case_when(vname %in% c("wt", "c00100", "e00200", "taxbc") ~ 100,
                                       fn %in% c("n.sum", "val.sum") ~ 100,
                                       TRUE ~ 1))  %>%
      left_join(puf.vnames %>% dplyr::select(vname, vdesc))
    return(list(recipe.use=recipe.use, recipe.flagged=recipe.flagged))
  }
  
  getinplist <- function(syn, recipe.use){
    inputs <- list()
    inputs$recipe <- recipe.use
    # inputs$synsub <- syn[, unique(inputs$recipe$vname)] %>% mutate(wt=1)
    inputs$synsub <- syn %>% dplyr::select(unique(inputs$recipe$vname)) %>% mutate(wt=1)
    synlong <- inputs$synsub %>%
      dplyr::mutate(wtnum=row_number()) %>%
      gather(vname, value, -wtnum)

    # create a data frame with one row for each weight and obj.element combination
    coeffs <- expand.grid(wtnum=1:nrow(inputs$synsub), 
                          obj.element=inputs$recipe$obj.element, stringsAsFactors = FALSE) %>%
      ungroup %>%
      left_join(inputs$recipe %>% dplyr::select(obj.element, vname, fn, scale, priority.weight, target)) %>%
      left_join(synlong) %>%
      mutate(coeff=case_when(fn=="val.sum" ~ value,
                            fn=="val.pos" ~ value*(value>0),
                            fn=="val.neg" ~ value*(value<0),
                            fn=="n.sum" ~ 1,
                            fn=="n.pos" ~ 1*(value>0),
                            fn=="n.neg" ~ 1*(value<0),
                            TRUE  ~ 0)) %>%
      dplyr::select(obj.element, vname, fn, wtnum, scale, priority.weight, value, coeff, target)

    inputs$coeffs <- coeffs
    return(inputs)
  }
  
  base <- left_join(idfile %>% filter(group==group.ind), mrgdf)
  puf <- base %>% filter(ftype=="puf")
  syn <- base %>% filter(ftype=="syn")
  recipes <- getrec(puf, syn, recipe, puf.vnames)
  recipe.use <- recipes$recipe.use
  inputs <- getinplist(syn, recipe.use)

  # bounds on the weights
  xlb <- rep(1, nrow(syn))
  xub <- rep(1.5*max(puf$wt), nrow(syn)) # FIX THIS djb

  # starting point:
  x0 <- (xlb + xub) / 2
  x0 <- x0 * sum(puf$wt / sum(x0))
  
  opts <- list("algorithm"="NLOPT_LD_MMA",
               "xtol_rel"=1.0e-8,
               "maxeval"=500)
  result <- nloptr(x0, 
               eval_f=eval_f_wtfs,
               eval_grad_f = eval_grad_f_wtfs,
               lb = xlb, ub = xub,
               opts = opts, inputs=inputs)

  # result <- mma(x0, fn=eval_f_wtfs, gr=eval_grad_f_wtfs,
  #               lower=xlb, upper=xub,
  #               nl.info = FALSE, inputs=inputs)
  
  #.. CAUTION: some pieces of result are corrupted by parallel processing ----
  # set them to null or the output files become huge!!
  # result$eval_f and result$nloptr_environment were corrupted
  result$eval_f <- NULL
  result$nloptr_environment <- NULL
  # END CAUTION ----
  
  optim <- list()
  optim$result <- result
  optim$puf <- puf
  optim$syn <- syn
  optim$inputs <- inputs
  optim$recipe.flagged <- recipes$recipe.flagged
  # 
  saveRDS(optim, paste0(globals$tc.dir, "weight_pieces/optim_group_", group.ind, ".rds"))
  
  b <- proc.time()
  print(b - a)
  return(inputs)
}



library("doParallel")
cl <- makeCluster(6)
registerDoParallel(cl)

# define recipe if new one desired
recipe <- get_recipe_long(get_weighting_recipe("recipe5")) %>%
  filter(vname %in% names(puf)) %>%
  dplyr::select(vname, vname, fn)

packages <- c("magrittr", "tidyverse", "dplyr", "nloptr")
# CAUTION:  DO NOT PASS large items as function arguments
# instead export them
xport <- c("globals", "idfile", "recipe", "puf.vnames", "mrgdf",
           "n.neg", "n.pos", "n.sum", "val.neg", "val.pos", "val.sum",
           "naz", "eval_f_wtfs", "eval_grad_f_wtfs") 
popts <- list(.packages=packages, .export=xport)
popts

# 1:max(idfile$group)
# Warning message:
#   Unknown or uninitialised column: 'target'.

a <- proc.time()
warn <- options()$warn
options(warn=-1)
d <- llply(1:max(idfile$group), .progress="text", .parallel=TRUE, .paropts=popts, .fun=rungroup)
# d1 <- llply(4:4, .progress="text", .parallel=FALSE, .paropts=popts, .fun=rungroup)
options(warn=warn)
b <- proc.time()
b - a

stopCluster(cl)
# d2 <- d
# 
# d <- c(d1, d2)



#******************************************************************************************************************
# ERROR checking ####
#******************************************************************************************************************

#.. CAUTION: RESAVE without pieces that were corrupted by parallel processing ----
resave <- function(group.ind){
  # optim$result$eval_f and optim$result$nloptr_environment were corrupted
  optim <- readRDS(paste0(globals$tc.dir, "weight_pieces/optim_group_", group.ind, ".rds"))
  optim$result$eval_f <- NULL
  optim$result$nloptr_environment <- NULL
  saveRDS(optim, paste0(globals$tc.dir, "weight_pieces/optim_group_", group.ind, ".rds"))
}
n <- 214
l_ply(1:n, resave, .progress="text")
#.. END CAUTION -----

tmp <- readRDS(paste0(globals$tc.dir, "weight_pieces/optim_group_", 4, ".rds"))

names(tmp)
object.size(tmp)
object.size(tmp$result)
object.size(tmp$puf)
object.size(tmp$syn)
object.size(tmp$inputs)
object.size(tmp$recipe.flagged)
glimpse(tmp$syn)
memory()
# saveRDS(tmp, "d:/tmp1.rds")

# hmmm...tmp$result is causing the problem
# [1] "x0"                     "eval_f"                 "lower_bounds"           "upper_bounds"           "num_constraints_ineq"  
# [6] "eval_g_ineq"            "num_constraints_eq"     "eval_g_eq"              "options"                "local_options"         
# [11] "nloptr_environment"     "call"                   "termination_conditions" "status"                 "message"               
# [16] "iterations"             "objective"              "solution"               "version"                "num.evals"   
tmp4 <- list()
tmp4$val <- tmp$result$num.evals
saveRDS(tmp4, "d:/tmp4.rds")

# tmp4.par <- tmp

gind <- 10
my1 <- readRDS(paste0("D:/tcdir/weight_pieces/synthpop8_500iter/optim_group_", gind, ".rds"))
max1 <- readRDS(paste0(globals$tc.dir, "weight_pieces/optim_group_", gind, ".rds"))

nrow(my1$puf)
nrow(max1$puf)

nrow(my1$syn)
nrow(max1$syn)

my1$result$objective
max1$result$objective

my1$result$solution %>% ht
max1$result$solution %>% ht

str(tmp$result$nloptr_environment) # tmp$result$eval_f and tmp$result$nloptr_environment got corrupted somehow

# tmp4$result <- tmp$result
# tmp4$puf <- tmp$puf
# tmp4$syn <- tmp$syn
# tmp4$inputs <- tmp$inputs
# tmp4$recipe.flagged <- tmp$recipe.flagged

n <- 214
mypiece <- function(group.ind){
  optim <- readRDS(paste0("D:/tcdir/weight_pieces/synthpop8_500iter/optim_group_", group.ind, ".rds"))
}
mylist <- llply(1:n, mypiece, .progress="text")
optlist <- llply(1:n, getpiece, .progress="text")

max.obj <- laply(1:n, function(i) optlist[[i]]$result$objective)
my.obj <- laply(1:n, function(i) mylist[[i]]$result$objective)

objcomp <- tibble(i=1:length(max.obj), synpuf17=max.obj, synthpop8=my.obj)

cor(objcomp[, c("synpuf17", "synthpop8")])

objcomp %>%
  gather(file, obj, -i) %>%
  group_by(file) %>%
  do(qtiledf(.$obj, probs=c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1))) %>%
  dplyr::select(-n) %>%
  gather(ptile, value, -file, -n.notNA) %>%
  mutate(ptile=factor(ptile, levels=unique(ptile))) %>%
  spread(file, value) %>%
  dplyr::select(ptile, n.notNA, synthpop8, synpuf17) %>%
  kable(digits=2, format.args=list(big.mark = ','))

objcomp %>%
  arrange(-pmax(synpuf17, synthpop8)) %>%
  head(25) %>%
  kable(digits=1, format.args=list(big.mark = ','))
    

probs <- c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1)
quantile(max.obj, probs)
quantile(my.obj, probs)

split.rules %>% 
  filter(group %in% c(1, 180, 67, 177, 190)) %>%
  dplyr::select(-syncount) %>%
  kable(digits=0, format.args=list(big.mark = ','))

split.rules %>% 
  filter(group %in% c(209, 210, 57, 58, 124)) %>%
  dplyr::select(-syncount) %>%
  kable(digits=0, format.args=list(big.mark = ','))

#..END error checking ----


#******************************************************************************************************************
# Construct and save full file ####
#******************************************************************************************************************
getpiece <- function(group.ind){
  optim <- readRDS(paste0(globals$tc.dir, "weight_pieces/optim_group_", group.ind, ".rds"))
}

n <- 214
optlist <- llply(1:n, getpiece, .progress="text")
memory()

# optlist <- llply(2:3, getpiece, .progress="text")
# i <- 1; object.size(optlist[[i]]) / 1e6

length(optlist)
i <- 1
# i <- 2
names(optlist[[i]])
names(optlist[[i]]$inputs)
optlist[[i]]$inputs$recipe
optlist[[i]]$recipe.flagged
unique(optlist[[i]]$recipe.flagged$vname)
names(optlist[[i]]$result)
names(optlist[[i]]$syn)
glimpse(optlist[[i]]$syn)
# optlist[[i]]$syn$RECID %>% ht
# anyDuplicated(c(optlist[[1]]$syn$RECID, optlist[[2]]$syn$RECID))


# analyze summary result
# obj.vals <- laply(1:n, function(i) optlist[[i]]$result$value)
obj.vals <- laply(1:n, function(i) optlist[[i]]$result$objective)
# obj.vals <- {lapply(1:n, function(i) optlist[[i]]$result$value)} %>% unlist
ht(obj.vals)
quantile(obj.vals, probs=c(0, .01, .05, .25, .5, .75, .9, .95, .99, .995, 1))
obj.vals %>% round(3)
obj.vals %>% sort %>% round(3) 

message <- laply(1:n, function(i) optlist[[i]]$result$message)
count(tibble(message), message)

table(laply(1:n, function(i) optlist[[i]]$result$iterations))
# table(laply(1:n, function(i) optlist[[i]]$result$convergence))

count(idfile, mgroup, ftype, group) %>% spread(ftype, n) %>% mutate(diff=syn - puf, sum=puf + syn)

# match the groups up with their obj value
objdf <- split.rules %>%
  mutate(obj=obj.vals)

objdf %>% 
  kable(digits=c(rep(0, 6), 2), format.args=list(big.mark = ','))

objdf %>% 
  arrange(-obj) %>%
  kable(digits=c(rep(0, 6), 2), format.args=list(big.mark = ',')) %>%
  kable_styling()

# 10 worst groups:
#   | group| mgroup|      imin|      imax| pufcount| syncount|      obj|
#   |-----:|------:|---------:|---------:|--------:|--------:|--------:|
#   |   147|      3| 5,400,000|       Inf|      169|      502| 7,889.50|
#   |     5|      1|     4,000|     6,000|    2,592|    7,525| 4,783.01|
#   |     4|      1|     2,000|     4,000|    2,635|    8,035| 4,066.38|
#   |    21|      1|    36,000|    38,000|      962|    2,617| 2,334.46|
#   |   129|      3|    12,000|    14,000|      983|    2,853|   940.82|
#   |     1|      1|      -Inf|    -5,000|      994|    3,914|   888.06|
#   |    45|      2|      -Inf|  -100,000|    2,007|    6,873|   329.25|
#   |    50|      2|    14,000|    18,000|    1,243|    3,640|   310.03|
#   |    47|      2|    -5,000|     4,000|      807|    4,081|   212.88|
#   |   131|      3|    16,000|    18,000|      978|    2,661|   108.57|
# about another 10 with obj in [1, 100]; rest are obj < 1
# esp. bad groups: mars 1 >$5.4m; mars 1 $4-6k, mars 1 $2-4k, mars 1, $36-38k


# aggregate the pieces of the puf and synthetic files, and attach new weights
# first puf
puf.agg <- ldply(1:n, function(i) optlist[[i]]$puf)
names(puf.agg)
ht(puf.agg[, c(1:5, (ncol(puf.agg)-5):ncol(puf.agg))]) # NOTE that RECIDs are not in order
min(puf.agg$RECID); max(puf.agg$RECID)
puf.agg <- puf.agg %>% arrange(RECID)


# now synthetic
# syn.agg <- ldply(1:n, function(i) {optlist[[i]]$syn %>% mutate(wt.wtfs=optlist[[i]]$result$par)})
syn.agg <- ldply(1:n, function(i) {optlist[[i]]$syn %>% mutate(wt.wtfs=optlist[[i]]$result$solution)})
names(syn.agg) %>% sort
ht(syn.agg[, c(1:5, ncol(syn.agg))]) # RECIDs not in order here, either
syn.agg <- syn.agg %>%
  mutate(wt.syn=wt, wt=wt.wtfs, ftype="syn")

# save the wtfs file; also save puf counterpart
#sfname <- "synthpop10"
sfname <- "synpuf20_lowmatch"
saveRDS(syn.agg, paste0(globals$tc.dir, sfname, "_wfs", ".rds"))
write_csv(syn.agg, paste0(globals$synd, "syntheses/", sfname, "_wfs", ".csv"))
# saveRDS(syn.agg, paste0(globals$tc.dir, "synthpop5_all_wtfs_new", ".rds"))

# now save the stacked output as csv
glimpse(puf.agg)
glimpse(syn.agg)
setdiff(names(puf.agg), names(syn.agg))
setdiff(names(syn.agg), names(puf.agg))
stack <- bind_rows(puf.agg, syn.agg)
count(stack, ftype)
summary(stack %>% dplyr::select(starts_with("wt")))
saveRDS(stack, paste0(globals$tc.dir, sfname, "_wfs_stack", ".rds"))
write_csv(stack, paste0(globals$synd, sfname, "_wfs_stack", ".csv"))
# synpuf17_wfs_stack.csv

#******************************************************************************************************************
# Explore results vs syn and vs puf counterpart ####
#******************************************************************************************************************
n <- 214
altpiece <- function(group.ind){
  optim <- readRDS(paste0("D:/tcdir/weight_pieces/synthpop8_500iter/optim_group_", group.ind, ".rds"))
}
altlist <- llply(1:n, altpiece, .progress="text")
puf.agg <- ldply(1:n, function(i) altlist[[i]]$puf)
syn.agg <- readRDS(paste0(globals$tc.dir, "synthpop8", "_wfs", ".rds"))
names(syn.agg)

nsamples <- max(syn.agg$m)
nsamples <- 5
stack <- bind_rows(puf.agg,
                   syn.agg %>% mutate(wt=wt.syn, ftype="syn"),
                   syn.agg %>% mutate(wt=wt.wtfs, ftype="wtfs")) %>%
  mutate(wt=ifelse(ftype=="syn", wt / 5, wt)) # needed for synthpop8

count(stack, ftype)
stack %>%
  group_by(ftype) %>%
  # mutate(wt=ifelse(ftype=="syn", w 5, wt)) %>%
  do(qtiledf(.$wt, probs=c(0, .01, .05, .1, .25, .5, .75, .8, .85, .9, .95, .99, 1)))

p <- stack %>%
  mutate(wt=ifelse(ftype=="puf", wt / 3, wt)) %>%
  ggplot(aes(x=wt, y = ..density..)) +
  geom_histogram(binwidth=25, fill="blue") +
  geom_vline(aes(xintercept = median(wt))) +
  scale_x_continuous(breaks=seq(0, 5000, 250), limits=c(0, 1500)) +
  theme(axis.text.x=element_text(size=8, angle=30)) +
  facet_wrap(~ftype, nrow=3) +
  ggtitle("Distribution of weights")
p

# summaries by income range
agiranges <- c(-Inf, 0, 25e3, 50e3, 75e3, 100e3, 200e3, 500e3, 1e6, 10e6, Inf)
vlist <- c("c00100", "e00200", "e00300", "e00600", "e01700", "p23250", "taxbc")
dfsums <- stack %>%
  mutate(agirange=cut(c00100, agiranges, right=FALSE),
         wtone=1e9) %>%
  dplyr::select(ftype, agirange, wt, wtone, vlist) %>%
  gather(vname, value, -ftype, -agirange, -wt) %>%
  group_by(ftype, agirange, vname) %>%
  summarise(n=n(), wtsum.m=sum(wt) / 1e6, valsum.b=sum(wt * value) / 1e9) %>%
  left_join(puf.vnames %>% dplyr::select(vname, vdesc))
dfsums

f <- function(vname.in, stat.in){
  dfsums %>%
    filter(vname==vname.in) %>%
    dplyr::select(ftype, agirange, stat=stat.in, vname, vdesc) %>%
    spread(ftype, stat) %>%
    janitor::adorn_totals(where="row") %>%
    mutate_at(vars(syn, wtfs), funs(diff=. - puf, pdiff=(. - puf) / puf * 100)) %>%
    dplyr::select(-c(vname, vdesc), everything(), vname, vdesc) %>%
    kable(digits=c(0, rep(1, 5), rep(1, 2)), format.args=list(big.mark = ','))
}

f("c00100", stat="valsum.b")
f("e00200", stat="valsum.b")
f("taxbc", stat="valsum.b")
f("e01700", stat="valsum.b") # way off

dfsums.m <- stack %>%
  mutate(agirange=cut(c00100, agiranges, right=FALSE),
         wtone=1e9) %>%
  dplyr::select(ftype, mgroup, agirange, wt, wtone, vlist) %>%
  gather(vname, value, -ftype, -mgroup, -agirange, -wt) %>%
  group_by(ftype, mgroup, agirange, vname) %>%
  summarise(n=n(), wtsum.m=sum(wt) / 1e6, valsum.b=sum(wt * value) / 1e9) %>%
  left_join(puf.vnames %>% dplyr::select(vname, vdesc))
dfsums.m

fm <- function(vname.in, stat.in, mgroup.in=1:3){
  dfsums.m %>%
    filter(vname==vname.in, mgroup %in% mgroup.in) %>%
    dplyr::select(ftype, mgroup, agirange, stat=stat.in, vname, vdesc) %>%
    spread(ftype, stat) %>%
    janitor::adorn_totals(where="row") %>%
    mutate_at(vars(syn, wtfs), funs(diff=. - puf, pdiff=(. - puf) / puf * 100)) %>%
    mutate(stat=stat.in) %>%
    dplyr::select(-c(vname, vdesc), everything(), stat, vname, vdesc) %>%
    kable(digits=c(0, 0, rep(2, 5), rep(1, 2)), format.args=list(big.mark = ','))
}

fm("c00100", stat="wtsum.m")
fm("c00100", stat="valsum.b")
fm("e00200", stat="valsum.b")
fm("taxbc", stat="valsum.b")
fm("e01700", stat="valsum.b")

fm("c00100", stat="wtsum.m", 1)
fm("c00100", stat="valsum.b", 3)
fm("e00200", stat="valsum.b")
fm("taxbc", stat="valsum.b", 2)
fm("e01700", stat="valsum.b")

# optlist[[1]]$inputs$recipe
# 53 target elements x 147 file groups ~7.8k targets
# esp. bad groups: mars 1 >$5.4m; mars 1 $4-6k, mars 1 $2-4k, mars 1, $36-38k
# but not many people or much $ in these groups

# quick check on stack ----
sfname <- "synthpop10"
df <- readRDS(paste0(globals$tc.dir, sfname, "_wfs_stack", ".rds"))
names(df) %>% sort

f <- function(x, wt) sum(x * wt) / 1e9
df %>%
  group_by(ftype, msname) %>%
  summarise(nret=sum(wt) / 1e6, 
            c00100=f(c00100, wt), 
            e00600=f(e00600, wt),
            p22250=f(p22250, wt),
            p23250=f(p23250, wt),
            taxbc=f(taxbc, wt)) %>%
  pivot_longer(cols=-c(ftype, msname), names_to = "variable", values_to = "value") %>%
  pivot_wider(names_from = ftype, values_from = value) %>%
  mutate(diff=syn - puf,
         pdiff=diff / puf * 100) %>%
  kable(digits=2)



#******************************************************************************************************************
# DEFUNCT: Temporary clunky approach to getting all 3 weights for the synfile ####
#******************************************************************************************************************

# tmp <- readRDS(paste0(globals$tc.dir, sfname, "_reweighted_stackedfiles.rds"))
# weights <- tibble(rownum=tmp$RECID[tmp$ftype=="puf.full"],
#                   puf.RECID=tmp$RECID[tmp$ftype=="puf.full"],
#                   puf.wt=tmp$wt[tmp$ftype=="puf.full"],
#                   syn.RECID=tmp$RECID[tmp$ftype=="synthpop3"],
#                   syn.wt=tmp$wt[tmp$ftype=="synthpop3"],
#                   syn.rwt=tmp$wt[tmp$ftype=="synthpop3.rwt"])
# ht(weights)






#******************************************************************************************************************
# DEFUNCT -- look at other nonlinear solvers ####
#******************************************************************************************************************

#******************************************************************************************************************
#  DEoptim ####
#******************************************************************************************************************
library("RcppDE")
Rosenbrock <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
## DEoptim searches for minima of the objective function between
## lower and upper bounds on each parameter to be optimized. Therefore
## in the call to DEoptim we specify vectors that comprise the
## lower and upper bounds; these vectors are the same length as the
## parameter vector.
lower <- c(-10,-10)
upper <- -lower
## run DEoptim and set a seed first for replicability
set.seed(1234)
DEoptim(Rosenbrock, lower, upper)
## increase the population size
DEoptim(Rosenbrock, lower, upper, DEoptim.control(NP = 100))

a <- proc.time()
tmp <- DEoptim(eval_f_wtfs, xlb, xub, control = DEoptim.control(trace = FALSE, NP = 1000), inputs=inputs)
b <- proc.time()
b - a

names(tmp)
names(tmp$optim)
tmp$optim$bestval
tmp$optim$nfeval
tmp$optim$iter

w.sol <- tmp$optim$bestmem
quantile(w.sol)


#******************************************************************************************************************
# trustOptim ####
#******************************************************************************************************************
# NO GOOD - cannot set bounds, gives negative weights
install.packages("trustOptim")
library("trustOptim")

val <- trust.optim(x0, fn=eval_f_full_scaled, gr=eval_grad_f_full_scaled, hs=NULL,
                   method = "SR1", control = list(report.precision=1L, function.scale.factor=-1),
                   inputs=inputs)


#******************************************************************************************************************
# optimx ####
#******************************************************************************************************************
# install.packages("optimx")
# install.packages("numDeriv")
# numDeriv
library("optimx")
# c("Nelder-Mead","BFGS")
# methods that allow box constraints
# Rcgmin bobyqa L-BFGS-B Rvmmmin maybe spg

grad.nd <- function(x, inputs) {
  require(numDeriv)
  grad.nd <- grad(eval_f_full_scaled, x, inputs=inputs)
  return(grad.nd)
}

opx <- optimx(x0, fn=eval_f_full_scaled, gr=grad.nd, hess=NULL,
              lower=xlb, upper=xub,
              method="bobyqa", itnmax=100, hessian=FALSE,
              control=list(trace=3),
              inputs=inputs)

opx <- optimx(x0, fn=eval_f_full_scaled, gr=eval_grad_f_full_scaled, hess=NULL,
              lower=xlb, upper=xub,
              method="bobyqa", itnmax=100, hessian=FALSE,
              control=list(trace=3),
              inputs=inputs)





