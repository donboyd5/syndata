


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
sfname <- "synthpop3"
synprep <- readRDS(paste0(globals$tc.dir, sfname, "_rwprep.rds"))
# synprep <- readRDS(paste0(globals$synd, sfname, "_rwprep.rds")) # this file is now in synpuf in Google Drive
names(synprep)

# now get the reforms and merge in taxbc

# merge and then split
tcvars <- c("c00100", "taxbc") # taxcalc vars
mrgdf <- left_join(synprep$tc.base, synprep$tc.output %>% dplyr::select(RECID, tcvars)) %>%
  # backwards compatibility for synthpop3
  mutate(ftype=case_when(ftype=="puf.full" ~ "puf",
                         ftype=="synthpop3" ~ "syn"))
glimpse(mrgdf)
count(mrgdf, ftype)


#******************************************************************************************************************
#  pick a subset ####
#******************************************************************************************************************
mrgdf2 <- mrgdf %>%
  filter(MARS==2, c00100>=0, c00100<=25e3)

puffile <- mrgdf2 %>% filter(ftype=="puf") # this has original puf weight
synfile <- mrgdf2 %>% filter(ftype=="syn") 

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
recipe <- get_recipe_long(get_weighting_recipe("recipe3")) %>%
  filter(vname %in% names(puffile)) %>%
  select(vname, vname, fn)


# start here to adjust a previously created recipe ----
tscale <- 1
recipe <- recipe %>%
  rowwise() %>%
  mutate(target=do.call(fn, list(puffile, vname, puffile$wt))) %>%
  ungroup %>%
  mutate(scale=ifelse(target!=0, abs(target / tscale), 1/ tscale),
         obj.element=paste0(vname, "_", fn)) %>%
  select(obj.element, vname, fn, scale, target) %>%
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
  mutate_at(vars(starts_with("flag")), funs(naz)) %>%
  ungroup %>%
  arrange(vname, fn)
recipe.flagged

# remove recipe elements where the target is zero
recipe.use <- recipe.flagged %>%
  filter(!(flag.dropneg | flag.dropdupsum | flag.dropdupn)) %>%
  filter(target!=0) %>%
  select(obj.element, vname, fn, scale, target)
recipe.use

# finally, add priority weights
recipe.use <- recipe.use %>%
  mutate(priority.weight=case_when(vname %in% c("wt", "c00100", "e00200", "taxbc") ~ 100,
                                   fn %in% c("n.sum", "val.sum") ~ 10,
                                   TRUE ~ 1))  %>% 
  left_join(puf.vnames %>% select(vname, vdesc))
recipe.use

recipe.use %>% arrange(-priority.weight)

# What would our objective function be if each targeted variable was off by a given % (as a decimal)?
pct <- .01
sum(recipe.use$priority.weight * (pct^2))
# what if they were off by that same pct on average but with a random variation?
pctv <- rnorm(nrow(recipe.use), pct, sd=.05)
sum(recipe.use$priority.weight * (pctv^2))
(pctv *100) %>% round(., 1)


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
  left_join(inputs$recipe %>% select(obj.element, vname, fn, scale, priority.weight, target)) %>%
  left_join(synlong) %>%
  mutate(coeff=case_when(fn=="val.sum" ~ value,
                         fn=="val.pos" ~ value*(value>0),
                         fn=="val.neg" ~ value*(value<0),
                         fn=="n.sum" ~ 1,
                         fn=="n.pos" ~ 1*(value>0),
                         fn=="n.neg" ~ 1*(value<0),
                         TRUE  ~ 0)) %>%
  select(obj.element, vname, fn, wtnum, scale, priority.weight, value, coeff, target)
# glimpse(coeffs)
# ht(coeffs)
inputs$coeffs <- coeffs


#******************************************************************************************************************
#  run ipoptr ####
#******************************************************************************************************************
# inputs$recipe

# bounds on the weights
xlb <- rep(1, nrow(synfile))
xub <- rep(1.5*max(puffile$wt), nrow(synfile))
#xub <- rep(max(puffile$wt), nrow(synfile))

# starting point:
# x0 <- (xlb + xub) / 2
x0 <- rep(mean(puffile$wt), nrow(synfile))

# alternatively:
set.seed(1234)
x0 <- rnorm(nrow(synfile), mean(puffile$wt), sd(puffile$wt)); x0 <- pmax(x0, xlb); x0 <- pmin(x0, xub); quantile(x0)


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


# control = list(xtol_rel = 1e-8)
nloptr.print.options()

a <- proc.time()
t2 <- mma(x0, fn=eval_f_wtfs, gr = eval_grad_f_wtfs,
          lower=xlb, upper=xub,
          nl.info = FALSE, inputs=inputs)
b <- proc.time()
b - a

names(t2)
# t2$par
t2$value
t2$iter
t2$convergence
t2$message
w.sol <- t2$par


#******************************************************************************************************************
#  Examine results ####
#******************************************************************************************************************
# retrieve a previous run or else use the results from above


# ------------------
w.sol <- result$solution
# w.sol <- val$solution

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
# Solve the problem for the full file, in pieces ####
#******************************************************************************************************************
sfname <- "synthpop5"

synprep <- readRDS(paste0(globals$tc.dir, sfname, "_all_rwprep.rds"))
# synprep <- readRDS(paste0(globals$synd, sfname, "_all_rwprep.rds")) # this file is now in synpuf in Google Drive
names(synprep)

# now get the reforms and merge in taxbc, then merge and then split
tcvars <- c("c00100", "taxbc") # taxcalc vars
mrgdf <- left_join(synprep$tc.base, synprep$tc.output %>% dplyr::select(RECID, tcvars))
glimpse(mrgdf)
count(mrgdf, ftype)


# get full versions of the two files
puf <- mrgdf %>% filter(ftype=="puf")
syn <- mrgdf %>% filter(ftype=="syn") %>% mutate(wt=wt / 3)
sum(puf$wt)
sum(syn$wt)
# note that RECID is sequential from 1 on puf to the highest value on syn and is unique

#.. examine the data to estimate where we should split it ----
# determine different numbers of cuts for each marital group
tmp <- puf %>% 
  select(MARS, wt, c00100) %>%
  mutate(mgroup=ifelse(MARS %in% 1:2, MARS, 3),
         inum=case_when(mgroup==1 ~ 12,
                        mgroup==2 ~ 20,
                        mgroup==3 ~ 10,
                        TRUE ~ 0)) %>%
  group_by(mgroup) %>%
  mutate(igroup=ntile(c00100, inum))

gbreaks <- tmp %>%
  group_by(mgroup, igroup) %>%
  summarise(imin=min(c00100), imax=max(c00100), n=n(), wt.sum=sum(wt), agib=sum(wt * c00100) / 1e9)
gbreaks

# how does synfile look if we choose these exact breaks?
i1 <- gbreaks %>% filter(mgroup==1) %>% .[["imin"]]
i1 <- c(-Inf, i1[-1], Inf)
i1

i2 <- gbreaks %>% filter(mgroup==2) %>% .[["imin"]]
i2 <- c(-Inf, i2[-1], Inf)
i2

i3 <- gbreaks %>% filter(mgroup==3) %>% .[["imin"]]
i3 <- c(-Inf, i3[-1], Inf)
i3

# adjust the breaks so that there is always a zero
adjb <- function(b) c(b[1], 0, b[2:length(b)])
adjb(i3)

i1a <- adjb(i1)
i2a <- adjb(i2)
i3a <- adjb(i3)

f <- function(c00100, mgroup){
  brks <- if(mgroup[1]==1) i1a else
    if(mgroup[1]==2) i2a else
      if(mgroup[1]==3) i3a
  igroup.element <- function(c00100) min(which(c00100 < brks)) - 1
  igroup <- sapply(c00100, igroup.element)
  return(igroup)
}

# look at two files together
btmp <- bind_rows(puf %>% mutate(type="puf"),
                  syn %>% mutate(type="syn")) %>%
  select(type, MARS, wt, c00100) %>%
  mutate(mgroup=ifelse(MARS %in% 1:2, MARS, 3)) %>%
  group_by(type, mgroup) %>%
  mutate(igroup=f(c00100, mgroup))
btmp %>% ungroup %>% filter(mgroup==1) %>% count(igroup)

btmp %>%
  group_by(type, mgroup, igroup) %>%
  summarise(n=n()) %>%
  spread(type, n) %>%
  mutate(diff=syn - puf) %>%
  ungroup %>%
  arrange(mgroup, igroup)

#.. create split rules ----
m1 <- tibble(imin=i1a[-length(i1a)], imax=i1a[-1], mgroup=1)
m2 <- tibble(imin=i2a[-length(i2a)], imax=i2a[-1], mgroup=2)
m3 <- tibble(imin=i3a[-length(i3a)], imax=i3a[-1], mgroup=3)
split.rules <- bind_rows(m1, m2, m3) %>% mutate(group=row_number()) %>% select(group, mgroup, imin, imax)
split.rules

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
  select(ftype, RECID, mgroup, group) %>%
  arrange(RECID)
ht(idfile)
count(idfile, mgroup, ftype, group) %>% spread(ftype, n) %>% mutate(diff=syn - puf, sum=puf + syn)

# now we are ready to run the file in pieces

recipe <- get_recipe_long(get_weighting_recipe("recipe3_sub")) %>%
  filter(vname %in% names(puf)) %>%
  dplyr::select(vname, vname, fn)


getrec <- function(puf, syn, recipe, puf.vnames){
  # adjust a previously created recipe
  tscale <- 1
  recipe <- recipe %>%
    rowwise() %>%
    mutate(target=do.call(fn, list(puf, vname, puf$wt))) %>%
    ungroup %>%
    mutate(scale=ifelse(target!=0, abs(target / tscale), 1/ tscale),
           obj.element=paste0(vname, "_", fn)) %>%
    dplyr::select(obj.element, vname, fn, scale, target) %>%
    arrange(vname, fn)
  
  #..weed out unnecessary elements of the recipe ----
  # if the target is 0 for negative values we can drop the neg versions AS LONG AS WE HAVE NO SYNTH NEG VERSIONS
  # if the val.pos and val.sum versions are identical then we can drop the val.sum version
  # can drop the "neg" and "pos" versions
  recipe.flagged <- recipe %>%
    rowwise() %>%
    mutate(syn.unwtd=do.call(fn, list(syn, vname, rep(1, nrow(syn))))) %>% # so we can check if negs!
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
                                     fn %in% c("n.sum", "val.sum") ~ 10,
                                     TRUE ~ 1))  %>% 
    left_join(puf.vnames %>% dplyr::select(vname, vdesc))
  return(list(recipe.use=recipe.use, recipe.flagged=recipe.flagged))
}


getinplist <- function(syn, recipe.use){
  inputs <- list()
  inputs$recipe <- recipe.use
  inputs$synsub <- syn[, unique(inputs$recipe$vname)] %>% mutate(wt=1)
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

  inputs$coeffs <- coeffs
  return(inputs)
}


rungroup <- function(group.ind){
  
  n.sum <- function(df, var, weight, condition=TRUE){
    # get the weighted number of records for which a logical condition related to the variable is met
    # var is a numeric column in df
    # weight is a vector of weights (could be a column in df, or could be external)
    # condition: boolean expression as text
    # returns a scalar which is the weighted number of records
    condition <- parse(text=condition)
    df %>%
      select(variable=!!var) %>%
      summarise(value=sum(weight * eval(condition))) %>%
      .[[1]]
  }
  
  
  val.sum <- function(df, var, weight, condition=TRUE){
    # get the weighted value of a variable for which a logical condition related to the variable is met
    # var is a numeric column in df
    # weight is a vector of weights (could be a column in df, or could be external)
    # condition: boolean expression as text
    # returns a scalar which is the weighted value of the variable
    condition <- parse(text=condition)
    df %>% 
      select(variable=!!var) %>%
      summarise(value=sum(variable * weight * eval(condition))) %>%
      .[[1]]
  }
  
  getrec <- function(puf, syn, recipe, puf.vnames){
    # adjust a previously created recipe
    tscale <- 1
    print(recipe)
    # recipe <- recipe %>%
    #   rowwise() %>%
    #   mutate(target=n.sum(puf, vname, puf$wt))
    recipe$target <- NA_real_
    for(i in 1:nrow(recipe)){
      recipe$target[i] <- get(recipe$fn[i])(df=puf, var=recipe$vname[i], puf$wt)
    }
      #mutate(target=do.call(fn, list(puf, vname, puf$wt))) # %>%
    #   ungroup %>%
    #   mutate(scale=ifelse(target!=0, abs(target / tscale), 1/ tscale),
    #          obj.element=paste0(vname, "_", fn)) %>%
    #   dplyr::select(obj.element, vname, fn, scale, target) %>%
    #   arrange(vname, fn)
    # 
    # #..weed out unnecessary elements of the recipe ----
    # # if the target is 0 for negative values we can drop the neg versions AS LONG AS WE HAVE NO SYNTH NEG VERSIONS
    # # if the val.pos and val.sum versions are identical then we can drop the val.sum version
    # # can drop the "neg" and "pos" versions
    # recipe.flagged <- recipe %>%
    #   rowwise() %>%
    #   mutate(syn.unwtd=do.call(fn, list(syn, vname, rep(1, nrow(syn))))) %>% # so we can check if negs!
    #   group_by(vname) %>%
    #   mutate(flag.dropneg=ifelse(str_detect(fn, "neg") & target==0 & syn.unwtd==0, 1, 0),
    #          flag.dropdupsum=ifelse(target==target[match("val.sum", fn)] & (fn=="val.pos"), 1, 0),
    #          flag.dropdupn=ifelse(target==target[match("n.sum", fn)] & (fn=="n.pos"), 1, 0)) %>%
    #   mutate_at(vars(starts_with("flag")), funs(naz)) %>%
    #   ungroup %>%
    #   arrange(vname, fn)
    # 
    # # remove recipe elements where the target is zero
    # recipe.use <- recipe.flagged %>%
    #   filter(!(flag.dropneg | flag.dropdupsum | flag.dropdupn)) %>%
    #   filter(target!=0) %>%
    #   dplyr::select(obj.element, vname, fn, scale, target)
    # 
    # # finally, add priority weights
    # recipe.use <- recipe.use %>%
    #   mutate(priority.weight=case_when(vname %in% c("wt", "c00100", "e00200", "taxbc") ~ 100,
    #                                    fn %in% c("n.sum", "val.sum") ~ 10,
    #                                    TRUE ~ 1))  %>% 
    #   left_join(puf.vnames %>% dplyr::select(vname, vdesc))
    return(list(recipe.use=recipe, recipe.flagged=recipe))
    #return(list(recipe.use=recipe.use, recipe.flagged=recipe.flagged))
  }
  
  getinplist <- function(syn, recipe.use){
    inputs <- list()
    inputs$recipe <- recipe.use
    inputs$synsub <- syn[, unique(inputs$recipe$vname)] %>% mutate(wt=1)
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
    
    inputs$coeffs <- coeffs
    return(inputs)
  }
  
  # a <- proc.time()
  base <- left_join(idfile %>% filter(group==group.ind), mrgdf)
  puf <- base %>% filter(ftype=="puf")
  syn <- base %>% filter(ftype=="syn")
  recipes <- getrec(puf, syn, recipe, puf.vnames)
  # recipe.use <- recipes$recipe.use
  # inputs <- getinplist(syn, recipe.use)
  # 
  # # bounds on the weights
  # xlb <- rep(1, nrow(syn))
  # xub <- rep(1.5*max(puf$wt), nrow(syn))
  # 
  # # starting point:
  # x0 <- (xlb + xub) / 2
  # x0 <- x0 * sum(puf$wt / sum(x0))
  # 
  # result <- mma(x0, fn=eval_f_wtfs, gr=eval_grad_f_wtfs,
  #               lower=xlb, upper=xub,
  #               nl.info = FALSE, inputs=inputs)
  # 
  # optim <- list()
  # optim$result <- result
  # optim$puf <- puf
  # optim$syn <- syn
  # optim$inputs <- inputs
  # optim$recipe.flagged <- recipes$recipe.flagged
  tmp <- tibble(a=group.ind)
  saveRDS(tmp, paste0(globals$tc.dir, "weight_pieces/optim_group_", group.ind, ".rds"))
  #write_csv(tmp, globals$tc.dir, "weight_pieces/optim_group_", group.ind, ".csv")
  # saveRDS(ab.dist, paste0(globals$tc.dir, "dist_pieces/", "ab.dist_m", mgroup.in, "_ntile", ntile.in, ".rds"))
  # b <- proc.time()
  #print(b - a)
  return(recipes)
}


# TODO: make these computations parallel

library("doParallel")
cl <- makeCluster(2)
registerDoParallel(cl)

packages <- c("magrittr", "tidyverse", "dplyr", "nloptr")
# CAUTION:  DO NOT PASS large items as function arguments
# instead export them
xport <- c("globals", "idfile", "recipe", "puf.vnames", "mrgdf", "n.sum", "val.sum")
xport <- c("globals", "idfile", "recipe", "puf.vnames", "mrgdf") 
popts <- list(.packages=packages, .export=xport)

# 1:max(idfile$group)

a <- proc.time()
# ldply(1:2, .fun=rungroup, idfile=idfile, recipe=recipe, puf.vnames=puf.vnames, mrgdf=mrgdf,
#       .progress="text", .parallel=TRUE, .paropts=popts)
llply(1:1, .progress="text", .parallel=TRUE, .paropts=popts, .fun=rungroup)
b <- proc.time()
b - a # 

# stopCluster(cl)
# stopImplicitCluster()

match.fun(funcList[[1]])
match.fun(list("n.sum", "val.sum")[[1]])
df <- tibble(a=c(1, .5, 2, 3, 4), b=6:10, fn=c("max", "mean", "sum", "min", "mean"))
df %>%
  rowwise() %>%
  mutate(c=get(fn)(c(a, b)))

df$c <- get(df$fn)(c(df$a, df$b))
df
for(i in 1:nrow(df)){
  df$c[i] <- get(df$fn[i])(c(df$a[i], df$b[i]))
}
df

lapply(df$fn, get, x=c(df$a, df$b))

a <- proc.time()
ldply(1:2, .fun=rungroup, idfile=idfile, recipe=recipe, puf.vnames=puf.vnames, mrgdf=mrgdf,
      .progress="text", .parallel=TRUE, .paropts=popts)
b <- proc.time()
b - a # 

rungroup <- function(group.ind){
  # tmp <- tibble(a=group.ind)
  tmp <- idfile[1, 1:2]
  return(tmp)
}

packages <- c("magrittr", "tidyverse", "dplyr", "nloptr")
# xport <- c("globals", "idfile", "mrgdf", "getrec", "getinplist")
# xport <- c("globals", "getrec", "getinplist")
xport <- c("globals", "idfile", "recipe", "puf.vnames", "mrgdf")
xport <- c("globals")
popts <- list(.packages=packages, .export=xport)

a <- proc.time()
ldply(1:2, .progress="text", .parallel=TRUE, .paropts=popts, .fun=rungroup)
# ldply(1:2, .fun=rungroup, idfile=idfile, recipe=recipe, puf.vnames=puf.vnames, mrgdf=mrgdf, .progress="text", .parallel=FALSE, .paropts=popts)
b <- proc.time()
b - a # 


#******************************************************************************************************************
# Construct and save full file ####
#******************************************************************************************************************


getpiece <- function(group.ind){
  optim <- readRDS(paste0(globals$tc.dir, "weight_pieces/optim_group_", group.ind, ".rds"))
}

n <- 45
optlist <- llply(1:n, getpiece, .progress="text")
memory()

# group.ind <- 2

length(optlist)
names(optlist[[1]])
names(optlist[[1]]$inputs)
optlist[[1]]$inputs$recipe
optlist[[1]]$recipe.flagged
unique(optlist[[1]]$recipe.flagged$vname)
names(optlist[[1]]$result)
names(optlist[[1]]$syn)

# analyze summary result
obj.vals <- laply(1:n, function(i) optlist[[i]]$result$par)
quantile(obj.vals) %>% round(2)

message <- laply(1:n, function(i) optlist[[i]]$result$message)
count(tibble(message), message)

table(laply(1:n, function(i) optlist[[i]]$result$value))
table(laply(1:n, function(i) optlist[[i]]$result$iter))
table(laply(1:n, function(i) optlist[[i]]$result$convergence))


# aggregate the pieces of the puf and synthetic files, and attach new weights
# first puf
puf.agg <- ldply(1:n, function(i) optlist[[i]]$puf)
names(puf.agg)
ht(puf.agg[, c(1:5, (ncol(puf.agg)-5):ncol(puf.agg))]) # NOTE that RECIDs are not in order
min(puf.agg$RECID); max(puf.agg$RECID)
puf.agg <- puf.agg %>% arrange(RECID)


# now synthetic
syn.agg <- ldply(1:n, function(i) {optlist[[i]]$syn %>% mutate(wt.wtfs=optlist[[i]]$result$par)})
names(syn.agg) %>% sort
ht(syn.agg[, c(1:5, ncol(syn.agg))]) # RECIDs not in order here, either
syn.agg <- syn.agg %>%
  mutate(wt.syn=wt, wt=wt.wtfs, ftype="syn.wtfs")

# save the wtfs file; also save puf counterpart
saveRDS(syn.agg, paste0(globals$tc.dir, "synthpop5_all_wtfs", ".rds"))



#******************************************************************************************************************
# Explore results vs syn and vs puf counterpart ####
#******************************************************************************************************************
stack <- bind_rows(puf.agg,
                   syn.agg %>% mutate(wt=wt.syn / 3, ftype="syn"),
                   syn.agg %>% mutate(wt=wt.wtfs, ftype="wtfs"))
count(stack, ftype)

p <- stack %>%
  ggplot(aes(x=wt, y = ..density..)) +
  geom_histogram(binwidth=25, fill="blue") +
  geom_vline(aes(xintercept = median(wt))) +
  scale_x_continuous(breaks=seq(0, 5000, 250), limits=c(0, 3000)) +
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
f("e01700", stat="valsum.b")


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





