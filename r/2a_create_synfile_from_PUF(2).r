


#****************************************************************************************************
#                Libraries ####
#****************************************************************************************************
library("magrittr")
library("plyr") # needed for ldply; must be loaded BEFORE dplyr
library("tidyverse")
options(tibble.print_max = 60, tibble.print_min = 60) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats

library("readxl") # readxl, for .xls and .xlsx files.

library("scales")

library("btools") # should not need this. If you do, it's at https://github.com/donboyd5/btools

library("knitr")

library("janitor") # for adorn_totals

library("synthpop") # note: masks select in dplyr


#****************************************************************************************************
#                Includes ####
#****************************************************************************************************
source("./r/includes/globals_system_specific_boyd.r") # use a different version of this file if changing systems
source("./r/includes/globals_other.r")

source("./r/includes/functions_general.r")
source("./r/includes/functions_synthesis.r")


#****************************************************************************************************
#                Globals ####
#****************************************************************************************************
sub1 <- expression(MARS==1)
sub2 <- expression(MARS==2)
sub34 <- expression(MARS %in% c(3, 4))


#****************************************************************************************************
#                Get puf.names, psum, and puf ####
#****************************************************************************************************
puf.vnames <- get_puf_vnames()
# puf.names

# really should do this on the full file to fine out what the big vars are. Var lists below reflect full file.
psum <- get_pufbase_sums()
psum

puf <- get_puf.base()

names(puf) %>% sort
count(puf, MARS)


#****************************************************************************************************
#                CAUTION: DO THIS OR... - Create puf.base, to be synthesized, from a file from Max, or... ####
#****************************************************************************************************

puf.from.max <- read_csv(paste0(globals$synd, "puf_samples/", "train50.csv"),
                         col_types = cols(.default= col_double()),
                         n_max=-1)

glimpse(puf.from.max)
names(puf.from.max) %>% sort


# CAUTION: get needed variables ----
# "e04600" "p04470" "e04800" "e62100" "e05800" "e08800" "e59560" "e26190"
needvars <- c("e04600", "p04470", "e04800", "e62100", "e05800", "e08800", "e59560", "e26190")

pfm2 <- puf.from.max %>% 
  setNames(change_case(names(.))) %>%
  left_join(puf %>% dplyr::select(RECID, needvars)) %>%
  # create variables not in Max's file but needed so that I don't have to modify synthesis routines
  mutate(ftype="puf",
         e00600=e00600_minus_e00650 + e00650,
         e01500=e01500_minus_e01700 + e01700,
         pufseqn=row_number(),
         wt=s006 / 100)

pfm2 <- puf.from.max %>%
  mutate(ftype="puf",
         pufseqn=row_number(),
         wt=s006 / 100)

glimpse(pfm2)
summary(pfm2)
names(pfm2) %>% sort
cor(pfm2 %>% dplyr::select(p08000, p23250, e09700))


# prepare file for synthesis
puf.base <- pfm2 %>%
  # filter(!RECID %in% 999996:999999) %>%
  mutate(e00600_minus_e00650=e00600 - e00650,
         e01500_minus_e01700=e01500 - e01700,
         divratio=e00650 / e00600,
         divratio=ifelse(is.na(divratio), 0, divratio),
         penratio=e01700 / e01500,
         penratio=ifelse(is.na(penratio), 0, penratio))
glimpse(puf.base)

#****************************************************************************************************
#                CAUTION: OR DO THIS... Create puf.base, to be synthesized, from puf ####
#****************************************************************************************************

# prepare file for synthesis
puf.base <- puf %>%
  filter(!RECID %in% 999996:999999) %>%
  mutate(e00600_minus_e00650=e00600 - e00650,
         e01500_minus_e01700=e01500 - e01700,
         divratio=e00650 / e00600,
         divratio=ifelse(is.na(divratio), 0, divratio),
         penratio=e01700 / e01500,
         penratio=ifelse(is.na(penratio), 0, penratio))


#****************************************************************************************************
#                Explore relationships ####
#****************************************************************************************************
df <- puf.base %>%
  filter(eval(sub1))
glimpse(df)

corvars <- names(df)[(names(df) %>% str_sub(., 1, 1) %in% c("e", "p")) | names(df)=="wt"]
corvars <- setdiff(corvars, c("efi", "elect", "prep", "penratio"))
corvars

corpuf <- df %>%
  dplyr::select(corvars) %>%
  do(cordf(.)) %>%
  do(trimcor(.)) %>%
  separate(combo, c("vname1", "vname2"), remove=TRUE) %>%
  ungroup %>%
  left_join(puf.vnames %>% dplyr::select(vname1=vname, vdesc1=vdesc)) %>%
  left_join(puf.vnames %>% dplyr::select(vname2=vname, vdesc2=vdesc))
ht(corpuf)

txt <- "pension"
corpuf %>%
  filter(str_detect(vdesc1, coll(txt, ignore.case=TRUE)) |
           str_detect(vdesc2, coll(txt, ignore.case=TRUE))) %>%
  arrange(-abs(value)) %>%
  kable(digits=3)

#.... ----- START OF SYNTHESIS FOR REAL ....----

#****************************************************************************************************
#                Get synthesis recipe ####
#****************************************************************************************************
# sheetname <- "recipe5"
# sheetname <- "recipe6"
# sheetname <- "recipe7"
sheetname <- "simple"
synrecipe <- read_excel("./data/synthesis_recipes.xlsx", sheet = sheetname) %>%
  mutate(value=as.numeric(value))
# synrecipe


synrecipe.use <- synrecipe %>%
  dplyr::select(idvar, xvar.keep, xvar.drop, synvar, vname, groupname, value, vdesc) %>%
  filter(idvar | xvar.keep | xvar.drop | synvar) %>%
  mutate(synorder=row_number()) %>%
  arrange(-idvar, -xvar.keep, -xvar.drop, synorder)
synrecipe.use

# CAUTION: Drop variables that are not in Max's file!! ----
# "e04600" "p04470" "e04800" "e62100" "e05800" "e08800" "e59560" "e26190"
badvars <- setdiff(unique(synrecipe.use$vname), names(puf.base))
psum %>% filter(vname %in% badvars)
# END OF CAUTION: ----


synrecipe.use %>%
  dplyr::select(idvar, xvar.keep, xvar.drop, synvar) %>%
  gather(variable, value) %>%
  group_by(variable) %>%
  summarise(n=sum(value))


#****************************************************************************************************
#                Prepare synthesis variables ####
#****************************************************************************************************

idvars <- synrecipe.use %>% filter(idvar==1) %>% .[["vname"]]
xvars.keep <- synrecipe.use %>% filter(xvar.keep==1) %>% .[["vname"]]
xvars.drop <- synrecipe.use %>% filter(xvar.drop==1) %>% .[["vname"]]
synvars <- synrecipe.use %>% filter(synvar==1) %>% .[["vname"]]

allvars <- c(idvars, xvars.keep, xvars.drop, synvars)
allvars

# move variables as needed to create the visit sequence in desired order
synvars.vs <- synvars
synvars.vs
# synvars.vs <- move(synvars.vs, "e01500", "e01500_minus_e01700")
# synvars.vs <- move(synvars.vs, "e00600", "e00600_minus_e00650")
# synvars.vs <- move(synvars.vs, "e01000", "e00650")
synvars.vs

# the full visit sequence includes all variables
visit.sequence <- c(idvars, xvars.keep, xvars.drop, synvars.vs)
visit.sequence

# setdiff(names(puf.base), visit.sequence)
# setdiff(visit.sequence, names(puf.base))

#****************************************************************************************************
#                Set up synthesis methods ####
#****************************************************************************************************
# note variables that must be created before running syn: divratio, penratio
# note that taxcalc requires e00650 (qualified div) <= e00600 (total div)
# I create and use divratio = e00600 / e00650 to ensure this
# similar reason, I create penratio=E01700 / E01500

methods <- rep("cart", length(visit.sequence)) # set default method

# override default method for some variables

not.synthesized <- c(idvars, xvars.keep, xvars.drop)
methods[1:length(not.synthesized)] <- "" # these will not be estimated

# methods[visit.sequence=="e01700"] <- "~I(e01500 * penratio)"
methods[visit.sequence=="e01500"] <- "~I(e01700 + e01500_minus_e01700)"

methods[visit.sequence=="e00650"] <- "~I(e00600 * divratio)"
# methods[visit.sequence=="e00600"] <- "~I(e00650 + e00600_minus_e00650)"

# methods

cbind(visit.sequence, methods)

# write_csv(tibble(vars=visit.sequence, method=methods), "d:/tcdir/synthpop6_vars.csv")

#****************************************************************************************************
#                TEST: Synthesize a single file part, or a file all at once ####
#****************************************************************************************************
# create synthpop from parts

pufsyn <- puf.base %>%
  filter(eval(sub1)) %>%
  #sample_frac(size=.2) %>%
  dplyr::select(visit.sequence)
glimpse(pufsyn)

# set the predictor matrix so that we do not use id variables in predictions
# a column of the pm tells which vars the colvar is a predictor for
test <- syn(pufsyn, visit.sequence=visit.sequence, method=methods, seed=1234, m=0, k=nrow(pufsyn), proper=FALSE) # m: 1=32, 2=64, 3=92
test$predictor.matrix
colnames(test$predictor.matrix)
test$predictor.matrix[, "p22250"]
test$predictor.matrix[, "e02100"]
test$predictor.matrix[, "pufseqn"]
pred.mat.corrected <- test$predictor.matrix
pred.mat.corrected[, "pufseqn"] <- 0
pred.mat.corrected # use this revised predictor matrix in the synthesis
pred.mat.corrected[, "p22250"]

smoothing.rules <- list(e00200="density")
synx <- syn(pufsyn, visit.sequence=visit.sequence, method=methods, smoothing=smoothing.rules, 
            seed=1, m=1, k=nrow(pufsyn), proper=FALSE, predictor.matrix = pred.mat.corrected)
tmp <- synx$syn
glimpse(tmp)

# examine the distribution of wages ----
# unsmoothed <- synx
# smoothed <- synx
compwages <- tibble(rownum=1:nrow(pufsyn), 
                    wages=pufsyn$e00200,
                    unsmoothed=unsmoothed$syn$e00200,
                    smoothed=smoothed$syn$e00200)
ht(compwages)
cor(compwages %>% dplyr::select(-rownum)) # correlation between unsmoothed and smoothed is 1

f <- function(x) length(unique(x))
compwages %>%
  summarise_at(vars(wages, smoothed, unsmoothed), ~f(.)) # we now have 8x as many unique values

# compwages %>%
#   gather(vname, value, -rownum, -wages) %>%
#   filter(wages<10000) %>%
#   ggplot(aes(wages, value, colour=vname)) +
#   geom_point()

# end ----


seed <- 1234
# synx <- syn(pufsyn, visit.sequence=visit.sequence, method=methods, seed=seed, m=0)
a <- proc.time()
# synx <- syn(pufsyn, visit.sequence=visit.sequence, method=methods, seed=seed, m=1, k=4*nrow(pufsyn)) # k: 1=32, 2=45, 3=57
synx <- syn(pufsyn, visit.sequence=visit.sequence, method=methods, seed=seed, m=1, k=nrow(pufsyn)) # m: 1=32, 2=64, 3=92
b <- proc.time()
b - a

# names(synx)
# str(synx)
# ht(synx$syn)


#****************************************************************************************************
#                FULL RUN: Loop through subsets and synthesize each ####
#****************************************************************************************************
#.. Functions for full run ----
get_synrecipe <- function(recipe_name){
  synrecipe <- read_excel("./data/synthesis_recipes.xlsx", sheet = recipe_name) %>%
    mutate(value=as.numeric(value))

  synrecipe.use <- synrecipe %>%
    dplyr::select(idvar, xvar.keep, xvar.drop, synvar, vname, groupname, value, vdesc) %>%
    filter(idvar | xvar.keep | xvar.drop | synvar) %>%
    mutate(synorder=row_number()) %>%
    arrange(-idvar, -xvar.keep, -xvar.drop, synorder)
  
  return(synrecipe.use)
}

get_vargroups <- function(synrecipe.use){
  idvars <- synrecipe.use %>% filter(idvar==1) %>% .[["vname"]]
  xvars.keep <- synrecipe.use %>% filter(xvar.keep==1) %>% .[["vname"]]
  xvars.drop <- synrecipe.use %>% filter(xvar.drop==1) %>% .[["vname"]]
  synvars <- synrecipe.use %>% filter(synvar==1) %>% .[["vname"]]
  
  allvars <- c(idvars, xvars.keep, xvars.drop, synvars)
  
  # move variables as needed to create the visit sequence in desired order
  synvars.vs <- synvars
  # synvars.vs <- move(synvars.vs, "e01500", "e01500_minus_e01700")
  
  # the full visit sequence includes all variables
  visit.sequence <- c(idvars, xvars.keep, xvars.drop, synvars.vs)
  
  vargroups <- list()
  vargroups$idvars <- idvars
  vargroups$xvars.keep <- xvars.keep
  vargroups$xvars.drop <- xvars.drop
  vargroups$synvars.vs <- synvars.vs
  vargroups$visit.sequence <- visit.sequence
  return(vargroups)
}


get_methods <- function(vargroups){
  methods <- rep("cart", length(vargroups$visit.sequence)) # set default method
  
  # override default method for some variables
  not.synthesized <- c(vargroups$idvars, vargroups$xvars.keep, vargroups$xvars.drop)
  methods[1:length(not.synthesized)] <- "" # these will not be estimated
  
  # methods[visit.sequence=="e01700"] <- "~I(e01500 * penratio)"
  methods[vargroups$visit.sequence=="e01500"] <- "~I(e01700 + e01500_minus_e01700)"
  
  methods[vargroups$visit.sequence=="e00650"] <- "~I(e00600 * divratio)"
  # methods[visit.sequence=="e00600"] <- "~I(e00650 + e00600_minus_e00650)"

  # cbind(visit.sequence, methods)
  return(methods)
}


runsyn <- function(i, grouplist, groupnames){
  # get easier names
  idvars <- inputs$vargroups$idvars
  xvars.drop <- inputs$vargroups$xvars.drop
  visit.sequence <- inputs$vargroups$visit.sequence
  methods <- inputs$methods
  smoothing.rules <- inputs$smoothing.rules
  synname <- inputs$synname
  m <- inputs$m
  sampfrac <- inputs$sampfrac
  cart.minbucket <- 5
  if(!is.null(inputs$cart.minbucket)) cart.minbucket <- inputs$cart.minbucket
  
  correct_predictor_matrix <- function(visit.sequence){
    pufsyn <- puf.base %>%
      dplyr::select(visit.sequence)
    test <- syn(pufsyn, visit.sequence=visit.sequence, method=methods, seed=1234, m=0)
    pred.mat.corrected <- test$predictor.matrix
    pred.mat.corrected[, "pufseqn"] <- 0
    return(pred.mat.corrected)
  }
  
  print(paste0("Synthesizing group ", i, " of ", length(grouplist)))
  group <- grouplist[[i]]
  groupname <- groupnames[[i]]; print(groupname)
  outfname <- paste0(globals$tc.dir, "pufsyn_stack/", synname, "_", groupname, ".rds")
  
  pufsyn <- puf.base %>%
    ungroup %>% # just to be safe
    dplyr::filter(eval(group)) %>%
    sample_frac(sampfrac) %>%
    dplyr::select(visit.sequence)
  
  pred.mat.corrected <- correct_predictor_matrix(visit.sequence)
  
  seed <- 1234
  synx <- syn(pufsyn, visit.sequence=visit.sequence, method=methods, smoothing=smoothing.rules, cart.minbucket = cart.minbucket, seed=seed, m=m, k=nrow(pufsyn), predictor.matrix = pred.mat.corrected, print.flag = TRUE)

  #.. combine results - works even when m==1 ----
  synout <- bind_rows(synx$syn) %>% mutate(m=rep(1:synx$m, each=synx$n))
  stack <- bind_rows(pufsyn %>% mutate(ftype="puf", m=1),
                     synout %>% mutate(ftype="syn")) %>%
    dplyr::select(-xvars.drop) %>%
    dplyr::select(-m, everything(), m) %>%
    mutate(msname=groupname)
  
  # save the run details -- ASSUME that the data file will be ok when stack is returned
  synout <- list()
  # synx$syn <- NULL # save space and time
  synout$synx <- synx
  saveRDS(synout, outfname)
  
  return(stack)
}


#.. preparation needed for a parallel run ----
library("doParallel")
cl <- makeCluster(3)
registerDoParallel(cl)
# NOTE: when all done, run: stopCluster(cl)
showConnections()
#.. end parallel prep ----


#.. Set up for full run ----
puf.base <- puf %>%
  filter(!RECID %in% 999996:999999) %>%
  mutate(e00600_minus_e00650=e00600 - e00650,
         e01500_minus_e01700=e01500 - e01700,
         divratio=e00650 / e00600,
         divratio=ifelse(is.na(divratio), 0, divratio),
         penratio=e01700 / e01500,
         penratio=ifelse(is.na(penratio), 0, penratio))

#.. start the run ----
sub1 <- expression(MARS==1)
sub2 <- expression(MARS==2)
sub34 <- expression(MARS %in% c(3, 4))
mslist <- list(sub1, sub2, sub34)
msnames <- c("MARS1", "MARS2", "MARS34")


# define variables that control the size of the run
inputs <- list()
inputs$synname <- "synthpop10"
inputs$sampfrac <- 1
inputs$m <- 5
inputs$syngroups <- 2
inputs$cart.minbucket <- 10
inputs$synrecipe.use <- get_synrecipe("recipe7")
inputs$vargroups <- get_vargroups(inputs$synrecipe.use)
inputs$methods <- get_methods(inputs$vargroups)
cbind(inputs$vargroups$visit.sequence, inputs$methods)
# choose one of the following
# inputs$smoothing.rules <- list(e00200="density", e00900="density")
inputs$smoothing.rules <- list(e00200="density") 
# inputs$smoothing.rules <- NULL
inputs

packages <- c("magrittr", "tidyverse", "synthpop")
xport <- c("globals", "puf.base", "inputs")
popts <- list(.packages=packages, .export=xport)
popts


a <- proc.time()
stack <- ldply(inputs$syngroups, runsyn, mslist, msnames, .progress="text", .parallel=FALSE, .paropts=popts)
output <- list()
output$stack <- stack %>% arrange(ftype, pufseqn) %>% mutate(rownum=row_number())
output$inputs <- inputs
saveRDS(output, paste0(globals$tc.dir, "pufsyn_stack/", inputs$synname, "_stack.rds"))
write_csv(output$stack, paste0(globals$tc.dir, "pufsyn_stack/", inputs$synname, "_stack.csv"))
b <- proc.time()
b - a # 3.5 hours

stopCluster(cl)

# CAUTION: Optional: code to build the stacked file from pieces if it was run in parts ----
# MAKE SURE WE STILL HAVE INPUTS AVAILABLE, or else don't write the rds file !!!
# f <- function(groupnum) {
#   group <- mslist[[groupnum]]
#   groupname <- msnames[groupnum]
#   
#   pufsyn <- puf.base %>%
#     ungroup %>% # just to be safe
#     dplyr::filter(eval(group)) %>%
#     sample_frac(inputs$sampfrac) %>%
#     dplyr::select(inputs$vargroups$visit.sequence)
#   
#   synx <- readRDS(paste0(globals$tc.dir, "pufsyn_stack/", inputs$synname,  "_", groupname, ".rds"))
#   synout <- bind_rows(synx$synx$syn) %>% mutate(m=rep(1:synx$synx$m, each=synx$synx$n))
#   stack <- bind_rows(pufsyn %>% mutate(ftype="puf", m=1),
#                      synout %>% mutate(ftype="syn")) %>%
#     dplyr::select(-c(inputs$vargroups$xvars.drop)) %>%
#     dplyr::select(-m, everything(), m) %>%
#     mutate(msname=groupname)
#   return(stack)
# }
# 
# stack <- ldply(1:3, f, .progress="text")
# count(stack, ftype, m)
# output <- list()
# output$stack <- stack %>% arrange(ftype, pufseqn) %>% mutate(rownum=row_number())
# output$inputs <- inputs
# saveRDS(output, paste0(globals$tc.dir, inputs$synname, "_stack.rds"))
# write_csv(output$stack, paste0(globals$tc.dir, inputs$synname, "_stack.csv"))
# CAUTION end build stacked file ----

names(synx$synx)

glimpse(stack)
count(stack, m)
count(stack, ftype)
ht(stack[, c("pufseqn", "MARS", "ftype")])

# NOTE: for synthpop10:
# - using recipe7
# - inputs$smoothing.rules <- list(e00200="density", e00900="density")
# - and ERRONEOUSLY ctree.minbucket <- 10 (I wanted cart.minbucket <- 10)
#   (NOTE: this means that minbucket was default)
# - it worked fine for MARS1 and MARS34, but for MARS2 I got the error
#     Error in do.ply(i) : task 2 failed - "sample is too sparse to find TD"
#   This is some sort of error in bw.SJ kernel density estimation on e00900 (I presume, because)
#   density worked fine when applied only to e00200) - see https://stat.ethz.ch/pipermail/r-help/2008-April/158701.html
#   My guess is that it occurs when the leaves of a terminal node all have one value (e.g., 0), as might happen
#   with e00900. Code: https://github.com/SurajGupta/r-source/blob/master/src/library/stats/R/bandwidths.R
# I reran on MARS2 (a) removing density on e00900, AND (b) setting cart.minbucket <- 10

# save a true-PUF unstacked version of the stacked file to tc.dir. copy MANUALLY to synpuf when ready ----
synname <- "synthpop10"
stack <- readRDS(paste0(globals$tc.dir, synname, "_stack.rds"))$stack
glimpse(stack)
names(stack) %>% sort
count(stack, ftype)

tpuf <- stack %>%
  filter(ftype=="syn") %>%
  dplyr::select(-c(divratio, e01500_minus_e01700, ftype, m, msname, pufseqn, rownum, wt)) %>%
  set_names(str_to_upper(names(.))) %>%
  mutate(ROWNUM=row_number())
glimpse(tpuf)
write_csv(tpuf, paste0(globals$tc.dir, synname, ".csv"))


#****************************************************************************************************
#                Save results ####
#****************************************************************************************************
# create list for the full file with idvars, xvars, synvars.vs, and stack

# use the output from the smallest file as a base
output <- readRDS(paste0(globals$tc.dir, "pufsyn_stack/", synname, "_", "MARS34", ".rds")) # a list

visit.sequence <- output$synx$visit.sequence
output$synx <- list()
output$synx$visit.sequence <- visit.sequence
output$stack <- stack
saveRDS(output, paste0(globals$tc.dir, "pufsyn_stack/", synname, ".rds"))


# save files for unweighted file evaluation - and file for Max
output <- readRDS(paste0(globals$tc.dir, "pufsyn_stack/", synname, ".rds"))
names(output)
names(output$stack)
output$idvars
output$xvars.keep
output$xvars.drop
output$synvars.vs
count(output$stack, ftype)

# sequential, parallel time for 3% 183, 122  = -33%
# sequential, parallel time for 10% 819, 571 = -30%

#****************************************************************************************************
#                Check saved results ####
#****************************************************************************************************
synname <- "synthpop9"
# synname <- "train50"
st2 <- readRDS(paste0(globals$tc.dir, "pufsyn_stack/", synname, "_stack.rds"))
names(synlist)
stack <- synlist$stack
synlist$stack %>% glimpse
names(synlist$stack) %>% sort

synlist$stack %>%
  count(ftype, m, msname)

synlist$stack %>%
  count(ftype, m)

# get wsamp for Dan...?

# write_csv(synlist$stack, paste0(globals$synd, synname, "_all_stack.csv"))

# write file for Max
# synname <- "train50"
# synlist <- readRDS(paste0(globals$tc.dir, "pufsyn_stack/", synname, ".rds"))
# names(synlist)
# glimpse(synlist$synx)
# synlist$idvars
# synlist$xvars.keep
# synlist$xvars.drop
# synlist$synvars.vs
# glimpse(synlist$stack)
# write_csv(synlist$stack %>%
#             filter(ftype=="syn") %>%
#             dplyr::select(-m, -msname),
#           paste0(globals$synd, "puf_samples/", synname, "_synthpop.csv"))



# 
# m1 <- readRDS(paste0(globals$tc.dir, "pufsyn_stack/", "train50_MARS1", ".rds"))
# names(m1)
# glimpse(m1$stack)
# names(m1$synx)
# pm <- m1$synx$predictor.matrix
# colnames(pm)
# pm[, "p22250"]
# pm[, "p22250"]
# pm %>% as.data.frame() %>%
#   gather(vname, value) %>%
#   group_by(vname) %>%
#   summarise(n=sum(value)) %>%
#   ht

# Look at good file and bad file and train
# gf <- readRDS(paste0(globals$tc.dir, "pufsyn_stack/", "train50", ".rds"))$stack %>% dplyr::select(-m, -msname)
# bf <- readRDS(paste0(globals$tc.dir, "pufsyn_stack/BAD/", "train50", ".rds"))$stack %>% dplyr::select(-m, -msname)
# tf <- puf.base %>% dplyr::select(one_of(names(gf)))
# names(tf) %>% sort
# 
# summary(tf)
# summary(gf)
# 
# af <- bind_rows(tf %>% mutate(check="tf"),
#                 gf %>% mutate(check="gf"),
#                 bf %>% mutate(check="bf") %>% dplyr::select(-p22250))
# glimpse(af)
# 
# tmp <- af %>% 
#   dplyr::select(check, ftype, e00200, e01700, e09600, p08000, e09700) %>%
#   gather(vname, value, -check, -ftype) %>%
#   group_by(check, ftype, vname) %>%
#   summarise(mean=mean(value), sd=sd(value))
# 
# tmp %>%
#   arrange(vname, ftype, check)
# 
# tmp %>%
#   dplyr::select(-sd) %>%
#   spread(check, mean) %>%
#   arrange(vname, ftype) %>%
#   dplyr::select(vname, ftype, everything())
# 
# tmp %>%
#   dplyr::select(-mean) %>%
#   spread(check, sd) %>%
#   arrange(vname, ftype) %>%
#   dplyr::select(vname, ftype, everything())
  


# Experimental multidplyr ----

runsyn_mdp <- function(df){
  correct_predictor_matrix <- function(visit.sequence){
    pufsyn <- puf.base %>%
      dplyr::select(visit.sequence)
    test <- syn(pufsyn, visit.sequence=visit.sequence, method=methods, seed=1234, m=0)
    pred.mat.corrected <- test$predictor.matrix
    pred.mat.corrected[, "pufseqn"] <- 0
    return(pred.mat.corrected)
  }
  
  # get easier names
  idvars <- inputs$vargroups$idvars
  xvars.drop <- inputs$vargroups$xvars.drop
  visit.sequence <- inputs$vargroups$visit.sequence
  methods <- inputs$methods
  smoothing.rules <- inputs$smoothing.rules
  m <- inputs$m
  
  # groupname <- groupnames[[i]]; print(groupname)
  # outfname <- paste0(globals$tc.dir, "pufsyn_stack/", synname, "_", groupname, ".rds")
  pufsyn <- df %>%
    dplyr::select(visit.sequence)
  
  pred.mat.corrected <- correct_predictor_matrix(visit.sequence)
  
  seed <- 1234
  synx <- syn(pufsyn, visit.sequence=visit.sequence, method=methods, smoothing=smoothing.rules, seed=seed, m=m, k=nrow(pufsyn), predictor.matrix = pred.mat.corrected, print.flag = TRUE)
  
  #.. combine results - works even when m==1 ----
  synout <- bind_rows(synx$syn) %>% mutate(m=rep(1:synx$m, each=synx$n))
  stack <- bind_rows(pufsyn %>% mutate(ftype="puf", m=1),
                     synout %>% mutate(ftype="syn")) %>%
    dplyr::select(-xvars.drop) %>%
    dplyr::select(-m, everything(), m)
  
  # output <- list()
  # output$inputs <- inputs
  # output$synx <- synx
  # output$stack <- stack
  # saveRDS(output, outfname)
  
  return(stack)
}

# IN PROGRESS: try out multidplyr
# library("multidplyr")
# cluster <- create_cluster(cores = 3)
# cluster_copy(cluster, runsyn_mdp) # Register function -- reregister if function changes
# cluster_copy(cluster, inputs)
# cluster_ls(cluster) # Check registered items
# cluster_assign_value(cluster, 'runsyn_mdp', runsyn_mdp)
# # cluster_get(cluster, 'inputs') # Get the item
# # cluster_rm(cluster, c('runsyn_mdp')) # Unregister function [or other items]
# 
# a <- proc.time()
# stack <- puf.base %>%
#   mutate(msgroup=case_when(eval(sub1) ~ msnames[1],
#                            eval(sub2) ~ msnames[2],
#                            eval(sub34) ~ msnames[3],
#                            TRUE ~ "error")) %>%
#   filter(msgroup %in% msnames[inputs$syngroups]) %>%
#   sample_frac(inputs$sampfrac) %>%
#   partition(msgroup) %>%
#   do(runsyn_mdp(.)) %>%
#   collect() %>% 
#   ungroup
# b <- proc.time()
# b - a

# count(stack, msgroup)
