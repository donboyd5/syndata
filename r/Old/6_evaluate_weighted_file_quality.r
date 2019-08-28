

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

library("btools") # library that I created (install from github)

# library("ipoptr")
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
# source("./r/includes/functions_target_setup_and_analysis.r")
# source("./r/includes/functions_ipopt.r")


#****************************************************************************************************
#  Define syn file and get base file we will need for much work ####
#****************************************************************************************************
# two major steps:
# 1. run tax calc on the 2 UNWEIGHTED files
# 2. get a set of weights to look at weighted results

# sfname <- "synthpop4"
# sfname <- "synpuf8"
# sfname <- "synthpop5_all"
sfname <- "synthpop6_all"

# goal:
# create stacked puf and syn files (from prep), and create unique RECID (MUST be called this as we use to merge with tc output)
# create weights file for the synfile with same unique RECID
# merge results and analyze

prep <- readRDS(paste0(globals$tc.dir, sfname, "_rwprep.rds"))
names(prep)
names(prep$tc.output) %>% sort
names(prep$tc.base) %>% sort
groups(prep$tc.base)
groups(prep$tc.output)


# anyDuplicated(prep$tc.base$RECID)
# prep$tc.base$RECID %>% sort %>% ht
# prep$tc.base %>%
#   group_by(ftype) %>%
#   summarise(rmin=min(RECID), rmax=max(RECID))



#****************************************************************************************************
#  Get weights so that we have them for post-reform analysis ####
#****************************************************************************************************
#.. 2. get the weights ----
# get the original weights from puf syn
# names(prep$tc.base) %>% sort
wt.puf <- prep$tc.base %>%
  filter(ftype=="puf") %>%
  dplyr::select(ftype, m, pufseqn, wt.puf=wt) %>%
  mutate(wtype="puf")

wt.syn <- prep$tc.base %>%
  filter(ftype=="syn") %>%
  dplyr::select(ftype, m, pufseqn, wt.syn=wt) %>%
  mutate(wtype="syn")
wt.syn$RECID %>% sort

wt.rwt <- readRDS(paste0(globals$tc.dir, sfname, "_rwt.rds")) %>%
  dplyr::select(ftype, m, pufseqn, wt.rwt) %>%
  mutate(wtype="rwt")

wt.wtfs <- readRDS(paste0(globals$tc.dir, sfname, "_wtfs.rds")) %>%
  dplyr::select(ftype, m, pufseqn, wt.wtfs) %>%
  mutate(wtype="wtfs")

# wt.wtfsv2 <- readRDS(paste0(globals$tc.dir, sfname, "_wtfs_v2.rds")) %>%
#   dplyr::select(ftype, m, pufseqn, wt.wtfsv2=wt.wtfs) %>%
#   mutate(wtype="wtfsv2")
# names(wt.wtfsv2) %>% sort

# weights <- bind_rows(
#   wt.puf %>% dplyr::select(m, pufseqn, wt=wt.puf, wtype),
#   wt.syn %>% dplyr::select(m, pufseqn, wt=wt.syn, wtype),
#   wt.rwt %>% dplyr::select(m, pufseqn, wt=wt.rwt, wtype),
#   wt.wtfs %>% dplyr::select(m, pufseqn, wt=wt.wtfs, wtype),
#   wt.wtfsv2 %>% dplyr::select(m, pufseqn, wt=wt.wtfsv2, wtype))

weights <- bind_rows(
  wt.puf %>% dplyr::select(m, pufseqn, wt=wt.puf, wtype),
  wt.syn %>% dplyr::select(m, pufseqn, wt=wt.syn, wtype),
  # wt.rwt %>% dplyr::select(m, pufseqn, wt=wt.rwt, wtype),
  wt.wtfs %>% dplyr::select(m, pufseqn, wt=wt.wtfs, wtype))

count(weights, wtype)

weights %>%
  filter(wtype!="puf") %>% 
  select(m, pufseqn, wtype, wt) %>% 
  spread(wtype, wt) %>% 
  select(-m, -pufseqn) %>%
  cor(.) %>% 
  kable(digits=4)

# look at distribution of weights
nsamples <- max(prep$tc.base$m)
wdist <- weights %>%
  mutate(wt=ifelse(wtype %in% c("puf", "syn"), wt / nsamples, wt)) %>% # to put them on even footing for distribution
  mutate(wtype=factor(wtype, levels=c("puf", "syn", "rwt", "wtfs", "wtfsv2")))
count(wdist, wtype)

wdist %>%
  group_by(wtype) %>%
  do(qtiledf(.$wt, probs=c(0, .005, .01, .05, .1, .25, .5, .75, .9, .95, .99, .995, 1)))

p <- wdist %>%
  ggplot(aes(x=wt, y = ..density..)) +
  geom_histogram(binwidth=25, fill="blue") +
  geom_vline(aes(xintercept = median(wt))) +
  scale_x_continuous(breaks=seq(0, 5000, 250), limits=c(0, 1200)) +
  theme(axis.text.x=element_text(size=8, angle=30)) +
  facet_wrap(~wtype, nrow=5) +
  ggtitle("Distribution of weights, adjusted for # of records (syn has 3x # of records that puf has)")
p


#****************************************************************************************************
#  Run tax calc on the 2 UNWEIGHTED files and save results ####
#****************************************************************************************************
#.. 1. run tax calc on the 2 UNWEIGHTED files ----

# IMPORTANT: make sure we already have a proper RECID

pufsyn <- prep$tc.base %>%
  mutate(c00100=prep$tc.output$c00100, taxbc=prep$tc.output$taxbc) %>%
  arrange(ftype, m, pufseqn)
count(pufsyn, ftype)

stack <- pufsyn %>% 
  arrange(RECID)
glimpse(stack)
count(stack, ftype)
names(stack) %>% sort
# note that we have taxbc and RECID on the file

# $ conda update -c PSLmodels taxcalc -- make sure cmd prompt is run as administrator
# we call get_pufsyn_taxplan_results(stacked.pufsyn, taxplan.fn=NULL)
# get_pufsyn_taxplan_results requires as inputs:
#    stacked.pufsyn -- a dataframe with stacked tax syn file and its corresponding puf (same variables) -- variable name case will be adjusted
#    taxplan.fn -- character variable holding the name (without directory) of the taxplan json file e.g., "rate_cut.json"

# for examples, see https://github.com/PSLmodels/Tax-Calculator/tree/master/taxcalc/reforms

# reform <- "rate_cut"
# reform <- "toprate"
# reform <- "EITC"
# reform <- "Trump2013"
# reform <- "rt5_2013"
# reform <- "brk4_2013"
# reform <- "brk4_1k_2013"
# reform <- "brk4_10k_2013"

a <- proc.time()
reform <- "rt5_2013"
reform.fn <- paste0(reform, ".json")
system.time(df1 <- get_pufsyn_taxplan_results(stack, reform.fn)) # takes about 5 minutes for # of recs in puf
saveRDS(df1, paste0(globals$tc.dir, sfname, "_", reform, ".rds"))
b <- proc.time()
b - a
# names(df1$tc.output) %>% sort

# tc.wincmd("tcbase.csv", globals$tc.dir, globals$tc.cli, taxyear=2013, taxplan.fn="brk4_2013.json", taxplans.dir=globals$taxplans.dir)


#****************************************************************************************************
#  Retrieve results for tax plans and examine weighted results ####
#****************************************************************************************************
# tcvars <- c("c00100", "standard", "c17000", "c18300", "c21060", "c04800", "taxbc", "c09600", "c05800")
# c00100 AGI
# c17000 Sch A: Medical expenses deducted
# c18300 Sch A: State and local taxes deducted
# c21060 Itemized deductions before phase-out (zero for non-itemizers)
# standard Standard deduction (zero for itemizers)
# c04800 Regular taxable income
# taxbc regular tax before credits
# c09600 Alternative Minimum Tax (AMT) liability
# c05800 taxbc plus AMT liability

# df1 <- readRDS(paste0(globals$tc.dir, sfname, "_", str_remove("rate_cut.json", ".json"), ".rds"))
# df1 <- readRDS(paste0(globals$tc.dir, sfname, "_", str_remove("rt5_2013.json", ".json"), ".rds"))
# df2 <- readRDS(paste0(globals$tc.dir, sfname, "_", str_remove("Trump2013.json", ".json"), ".rds"))
# count(df1$tc.base, ftype, m)
# count(df2$tc.base, ftype, m)

# synthpop5_all_rate_cut.rds
# synthpop5_all_Trump2013.rds
# synthpop5_all_rt5_2013.rds
# synthpop5_all_brk4_2013.rds
# synthpop5_all_brk4_1k_2013.rds
# synthpop5_all_brk4_10k_2013.rds
reforms <- c("rate_cut", "Trump2013", "rt5_2013", "brk4_2013", "brk4_1k_2013", "brk4_10k_2013")

getref <- function(refname, sfname="synthpop5_all"){
  readRDS(paste0(globals$tc.dir, sfname, "_", refname, ".rds"))
}

df1 <- getref(reforms[1])
df2 <- getref(reforms[2])
df3 <- getref(reforms[3])
df4 <- getref(reforms[4])
df5 <- getref(reforms[5])
df6 <- getref(reforms[6])

# names(df$tc.output) %>% sort

df <- df6; reform.fn <- reforms[6]
df <- df3; reform.fn <- reforms[3]
df <- df2; reform.fn <- reforms[2]

mrgdf2 <- left_join(df$tc.base, df$tc.output %>% dplyr::select(RECID, taxbc.reform=taxbc), by="RECID")
glimpse(mrgdf2)
count(mrgdf2, ftype)
names(mrgdf2) %>% sort

# highlight all of this to run!
f1 <- mrgdf2 %>% 
  filter(ftype=="puf") %>% 
  mutate(wtype="puf")
ht(f1)

f2 <- mrgdf2 %>% 
  filter(ftype=="syn") %>% 
  dplyr::select(-wt) %>% 
  left_join(weights %>% filter(wtype=="syn") %>% dplyr::select(wtype, m, pufseqn, wt))
count(f2, wtype)

f3 <- mrgdf2 %>% 
  filter(ftype=="syn") %>% 
  dplyr::select(-wt) %>% 
  left_join(weights %>% filter(wtype=="rwt") %>% dplyr::select(wtype, m, pufseqn, wt))

f4 <- mrgdf2 %>% 
  filter(ftype=="syn") %>% 
  dplyr::select(-wt) %>% 
  left_join(weights %>% filter(wtype=="wtfs") %>% dplyr::select(wtype, m, pufseqn, wt))

f5 <- mrgdf2 %>% 
  filter(ftype=="syn") %>% 
  dplyr::select(-wt) %>% 
  left_join(weights %>% filter(wtype=="wtfsv2") %>% dplyr::select(wtype, m, pufseqn, wt))

stack2 <- bind_rows(f1, f2, f3, f4, f5)
count(stack2, wtype)
glimpse(stack2)


f <- function(var, wt){sum(var * wt) / 1e9}
stack2 %>%
  group_by(wtype) %>%
  summarise(wtm=sum(wt / 1e6), 
            agib=f(c00100, wt), 
            wageb=f(e00200, wt), 
            taxbcb=f(taxbc, wt),
            taxbc.reform=f(taxbc.reform, wt)) %>%
  gather(vname, value, -wtype) %>%
  spread(wtype, value) %>%
  mutate_at(vars(syn, rwt, wtfs, wtfsv2), funs(diff=. - puf, pdiff=(. - puf) / puf * 100)) %>%
  kable(digits=1)

# get the % change in tax by file
dtot <- function(df){
  dsums <- df %>% 
    summarise_at(vars(-wtype, -agirange), funs(sum)) %>%
    mutate(wtype=df$wtype[1], agirange="Total")
  dfout <- bind_rows(df, dsums)
  return(dfout)
}


agiranges <- c(-Inf, 0, 25e3, 50e3, 75e3, 100e3, 200e3, 500e3, 1e6, 10e6, Inf)

diffs <- stack2 %>%
  mutate(agirange=cut(c00100, agiranges, right=FALSE)) %>%
  group_by(wtype, agirange) %>%
  summarise_at(vars(taxbc, taxbc.reform), funs(sum(. * wt) / 1e9)) %>%
  do(dtot(.)) %>%
  mutate(agirange=factor(agirange, levels=unique(agirange), ordered=TRUE),
         diff=taxbc.reform - taxbc,
         pdiff=diff / taxbc * 100)

diffs %>%
  dplyr::select(wtype, agirange, taxbc) %>%
  spread(wtype, taxbc) %>%
  mutate(var="taxbc$b", units="$ billions", reform=reform.fn) %>%
  kable(digits=1)

diffs %>%
  dplyr::select(wtype, agirange, taxbc) %>%
  spread(wtype, taxbc) %>%
  mutate_at(vars(syn, rwt, wtfs, wtfsv2), funs(. / puf * 100 - 100)) %>%
  mutate(var="baseline taxbc: synfiles % diff from puf", reform=reform.fn) %>%
  kable(digits=1)

diffs %>%
  dplyr::select(wtype, agirange, taxbc) %>%
  spread(wtype, taxbc) %>%
  mutate_at(vars(syn, rwt, wtfs, wtfsv2), funs(. - puf)) %>%
  mutate(var="baseline taxbc: synfiles minus puf", reform=reform.fn) %>%
  kable(digits=1)

diffs %>%
  dplyr::select(wtype, agirange, taxbc.reform) %>%
  spread(wtype, taxbc.reform) %>%
  mutate(var="taxbc.reform$b", units="$ billions", reform=reform.fn) %>%
  kable(digits=1)

diffs %>%
  dplyr::select(wtype, agirange, taxbc.reform) %>%
  spread(wtype, taxbc.reform) %>%
  mutate_at(vars(syn, rwt, wtfs, wtfsv2), funs(. / puf * 100 - 100)) %>%
  mutate(var="tax reform taxbc: synfiles % diff from puf", reform=reform.fn) %>%
  kable(digits=1)

diffs %>%
  dplyr::select(wtype, agirange, diff) %>%
  spread(wtype, diff) %>%
  mutate(units="$ billions change from baseline", reform=reform.fn) %>%
  kable(digits=1)

diffs %>%
  dplyr::select(wtype, agirange, pdiff) %>%
  spread(wtype, pdiff) %>%
  mutate(units="% change from baseline", reform=reform.fn) %>%
  kable(digits=1)

diffs %>%
  dplyr::select(wtype, agirange, pdiff) %>%
  spread(wtype, pdiff) %>%
  mutate_at(vars(syn, rwt, wtfs, wtfsv2), funs(. - puf)) %>%
  mutate(units="puf=% change from baseline, syn=synfile % change minus puf % change", reform=reform.fn) %>%
  kable(digits=1)

diffs %>%
  dplyr::select(wtype, agirange, diff) %>%
  spread(wtype, diff) %>%
  mutate_at(vars(syn, rwt, wtfs, wtfsv2), funs(. - puf)) %>%
  mutate(units="puf=$ billions reform impact, syn=synfile $b impact minus puf impact", reform=reform.fn) %>%
  kable(digits=1)

diffs %>%
  dplyr::select(wtype, agirange, diff) %>%
  spread(wtype, diff) %>%
  mutate_at(vars(syn, rwt, wtfs, wtfsv2), funs(. / puf * 100 - 100)) %>%
  mutate(var="reform impact: synfiles % diff from puf", reform=reform.fn) %>%
  kable(digits=1)


#****************************************************************************************************
#  ONETIME: create a file with wsamp for Dan, to write as csv to synpuf ####
#****************************************************************************************************
#.. create a file to write as csv to synpuf ----
# get wsamp for Dan
wsamp <- readRDS(paste0(globals$tc.dir, "puf_lc.rds")) %>%
  dplyr::select(pufseqn, wsamp) %>%
  mutate(ftype="puf")
ht(wsamp) # includes aggregate records

file.out <- pufsyn %>%
  mutate(wt.puf=ifelse(ftype=="puf", wt, NA_real_)) %>%
  left_join(weights) %>% 
  dplyr::select(-wt, -divratio, -e01500_minus_e01700) %>%
  left_join(wsamp)
file.out$RECID %>% sort %>% ht
summary(file.out %>% select(starts_with("w")))
names(file.out) %>% sort
anyDuplicated(file.out$RECID)
anyDuplicated(file.out$pufseqn)
count(file.out, ftype, m)

file.out %>%
  write_csv(paste0(globals$synd, "synthpop5_all_for_revenuescores.csv"))

#.. end write to csv ----
