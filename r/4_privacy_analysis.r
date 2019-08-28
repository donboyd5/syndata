

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


#****************************************************************************************************
#                Includes ####
#****************************************************************************************************
source("./r/includes/globals_system_specific_boyd.r") # use a different version of this file if changing systems
source("./r/includes/globals_other.r")

source("./r/includes/functions_general.r")
source("./r/includes/functions_synthesis.r")


#****************************************************************************************************
#                Functions ####
#****************************************************************************************************


#****************************************************************************************************
#                Get puf.names, and puf or other files ####
#****************************************************************************************************
puf.vnames <- get_puf_vnames()
# puf.vnames

psum <- get_pufbase_sums()
psum

#****************************************************************************************************
#                Get data ####
#****************************************************************************************************
sfname <- "synthpop5_all.rds"
synprep <- readRDS(paste0(globals$tc.dir, str_remove(sfname, ".rds"), "_rwprep.rds"))

names(synprep)
# merge and then split
tcvars <- c("c00100", "taxbc") # taxcalc vars
mrgdf <- left_join(synprep$tc.base, synprep$tc.output %>% dplyr::select(RECID, tcvars))
glimpse(mrgdf)
count(mrgdf, ftype)


#******************************************************************************************************************
#  pick a subset ####
#******************************************************************************************************************


#******************************************************************************************************************
#  FULL FILE distance calculations ####
#******************************************************************************************************************
glimpse(mrgdf)
count(mrgdf, ftype)

#.. preparation BEFORE distance calcs ----
stack <- mrgdf %>%
  filter(m==1) %>%
  # sample_frac(size=.1) %>%
  mutate(mgroup=case_when(MARS==1 ~ 1,
                          MARS==2 ~2,
                          MARS %in% 3:4 ~3,
                          TRUE ~ -9))
count(stack, ftype, mgroup)

#.. define variables to use in distance calculations ----
# distvars <- c("c00100", "e00200", "e00300", "e00400", "XTOT") # for testing

# distance vars for real
stacknames <- names(stack)
stacknames %>% sort
keep1 <- stacknames[str_sub(stacknames, 1, 1) %in% c("c", "e", "p", "n", "t")]
keep.superset <- c(keep1, "XTOT")
dropvars <- c(stacknames[str_sub(stacknames, -1) %in% c("p", "s")], "e01500_minus_e01700", "pufseqn")
distvars <- setdiff(keep.superset, dropvars) %>% sort
distvars
nvars <- length(distvars)

# e30500 is a problem variable...several groups have sd=0 (can't divide by that when scaling)
stack %>% 
  rename(val=e30500) %>%
  group_by(ftype, mgroup) %>%
  summarise(min=min(val), mean=mean(val), max=max(val), sd=sd(val),
            p10=pany(val, .1), p50=pany(val, .5), p90=pany(val, .9))

stack.scaled <- stack %>%
  select(RECID, ftype, mgroup, distvars) %>%
  # scale and do ntiles within each ftype mgroup combination
  group_by(ftype, mgroup) %>%
  # within file type, mgroup, scale distances to [0, 1] -- subtract mean then divide by sd
  do(data.frame(RECID=.$RECID, scale(.[, distvars]))) %>%
  mutate_at(vars(distvars), funs(ifelse(is.nan(.), NA_real_, .))) %>% 
  mutate(ntile=ntile(c00100, 100)) %>%
  ungroup %>%
  arrange(ftype, mgroup, RECID)
glimpse(stack.scaled)
ht(stack.scaled)
summary(stack.scaled)
stack.scaled %>%
  group_by(ftype, mgroup) %>%
  summarise(e00200.mean=mean(e00200), e00200.sd=sd(e00200))

# check to see if any remaining bad data
# stack.scaled %>%
#   select(ftype, mgroup, e30500) %>%
#   mutate(value=case_when(is.nan(e30500) ~ "nan",
#                          is.na(e30500) ~ "na",
#                          is.infinite(e30500) ~ "inf",
#                          is.finite(e30500) ~ "finite",
#                          TRUE ~ "error")) %>%
#   group_by(ftype, mgroup, value) %>%
#   summarise(n=n())

a <- stack.scaled %>% 
  filter(ftype=="puf")

b <- stack.scaled %>% 
  filter(ftype=="syn")

glimpse(a)
glimpse(b)



get_nd <- function(ntile.in, mgroup.in, a, b, distvars, globals){
  # get neighbor and distance for a particular mgroup and ntile
  # define the needed functions in here to make parallelization easy
  edist <- function(vec, mat) {
    # euclidean distance between a vector and each row of a matrix
    # returns a vector of distances with same length as vec
    mat <- as.matrix(mat)
    matdiff <- sweep(mat, 2, vec) # subtract vec from each row of mat 
    # convert NA values to 0, which result from scaling of variables with undefined SD, so they don't affect distance
    matdiff[is.na(matdiff)] <- 0
    return(sqrt(rowSums(matdiff ^2)))
  }
  
  edist.ab <- function(amat, bmat){
    # euclidean distance between each row of an a matrix and each row of a b matrix
    # returns a long vector with all of the distances
    # the first subvector is distances between row 1 of amat and all rows of bmat
    # the next subvector is distances between row 2 of amat and all rows of bmat
    # and so on
    ab <- apply(amat, 1, edist, bmat)
    return(as.vector(ab))
  }
  
  # set up the range over which we will search
  nbound <- 10
  b.min <- max(ntile.in - nbound, 1)
  b.max <- min(ntile.in + nbound, 100)
  b.ntiles <- b.min:b.max
  
  # get an a record for this mgroup.in and ntile.in
  # get ALL b records for this mgroup.in that are in the b.ntiles
  a.ntile.df <- a %>% filter(mgroup==mgroup.in, ntile==ntile.in) %>% arrange(mgroup, ntile, RECID)
  b.ntiles.df <- b %>% filter(mgroup==mgroup.in, ntile %in% b.ntiles) %>% arrange(mgroup, ntile, RECID)
  
  ab.ids <- full_join(a.ntile.df %>% select(a.RECID=RECID, a.mgroup=mgroup, a.ntile=ntile) %>% mutate(fake=1),
                      b.ntiles.df %>% select(b.RECID=RECID, b.mgroup=mgroup, b.ntile=ntile) %>% mutate(fake=1),
                      by="fake") %>%
    select(a.mgroup, a.ntile, a.RECID, b.RECID, b.ntile, b.mgroup) %>%
    arrange(a.mgroup, a.ntile, a.RECID, b.mgroup, b.ntile, b.RECID)
  
  ab.dist <- ab.ids %>%
    mutate(dist=edist.ab(a.ntile.df[, distvars], b.ntiles.df[, distvars]))
  # write this all-distances file to a permanent directory, but only return the nearest neighbor
  # write_csv(ab.dist, paste0(globals$tc.dir, "dist_pieces/", "ab.dist_m", mgroup.in, "_ntile", ntile.in, ".rds"))
  saveRDS(ab.dist, paste0(globals$tc.dir, "dist_pieces/", "ab.dist_m", mgroup.in, "_ntile", ntile.in, ".rds"))
  
  nneighbor <- ab.dist %>%
    group_by(a.mgroup, a.ntile, a.RECID) %>%
    arrange(dist) %>%
    filter(row_number()==1) %>% # don't allow ties
    ungroup %>%
    arrange(a.mgroup, a.ntile, a.RECID, dist)
  return(nneighbor)
}


library("doParallel")
cl <- makeCluster(6)
registerDoParallel(cl)

packages <- c("magrittr", "tidyverse")
xport <- c("globals")
popts <- list(.packages=packages, .export=xport)

ntiles.todo <- 100

time.a <- proc.time()
nn.df1 <- ldply(1:ntiles.todo, get_nd, mgroup.in=1, a=a, b=b, distvars=distvars, globals=globals, .progress="none", .parallel=TRUE, .paropts=popts)
nn.df2 <- ldply(1:ntiles.todo, get_nd, mgroup.in=2, a=a, b=b, distvars=distvars, globals=globals, .progress="none", .parallel=TRUE, .paropts=popts)
nn.df3 <- ldply(1:ntiles.todo, get_nd, mgroup.in=3, a=a, b=b, distvars=distvars, globals=globals, .progress="none", .parallel=TRUE, .paropts=popts)
time.b <- proc.time()
time.b - time.a # a bit less than 2 hours on full puf, all vars, vs. 1 sample from full syn

nn.df <- bind_rows(nn.df1, nn.df2, nn.df3)
saveRDS(nn.df, paste0(globals$tc.dir, "synthpop5_nn.rds"))

# nn.df1
# nn.df2
# nn.df3


#******************************************************************************************************************
#  Examine results -- FULL FILE distance calculations ####
#******************************************************************************************************************

badvars <- c("divratio", "e01500_minus_e01700", "MARS", "m", "pufseqn", "wt")



glimpse(mrgdf)
count(mrgdf, ftype)

#.. get the needed a and b files ----
tcvars <- c("c00100", "taxbc") # taxcalc vars
mrgdf <- left_join(synprep$tc.base, synprep$tc.output %>% dplyr::select(RECID, tcvars))

stack <- mrgdf %>%
  filter(m==1) %>%
  mutate(mgroup=case_when(MARS==1 ~ 1,
                          MARS==2 ~2,
                          MARS %in% 3:4 ~3,
                          TRUE ~ -9))
count(stack, ftype, mgroup)

a.raw <- stack %>% 
  filter(ftype=="puf")

b.raw <- stack %>% 
  filter(ftype=="syn")
glimpse(a.raw)


#.. get nearest neighbor indexes ----
nn.df <- readRDS(paste0(globals$tc.dir, "synthpop5_nn.rds"))
ht(nn.df)
nn.df %>%
  group_by(a.mgroup, b.mgroup) %>%
  do(qtiledf(.$dist))

# get average distance
nn.df %>%
  summarise(mean=mean(dist), median=median(dist), min=min(dist), max=max(dist), sd=sd(dist)) %>%
  gather(dist.stat, stat.value) %>%
  mutate(stat.avg.per.var=stat.value / nvars) %>%
  select(dist.stat, stat.value) %>%
  kable(digits=3)


# now get the benedetto privacy measure:
#  for each variable 
#  - get rmse of puf to nearest syn record
#  - standardize as rmse / (sample sd of var)
# average across all variables
# worst protection is rmse=0 -- all syn vars equal puf record, for all vars, all records
# best protection is rmse=1 (actually, can be greater than 1)
# their table 1 shows rmse=.95 for net worth and 2.21 for fica earnings
# their figure 1 shows about 0.2 for privacy protection
# I get 0.48 for all vars, which is better protection than their 0.2
#   virtually every variable is above 0.2


nndata <- nn.df %>%
  left_join(a.raw %>% setNames(paste0("a.", names(a.raw)))) %>%
  left_join(b.raw %>% setNames(paste0("b.", names(b.raw))))
glimpse(nndata)
nndata %>% dplyr::select(a.RECID, b.RECID, b.ntile, a.e00800, b.e00800) %>% ht
nndata %>% dplyr::select(a.RECID, b.RECID, b.ntile, a.e30500, b.e30500) %>% ht

t1 <- proc.time()
nn2 <- nndata %>%
  dplyr::select(-dist, -a.ftype, -b.ftype, -a.msname, -b.msname) %>%
  gather(var, value, -a.mgroup, -b.mgroup, -a.ntile, -b.ntile, -a.RECID, -b.RECID) %>%
  # separate(var, c("file", "vname")) # separate takes MUCH longer than mutate in this case
  mutate(file=str_sub(var, 1, 1),
         vname=str_sub(var, 3, -1)) %>%
  dplyr::select(-var) %>%
  spread(file, value)
proc.time() - t1 # 2+ mins
glimpse(nn2) # 17+m records
ht(nn2)


nnstats <- nn2 %>%
  filter(!vname %in% badvars) %>%
  group_by(vname) %>%
  summarise(rmse=rmse(a, b), sd=sd(a)) %>%
  mutate(rmse.std=rmse / sd)
nnstats$vname %>% sort

nnstats %>% arrange(vname)
nnstats %>% arrange(rmse.std) %>% ht


nnstats %>%
  summarise(rmse.std.mean=mean(rmse.std, na.rm=TRUE), rmse.std.mdn=median(rmse.std, na.rm=TRUE))



#******************************************************************************************************************
#  test analysis - distance calculations ####
#******************************************************************************************************************
# stack <- mrgdf %>%
#   filter(MARS==2, c00100>=0, c00100<=100e3, m==1)

stack <- mrgdf %>%
  filter(MARS==2, m==1)


stack %>% group_by(ftype) %>% summarise(RECID.min=min(RECID), RECID.max=max(RECID))
anyDuplicated(stack$RECID)
# good, we can use RECID as the id variable


# define variables to use in distance calculations
# distvars <- c("c00100", "e00200", "e00300", "e00400", "XTOT")
stacknames <- names(stack)
stacknames %>% sort
str_subset(stacknames, )
keep1 <- stacknames[str_sub(stacknames, 1, 1) %in% c("c", "e", "p", "n", "t")]
keep.superset <- c(keep1, "XTOT")
dropvars <- c(stacknames[str_sub(stacknames, -1) %in% c("p", "s")], "e01500_minus_e01700", "pufseqn")
distvars <- setdiff(keep.superset, dropvars) %>% sort
distvars

df <- stack %>%
  select(RECID, ftype, distvars) %>%
  group_by(ftype) %>%
  do(data.frame(RECID=.$RECID, scale(.[, distvars]))) %>% # scale the distances to [0, 1], within file type
  mutate(ntile=ntile(c00100, 100)) %>%
  ungroup %>%
  arrange(ftype, RECID)
glimpse(df)
ht(df)


a <- df %>% 
  filter(ftype=="puf")

b <- df %>% 
  filter(ftype=="syn")

glimpse(a)
glimpse(b)
ht(a)

nbound <- 10
a.ntile <- 1
b.min <- max(a.ntile - nbound, 1)
b.max <- min(a.ntile + nbound, 100)
b.ntiles <- b.min:b.max
b.ntiles
an.df <- a %>% filter(ntile==a.ntile) %>% arrange(ntile, RECID)
bns.df <- b %>% filter(ntile %in% b.ntiles) %>% arrange(ntile, RECID)

# 
system.time(d <- edist.ab(an.df[, distvars], bns.df[, distvars])) # 12 secs with 90 distance variables
# for MARS==2, m==1
# it is approx 80 secs
# we need to do this ~ 300 times
# so total runtime ~ 6.7 hours
# or 
length(d)
d
min(d)
quantile(d)

# 100 x 12 x 3 ~ 1 hour

ab.ids <- full_join(an.df %>% select(a.RECID=RECID, a.ntile=ntile) %>% mutate(fake=1),
                    bns.df %>% select(b.RECID=RECID, b.ntile=ntile) %>% mutate(fake=1)) %>%
  select(a.ntile, a.RECID, b.RECID, b.ntile) %>%
  arrange(a.ntile, a.RECID, b.ntile, b.RECID)

ab.dist <- ab.ids %>%
  mutate(dist=d)
count(ab.ids, a.ntile)
length(unique(ab.ids$a.RECID))

nneighbor <- ab.dist %>%
  group_by(a.ntile, a.RECID) %>%
  arrange(dist) %>%
  filter(row_number()==1) %>% # don't allow ties
  ungroup %>%
  arrange(a.ntile, a.RECID, dist)
count(nneighbor, b.ntile)

# now get the benedetto privacy measure:
#  for each variable 
#  - get rmse of puf to nearest syn record
#  - standardize as rmse / (sample sd of var)
# average across all variables
# worst protection is rmse=0 -- all syn vars equal puf record, for all vars, all records
# best protection is rmse=1 (actually, can be greater than 1)
# their table 1 shows rmse=.95 for net worth and 2.21 for fica earnings
# their figure 1 shows about 0.2 for privacy protection
# I get 2.09 for all vars

nndata <- nneighbor %>%
  left_join(a %>% setNames(paste0("a.", names(a)))) %>%
  left_join(b %>% setNames(paste0("b.", names(b))))
glimpse(nndata) 
nndata %>% select(a.RECID, b.RECID, b.ntile, a.e00800, b.e00800)


nn2 <- nndata %>%
  select(-a.ftype, -b.ftype, -b.ntile, -dist) %>%
  gather(var, value, -a.ntile, -a.RECID, -b.RECID) %>%
  separate(var, c("file", "vname")) %>%
  spread(file, value)
glimpse(nn2)

nnstats <- nn2 %>%
  group_by(a.ntile, vname) %>%
  summarise(rmse=rmse(a, b), sd=sd(a)) %>%
  mutate(rmse.std=rmse / sd)
nnstats %>% arrange(vname)
nnstats %>% arrange(rmse.std)

nnstats %>% arrange(rmse.std) %>% ht(10)

nnstats %>%
  filter(!is.infinite(rmse.std)) %>%
  group_by(a.ntile) %>%
  summarise(rmse.std.mean=mean(rmse.std, na.rm=TRUE), rmse.std.mdn=median(rmse.std, na.rm=TRUE))


#******************************************************************************************************************
#  Categorical variables and sample uniques ####
#******************************************************************************************************************
# find the categoricals
names(mrgdf) %>% sort
cats <- c("n24", "MARS", "xfpt", "xfst", "xocah", "XTOT")

cats.e <- expression(n24, MARS, xfpt, xfst, xocah, XTOT)

catcounts <- mrgdf %>% 
  group_by(ftype, n24, MARS, xfpt, xfst, xocah, XTOT) %>%
  summarise(n=n())
catcounts %>% arrange(n)

catcounts %>%
  filter(n==1) %>%
  arrange(n24, MARS, xfpt, xfst, xocah, XTOT)

matches <- mrgdf %>%
  filter(n24==0, MARS==2, xfpt==1, xfst==0, xocah==0, XTOT==1)
matches # no one would mistake these for the same returns

# closeness of sample unique to population unique
matches %>%
  select(-m, -divratio, -e01500_minus_e01700, -RECID, -pufseqn, -msname) %>%
  select(-ends_with("p"), -ends_with("s")) %>%
  gather(vname, value, -ftype) %>%
  spread(ftype, value) %>%
  left_join(psum %>% select(vname, vdesc, puf.billions=value)) %>%
  mutate(sortgroup=case_when(vname %in% c("n24", "MARS", "xfpt", "xfst", "xocah", "XTOT", "wt") ~ 1,
                             vname %in% c("c00100", "taxbc") ~ 2,
                             TRUE ~ 3)) %>%
  group_by(sortgroup) %>%
  arrange(sortgroup, -puf.billions) %>%
  kable(digits=c(0, 0, 0, 0, 1, 0), format.args=list(big.mark = ','))


# k-anonymity ----
# This is simply a count of the number of individuals with a sample frequency of their key lower than k. 
# In a safe dataset, the number of individuals sharing the same combination of values (keys) of
# categorical quasi-identifiers should be higher than a specified threshold k

catcounts %>%
  filter(ftype=="syn", n<=20) %>%
  arrange(n) %>%
  mutate(k=n) %>%
  group_by(k) %>%
  summarise(n=sum(n)) %>%
  mutate(n_violated=lag(cumsum(n))) %>%
  select(-n) %>%
  filter(k > 1)

library("sdcMicro")
selectedKeyVars <- c("n24", "MARS", "xfpt", "xfst", "xocah", "XTOT")

# selected linked variables (ghost variables)
#selectedGhostVars = c('urbrur')

# selected categorical numerical variables
#selectedNumVar = c('wage', 'savings')

# weight variable
#selectedWeightVar = c('wgt')

# selected pram variables
#selectedPramVars = c('roof', 'wall')

# household id variable (cluster)
#selectedHouseholdID = c('idh')

# stratification variable
#selectedStrataVar = c('strata')

# sensitive variables for l-diversity computation
#selectedSensibleVar = c('health')

# creating the sdcMicro object with the assigned variables
# sdcInitial <- createSdcObj(dat         = file,
#                            keyVars     = selectedKeyVars,
#                            ghostVars   = selectedGhostVars,
#                            numVar      = selectedNumVar,
#                            weightVar   = selectedWeightVar,
#                            pramVars    = selectedPramVars,
#                            hhId        = selectedHouseholdID,
#                            strataVar   = selectedStrataVar,
#                            sensibleVar = selectedSensibleVar)

sdcdat <- mrgdf %>% filter(ftype=="syn" ) %>% select(selectedKeyVars)
sdcInitial <- createSdcObj(dat = sdcdat,
                           keyVars = selectedKeyVars)


# Summary of object
sdcInitial
# this gives the same answer for k-anonymity as I calculated

# the sdc risk measure
tmp <- sdcInitial@risk$individual
str(tmp)
max(tmp[, 1])
tdf <- tibble(rn=1:nrow(tmp), risk=tmp[, 1])
glimpse(tdf)
tdf %>%
  arrange(-risk) %>%
  head(40)
