

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

library("Hmisc")

library("janitor")


#****************************************************************************************************
#                Includes ####
#****************************************************************************************************
source("./r/includes/globals_system_specific_boyd.r") # use a different version of this file if changing systems
source("./r/includes/globals_other.r")

source("./r/includes/functions_general.r")
source("./r/includes/functions_synthesis.r")

puf.vnames <- get_puf_vnames()

seedvars <- c("MARS", "wt", "e00100", "e04600", "p04470", "e04800", "e62100", "e05800", "e08800", "e59560", "e26190")
puf.vnames %>% filter(vname %in% seedvars)


pufbase <- get_puf.base()

#****************************************************************************************************
#               Get data ####
#****************************************************************************************************

synname <- "synthpop10" 
pufsyn <- read_csv(paste0(globals$tc.dir, synname, "_stack.csv"), 
                   col_types = cols(ftype=col_character(),
                                    msname=col_character(),
                                    .default= col_double()),
                   n_max=-1)
names(pufsyn) %>% sort
count(pufsyn, ftype)

pufsyn %>%
  group_by(ftype) %>%
  summarise(rnmin=min(rownum), rnmax=max(rownum))

# conform variable names, then group and check for dups
dropvars <- c("divratio", "e01500_minus_e01700", "m", "msname")
stack <- pufsyn %>%
  dplyr::select(-dropvars) %>%
  mutate(rownum=row_number())
summary(stack) # no missing values
names(stack) %>% sort

# tmp <- stack %>%
#   filter(e00300=="14930000", MARS==4) %>%
#   arrange(ftype)
# tmp %>% dplyr::select(ftype, pufseqn, rownum, MARS, e00300)
# glimpse(tmp)
# 
# tmp %>%
#   unite(ft_rn, ftype, rownum) %>%
#   gather(vname, value, -ft_rn) %>%
#   spread(ft_rn, value)

a <- proc.time()
(grouping_vars <- setdiff(names(stack), c("ftype", "wt", "pufseqn", "rownum")))
mcounts <- stack %>%
  # filter(MARS==1) %>%
  group_by_at(vars(one_of(grouping_vars))) %>%
  summarise(n=n(), 
            npuf=sum(ftype=="puf"), 
            nsyn=sum(ftype=="syn"),
            wt.pufmin=min(ifelse(ftype=="puf", wt, 1e9)),
            wt.pufmdn=median(ifelse(ftype=="puf", wt, NA), na.rm=TRUE),
            wt.pufmax=max(ifelse(ftype=="puf", wt, -1e9)),
            rnmin=min(rownum),
            rnmax=max(rownum)) %>%
  ungroup
b <- proc.time()
b - a # 10 minutes
saveRDS(mcounts, paste0(globals$tc.dir, synname, "_mcounts.rds"))

synname <- "synthpop10"
mcounts <- readRDS(paste0(globals$tc.dir, synname, "_mcounts.rds"))

glimpse(mcounts)
names(mcounts) %>% sort
sum(mcounts$n)
sum(mcounts$npuf)
sum(mcounts$nsyn)

dups <- mcounts %>%
  filter(n>1, npuf>=1, nsyn>=1) %>%
  dplyr::select(n, npuf, nsyn, wt.pufmin, wt.pufmdn, wt.pufmax, rnmin, rnmax, MARS, everything())
glimpse(dups)

# find number of nonzero columns
nvals <- apply(as.matrix(dups[, -c(1:8)]), 1, function(x) sum(x != 0))
dups2 <- dups %>%
  ungroup %>%
  mutate(nnz=nvals) %>%
  arrange(-nnz, -n) %>%
  dplyr::select(n, nnz, everything()) %>%
  mutate(rownum=row_number())

# summarize nonzeros ----
dups2 %>% 
  group_by(nnz) %>%
  summarise(npairs=n(), 
            npufrecs=sum(npuf),
            nsynrecs=sum(nsyn),
            min_pufwt=min(wt.pufmin),
            wt.pufmdn=median(wt.pufmdn),
            max_pufwt=max(wt.pufmax)) %>%
  janitor::adorn_totals() %>%
  kable(digits=0, format.arg=list(big.mark=","))


dlong <- dups2 %>%
  gather(vname, value, -c(rownum, nnz, n, npuf, nsyn, wt.pufmin, wt.pufmdn, wt.pufmax, rnmin, rnmax)) %>%
  filter(value!=0) %>%
  dplyr::select(rownum, everything()) %>%
  left_join(puf.vnames %>% dplyr::select(vname, vdesc))
vnames <- unique(dlong$vname)
vfirst <- c("MARS", "XTOT", "FDED", "MIDR", "n24", "EIC", "f6251")
vorder <- c(vfirst, setdiff(vnames, vfirst) %>% sort)           
vorder

dlong %>% filter(rownum==2)
dlong %>% filter(nnz==8) %>% arrange(rownum)
dlong %>% filter(nnz==7) %>% arrange(rownum)


dlong %>%
  mutate(vname=factor(vname, levels=vorder)) %>%
  filter(nnz==3) %>%
  # filter(rownum==294) %>%
  #filter(wt.pufmin <=20) %>%
  # filter(!vname %in% c("MARS", "XTOT", "FDED", "e00200"))
  arrange(rownum, vname) %>%
  kable(digits=0, format.args = list(big.mark=","))

# which variables?
dlong %>%
  group_by(nnz, vname, vdesc) %>%
  summarise(n=n()) %>%
  spread(nnz, n) %>%  
  mutate_at(vars(-vname, -vdesc), ~naz(.)) %>%
  adorn_totals("col") %>%
  arrange(-Total)

count(dlong, vname)
vfirst2 <- c("MARS", "FDED", "XTOT", "n24", "EIC", "MIDR", "DSI", "f2441", "f6251")
vorder <- c(vfirst2, setdiff(vname, vfirst2) %>% sort)           
vorder

# a <- c("v1", "XTOT", "MARS", "v2", "FDED", "abc")
# intersect(vfirst2, a)
f <- function(vname){
  vfirst2 <- c("MARS", "FDED", "XTOT", "n24", "EIC", "MIDR", "DSI", "f2441", "f6251")
  vfirst2 <- intersect(vfirst2, vname)
  vorder <- c(vfirst2, setdiff(vname, vfirst2) %>% sort)           
  paste(vorder, collapse="-")
}
tmp <- dlong %>%
  group_by(rownum) %>%
  summarise(nnz=max(nnz), nzmatchvars=f(vname))

tmp2 <- dlong %>%
  filter(!vname %in% vfirst2) %>%
  group_by(rownum) %>%
  summarise(nnz=max(nnz), nzmatchvars=f(vname))

ht(dlong)

catvars <- c("MARS", "FDED", "XTOT", "n24", "EIC", "MIDR", "DSI", "f2441", "f6251")
tmp2 <- dlong %>%
  mutate(contvar=ifelse(vname %in% catvars, 0, 1)) %>%
  group_by(nnz, rownum) %>%
  summarise(ncontvars=sum(contvar)) %>%
  group_by(nnz, ncontvars) %>%
  summarise(ncombos=n()) %>%
  spread(ncontvars, ncombos) %>%
  mutate_at(vars(-nnz), ~naz(.)) %>%
  adorn_totals() %>%
  adorn_totals(where="col")
tmp2 %>%
  kable(digits=0, format.args=list(big.mark=","))

tmp3 <- dlong %>%
  mutate(contvar=ifelse(vname %in% catvars, 0, 1)) %>%
  filter(contvar==1) %>%
  group_by(nnz, vname, vdesc) %>%
  summarise(n=n()) %>%
  group_by(nnz) %>%
  arrange(desc(n)) %>%
  spread(nnz, n) %>%
  mutate_at(vars(-c(1:2)), ~naz(.)) %>%
  adorn_totals("col") %>%
  adorn_totals() %>%
  arrange(-Total)
tmp3 %>%
  kable(digits=0, format.args=list(big.mark=","), caption="Counts of nonzero continuous variables in matched pairs with nnz (column)")


# wages still involved in some matches
# unemployment comp involved in some matches
# most other are trivial


# interlude -- find nonmatches ----
goodrecs <- mcounts %>%
  filter(n>1, npuf==0, nsyn>=1) %>%
  select(-starts_with("wt")) %>%
  select(rnmin, rnmax, n, npuf, nsyn, MARS, everything())

# fine a few records for max
maxrecs <- goodrecs %>% 
  filter(n %in% c(10, 62))
names(maxrecs) %>% sort

maxrecs.syn <-
  maxrecs %>%
  select(-starts_with("wt")) %>%
  left_join(stack)
names(maxrecs.syn) %>% sort

nonmatches <- maxrecs.syn %>%
  select(-rnmin, -rnmax, -pufseqn) %>%
  select(ftype, rownum, n, npuf, nsyn, wt, MARS, everything()) %>%
  arrange(n, rownum)

# end interlude ----



