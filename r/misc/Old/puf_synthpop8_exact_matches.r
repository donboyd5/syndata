

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


#****************************************************************************************************
#               Max's question ####
#****************************************************************************************************

# dir <- "D:/Google Drive/synpuf/syntheses/"
# s8bad <- read_csv(paste0(dir, "synthpop8.csv"),
#                col_types = cols(ftype=col_character(), .default= col_double()))
# names(s8bad) %>% sort
# count(s8bad, ftype)
# ftype      n
# <dbl>  <int>
#   1    NA 163786
# oops, synthpop8 should have 5x the number of records as puf

# try again -- this looks good
s8 <- read_csv(paste0(globals$synd, "synthpop8_stack.csv"), 
                    col_types = cols(ftype=col_character(),
                                     msname=col_character(),
                                     .default= col_double()),
                    n_max=-1)
names(s8) %>% sort
count(s8, ftype)

# conform variable names, then group and check for dups
dropvars <- c("divratio", "e01500_minus_e01700", "m", "msname")
stack <- s8 %>%
  select(-dropvars) %>%
  mutate(rownum=row_number())
summary(stack) # no missing values
names(stack) %>% sort

(grouping_vars <- setdiff(names(stack), c("ftype", "wt", "pufseqn", "rownum")))
mcounts <- stack %>%
  filter(MARS==1) %>%
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

glimpse(mcounts)
names(mcounts) %>% sort
sum(mcounts$n)
sum(mcounts$npuf)
sum(mcounts$nsyn)

dups <- mcounts %>%
  filter(n>1, npuf>=1, nsyn>=1) %>%
  select(n, npuf, nsyn, wt.pufmin, wt.pufmdn, wt.pufmax, MARS, everything())

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

write_csv(nonmatches, paste0(globals$synd, "synthpop8_selected_nonmatches.csv"))

rn <- 974622
stack %>%
  filter(rownum==rn) %>%
  select(ftype, rownum, wt)
# end interlude ----



#****************************************************************************************************
#               unworkable versions of synthpop8 ####
#****************************************************************************************************
# dir2 <- "D:/Google Drive/synpuf/"
# s8 <- read_csv(paste0(dir2, "synthpop8_wfs_stack.csv"), # this is syn and puf but ftype has all NA!
#                col_types = cols(.default= col_double()))
# names(s8) %>% sort
# count(s8, ftype)
# 
# s8 <- readRDS(paste0(globals$tc.dir, "synthpop8_wfs", ".rds")) # this is syn only
# names(s8) %>% sort
# count(s8, ftype)
# 
# s8 <- readRDS(paste0(globals$tc.dir, "synthpop8_rwprep.rds")) # this is a list with tc.base and tc.output
# names(s8) %>% sort
# count(s8, ftype)




#****************************************************************************************************
#               synthpuf8 ####
#****************************************************************************************************
# vars <- c("MARS", "wt", "e00100", "e04600", "p04470", "e04800", "e62100", "e05800", "e08800", "e59560", "e26190")

# stack < read_csv(paste0(globals$synd, "synthpop8_stack.csv"),
#                   col_types = cols(.default= col_double()))
# names(stack) %>% sort
# s8 <- readRDS(paste0(globals$tc.dir, "synthpop8_wfs_stack.rds"))

s8 <- readRDS(paste0(globals$tc.dir, "synthpop8_rwprep.rds"))
# names(s8) %>% sort
names(s8$tc.base) %>% sort
names(s8$tc.output) %>% sort
str_subset(names(s8$tc.output), "26190")
# keep c00100 c04600 c04470 c04800 c62100 c05800 THERE IS NO: 08800, 59560

s8a <- s8$tc.base %>%
  left_join(s8$tc.output %>% select(RECID, c00100, c04600, c04470, c04800, c62100, c05800))
s8b <- s8a %>%
  select(-c(divratio, e01500_minus_e01700, m, msname, pufseqn, RECID))
names(s8b)

puf <- s8b %>% filter(ftype=="puf")
syn <- s8b %>% filter(ftype=="syn")

pufu <- puf %>% unique
synu <- syn %>% unique

# find the number of exact matches between pufu and synu
# potential blocks MARS c00100 c04600 c04470 c04800 c62100 c05800 THERE IS NO: 08800, 59560
cseeds <- c("c00100", "c04600", "c04470", "c04800", "c62100", "c05800")
count(s8b, ftype, MARS) %>% spread(ftype, n) %>% arrange(-syn) # 4
count(s8b, ftype, wt) %>% spread(ftype, n) %>% arrange(-syn) # 92
count(s8b, ftype, c00100) %>% spread(ftype, n) %>% arrange(-syn) # 420k rows
count(s8b, ftype, c04600) %>% spread(ftype, n) %>% arrange(-syn) # 29k rows
# c04600   puf    syn
# <dbl> <int>  <int>
#   1     0  44732 222697
# 2  3900  43592 218447
# 3  7800  35680 179302
# 4 11700  14599  73253
# 5 15600  12456  62467
# 6 19500   6518  32732
# 7  2258.     1      7


names(stack2)
tmp <- s8$tc.base %>% 
  filter(e01700==32900, e01500==32900, e09900==810, MARS==3, FDED==2, XTOT==1) %>%
  select(ftype, pufseqn, everything())
names(tmp) %>% sort


# look for exact matches on several variables
# puf3 <- puf %>% filter(MARS==3) %>% select(-cseeds, -ends_with("p"), -ends_with("s"), -wt)
puf3 <- puf %>% filter(MARS==3) %>% select(-cseeds, -ends_with("p"), -ends_with("s"))
syn3 <- syn %>% filter(MARS==3) %>% select(-cseeds, -ends_with("p"), -ends_with("s"))
names(puf3) %>% sort

# look for exact matches
stack2 <- bind_rows(puf3, syn3)
mcounts <- stack2 %>%
  mutate(MARS=3) %>%
  group_by_at(vars(one_of(setdiff(names(.), c("ftype", "wt"))))) %>%
  summarise(n=n(), npuf=sum(ftype=="puf"), nsyn=sum(ftype=="syn"),
            wt.pufmin=min(ifelse(ftype=="puf", wt, 1e9)),
            wt.pufmdn=median(ifelse(ftype=="puf", wt, NA), na.rm=TRUE),
            wt.pufmax=max(ifelse(ftype=="puf", wt, -1e9)))
glimpse(mcounts)
names(mcounts) %>% sort
sum(mcounts$n)
sum(mcounts$npuf)
sum(mcounts$nsyn)

dups <- mcounts %>%
  filter(n>1, npuf>=1, nsyn>=1) %>%
  select(n, npuf, nsyn, wt.pufmin, wt.pufmdn, wt.pufmax, MARS, everything())



# find number of nonzero columns
nvals <- apply(as.matrix(dups[, -c(1:6)]), 1, function(x) sum(x != 0))
dups2 <- dups %>%
  ungroup %>%
  mutate(nnz=nvals) %>%
  arrange(-nnz, -n) %>%
  select(n, nnz, everything()) %>%
  mutate(rownum=row_number())

write_csv(dups2, "./results/dups2.csv")

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
  kable(digits=0)


dlong <- dups2 %>%
  gather(vname, value, -c(rownum, nnz, n, npuf, nsyn, wt.pufmin, wt.pufmdn, wt.pufmax)) %>%
  filter(value!=0) %>%
  select(rownum, everything()) %>%
  left_join(puf.vnames %>% select(vname, vdesc))
vnames <- unique(dlong$vname)
vfirst <- c("MARS", "XTOT", "FDED", "MIDR", "n24", "f6251")
vorder <- c(vfirst, setdiff(vnames, vfirst) %>% sort)           
vorder

dlong %>% filter(rownum==2)


dlong %>%
  mutate(vname=factor(vname, levels=vorder)) %>%
  filter(nnz==3) %>%
  # filter(rownum==294) %>%
  #filter(wt.pufmin <=20) %>%
  # filter(!vname %in% c("MARS", "XTOT", "FDED", "e00200"))
  arrange(rownum, vname) %>%
  kable(digits=0, format.args = list(big.mark=","))

  # 157 294

# which variables?
dlong %>%
  group_by(nnz, vname, vdesc) %>%
  summarise(n=n()) %>%
  spread(nnz, n) %>%
  mutate_at(vars(-vname, -vdesc), ~naz(.)) %>%
  mutate(tot=`3` + `4` + `5` + `6`) %>%
  arrange(-tot)



