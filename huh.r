

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

vars <- c("MARS", "wt", "e00100", "e04600", "p04470", "e04800", "e62100", "e05800", "e08800", "e59560", "e26190")
puf.vnames %>% filter(vname %in% vars)


df <- get_puf.base()

df3 <- df %>%
  select(vars) %>%
  gather(vname, value) %>%
  group_by(vname) %>%
  mutate(nunique=length(unique(value))) %>%
  group_by(vname, nunique, value) %>%
  summarise(value_count=n())
glimpse(df3)
ht(df3)

df3 %>%
  group_by(vname) %>%
  arrange(desc(value_count)) %>%
  top_n(5) %>%
  ungroup %>%
  left_join(puf.vnames %>% select(vname, vdesc)) %>%
  arrange(nunique, desc(value_count)) %>%
  select(vname, nunique, value_count, value, vdesc)

drops <- c("wt", "p04470", "e04800")
df2 <- df %>%
  group_by_at(vars(one_of(setdiff(vars, drops)))) %>%
  summarise(n=n()) %>%
  ungroup
nrow(df2) / nrow(df)
ht(df2)


pufdups <- left_join(df2 %>% mutate(indf2=1),
                     df) %>%
  filter(indf2==1)

count(df2, n) %>%
  mutate(n.cum=cumsum(nn), n.pct=n.cum / sum(nn) * 100)
  

ht(df2)
quantile(df2$n, probs=c(0, .25, .5, .75, seq(.9, 1, .01)))
df3 <- df2 %>%
  arrange(-n)

df3 <- df %>%
  group_by(MARS, e00100) %>%
  summarise(n=n())
nrow(df3) / nrow(df)

lsyn <- readRDS(paste0("d:/tcdir/", "synthpop8_rwprep.rds"))
names(lsyn)
v2 <- names(lsyn$tc.output) %>% sort

# find variables that are good candidates for density calcs
synvars <- setdiff(names(lsyn$tc.base)[12:69], c("e01500_minus_e01700", "divratio"))
synvars

uvals <- lsyn$tc.base %>%
  filter(ftype=="puf") %>%
  summarise_at(vars(synvars), ~length(unique(.)))

uvals %>%
  gather(vname, count) %>%
  arrange(count)

count(lsyn$tc.base, e09700)


setdiff(vars, v2)
# "e08800", "e59560", 
vcalc <- c("c00100", "c04600", "c04470", "c04800", "c62100", "c05800")

f2 <- lsyn$tc.base %>%
  left_join(lsyn$tc.output %>% select(RECID, vcalc))
names(f2) %>% sort

f2a <- f2 %>%
  filter(!(ftype=="syn" & m>1))

f3 <- f2a %>%
  group_by_at(vars(one_of(c("MARS", vcalc, "ftype")))) %>%
  summarise(n=n()) %>%
  ungroup
ht(f3)
count(f3, n)

fwide <- f3 %>%
  filter(ftype=="puf") %>%
  rename(ftype.puf=ftype, n.puf=n) %>%
  left_join(f3 %>% filter(ftype=="syn") %>% rename(ftype.syn=ftype, n.syn=n))
ht(fwide)

dups <- fwide %>%
  filter(!is.na(n.syn))
ht(dups)

# which NON-seed vars have the fewest unique vals and 
names(dups)
# get puf recs that correspond to these values


dups %>%
  summarise_at(vars(vcalc), ~sum(.==0))

count(dups, n.puf, n.syn)
count(dups, n.puf) %>% mutate(n.cum=cumsum(n))

# ok, get all the records in f2a that have the same vars as dups with n.=1
dups2 <- f2a %>%
  left_join(dups %>% filter(n.puf==1, n.syn==1) %>% select(-ftype.puf, -ftype.syn)) %>%
  filter(n.puf==1) %>%
  arrange(MARS, c00100, c04600, c04470, c04800, c62100, c05800, ftype) %>%
  select(ftype, everything())
count(dups2, ftype)

start <- 81
dups2 %>% 
  dplyr::select(-vcalc, -m, -msname, -pufseqn, -RECID, -divratio, -ends_with("p"), -ends_with("s"), -n.puf, -n.syn) %>%
  slice(start:(start + 9))


#****************************************************************************************************
#               synthpuf8 distances ####
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
  mutate(MARS=1) %>%
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



