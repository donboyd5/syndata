

#****************************************************************************************************
#                Includes ####
#****************************************************************************************************
source(file.path(PROJHOME, "r/includes", "libraries.r"))
# library("readxl")
library("knitr")

# search() %>% sort # how did precis get loaded?

source(file.path(PROJHOME, "r/includes", "globals_system_specific_boyd.r"))
source(file.path(PROJHOME, "r/includes", "globals_other.r"))
source(file.path(PROJHOME, "r/includes", "functions_general.r"))
# source("./r/includes/functions_synthesis.r")

search() %>% sort

ns <- function(df) {names(df) %>% sort}

#****************************************************************************************************
#                Get data ####
#****************************************************************************************************
puf.vnames <- get_puf_vnames()

puf <- read_csv(paste0("D:/Google Drive/synpuf/puf/", "puf2011.csv"), 
                col_types = cols(.default= col_double()), 
                n_max=-1) %>%
  mutate(ftype="puf", id=row_number() - 1)
glimpse(puf)
puf <- puf %>% filter(!RECID %in% 999996:999999)

rf <- read_csv(paste0("D:/Google Drive/synpuf/syntheses/", "synpuf20.csv"), 
               col_types = cols(.default= col_double()), 
               n_max=-1) %>%
  mutate(ftype="rf", id=row_number() - 1)
glimpse(rf)
anyDuplicated(rf$RECID)
count(rf, ftype)


cart <- read_csv(paste0("D:/tcdir/", "synthpop10.csv"), 
                  col_types = cols(.default= col_double()), 
                  n_max=-1) %>%
  mutate(ftype="cart", id=row_number() - 1) %>%
  select(-ROWNUM) # ROWNUM is not in rf
count(cart, ftype)

ns(cart)
sum(cart$E00200)


rfnearest <- read_csv(paste0("D:/Google Drive/synpuf/nearest/", "nearest20.csv"), 
                      col_types = cols(.default= col_double()), 
                      n_max=-1) %>% 
  select(rfid=id_A, pufid=id_B1, dist1)

cartnearest <- read_csv(paste0("D:/Google Drive/synpuf/nearest/", "nearest19.csv"), 
                      col_types = cols(.default= col_double()), 
                      n_max=-1) %>% 
  select(cartid=id_A, pufid=id_B1, dist1)

# quick check on puf-cart dups of a certain kind
dcheck <- bind_rows(puf %>% select(names(cart)),
                    cart,
                    rf %>% select(names(cart))) %>%
  filter(MARS==1, XTOT==1, FDED==2) %>%
  select(-MARS, -XTOT, -FDED) %>%
  gather(vname, value, -id, -ftype) %>%
  group_by(ftype, id) %>%
  summarise(n=n(), nnz=sum(value==0)) %>%
  filter(n==nnz)
dcheck %>% group_by(ftype) %>% summarise(n=n())
# ns(dcheck)


#****************************************************************************************************
#                Compare cart matches vs nearest ####
#****************************************************************************************************
glimpse(puf)
glimpse(cart)
glimpse(cartnearest)

ns(puf)
ns(cart)

# let's do this like they did for sipp - for each true rec, get the corresponding syn rec
# if a puf rec is matched to 2 or more syn recs, get just the first syn rec
# cartnear2 <- cartnearest %>%
#   filter(dist1==0)
# glimpse(cartnear2) # 72,674 records match exactly
# anyDuplicated(cartnear2) # no duplicates
# length(unique(cartnear2$pufid)) # 10923 unique
# length(unique(cartnear2$cartid)) # 72674 unique
# # cartnear2a <- cartnear2 %>% distinct() # same

# now do it my way - what are exact matches - stack files
ns(cart)
cartxvars <- c("ftype", "id") # cart vars not to match on
cartmvars <- setdiff(names(cart), cartxvars) # 64 vars
cartmvars %>% sort # good, 64
setdiff(names(cart), names(puf)) # good, all syn vars are in the puf
cartstack <- bind_rows(puf %>% select(names(cart)), cart)
ns(cartstack)

a <- proc.time()
cartdups_marked <- cartstack %>%
  mutate(dup1=duplicated(select(., -ftype, -id)),
         dup2=duplicated(select(., -ftype, -id), fromLast=TRUE),
         dup=(dup1==TRUE) | (dup2==TRUE))
b <- proc.time()
b - a # 99 secs

a <- proc.time()
cartdups <- cartdups_marked %>%
  filter(dup==TRUE) %>%
  select(-dup, -dup1, -dup2) %>%
  mutate(group_id = seplyr::group_indices_se(., cartmvars)) %>%
  arrange(group_id, ftype) %>%
  group_by(group_id) %>%
  mutate(ndups=n(), npuf=sum(ftype=="puf"), ncart=sum(ftype=="cart")) %>%
  ungroup %>%
  filter(npuf>0, ncart>0)
b <- proc.time()
b - a # 16 secs
glimpse(cartdups) # 18521 records -- 6,041 cart, 12,280 puf
count(cartdups, ftype) # 6,041 cart records and 12,280 puf records where a combo of values is duplicated
length(unique(cartdups$group_id)) # 2,740 combos of values that are in both files
count(cartdups, group_id)
cartdups %>% filter(group_id==7) %>% select(1:10, E01200:ncart)

# hmmm. cartnear2 (max) suggests 72,674 records match, this suggests 18,521 match
cartids <- cartdups %>% filter(ftype=="cart") %>% select(group_id, cartid=id)
pufids <- cartdups %>% filter(ftype=="puf") %>% select(group_id, pufid=id)
cartpairs <- full_join(pufids, cartids) # 72,715 pairs
ht(cartpairs)
length(unique(cartpairs$pufid)) # 12,480 unique puf records match, Max: 34,875 unique
12480 / 163786 # = 7.6% of puf recs found in syn
length(unique(cartpairs$cartid)) # 6,041 unique cart records match, Max:  10293 unique
6041 / 818930 # 0.7% of CART records found in PUF

cart_distinct <- cartdups %>%
  select(-id) %>%
  group_by(ftype) %>%
  distinct # these are the 2,740 unique combinations of values
count(cart_distinct, ftype)
max(cart_distinct$npuf)
max(cart_distinct$ncart)
cart_distinct %>%
  filter(ncart==max(ncart)) %>%
  as.data.frame

# count continuous variables
# DSI EIC FDED F2441 F6251 MARS MIDR N24 XTOT
catvars <- c("MARS", "XTOT", "EIC", "F6251", "F2441", "DSI", "MIDR", "N24", "FDED")
contvars <- setdiff(names(cart_distinct), c(catvars, c("ftype", "group_id", "ndups", "npuf", "ncart")))
contvars
cart_contvars <- cart_distinct %>%
  ungroup %>%
  filter(ftype=="puf") %>%
  select(c(group_id, npuf, ncart, contvars)) %>%
  gather(variable, value, contvars) %>%
  filter(value != 0) %>%
  group_by(group_id, npuf, ncart) %>%
  arrange(variable) %>%
  summarise(nnz=n(), matchvars=paste(variable, collapse="-")) %>%
  ungroup
count(cart_contvars, nnz)

glimpse(cart_contvars)
cart_contvars %>%
  group_by(nnz, matchvars) %>%
  summarise(nm=n()) %>%
  arrange(nnz, -nm)

cart_contvars %>%
  group_by(nnz) %>%
  summarise(npuf=sum(npuf))


# glimpse(rfdups)
# count(rfdups, ftype) # 134203 rf, 49807 puf
# rfdups %>%
#   group_by(ftype) %>%
#   summarise(n=n(), nunique=length(unique(id))) # uniques: 134203 rf, 49807 puf
# # I find more unique puf records than Max because he searched from rf to puf; he did not
# # search from every puf - so some rf records could have matched the same puf
# # save the rfdups so that I can create a lowmatch file
# write_csv(rfdups %>% 
#             select(ftype, RECID, id, group_id, ndups, npuf, nrf),
#           paste0(globals$tc.dir, "rfdups.csv"))
# 
# # what are the group sizes?
# rfdups %>% select(group_id, ndups, npuf, nrf) %>% distinct() %>% arrange(-ndups)
# count(rfdups, ndups) %>% arrange(-n)
# 
# # do we find every one of Max's puf-syn matches in my data?
# rfids1 <- rfdups %>%
#   select(group_id, ftype, id, ndups, npuf, nrf)
# 
# # write some groups
# # group_id 318 has 37 puf, 306 rf records
# # 15848 has 257 dups
# gid <- 15848
# gid <- 6892 # 22 variables match
# rfdups %>% 
#   filter(group_id==gid) %>% 
#   select(RECID, ftype, id, group_id, ndups, npuf, nrf, everything()) %>%
#   write_csv(paste0(globals$tc.dir, "group", gid, ".csv"))
# 
# rfids1 %>% filter(row_number() %in% c(373, 374))
# rfids1 %>% filter(group_id==133)
# d <- puf %>% filter(id==153879)
# 
# glimpse(rfdups)
# rfdups_unique <- rfdups %>%
#   select(-xvars) %>%
#   group_by(group_id) %>%
#   distinct() %>%
#   ungroup # 34875 unique combinations
# 
# rfdups_unique %>%
#   summarise(n=n(), ndups=sum(ndups), npuf=sum(npuf), nrf=sum(nrf))
# 





#****************************************************************************************************
#                Compare RF matches vs nearest ####
#****************************************************************************************************
glimpse(puf)
glimpse(rf)
glimpse(rfnearest)

names(puf) %>% sort
ns(rf)

# let's do this like they did for sipp - for each true rec, get the corresponding syn rec
# if a puf rec is matched to 2 or more syn recs, get just the first syn rec
rfnear2 <- rfnearest %>%
  filter(dist1==0)
glimpse(near2) # 134203 records match
anyDuplicated(rfnear2) # no duplicates
length(unique(rfnear2$pufid)) # 34875 unique
length(unique(rfnear2$rfid)) # 134203 unique
# near2a <- near2 %>% distinct()

# puf records that have exact match in syn
# puf2 <- near2 %>%
#   select(pufid, synid) %>%
#   left_join(puf %>% select(pufid=id, vars))

# now do it my way - what are exact matches - stack files
xvars <- c("S006", "ftype", "RECID", "id")
mvars <- setdiff(names(rf), xvars) # 64 vars
mvars %>% sort
setdiff(names(rf), names(puf)) # good, all syn vars are in the puf
rfstack <- bind_rows(puf %>% select(names(rf)), rf)

a <- proc.time()
rfdups_marked <- rfstack %>%
  mutate(dup1=duplicated(select(., -ftype, -RECID, -id, -S006)),
         dup2=duplicated(select(., -ftype, -RECID, -id, -S006), fromLast=TRUE),
         dup=(dup1==TRUE) | (dup2==TRUE))
b <- proc.time()
b - a # 90 secs

a <- proc.time()
rfdups <- rfdups_marked %>%
  filter(dup==TRUE) %>%
  select(-dup, -dup1, -dup2) %>%
  mutate(group_id = seplyr::group_indices_se(., mvars)) %>%
  arrange(group_id, ftype) %>%
  group_by(group_id) %>%
  mutate(ndups=n(), npuf=sum(ftype=="puf"), nrf=sum(ftype=="rf")) %>%
  ungroup %>%
  filter(npuf>0, nrf>0)
b <- proc.time()
b - a # 16 secs

glimpse(rfdups)
count(rfdups, ftype) # 134203 rf, 49807 puf
rfdups %>%
  group_by(ftype) %>%
  summarise(n=n(), nunique=length(unique(id))) # uniques: 134203 rf, 49807 puf
# I find more unique puf records than Max because he searched from rf to puf; he did not
# search from every puf - so some rf records could have matched the same puf
# save the rfdups so that I can create a lowmatch file
write_csv(rfdups %>% 
            select(ftype, RECID, id, group_id, ndups, npuf, nrf),
          paste0(globals$tc.dir, "rfdups.csv"))

# what are the group sizes?
rfdups %>% select(group_id, ndups, npuf, nrf) %>% distinct() %>% arrange(-ndups)
count(rfdups, ndups) %>% arrange(-n)

# do we find every one of Max's puf-syn matches in my data?
rfids1 <- rfdups %>%
  select(group_id, ftype, id, ndups, npuf, nrf)

# write some groups
# group_id 318 has 37 puf, 306 rf records
# 15848 has 257 dups
gid <- 15848
gid <- 6892 # 22 variables match
rfdups %>% 
  filter(group_id==gid) %>% 
  select(RECID, ftype, id, group_id, ndups, npuf, nrf, everything()) %>%
  write_csv(paste0(globals$tc.dir, "group", gid, ".csv"))

rfids1 %>% filter(row_number() %in% c(373, 374))
rfids1 %>% filter(group_id==133)
d <- puf %>% filter(id==153879)

glimpse(rfdups)
rfdups_unique <- rfdups %>%
  select(-xvars) %>%
  group_by(group_id) %>%
  distinct() %>%
  ungroup # 34875 unique combinations

rfdups_unique %>%
  summarise(n=n(), ndups=sum(ndups), npuf=sum(npuf), nrf=sum(nrf))


#****************************************************************************************************
#                Analyze the unique rf combinations ####
#****************************************************************************************************
# first, what unique categorical combos do we have in the puf?
rfcatvars <- c("DSI", "EIC", "F2441", "F6251", "FDED", "MARS", "MIDR", "N24", "XTOT")
catgroups <- puf2 %>% 
  select(rfcatvars) %>%
  group_by_at(vars(one_of(rfcatvars))) %>%
  summarise(n=n()) %>%
  ungroup
ht(catgroups)
catgroups %>% arrange(n)
count(catgroups, n) # 157 combos of categoricals have only 1 record! 2 have 67
sum(catgroups$n)
catgroups %>% filter(n==1) %>% arrange(MARS, XTOT)

# now let's look at the continuous variables
rfcontvars <- setdiff(mvars, rfcatvars)
glimpse(rfdups_unique)

# which continuous vars are nz?
rfdu_long <- rfdups_unique %>%
  select(group_id, ndups, npuf, nrf, rfcontvars) %>%
  gather(vname, value, rfcontvars) %>%
  filter(value!=0)
glimpse(rfdu_long)
count(rfdu_long, vname) %>%
  left_join(puf.vnames %>% mutate(vname=str_to_upper(vname)) %>% select(vname, vdesc)) %>%
  arrange(-n)

matchtypes <- rfdu_long %>%
  arrange(group_id, vname) %>%
  group_by(group_id, ndups, npuf, nrf) %>%
  summarise(nnz=n(), nzmatches=paste0(vname, collapse="-")) %>%
  ungroup
glimpse(matchtypes)
count(matchtypes, nnz)
write_csv(matchtypes, paste0(globals$tc.dir, "rfmatchtypes.csv"))

matchtypes %>%
  group_by(nnz) %>%
  summarise(npuf=sum(npuf), nrf=sum(nrf)) %>%
  arrange(-nnz) %>%
  mutate(cnrf=cumsum(nrf))
# if we wanted to eliminate matches 
# row  nnz  npuf   nrf   cnrf
# <int> <int> <int>  <int>
# 1    22     1     1      1
# 2    21     1     1      2
# 3    20     2     2      4
# 4    19     1     1      5
# 5    18     2     2      7
# 6    17     2     2      9
# 7    16     9     9     18
# 8    15    12    12     30
# 9    14    18    18     48
# 10    13    31    31     79
# 11    12    45    46    125
# 12    11    72    72    197
# 13    10   132   135    332
# 14     9   214   219    551
# 15     8   348   363    914
# 16     7   557   613   1527
# 17     6   923  1054   2581
# 18     5  1680  2114   4695
# 19     4  3021  4319   9014
# 20     3  5009  8433  17447
# 21     2 12099 28598  46045
# 22     1 25461 87094 133139

matchtypes %>%
  group_by(nnz) %>%
  summarise(npuf=sum(npuf), nrf=sum(nrf)) %>%
  arrange(-nnz) %>%
  mutate(cnrf=cumsum(nrf))

write_csv(rfdups %>% 
            select(ftype, RECID, id, group_id, ndups, npuf, nrf),
          paste0(globals$tc.dir, "rfdups.csv"))

matchtypes %>%
  filter(nnz < 5) %>%
  summarise(npuf_sum=sum(npuf), nrf=sum(nrf), npuf_min=min(npuf))

matchtypes %>%
  group_by(npuf) %>%
  summarise(nrf=sum(nrf)) %>%
  arrange(-nnz) %>%
  mutate(cnrf=cumsum(nrf))


#****************************************************************************************************
#                Create and save a low-match synpuf20 ####
#****************************************************************************************************
rf_raw <- read_csv(paste0("D:/Google Drive/synpuf/syntheses/", "synpuf20.csv"), 
               col_types = cols(.default= col_double()), 
               n_max=-1)

rfmatchtypes <- read_csv(paste0(globals$tc.dir, "rfmatchtypes.csv"))
rfdups <- read_csv(paste0(globals$tc.dir, "rfdups.csv"))

glimpse(rf_raw)
glimpse(rfdups)
glimpse(rfmatchtypes)

rfdups2 <- rfdups %>%
  filter(ftype=="rf") %>%
  select(RECID, group_id, npuf, nrf) %>%
  left_join(rfmatchtypes %>% select(group_id, nnz, nzmatches))

rfxrecs <- rfdups2 %>%
  filter(nnz >= 5) # ge 5 4,695; ge 4 9,014

rf_lowmatch <- rf_raw %>%
  filter(!RECID %in% rfxrecs$RECID)
glimpse(rf_lowmatch)

# write the low-match file
write_csv(rf_lowmatch,
          paste0("D:/Google Drive/synpuf/syntheses/", "synpuf20_lowmatch.csv"))

rfdups %>% filter(group_id==2)
d <- rf_raw %>% filter(RECID==668453)




#****************************************************************************************************
#                Could the problem be the data I am using? NO. ####
#****************************************************************************************************
pufsyn <- readRDS(paste0(globals$tc.dir, "synth10syn20_stack.rds"))
glimpse(pufsyn)
# count(pufsyn, ftype)
# names(pufsyn) %>% sort
# pufsyn %>%
#   group_by(ftype) %>%
#   summarise(RECID=anyDuplicated(RECID), recna=sum(is.na(RECID)))

# define matching vars
exclvars <- c(c("ftype", "c00100", "divratio", "group", "m", "mgroup", "msname",
                "pufseqn", "RECID", "rownum", "taxbc", "wfsrecid"), 
              str_subset(names(pufsyn), "minus"),
              str_subset(names(pufsyn), "0p"),
              str_subset(names(pufsyn), "0s"),
              str_subset(names(pufsyn), "wt"))
exclvars
(matchvars <- setdiff(names(pufsyn), exclvars))
# how do my match vars compare to Max's 
matchvars %>% str_to_upper() %>% sort
mvars %>% str_to_upper() %>% sort
identical(matchvars %>% str_to_upper() %>% sort,
          mvars %>% str_to_upper() %>% sort)
# perfect

djb <- pufsyn %>%
  filter(ftype %in% c("puf", "rf")) %>%
  select(ftype, RECID, matchvars)
glimpse(djb)

# now do it my way - what are exact matches - stack files
xvarsdjb <- c("ftype", "RECID")
mvarsdjb <- setdiff(names(djb), xvarsdjb) # 64 vars
mvarsdjb %>% sort

a <- proc.time()
djb_marked <- djb %>%
  mutate(dup1=duplicated(select(., -ftype, -RECID)),
         dup2=duplicated(select(., -ftype, -RECID), fromLast=TRUE),
         dup=(dup1==TRUE) | (dup2==TRUE))
b <- proc.time()
b - a # 90 secs

a <- proc.time()
djbdups <- djb_marked %>%
  filter(dup==TRUE) %>%
  select(-dup, -dup1, -dup2) %>%
  mutate(group_id = seplyr::group_indices_se(., mvarsdjb)) %>%
  arrange(group_id, ftype) %>%
  group_by(group_id) %>%
  mutate(ndups=n(), npuf=sum(ftype=="puf"), nrf=sum(ftype=="rf")) %>%
  ungroup %>%
  filter(npuf>0, nrf>0)
b <- proc.time()
b - a # 16 secs

glimpse(djbdups)
count(djbdups, ftype) # 134203 rf, 49807 puf - great same as Max gets



#****************************************************************************************************
#                check ####
#****************************************************************************************************
glimpse(df)
names(df) %>% sort
count(df, ftype)
# check on number of exact matches
# collapse
# dfcol <- df %>%
#   group_by_at(matchvars)
# create unique records of puf and syn
upuf <- df %>%
  filter(ftype=="puf") %>%
  distinct()

ucart <- df %>%
  filter(ftype=="cart") %>%
  distinct()

uboth <- bind_rows(upuf, ucart)

ubdup <- uboth %>%
  group_by_at(vars(-ftype))

tmp <- ubdup %>%
  summarise(n=n())
names(tmp) %>% sort
ht(tmp %>% ungroup %>% select(n))

u2 <- uboth %>%
  mutate(group_id = seplyr::group_indices_se(., matchvars))
tmp <- uboth %>%
  group_by(group_id) %>%
  summarise(ndups=n(), npuf=sum(ftype=="puf"), ncart=sum(ftype=="cart"))



u2 <- uboth %>%
  mutate(ftype=uboth$ftype)

count(u2, ftype)


uboth <- bind_rows(upuf, ucart)
(grouping_vars <- setdiff(names(uboth), "ftype"))

a <- proc.time()
mcounts <- uboth %>%
  # filter(MARS==1) %>%
  group_by_at(vars(one_of(grouping_vars))) %>%
  summarise(ndups=n(), 
            npuf=sum(ftype=="puf"), 
            ncart=sum(ftype=="cart")
            rnmin=min(rownum),
            rnmax=max(rownum)) %>%
  ungroup
b <- proc.time()
b - a # 10 minutes
saveRDS(mcounts, paste0(globals$tc.dir, "mcounts_new.rds"))


#****************************************************************************************************
#                Prep ####
#****************************************************************************************************

puf.vnames <- get_puf_vnames()

seedvars <- c("MARS", "wt", "e00100", "e04600", "p04470", "e04800", "e62100", "e05800", "e08800", "e59560", "e26190")
puf.vnames %>% filter(vname %in% seedvars)


pufbase <- get_puf.base()


pufsyn <- readRDS(paste0(globals$tc.dir, "synth10syn20_stack.rds"))
glimpse(pufsyn)
count(pufsyn, ftype)
names(pufsyn) %>% sort
pufsyn %>%
  group_by(ftype) %>%
  summarise(RECID=anyDuplicated(RECID), recna=sum(is.na(RECID)))

# define matching vars
exclvars <- c(c("ftype", "c00100", "divratio", "group", "m", "mgroup", "msname",
                "pufseqn", "RECID", "rownum", "taxbc", "wfsrecid"), 
              str_subset(names(pufsyn), "minus"),
              str_subset(names(pufsyn), "0p"),
              str_subset(names(pufsyn), "0s"),
              str_subset(names(pufsyn), "wt"))
exclvars
(matchvars <- setdiff(names(pufsyn), exclvars))

df <- pufsyn %>%
  filter(ftype!="rf") %>%
  select(ftype, RECID, matchvars)
glimpse(df)

catvars <- c("DSI", "EIC", "f2441", "f6251", "FDED", "MARS", "MIDR", "n24", "XTOT")
catgroups <- df %>% 
  filter(ftype=="puf") %>%
  select(catvars) %>%
  group_by_at(vars(one_of(catvars))) %>%
  summarise(n=n()) %>%
  ungroup
ht(catgroups)
catgroups %>% arrange(n)
count(catgroups, n)
sum(catgroups$n)

count(df, MARS)

# either use the full df, or make a small version to look at
# mv <- matchvars[1:15]
# df2 <- df %>%
#   filter(MARS==3) %>%
#   select(ftype, mv)


# check on number of exact matches
# collapse
# dfcol <- df %>%
#   group_by_at(matchvars)
  
# use df or df2
a <- proc.time()
dups <- df %>%
  mutate(dup1=duplicated(select(., -ftype, -RECID)),
         dup2=duplicated(select(., -ftype, -RECID), fromLast=TRUE),
         dup=(dup1==TRUE) | (dup2==TRUE))
b <- proc.time()
b - a # 90 secs

ht(dups)
count(dups, dup1, dup2)
count(dups, ftype, dup1, dup2)
count(dups, ftype, dup)
count(dups %>% mutate(dup=(dup1==TRUE) | (dup2==TRUE)), ftype, dup)
saveRDS(dups, paste0(globals$tc.dir, "dupsmarked.rds"))

#dups <- readRDS(paste0(globals$tc.dir, "dupsmarked.rds"))
count(dups, ftype)
count(dups, ftype, dup1, dup2)

a <- proc.time()
dups2 <- dups %>%
  filter(dup==TRUE) %>%
  select(-dup, -dup1, -dup2) %>%
  mutate(group_id = seplyr::group_indices_se(., matchvars)) %>%
  arrange(group_id, ftype) %>%
  group_by(group_id) %>%
  mutate(ndups=n(), npuf=sum(ftype=="puf"), nsyn=sum(ftype=="cart")) %>%
  ungroup %>%
  filter(npuf>0, nsyn>0)
b <- proc.time()
b - a # 16 secs
saveRDS(dups2, paste0(globals$tc.dir, "dupstrue.rds"))
count(dups2, ftype)
glimpse(dups2)
tmp <- dups2 %>%
  filter(ftype=="puf", group_id==2)
tmp2 <- df %>%
  filter(ftype=="puf", RECID %in% tmp$RECID)
write_csv(tmp2, "d:/tcdir/tmp2.csv")

dups3 <- dups2 %>%
  select(-RECID) %>%
  group_by(group_id, ftype) %>%
  distinct() %>%
  ungroup
count(dups3, ftype)
dups3 %>%
  group_by(ftype) %>%
  summarise(npuf=sum(npuf), nsyn=sum(nsyn))
glimpse(dups3)

# one record for every unique combination
dups4 <- dups3 %>%
  select(-ftype) %>%
  distinct()
dups4 %>%
  summarise(n=n(), ndups=sum(ndups), npuf=sum(npuf), nsyn=sum(nsyn))
# n    ndups  npuf nsyn
# 2740 18521 12480  6041
# so we have 2,740 unique kinds of exact matches
12480 / 163786
6041

# what are they like?
names(dups4) %>% sort # 68 vars, incl group_id, ndups, npuf, nsyn
idvars <- c("group_id", "ndups", "npuf", "nsyn")
# EIC Earned Income Credit Code # children claimed
# f2441 Number of qualifying individuals child care credit
# f6251 AMT return?
# fded form of deduction code
# MIDR Married Filing Separately Itemized Deductions Requirement Indicator
# n24 Number of Children for Child Tax Credit
catvars <- c("DSI", "EIC", "f2441", "f6251", "FDED", "MARS", "MIDR", "n24", "XTOT")

duptypes <- dups4 %>%
  gather(vname, value, -idvars) %>%
  mutate(vtype=case_when(vname %in% catvars ~ "categorical",
                         !vname %in% catvars ~ "continuous",
                         TRUE ~ "other"))
count(duptypes, vtype)
count(duptypes, vtype, vname)
ht(duptypes)

# which continuous vars are nz?
contvars <- duptypes %>%
  filter(vtype=="continuous", value!=0) %>%
  arrange(group_id, vname) %>%
  group_by(group_id, ndups, npuf, nsyn) %>%
  summarise(nnz=n(), vname=paste0(vname, collapse="-")) %>%
  ungroup
count(contvars %>% filter(nnz==1), vname) %>%
  left_join(puf.vnames %>% select(vname, vdesc)) %>%
  arrange(-n)
count(contvars, nnz)

count(contvars %>% filter(nnz==2), vname) %>%
  arrange(-n)

count(contvars %>% filter(nnz==3), vname) %>%
  arrange(-n)

duptypes %>%
  filter(vtype=="continuous", value!=0) %>%
  select(-vtype) %>%
  group_by(group_id) %>%
  mutate(nnz=n()) %>%
  filter(nnz==4) %>%
  left_join(puf.vnames %>% select(vname, vdesc)) %>%
  arrange(group_id, vname)



typecounts <- duptypes %>%
  group_by(group_id, ndups, npuf, nsyn, vtype) %>%
  summarise(n=n(), nnz=sum(value!=0))
typecounts

typecounts %>%
  filter(vtype=="continuous") %>%
  group_by(nnz) %>%
  summarise(n=n())


typecounts %>%
  filter(vtype=="continuous") %>%
  arrange(-nnz)

# group_id 5032 has 5 nnz continuous matches
duptypes %>%
  filter(group_id==5032) %>%
  left_join(puf.vnames %>% select(vname, vdesc)) %>%
  arrange(desc(value))
# group_id ndups  npuf  nsyn vname  value vtype       vdesc                                                        
# 1     5032     2     1     1 e02400 17100 continuous  Gross Social Security benefits                               
# 2     5032     2     1     1 e01700  6510 continuous  Pensions and annuities included in AGI                       
# 3     5032     2     1     1 e01500  6510 continuous  Total pensions and annuities received                        
# 4     5032     2     1     1 e01400  3950 continuous  Taxable IRA distribution                                     
# 5     5032     2     1     1 e00300   990 continuous  Interest received                                            
# 6     5032     2     1     1 FDED       2 categorical NA                                                           
# 7     5032     2     1     1 MARS       1 categorical Marital status                                               
# 8     5032     2     1     1 XTOT       1 categorical Exemptions                                                   
# 9     5032     2     1     1 EIC        0 categorical NA                                                           
# 10     5032     2     1     1 f6251      0 categorical Form 6251, Alternative Minimum Tax                           
# 11     5032     2     1     1 f2441      0 categorical Form 2441, Child Care Credit Qualified Individual        
  

check <- dups3 %>%
  ungroup %>%
  filter(row_number()==1) %>%
  select(-ftype, -group_id, -ndups, -npuf, -nsyn)

tmp <- check %>% left_join(df)
count(tmp, ftype)

sum(dups3$npuf)

# dups2 <- readRDS(paste0(globals$tc.dir, "dupstrue.rds"))
count(dups2, ftype)
# 1 puf   12480
# 2 cart   6041
length(unique(dups2$group_id)) # 2,740 unique groups
check <- dups2 %>% filter(ndups==max(ndups))
write_csv(check, "d:/tcdir/check.csv")

# get unduplicated count of number of dup records
dupcounts <- dups2 %>%
  group_by(group_id) %>%
  summarise(ndups=first(ndups),
            npuf=first(npuf),
            nsyn=first(nsyn))
dupcounts %>%
  summarise_at(vars(ndups, npuf, nsyn), ~sum(.))



dups2 %>%
  filter(group_id==6)

# alternative approach - grouping
dupsalt <- df %>%
  group_by_at(matchvars) %>%
  mutate(ndups=n(), npuf=sum(ftype=="puf"), nsyn=sum(ftype=="cart")) %>%
  ungroup %>%
  filter(npuf>0, nsyn>0)
glimpse(dupsalt)
count(dupsalt, ftype)


# % of recs that are duplicated in a puf-syn group
+18521 / 163786 # puf 11.3%
+12480 / 163786 # puf 7.6%
+6041 / 818930  # syn 0.7%

# Max says:
# 16.4% of RF synthetic records exactly match one or more true records
# 8.9% of CART synthetic records exactly match one or more true records




ht(dups2[, c(1:5, (ncol(dups3)-4):ncol(dups3))], 10)
length(unique(dups2$group_id))
sum(dups2$npuf)
sum(dups2$nsyn)           
           
# group_by_at(vars(one_of(matchvars))) %>%
group_indices(.data, ...)
names(dups3) %>% sort
ht(dups3)

mtcars %>% mutate(group = group_indices(.dots="cyl"))
mtcars %>% mutate(group = group_indices(., .dots=vars(select(c("cyl", "carb")))))


vars <- c("cyl", "carb")
mtcars %>% mutate(group = group_indices(!!enquo(vars)))
mtcars %>% mutate(group = group_indices(., !!!vars))

bygr <- mtcars %>% group_by(cyl, carb)
n_groups(bygr)
group_size(bygr)




mtcars %>% mutate(group = seplyr::group_indices_se(., vars))
group_indices_se(.data, groupingVars, add = FALSE)


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

# old stuff ----
tmp <- cartnear2 %>% 
  select(pufid, cartid) %>% mutate(pairtype="max")
# find pairs that max has that i don't
nonmatchpairs <- bind_rows(cartnear2 %>% 
                             select(cid=pufid, pid=cartid) %>% 
                             mutate(pairtype="max") %>%
                             rename(cartid=cid, pufid=pid),
                           cartpairs %>% 
                             select(pufid, cartid) %>% mutate(pairtype="don")) %>%
  arrange(pufid, cartid, pairtype) %>%
  mutate(dup1=duplicated(select(., -pairtype)),
         dup2=duplicated(select(., -pairtype), fromLast=TRUE),
         dup=(dup1==TRUE) | (dup2==TRUE)) %>%
  filter(dup==FALSE)

# 3 760862
nmmax1 <- bind_rows(puf %>% filter(id==3) %>% select(names(cart)) %>% mutate(ftype="puf"),
                    cart %>% filter(id==760862) %>% mutate(ftype="cart"))

pufrec <- 120794
check <- bind_rows(cart %>% filter(id==2) %>% mutate(ftype="cart"),
                   puf %>% filter(id %in% (pufrec-5):(pufrec+5)) %>% 
                     select(names(cart)) %>% 
                     mutate(ftype="puf"),
) %>%
  select(ftype, id, everything()) %>%
  arrange(ftype, id)

d <- puf %>% filter(E00200==3527)

cartrec <- 120794
check2 <- bind_rows(puf %>% filter(id==2) %>% 
                      mutate(ftype="puf") %>%
                      select(names(cart)),
                    cart %>% filter(id %in% (cartrec-5):(cartrec+5)) %>% 
                      mutate(ftype="cart")) %>%
  select(ftype, id, everything()) %>%
  arrange(ftype, id)

check3 <- cart %>% 
  filter(MARS==3, XTOT==1, FDED==2, E00200==17500)

matchids <- c(120794, 64622, 99507)
check2 <- bind_rows(cart %>% filter(id %in% matchids) %>% 
                      mutate(ftype="cart"),
                    puf %>% filter(id %in% matchids) %>% 
                      mutate(ftype="puf") %>%
                      select(names(cart))) %>%
  select(ftype, id, everything()) %>%
  arrange(ftype, id) %>%
  filter(ftype=="puf")

ns(check2)

check3 <- cart %>% 
  filter(MARS==1, XTOT==1, FDED==3, F6251==0,
         E00200==0, E00300==0, E00400==0, E00600==0, E00650==0, 
         E00700==0, E00800==0, E00900==0, E01100==0, E01200==0,
         E01400==0, E01500==0, E01700==0, E02000==0, E02100==0, 
         E02300==0, E02400==0, E03210==0, E03230==0, E03240==0,
         E03270==0, E07300==0, E62900==0,
         E20100==0, E26270==0, E87521==0,
         P22250==0, P23250==0)
check3 %>%
  select(-ftype, -id) %>%
  distinct() %>% 
  as.data.frame
check3$id



