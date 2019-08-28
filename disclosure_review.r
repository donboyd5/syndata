

#****************************************************************************************************
#                Libraries ####
#****************************************************************************************************
library("magrittr")
library("plyr") # needed for ldply; must be loaded BEFORE dplyr
library("tidyverse")
options(tibble.print_max = 60, tibble.print_min = 60) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats

# devtools::install_github("donboyd5/btools")
library("btools") # library that I created (install from github)


#****************************************************************************************************
#                Includes ####
#****************************************************************************************************
source("./r/includes/globals_system_specific_boyd.r") # use a different version of this file if changing systems
source("./r/includes/globals_other.r")

source("./r/includes/functions_general.r")


#****************************************************************************************************
#                globals ####
#****************************************************************************************************
gd <- "c:/Users/donbo/Google Drive/"

pufd <- paste0(gd, "synpuf/puf/")
synd <- paste0(gd, "synpuf/syntheses/")

puf_fn <- "puf2011.csv"
synrf_fn <- "synpuf20.csv"


#****************************************************************************************************
#                functions ####
#****************************************************************************************************
ns <- function(df) {names(df) %>% sort } # names sort 


#****************************************************************************************************
#                define disclosure-review variables ####
#****************************************************************************************************
df1 <- read_csv(paste0(synd, synrf_fn), n_max=10)
glimpse(df1)
ns(df1)
nodup_vars <- setdiff(names(df1), c("RECID", "S006"))
nodup_vars # 64 variables


#****************************************************************************************************
#                check duplicates methodology ####
#****************************************************************************************************
tmp <- tribble(
  ~a, ~b, ~c,
  #--|--|----
  10, 2, 3.6,
  20, 1, 8.5,
  10, 2, 3.6,
  10, 2, 3.6,
  10, 2, 3.6
)
tmp
unique(tmp) # unique will include 1 copy of each kind of record, including duplicates

duplicated(tmp) # doesn't catch that the first record is duplicated
flag_dups <- function(df) {duplicated(df) | duplicated(df, fromLast=TRUE)} # this gets all dups
flag_dups(tmp) # this is correct
sum(flag_dups(tmp))


#****************************************************************************************************
#                create a no-dups-allowed version of the puf ####
#****************************************************************************************************
puf <- read_csv(paste0(pufd, puf_fn))
ns(puf)

puf2 <- puf[, nodup_vars]
glimpse(puf2)
ns(puf2)

system.time(idups <- flag_dups(puf2)) # get a logical vector indicating which puf records are duplicates
sum(idups) # 19,527 dups
idups[1:100]

puf_undup <- puf2[!idups, ] # get the puf records that are unduplicated as these are the only ones we worry about


#****************************************************************************************************
#                create a no-dups version of the syn rf file ####
#****************************************************************************************************
syn <- read_csv(paste0(synd, synrf_fn))
glimpse(syn)
ns(syn)

system.time(idups_syn <- flag_dups(syn[, nodup_vars])) # 85 secs
sum(idups_syn) # 151767

syn2 <- syn %>%
  mutate(dup=ifelse(idups_syn==TRUE, 1, 0))
sum(syn2$dup)

system.time(usyn <- unique(syn[, nodup_vars])) # we only need 1 copy of each record - this will make things faster


#****************************************************************************************************
#                stack the no-dups puf and the no-dups version of the syn rf file ####
#****************************************************************************************************
stack <- bind_rows(puf_undup %>% mutate(ftype="puf"),
                   usyn %>% mutate(ftype="syn"))

# any dups in stack will be bad dups!
bad_dups <- flag_dups(stack[, nodup_vars])
sum(bad_dups) # counts them twice, once in puf, once in syn
# 60946 dups!
stack2 <- stack %>%
  mutate(dup=ifelse(bad_dups==TRUE, 1, 0))
sum(stack2$dup)

pufdups <- stack2 %>%
  filter(ftype=="puf", dup==1)
glimpse(pufdups)
summary(pufdups)
ns(pufdups) # 66 vars - drop dup and ftype


# find the bad dups records and delete them!
# synstack <- bind_rows(usyn %>% select(nodup_vars) %>% mutate(ftype="syn"),
#                       pufdups %>% select(nodup_vars) %>% mutate(ftype="pufdup"))
# ns(synstack)
# 
# synbad <- flag_dups(synstack[, nodup_vars])
# 
# synbad <- puf_dupssyn[, nodup_vars] %>%
#   left_join()

ns(syn)
synbad <- syn %>%
  left_join(pufdups[, nodup_vars] %>% mutate(ftype="puf"),
            by=nodup_vars) %>%
  mutate(ftype=ifelse(is.na(ftype), "nomatch", "match"))
count(synbad, ftype)
68583 / 818930

pct <- function(x) x / sum(x) * 100
synbad %>%
  group_by(ftype) %>%
  summarise(n=n(), wt=sum(S006) / 100, wagesb=sum(E00200 * S006/100) / 1e9) %>%
  mutate_at(vars(n, wt, wagesb), list(pct=pct))

# A tibble: 2 x 7
# ftype      n         wt wagesb n_pct wt_pct wagesb_pct
# <chr>  <int>      <dbl>  <dbl> <dbl>  <dbl>      <dbl>
#   1 puf    68583  95553639.  2281.  8.37   13.2       7.54
# 2 NA    750347 630740850. 27953. 91.6    86.8      92.5 
# Max has 750,347 also
















