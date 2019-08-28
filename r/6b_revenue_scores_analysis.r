
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


#****************************************************************************************************
#                Get Dan's revenue scores data ####
#****************************************************************************************************
rsdir <- paste0(globals$synd, "revenue_scores/")

base <- read_csv(paste0(globals$synd, "synthpop5_all_for_revenuescores.csv"), 
                 col_types = cols(ftype=col_character(), MARS=col_character(), msname=col_character(),
                                  .default= col_double()), 
                 n_max=-1)
glimpse(base)
summary(base)
base %>%
  group_by(ftype) %>%
  summarise(rmin=min(RECID), rmax=max(RECID))

rs.puf <- read_csv(paste0(rsdir, "pufscores.csv"))
glimpse(rs.puf)

rs.syn <- read_csv(paste0(rsdir, "synscores.csv"))
glimpse(rs.syn)


#****************************************************************************************************
#                Get subset we like and link to other data ####
#****************************************************************************************************
count(base, ftype)
count(rs.puf, name)

rsp2 <- rs.puf %>% 
  filter(name=="_rt5")

rss2 <- rs.syn %>% 
  filter(name=="_rt5")

ht(rsp2)
ht(rss2)

link <- base %>%
  left_join(rsp2 %>% select(-name)) %>%
  left_join(rss2 %>% select(-name))

tmp <- link %>%
  select(ftype, RECID, c00100, taxbc, starts_with("wt"), pbase, pscore, sbase, sscore)
cor(tmp$pbase, tmp$taxbc, use="pairwise.complete.obs")
cor(tmp$sbase, tmp$taxbc, use="pairwise.complete.obs")

agiranges <- c(-Inf, 0, 25e3, 50e3, 75e3, 100e3, 200e3, 500e3, 1e6, 10e6, Inf)
tsums <- tmp %>%
  mutate(wt=ifelse(ftype=="puf", wt.puf, wt.wtfs),
         base=ifelse(ftype=="puf", pbase, sbase),
         score=ifelse(ftype=="puf", pscore, sscore),
         agirange=cut(c00100, agiranges, right=FALSE),
         wtone=1e9) %>%
  group_by(ftype, agirange) %>%
  summarise(n=n(),
            wtsum=sum(wt),
            wtnnz=sum(wt * (score!=0)),
            taxbc=sum(taxbc * wt) / 1e9,
            base=sum(base * wt) / 1e9,
            score=sum(score * wt) / 1e9)

var <- "wtsum"
var <- "taxbc"
var <- "base"
var <- "score"
var <- "wtnnz"
tsums %>%
  select(ftype, agirange, value=var) %>%
  spread(ftype, value) %>%
  janitor::adorn_totals(.) %>%
  mutate(diff=syn - puf,
         pdiff=diff / puf * 100,
         var=var)


