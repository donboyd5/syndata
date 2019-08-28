

#****************************************************************************************************
#                Libraries ####
#****************************************************************************************************
library("magrittr")
library("plyr") # needed for ldply; must be loaded BEFORE dplyr
library("tidyverse")
options(tibble.print_max = 60, tibble.print_min = 60) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats
library("scales")

library("btools") # should not need this. If you do, it's at https://github.com/donboyd5/btools

library("knitr")

library("synthpop") # note: masks select in dplyr


#****************************************************************************************************
#                Includes ####
#****************************************************************************************************
source("./r/includes/globals_system_specific_boyd.r") # use a different version of this file if changing systems
source("./r/includes/globals_other.r")

source("./r/includes/functions_general.r")

puf.vnames <- get_puf_vnames()
# puf.names

puf <- get_puf.base()

names(puf) %>% sort
count(puf, MARS)


#****************************************************************************************************
#                Prepare checklist ####
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

psum <- puf.base %>% 
  select_if(is.numeric) %>%
  select(-c(RECID, s009, s006)) %>%
  summarise_at(vars(-wt), funs(sum(. * wt) / 1e9)) %>%
  gather(vname, value) %>%
  left_join(puf.vnames %>% dplyr::select(vname, vdesc)) %>% 
  arrange(-abs(value))
psum

write_csv(psum, "./data/checklist_base.csv")




