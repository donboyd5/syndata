
# NOTE: create version of puf with lower case variable names and with wt=S006 / 100


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
#                Create file ####
#****************************************************************************************************
puf.base.sums <- get_puf.base() %>%
  filter(!RECID %in% 999996:999999) %>%
  mutate(e00600_minus_e00650=e00600 - e00650,
         e01500_minus_e01700=e01500 - e01700,
         divratio=e00650 / e00600,
         divratio=ifelse(is.na(divratio), 0, divratio),
         penratio=e01700 / e01500,
         penratio=ifelse(is.na(penratio), 0, penratio)) %>% 
  select_if(is.numeric) %>%
  dplyr::select(-c(RECID, s009, s006)) %>%
  summarise_at(vars(-wt), funs(sum(. * wt) / 1e9)) %>%
  gather(vname, value) %>%
  left_join(puf.vnames %>% dplyr::select(vname, vdesc)) %>% 
  arrange(-abs(value))

saveRDS(puf.base.sums, "./data/puf.base.sums.rds")


