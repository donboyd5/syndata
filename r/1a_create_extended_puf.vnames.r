
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
#                Get puf variable names ####
#****************************************************************************************************
puf.vnames <- read_csv("./data/pufvars.csv")
addnames <- read_csv("vnum, vname, vdesc, category
                     0, c00100, Adjusted gross income (calculated), income
                     9000, e00200p, Wages-prime (calculated), income
                     9001, e00200s, Wages-spouse (calculated)
                     9010, taxbc, Tax before credit (calculated)
                     9020, MARS, Marital status
                     9030, XTOT, Exemptions")
puf.vnames <- bind_rows(addnames, puf.vnames) %>%
  arrange(vnum, vname)
ht(puf.vnames)
saveRDS(puf.vnames, "./data/puf.vnames.rds")


