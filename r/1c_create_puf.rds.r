
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
puf.fn <- "puf2011.csv"

# variables that I believe Tax-Calculator expects to be upper case
# note that FLPDYR and RECID are not synthesized
upper_case_vars <- globals$upper_case_vars #  c("RECID", "MARS", "XTOT", "DSI", "EIC", "FLPDYR", "MIDR") 

puf <- read_csv(paste0(globals$pufdir, puf.fn), 
                col_types = cols(.default= col_double()), 
                n_max=-1)
# quick exploration
glimpse(puf)
puf$RECID %>% sort # wow, RECID is not consecutive!
puf$RECID %>% sort %>% tail(20)

# change selected names to lower case
names(puf)
setdiff(upper_case_vars, names(puf)) # are all ucv's in puf?
# indexes for upper and lower case names
iuc <- which(names(puf) %in% upper_case_vars)
ilc <- which(!names(puf) %in% upper_case_vars)
adj.names <- names(puf)
adj.names[ilc] <- str_to_lower(names(puf)[ilc])
adj.names
names(puf)

puf2 <- puf
names(puf2) <- adj.names

puf2 <- puf2 %>% mutate(wt=s006 / 100, pufseqn=row_number())
glimpse(puf2)
saveRDS(puf2, paste0(globals$tc.dir, "puf_lc.rds"))


