
#****************************************************************************************************
#                Libraries ####
#****************************************************************************************************
library("magrittr")
library("plyr") # needed for ldply; must be loaded BEFORE dplyr
library("tidyverse")
options(tibble.print_max = 60, tibble.print_min = 60) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats

library("scales")
library("hms") # hms, for times.
library("lubridate") # lubridate, for date/times.
library("readxl") # readxl, for .xls and .xlsx files.
library("haven") # haven, for SPSS, SAS and Stata files.
library("vctrs")
# library("precis") # use precis2 in btools

library("tibbletime") # https://business-science.github.io/tibbletime/

library("grDevices")
library("knitr")

library("zoo") # for rollapply

library("btools") # library that I created (install from github)

# library("ipoptr")
# library("synthpop") # note: masks select in dplyr


#****************************************************************************************************
#                Globals ####
#****************************************************************************************************


#****************************************************************************************************
#                Includes ####
#****************************************************************************************************
source("./r/includes/globals_system_specific_boyd.r") # use a different version of this file if changing systems
source("./r/includes/globals_other.r")

source("./r/includes/functions_general.r")
# source("./r/includes/functions_target_setup_and_analysis.r")
# source("./r/includes/functions_ipopt.r")



# Steps in running tax plans:
# 1. Get the synthetic file and the puf
# 2. Create a puf counterpart and stack the files
# 3. Prepare the stacked file - split prime/spouse, etc.
# 4. Save the prepared stacked file
# 5. Run the prepared stacked file through each tax plan of interest
# 6. Create a list with the stacked file, and each of the reform results
# 7. Save the list

#****************************************************************************************************
#  0. ONETIME -- prepare a file from Max and its PUF counterpart ####
#****************************************************************************************************

# sfname <- "synpuf17"
# sfname <- "synpuf20"
sfname <- "synpuf20_lowmatch"

syn.from.max <- read_csv(paste0(globals$synd, "syntheses/", sfname, ".csv"), 
                         col_types = cols(.default= col_double()), 
                         n_max=-1)
glimpse(syn.from.max)

# variable diffs vs. my files:
# - sfm has:
#   -- s006 -- will drop and replace with wt.syn
#   -- e00600_minus_e00650 -- not targeted; however, keep in case desired in future
# - my files have:
#   -- pufseqn -- don't think we need
#   -- ftype -- will compute when I add puf records
#   -- m -- we need this, so compute m=1 for all
#   -- msname -- don't need
#   -- wt -- will compute 
#   -- divratio -- not needed
#   -- e09600 -- not needed, I synthesized by mistake

max.syn <- syn.from.max %>%
  setNames(change_case(names(.))) %>%
  mutate(e01500_minus_e01700=e01500 - e01700,
         wt.syn=s006 / 100,
         ftype="syn") %>%
  dplyr::select(-s006)

names(max.syn) %>% sort

max.puf <- get_puf.base() %>%
  dplyr::select(one_of(names(max.syn)), wt.puf=wt, ftype, pufseqn, RECID) %>%
  mutate(e00600_minus_e00650=e00600 - e00650,
         e01500_minus_e01700=e01500 - e01700)

setdiff(names(max.syn), names(max.puf))
setdiff(names(max.puf), names(max.syn))
  
max.stack <- bind_rows(max.puf, 
                       max.syn) %>% 
  mutate(m=1)
count(max.stack, ftype)
glimpse(max.stack)
summary(max.stack)

saveRDS(max.stack, paste0(globals$tc.dir, "pufsyn_stack/", sfname, ".rds"))

# nstack <- names(stack)
# nstack %>% sort
# 
# setdiff(names(max.syn), nstack)
# setdiff(nstack, names(max.syn))



#****************************************************************************************************
#  1. ONETIME PER SYNFILE - Optionally prepare a file and its PUF counterpart ####
#****************************************************************************************************
#.. ONLY RUN THIS IF NOT PREVIOUSLY DONE FOR A PARTICULAR SYNTHETIC FILE!! ----

# we call get_pufsyn_taxplan_results(stacked.pufsyn, taxplan.fn=NULL)
# get_pufsyn_taxplan_results requires as inputs:
#    stacked.pufsyn -- a dataframe with stacked tax syn file and its corresponding puf (same variables) -- variable name case will be adjusted
#    taxplan.fn -- character variable holding the name (without directory) of the taxplan json file e.g., "rate_cut.json"
#    BUT taxplan.fn SHOULD BE NULL BELOW


#.. My files ----
# sfname <- "synthpop5_all"
# sfname <- "synthpop6_all"
# sfname <- "synthpop7_all"
# sfname <- "synthpop8"
sfname <- "synthpop10"
stack <- readRDS(paste0(globals$tc.dir, "pufsyn_stack/", sfname, ".rds"))$stack

# beginning with synthpop10, _stack is in the file name
stack <- readRDS(paste0(globals$tc.dir, "pufsyn_stack/", sfname, "_stack.rds"))$stack

#.. Max's files start with synpuf, not synthpop, and are stored as just a file, not part of a list ----
# sfname <- "synpuf17"
# sfname <- "synpuf20"
sfname <- "synpuf20_lowmatch"
stack <- readRDS(paste0(globals$tc.dir, "pufsyn_stack/", sfname, ".rds"))

# Now that we have stack, we are ready to prepare a file for weighting ----
names(stack) %>% sort
glimpse(stack)
count(stack, ftype, m)

# stack$pufseqn[1:1000]

system.time(synprep <- get_pufsyn_taxplan_results(stack)) # this can take x minutes

quantile(stack$e00200)

# take a quick look at the results
names(synprep)
glimpse(synprep$tc.base)
glimpse(synprep$tc.output)
count(synprep$tc.base, ftype)
summary(synprep$tc.output)

saveRDS(synprep, paste0(globals$tc.dir, sfname, "_rwprep.rds"))


#****************************************************************************************************
#  take a look ####
#****************************************************************************************************

# sfname <- "synthpop6_all"
# sfname <- "synthpop7_all"
# sfname <- "synthpop8"
sfname <- "synthpop10"

# sfname <- "synpuf17"

synprep <- readRDS(paste0(globals$tc.dir, sfname, "_rwprep.rds"))
names(synprep)
glimpse(synprep)
names(synprep$tc.base) %>% sort

tcvars <- c("c00100", "taxbc", "c09600", "c05800")
mrgdf <- left_join(synprep$tc.base, synprep$tc.output %>% dplyr::select(RECID, tcvars))
glimpse(mrgdf)
count(mrgdf, ftype, m)

#.. OPTIONALLY write csv files to synpuf ----
globals
mrgdf %>% 
  filter(ftype=="puf") %>%
  write_csv(paste0(globals$synd, sfname, "_puf.csv"))

mrgdf %>% 
  filter(ftype=="syn") %>%
  write_csv(paste0(globals$synd, sfname, ".csv"))

