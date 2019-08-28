
source(file.path(PROJHOME, "r/includes", "libraries.r"))
library("readxl")
library("knitr")
library("kableExtra")
# search()
# getwd()

# source("./r/includes/globals_system_specific_boyd.r")
source(file.path(PROJHOME, "r/includes", "globals_system_specific_boyd.r"))
source(file.path(PROJHOME, "r/includes", "globals_other.r"))

puf <- read_csv(paste0("D:/Dropbox/OSPC - Shared/IRS_pubuse_2011/", "puf2011.csv"), 
                col_types = cols(.default= col_double()), 
                n_max=-1)
names(puf)
glimpse(puf)


puf %>% 
  filter(!RECID %in% 999996:999999) %>%
  summarise_at(vars(E00200, E01500), list(~n(), ~mean(.)))


puftaxdata <- read_csv(paste0("D:/Dropbox/OSPC - Shared/IRS_pubuse_2011/", "puf.csv"), 
                col_types = cols(.default= col_double()), 
                n_max=-1)
glimpse(puftaxdata)
count(puftaxdata, data_source)

puftaxdata %>% 
  filter(data_source==1, !RECID %in% 999996:999999) %>%
  summarise_at(vars(e00200, e01500), list(~n(), ~mean(.)))



