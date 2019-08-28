
source(file.path(PROJHOME, "r/includes", "libraries.r"))
library("readxl")
library("knitr")
library("kableExtra")
# search()
# getwd()

# source("./r/includes/globals_system_specific_boyd.r")
source(file.path(PROJHOME, "r/includes", "globals_system_specific_boyd.r"))
source(file.path(PROJHOME, "r/includes", "globals_other.r"))

dir <- "D:/tcdir/puf_vs_csv/"

puf <- read_csv(paste0(dir, "puf.csv"), 
                col_types = cols(.default= col_double()), 
                n_max=-1)
names(puf) %>% sort
glimpse(puf)
count(puf, data_source)

cps <- read_csv(paste0(dir, "cps.csv"), 
                col_types = cols(.default= col_double()), 
                n_max=-1)
names(cps) %>% sort
glimpse(cps)



intersect(names(puf), names(cps)) %>% sort
setdiff(names(puf), names(cps)) %>% sort
setdiff(names(cps), names(puf)) %>% sort


# baseline files
pufbl <- read_csv(paste0(dir, "puf-14-#-#-#.csv"), 
                  col_types = cols(.default= col_double()), 
                  n_max=-1)
names(pufbl) %>% sort


cpsbl <- read_csv(paste0(dir, "cps-14-#-#-#.csv"), 
                  col_types = cols(.default= col_double()), 
                  n_max=-1)


# quick checks
intersect(names(puf), names(pufbl)) %>% sort
intersect(names(cps), names(cpsbl)) %>% sort

sum(puf$s006)
sum(pufbl$s006)

sum(cps$s006)
sum(cpsbl$s006)


pufw <- puf %>% rename(s006_source=s006) %>%
  left_join(pufbl %>% select(RECID, s006, setdiff(names(pufbl), names(puf)))) %>%
  mutate(ftype="puf")
pufw$c05800
sum(pufw$c05800)

cpsw <- cps %>% rename(s006_source=s006) %>%
  left_join(cpsbl %>% select(RECID, s006, setdiff(names(cpsbl), names(cps)))) %>%
  mutate(ftype="cps")

pc <- bind_rows(pufw, cpsw)

pc %>%
  group_by(ftype) %>%
  summarise(n=n(), wts=sum(s006), 
            c00100=sum(s006 * c00100, na.rm=TRUE) / 1e9,
            taxbc=sum(s006 * taxbc, na.rm=TRUE) / 1e9, 
            c05800=sum(s006 * c05800, na.rm=TRUE) / 1e9)
saveRDS(pc, paste0(dir, "pufcps_2014.rds"))

# now run a basic ratecut reform
# puf-14-#-ratecut2014-#.csv
# cps-14-#-ratecut2014-#.csv
pufrc <- read_csv(paste0(dir, "puf-14-#-ratecut2014-#.csv"), 
                  col_types = cols(.default= col_double()), 
                  n_max=-1)
names(pufrc) %>% sort

cpsrc <- read_csv(paste0(dir, "cps-14-#-ratecut2014-#.csv"), 
                  col_types = cols(.default= col_double()), 
                  n_max=-1)
names(cpsrc) %>% sort


pufrcw <- puf %>% rename(s006_source=s006) %>%
  left_join(pufrc %>% select(RECID, s006, setdiff(names(pufrc), names(puf)))) %>%
  mutate(ftype="pufrc")

cpsrcw <- cps %>% rename(s006_source=s006) %>%
  left_join(cpsrc %>% select(RECID, s006, setdiff(names(cpsrc), names(cps)))) %>%
  mutate(ftype="cpsrc")

# simple checks
sum(pufrcw$c05800 * pufrcw$s006) / 1e9
sum(pufw$c05800 * pufw$s006) / 1e9

sum(cpsrcw$c05800 * cpsrcw$s006) / 1e9
sum(cpsw$c05800 * cpsw$s006) / 1e9

# complete the ratecut file
pcrc <- bind_rows(pufrcw, cpsrcw)
pcrc %>%
  group_by(ftype) %>%
  summarise(n=n(), wts=sum(s006), 
            c00100=sum(s006 * c00100, na.rm=TRUE) / 1e9,
            taxbc=sum(s006 * taxbc, na.rm=TRUE) / 1e9, 
            c05800=sum(s006 * c05800, na.rm=TRUE) / 1e9)
saveRDS(pcrc, paste0(dir, "pufcps_2014ratecut.rds"))

# now the complex reform
pufr4 <- read_csv(paste0(dir, "puf-14-#-reform42014-#.csv"), 
                  col_types = cols(.default= col_double()), 
                  n_max=-1)
names(pufr4) %>% sort

cpsr4 <- read_csv(paste0(dir, "cps-14-#-reform42014-#.csv"), 
                  col_types = cols(.default= col_double()), 
                  n_max=-1)
names(cpsr4) %>% sort


pufr4w <- puf %>% rename(s006_source=s006) %>%
  left_join(pufr4 %>% select(RECID, s006, setdiff(names(pufr4), names(puf)))) %>%
  mutate(ftype="pufr4")

cpsr4w <- cps %>% rename(s006_source=s006) %>%
  left_join(cpsr4 %>% select(RECID, s006, setdiff(names(cpsr4), names(cps)))) %>%
  mutate(ftype="cpsr4")

pcr4 <- bind_rows(pufr4w, cpsr4w)
pcr4 %>%
  group_by(ftype) %>%
  summarise(n=n(), wts=sum(s006), 
            c00100=sum(s006 * c00100, na.rm=TRUE) / 1e9,
            taxbc=sum(s006 * taxbc, na.rm=TRUE) / 1e9, 
            c05800=sum(s006 * c05800, na.rm=TRUE) / 1e9)
saveRDS(pcr4, paste0(dir, "pufcps_2014reform4.rds"))


# stack the baseline and reform files and analyze differences
stack <- bind_rows(pc %>% mutate(reform="base"), 
                   pcr4 %>% mutate(reform="reform")) %>%
         mutate(ftype=str_sub(ftype, 1, 3))
count(stack, ftype, reform)

# summaries by income range
agiranges <- c(-Inf, 0, 25e3, 50e3, 75e3, 100e3, 200e3, 500e3, 1e6, 10e6, Inf)

chg <- stack %>%
  mutate(agirange=cut(c00100, agiranges, right=FALSE),
         agirange=fct_expand(agirange, "total"),
         msgroup=case_when(MARS==1 ~ "married",
                           MARS==2 ~ "single",
                           TRUE ~ "other"),
         msgroup=fct_expand(msgroup, "total")) %>%
  group_by(ftype, reform, msgroup, agirange) %>%
  summarise(tax=sum(c05800 * s006) / 1e9)
# add in totals across income ranges before presenting
chg2 <- chg %>%
  bind_rows(chg %>% 
              group_by(ftype, reform, msgroup) %>%
              summarise(tax=sum(tax)) %>%
              mutate(agirange="total",
                     agirange=factor(agirange, levels=levels(chg$agirange))))
chg3 <- chg2 %>%
  bind_rows(chg2 %>% 
              group_by(ftype, reform, agirange) %>%
              summarise(tax=sum(tax)) %>%
              mutate(msgroup="total",
                     msgroup=factor(msgroup, levels=levels(chg$msgroup))))

chg4 <- chg3 %>%
  spread(reform, tax) %>%
  ungroup %>% # crucial so that column numbers are right if column indexing instead of names
  mutate(diff=reform - base,
         pdiff=diff / base * 100)

stat <- "pdiff"
chg4 %>% 
  select(ftype, msgroup, agirange, value=stat) %>%
  spread(ftype, value) %>%
  mutate(stat=stat, err=cps - puf) %>%
  select(msgroup, agirange, stat, puf, cps, err) %>%
  arrange(msgroup, agirange) %>%
  kable(digits=2, output=FALSE) %>%
  kableExtra::kable_styling(full_width=FALSE)

