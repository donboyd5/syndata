

#****************************************************************************************************
#                Includes ####
#****************************************************************************************************
source(file.path(PROJHOME, "r/includes", "libraries.r"))
# library("readxl")
library("knitr")

source(file.path(PROJHOME, "r/includes", "globals_system_specific_boyd.r"))
source(file.path(PROJHOME, "r/includes", "globals_other.r"))
source(file.path(PROJHOME, "r/includes", "functions_general.r"))

search() %>% sort

#****************************************************************************************************
#                functions ####
#****************************************************************************************************
mse <- function(y, yhat){
  sqerr <- (yhat - y)^2
  return(mean(sqerr))
}

rmse <- function(y, yhat){
  rmse <- sqrt(mse(y, yhat))
  return(rmse)
}

ns <- function(df) {names(df) %>% sort}


#****************************************************************************************************
#                Get data ####
#****************************************************************************************************
puf <- read_csv(paste0("D:/Google Drive/synpuf/puf/", "puf2011.csv"), 
                col_types = cols(.default= col_double()), 
                n_max=-1) %>%
  mutate(ftype="puf", id=row_number() - 1)
glimpse(puf)

syn <- read_csv(paste0("D:/Google Drive/synpuf/syntheses/", "synpuf20.csv"), 
                col_types = cols(.default= col_double()), 
                n_max=-1) %>%
  mutate(ftype="syn", id=row_number() - 1)
glimpse(syn)

nearest <- read_csv(paste0("D:/Google Drive/synpuf/nearest/", "nearest20.csv"), 
                    col_types = cols(.default= col_double()), 
                    n_max=-1) %>% 
  dplyr::select(synid=id_A, pufid=id_B1, dist1)


glimpse(puf)
glimpse(syn)
glimpse(nearest)

ns(syn)
ns(puf)

#****************************************************************************************************
#                define vars ####
#****************************************************************************************************
(keepvars <- setdiff(names(syn), c("RECID", "S006")))
idvars <- c("ftype", "id")
catvars <- c("DSI", "EIC", "F2441", "F6251", "FDED", "MARS", "MIDR", "N24", "XTOT") # categorical variables
(contvars <- setdiff(keepvars, c(idvars, catvars))) # 55 vars


#****************************************************************************************************
#                prep data ####
#****************************************************************************************************
puf2 <- puf %>% filter(!RECID %in% 999996:999999) %>% dplyr::select(keepvars)
syn2 <- syn %>% dplyr::select(keepvars)

# let's do this like they did for sipp - for each true rec, get the corresponding syn rec
# if a puf rec is matched to 2 or more syn recs, get just the first syn rec
near2 <- nearest %>%
  arrange(pufid, synid) %>%
  group_by(pufid) %>%
  arrange(synid) %>%
  filter(row_number()==1) %>%
  ungroup
anyDuplicated(near2$pufid)
anyDuplicated(near2$synid)
# good no dups

#****************************************************************************************************
#                quick look ####
#****************************************************************************************************
# Skewness -- lopsidedness:
#  + many low or zero values, distribution has hump on left and tails to right
#  - distribution tails off to the left
# Kurtosis -- heaviness of tails measured relative to normal distribution same variance
#  + relatively more outliers than normal distribution
#  - relatively fewer outliers
stats <- puf2 %>%
  dplyr::select(id, contvars) %>%
  gather(vname, value, -id) %>%
  group_by(vname) %>%
  summarise(n=n(),
            n.NA=sum(is.na(value)),
            mean=mean(value, na.rm=TRUE),
            median=median(value, na.rm=TRUE), 
            sd=sd(value, na.rm=TRUE),
            kurtosis=e1071::kurtosis(value, type=2, na.rm=TRUE),
            skewness=e1071::skewness(value, type=2, na.rm=TRUE))

# most variables highly skewed with hump left
stats %>% 
  as.data.frame %>% 
  arrange(skewness) %>%
  kable(digits=2, format.args = list(big.mark=",")) %>%
  kable_styling(full_width = FALSE)

# most variables heavier tails than normal
stats %>% 
  as.data.frame %>% 
  arrange(kurtosis) %>%
  kable(digits=2, format.args = list(big.mark=",")) %>%
  kable_styling(full_width = FALSE)


#****************************************************************************************************
#                create long file ####
#****************************************************************************************************
# get lowest distance syn for each puf record
# create long file with pufid, synid, puvval, synval, for each value
# then compute measures

puflong <- near2 %>%
  dplyr::select(pufid, synid, dist1) %>%
  left_join(puf2 %>% dplyr::select(pufid=id, contvars)) %>%
  gather(vname, puf, -pufid, -synid, -dist1) %>%
  arrange(pufid, synid, vname)
glimpse(puflong)

synlong <- near2 %>%
  dplyr::select(pufid, synid) %>%
  left_join(syn2 %>% dplyr::select(synid=id, contvars)) %>%
  gather(vname, syn,  -pufid, -synid) %>%
  arrange(pufid, synid, vname)
glimpse(synlong)

matched <- left_join(puflong, synlong)
head(matched, length(contvars) * 2) %>% as.data.frame

#****************************************************************************************************
#                calculate sipp value ####
#****************************************************************************************************
tcalc <- function(x){
  # the SIPP transformation measure per Gary Benedetto email
  # the log(total net worth) when total net worth>1,
  # the log(absolute value(total net worth)) when total net worth<-1
  # 0 otherwise
  ifelse((x < -1) | (x > 1), log(abs(x)), 0)
}
# tcalc(c(-10, -1.1, -1, -.5, 0, .5, 1, 1.1, 10))

# this takes a while
matched2 <- matched %>%
  gather(ftype, rawvalue, puf, syn) %>%
  mutate(tvalue=tcalc(rawvalue)) %>%
  gather(vartype, value, rawvalue, tvalue) %>%
  spread(ftype, value)
ht(matched2)

varstats <- matched2 %>%
  group_by(vname, vartype) %>%
  summarise(sd=sd(puf), rmse=rmse(puf, syn), ratio=rmse / sd)
varstats

# now get the average across all (included) variables
varstats %>%
  group_by(vartype) %>%
  summarise(avg=mean(ratio))


#****************************************************************************************************
#                create a table like their Table 1 ####
#****************************************************************************************************
# P23250 is LT gains less losses, pretty big (~$51b)
# do something like their table 1, but with wages
vars <- c("E00200", "P23250")
varstats %>%
  filter(vname %in% vars)


probs <- c(0, .25, .5, .75, .9, .95, .99, 1)
vprobs <- matched2 %>%
  filter(vname %in% vars) %>%
  group_by(vname, vartype) %>%
  do(qtiledf(.$puf, probs=probs))
vprobs

tab1dat <- matched2 %>%
  filter(vname %in% vars, vartype=="rawvalue", vname=="P23250", puf>=0) %>%
  group_by(vname) %>%
  arrange(puf) %>%
  mutate(ptile=percent_rank(puf),
         qtile=cut(ptile, c(-Inf, .75, .9, .95, .99, 1)))
ht(tab1dat)
count(tab1dat, qtile)

quantile(matched2 %>% filter(vartype=="rawvalue", vname=="P23250") %>% .[["puf"]],
         probs=c(0, .25, .5, .60, .70, .72, .73, .74, .75, .9, .95, .99, 1))

tab1dat %>%
  group_by(qtile) %>%
  summarise(n=n(), minp=min(ptile), maxp=max(ptile))

tab1dat %>%
  group_by(vname, qtile) %>%
  summarise(n=n(), 
            mean=mean(puf), 
            sd=sd(puf), 
            rmse=rmse(puf, syn),
            ratio=rmse / sd)




#****************************************************************************************************
#                analyze ####
#****************************************************************************************************

# do something like their table 1, but with wages
(ptiles <- quantile(matched$puf[matched$vname=="E00200"], probs=c(0, .25, .5, .75, .9, .95, .99, 1)))

tab1dat <- matched %>%
  filter(vname=="E00200") %>%
  mutate(ptile=percent_rank(puf),
         qtile=cut(ptile, c(-Inf, 0, .25, .5, .75, .9, .95, .99, 1)))
count(tab1dat, qtile)

tab1dat %>%
  group_by(qtile) %>%
  summarise(n=n(), 
            mean=mean(puf), 
            sd=sd(puf), 
            rmse=rmse(puf, syn),
            ratio=rmse / sd)


