
source(file.path(PROJHOME, "r/includes", "libraries.r"))
source(file.path(PROJHOME, "r/includes", "globals_system_specific_boyd.r"))
source(file.path(PROJHOME, "r/includes", "globals_other.r"))

mse <- function(y, yhat){
  sqerr <- (yhat - y)^2
  return(mean(sqerr))
}

rmse <- function(y, yhat){
  rmse <- sqrt(mse(y, yhat))
  return(rmse)
}


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
  select(synid=id_A, pufid=id_B1, dist1)


glimpse(puf)
glimpse(syn)
glimpse(nearest)

names(puf) %>% sort

# let's do this like they did for sipp - for each true rec, get the corresponding syn rec
# if a puf rec is matched to 2 or more syn recs, get just the first syn rec
near2 <- nearest %>%
  arrange(pufid, synid) %>%
  group_by(pufid) %>%
  arrange(synid) %>%
  filter(row_number()==1) %>%
  ungroup


vars <- c("E00200", "E00300", "E00600", "E01500", "E18400", "P23250")
puflong <- near2 %>%
  select(pufid, synid, dist1) %>%
  left_join(puf %>% select(pufid=id, vars)) %>%
  gather(vname, puf, vars) %>%
  arrange(pufid, synid, vname)

synlong <- near2 %>%
  select(pufid, synid) %>%
  left_join(syn %>% select(synid=id, vars)) %>%
  gather(vname, syn, vars) %>%
  arrange(pufid, synid, vname)
glimpse(synlong)

matched <- left_join(puflong, synlong)
head(matched, length(vars) * 2)

varstats <- matched %>%
  group_by(vname) %>%
  summarise(sd=sd(puf), rmse=rmse(puf, syn), ratio=rmse / sd)
varstats

# now get the average across all (included) variables
varstats %>%
  summarise(avg=mean(ratio))


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


