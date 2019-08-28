

# Shirakawa, Kiyomi, and Shinsuke Ito. "Creating Synthetic Microdata from Official Statistics: Random Number Generation in Consideration of Anscombe's Quartet," 2015, 14.



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
library("precis")

library("tibbletime") # https://business-science.github.io/tibbletime/

library("grDevices")
library("knitr")

library("zoo") # for rollapply

library("btools") # library that I created (install from github)
library("bdata")

library("mvtnorm")
library("portopt")

library("synthpop")

vignette(package = "synthpop")

#****************************************************************************************************
#                Get data ####
#****************************************************************************************************
fn <- "shirakawa.xlsx"

samp <- read_excel(paste0("./data/", fn),
                  sheet = "shirakawa_ito_tab5",
                  range = "A3:K23")
samp

basic <- read_excel(paste0("./data/", fn),
                    sheet = "shirakawa_ito_tab6",
                    range = "b2:e10")
basic

dtl <- read_excel(paste0("./data/", fn),
                  sheet = "shirakawa_ito_tab7",
                  range = "b3:h9")
dtl



#****************************************************************************************************
#                Take a look ####
#****************************************************************************************************
summary(samp)

# reproduce table 6, basic table
meas.order <- c("mean", "sd", "kurtosis", "skewness", "freq")
vars.order <- c("livexp", "food", "housing")
h1 <- samp %>%
  summarise_at(vars(livexp, food, housing),
               funs(mean(., na.rm=TRUE), 
                    sd(., na.rm=TRUE), 
                    kurtosis=e1071::kurtosis(., type=2, na.rm=TRUE),
                    skewness=e1071::skewness(., type=2, na.rm=TRUE),
                    freq=sum(!is.na(.)))) %>%
  gather(variable, value) %>%
  separate(variable, into=c("variable", "col1")) %>%
  mutate(col1=factor(col1, levels=meas.order),
         variable=factor(variable, levels=vars.order)) %>%
  spread(variable, value) %>%
  arrange(col1) %>%
  mutate(col1=as.character(col1))
h1

h2 <- samp %>% 
  select(livexp, food, housing) %>%
  cor(use="pairwise.complete.obs")
h2a <- data.frame(col1=paste0(rownames(h2), "_cor"), h2) %>%
  as.tibble %>%
  mutate(col1=as.character(col1))
bind_rows(h1, h2a) %>% kable(digits=4)
basic %>% kable(digits=4)

# reproduce table 7, details table
samp %>%
  group_by(groupno) %>%
  summarise_at(vars(livexp, food), funs(n(), mean, sd)) %>%
  select(groupno, starts_with("live"), starts_with("food"))
dtl


#****************************************************************************************************
#                Box Cox ####
#****************************************************************************************************

# - generate random mvt with mean sd from basic table

# sigma <- matrix(c(4,2,2,3), ncol=2)
# sigma
# x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma)
# x[1:10, ]
# colMeans(x)
# var(x)

# get correlation matrix and sd, create covariance matrix
mean.target <- basic  %>%
  filter(measure=="mean") %>%
  select(-measure) %>%
  unlist(., use.names=FALSE)

cormatdf <- basic %>%
  filter(str_detect(measure, coll("_cor")))

cormat <- cormatdf %>%
  select(-measure) %>%
  as.matrix
cormat
cormat[upper.tri(cormat)] <- cormat[lower.tri(cormat)]
cormat

sdvec <- basic %>%
  filter(measure=="sd") %>%
  select(-measure) %>%
  unlist(., use.names=FALSE)

sigma.target <- covmat(cormat, sdvec)

sigma.target
mean.target

x <- rmvnorm(n=20, mean=mean.target, sigma=sigma.target)
x

#****************************************************************************************************
#                Get puf data, create an excerpt ####
#****************************************************************************************************
pufraw <- readRDS("d:/Dropbox/RPrograms PC/OSPC/cps_puf_compare/data/rds/pufraw.rds")
names(pufraw)
glimpse(pufraw)

pufx <- pufraw %>%
  setNames(str_to_lower(names(.))) %>%
  select(mars, xtot, wt=s006, agi=e00100, wages=e00200, interest=e00300, div=e00600) %>%
  mutate(wt=wt / 100)
summary(pufx) # no NAs
# 1 agi 2 wages 3 interest 6 dividends; agi has neg values but 2, 3, 6 do not

quantile(pufx$wt, probs = seq())

pufx2 <- pufx %>%
  mutate(otherinc=agi - wages - interest - div)
summary(pufx2)

count(pufx2, mars, xtot)
# mars
# 0 aggregated return -- only 4 of these
# 1 single, 2 mfj or widower/child, 3 mfs, 4 hoh
# xtot ranges 0-5

# let's pick mars 2 xtot2 41.6k recs
pufx2 %>%
  filter(mars==2, xtot==2) %>%
  mutate(igroup=cut(agi, 4)) %>%
  group_by(igroup) %>%
  summarise(n=n(), agimin=min(agi), agimdn=median(agi), agimax=max(agi))

pufx3 <- pufx2 %>%
  filter(mars==2, xtot==2)
quantile(pufx3$agi, probs=seq(0, 1, .1))
quantile(pufx3$wt, probs=seq(0, 1, .1))

brks <- c(-1e9, 0, 25e3, 50e3, 75e3, 100e3, 200e3, 500e3, 1e9)
pufx3 <- pufx3 %>%
  mutate(igroup=cut(agi, brks))
count(pufx3, igroup)

pufx3 %>%
  group_by(igroup) %>%
  do(qtiledf(.$wt, c(0, .1, .25, .5, .75, .9, 1)))
# only a few recs have wts less than 1400 in the middle area


pufx3 %>%
  filter(wt>=1000) %>%
  group_by(igroup) %>%
  summarise_at(vars(agi, wages),
               funs(mean(., na.rm=TRUE), 
                    sd(., na.rm=TRUE), 
                    kurtosis=e1071::kurtosis(., type=1, na.rm=TRUE),
                    skewness=e1071::skewness(., type=1, na.rm=TRUE),
                    freq=sum(!is.na(.))))

pufx3 %>%
  group_by(igroup) %>%
  summarise_at(vars(agi, wages, interest, div, otherinc),
               funs(kurtosis=e1071::kurtosis(., type=1, na.rm=TRUE)))
# kurtosis really falls from >50-75k to >75-100k
levels(pufx3$igroup)[1]

igrp <- c("(5e+04,7.5e+04]", "(7.5e+04,1e+05]")
gdist <- function(levgrp) {
  pufx3 %>%
    filter(wt>=1000) %>%
    filter(igroup %in% levels(igroup)[levgrp]) %>%
    ggplot(aes(wages, colour=igroup)) + geom_density()
}
# agi group and wages dist
gdist(1) # neg agi: wages near 0
gdist(2) # 0-25k a bit bimodal, wages near 0 and 15k
gdist(3) # 25-50k bimodal 0, 35k
gdist(4) # 50-75k bi 0, 65k
gdist(5) # 75-100k, bi 0, 80k
gdist(6) # 100-200k bi 0 110k
gdist(7) # 200-500k bi 


#****************************************************************************************************
#                pufx3 from synthpop ####
#****************************************************************************************************
library("synthpop")
glimpse(pufx3)
pufx3.base <- pufx3 %>% dplyr::select(-igroup)
seed <- 1234
system.time(synx3 <- syn(pufx3.base, seed = seed)) # 22 secs!
system.time(synx3.big <- syn(pufx3.base, seed = seed, k=nrow(pufx3.base) * 2)) # 32 secs!
system.time(synx3.big <- syn(pufx3.base, seed = seed, k=nrow(pufx3.base) * 4)) # 51 secs!
system.time(synx3.big <- syn(pufx3.base, seed = seed, k=nrow(pufx3.base) * 10)) # 115 secs!
system.time(synx3.big <- syn(pufx3.base, seed = seed, m = 4)) # 86 secs!


pufx3.syn <- synx3$syn
glimpse(pufx3.syn)

summary(pufx3.syn)
summary(pufx3)
pufx3.syn <- pufx3.syn %>%
  mutate(igroup=cut(agi, brks))

gdist.s <- function(levgrp) {
  pufx3.syn %>%
    filter(wt>=1000) %>%
    filter(igroup %in% levels(igroup)[levgrp]) %>%
    ggplot(aes(wages, colour=igroup)) + geom_density()
}
gdist.s(1)
gdist.s(2)
gdist.s(3)
gdist.s(4)
gdist(4)
gdist.s(5)

taxmod <- function(df){

  # tibble(ti=seq(-1e3, 600e3, 1e3), tax=taxcalc(ti), etr=ifelse(ti>0, tax / ti, 0)) %>%
  #   ggplot(aes(ti, etr)) + geom_line() + geom_point()
  df <- df %>%
    mutate(tax=taxcalc(agi))
  return(df)
}

taxcalc <- function(ti){
  tax <- case_when(
    ti <= 0 ~ 0,
    ti <= 50e3 ~ .02 * ti,
    ti <= 100e3 ~ 1000 + .04 * (ti - 50e3),
    ti <= 200e3 ~ 3000 + .06 * (ti - 100e3),
    ti <= 500e3 ~ 9000 + .08 * (ti - 200e3),
    TRUE ~ 33000 + .10 * (ti - 500e3))
  return(tax)
}
brks <- c(-1e9, 0, 25e3, 50e3, 75e3, 100e3, 200e3, 500e3, 1e9)
pufx3.all <- bind_rows(pufx3.syn %>% mutate(type="syn", wt=wt / 4),
                       pufx3.base %>% mutate(type="base")) %>%
  mutate(agi2=wages + interest + div + otherinc,
         igroup=cut(agi2, brks),
         tax=taxcalc(agi2))

pufx3.sum <- pufx3.all %>%
  group_by(type, igroup) %>%
  summarise(n=n(), wt.sum=sum(wt), agi2.sum=sum(agi2), tax.sum=sum(tax), 
            agi2m.wtd=sum(wt * agi2) / 1e6,
            taxm.wtd=sum(wt * tax) / 1e6)

# var <- "n"
comp <- function(var){
  pufx3.sum %>%
    dplyr::select(type, igroup, !!var) %>%
    spread(type, !!var) %>%
    mutate(pdiff=syn / base * 100 - 100)
}
comp("n")
comp("wt.sum")  
comp("agi2.sum")
comp("tax.sum")
comp("agi2m.wtd")
comp("taxm.wtd")


#****************************************************************************************************
#                Kurtosis checks ####
#****************************************************************************************************
# how best to calculate kurtosis and skewness?
set.seed(1234)
n <- 100e3
x <- rnorm(n)

# e1071 seems more sophisticated in that it allows alternative calculations and provides better technical documentation

#..kurtosis ----
# formula for population kurtosis
k.p <- {n * sum((x - mean(x)) ^4)} /
  {sum((x - mean(x))^2) ^2}
k.p - 3 # excess
# same as moments, and e1071 type=1

# sample kurtosis
k.s <- {n * (n+1) * (n-1) * sum((x - mean(x)) ^4)} /
  {(n-2) * (n-3) * sum((x - mean(x))^2) ^2}
k.s - 3
# same as PerformanceAnalytics method="sample"

# moments calculates kurtosis, e1071 calculates excess kurtosis vs normal (subtracts 3)
moments::kurtosis(x) - 3 # same as population kurtosis

e1071::kurtosis(x)
e1071::kurtosis(x, type=1) # matches moments::kurtosis; e1071 doc says this is typical definition in older texts
e1071::kurtosis(x, type=2) # SAS and SPSS method; unbiased under normality
e1071::kurtosis(x, type=3) # default - MINITAB and BMDP

timeDate::kurtosis(x) # default is excess, same as e1071 type=3; options are method = c("excess", "moment", "fisher")
timeDate::kurtosis(x, method="excess") # same as default (e1071 type=3)
timeDate::kurtosis(x, method="moment") -3 # same result as excess after subtracting 3
timeDate::kurtosis(x, method="fisher") # doesn't match any e1071 calcs, matches PerformanceAnalytics

PerformanceAnalytics::kurtosis(x) # default is same as e1071 type=1
PerformanceAnalytics::kurtosis(x, method="excess") # same as default
PerformanceAnalytics::kurtosis(x, method="moment") - 3  # same as default after subtracting 3
PerformanceAnalytics::kurtosis(x, method="fisher") # matches timeDate fisher
PerformanceAnalytics::kurtosis(x, method="sample") - 3  # same as sample formula
PerformanceAnalytics::kurtosis(x, method="sample_excess") # same as e1071 type=2

??moments::kurtosis
??e1071::kurtosis
??timeDate::kurtosis
??PerformanceAnalytics::kurtosis


#..skewness ----
moments::skewness(x)
e1071::skewness(x)
e1071::skewness(x, type=1) # matches moments::skewness; e1071 doc says this is typical definition in older texts
e1071::skewness(x, type=2) # SAS and SPSS method; unbiased under normality
e1071::skewness(x, type=3) # default - MINITAB and BMDP
??moments::skewness
??e1071::skewness



#****************************************************************************************************
#                Try some puf things ####
#****************************************************************************************************
ns <- function(df) {names(df) %>% sort}
df <- read_csv(paste0("D:/Google Drive/synpuf/puf/", "puf2011.csv"), 
               col_types = cols(.default= col_double()), 
               n_max=-1)
ns(df)
vars <- c("RECID", "E00100", "E00200", "E00300", "E00400", "E00600", "E01500")
hipuf <- df %>%  
  filter(MARS==1, E00100 >= 200e3, !RECID %in% 999996:999999) %>%
  dplyr::select(vars)
ht(hipuf)

stats <- hipuf %>%
  gather(vname, value, -RECID) %>%
  group_by(vname) %>%
  summarise(n=n(),
            nz=sum(value==0),
            nnz=n - nz,
            mean=mean(value), 
            sd=sd(value), 
            kurt=e1071::kurtosis(value, type=2),
            skew=e1071::skewness(value, type=2))
stats

var <- "E00200"
hipuf %>%
  dplyr::select(RECID, value=var) %>%
  ggplot(aes(value)) +
  geom_density(size=1.5, colour="blue") +
  theme_bw()


tur <- c(1.0, 1.2, 1.1, 1.1, 2.4, 2.2, 2.6, 4.1, 5.0, 10.0, 4.0, 4.1, 4.2, 4.1, 5.1, 4.5, 5.0, 15.2, 10.0, 20.0, 1.1, 1.1, 1.2, 1.6, 2.2, 3.0, 4.0, 10.5)
box = boxcox(tur ~ 1,              # Transform Turbidity as a single vector
             lambda = seq(-6,6,0.1))      # Try values -6 to 6 by 0.1
str(box)
lambda <- tibble(lambda=box$x, loglik=box$y) %>%
  top_n(1, loglik) %>% # get the largest log likelihood
  .[["lambda"]]

data.frame(x = c(10, 4, 1, 6, 3, 1, 1)) %>% top_n(2)

bc <- tibble(rn=1:length(tur), tur, bctur=(tur ^ lambda - 1)/lambda)
bc  %>%
  gather(variable, value, -rn) %>%
  ggplot(aes(value, colour=variable)) +
  geom_density(size=1.5) +
  theme_bw()

getlam <- function(x){
  box <- boxcox(x ~ 1, lambda = seq(-6, 6, 0.1))
  lambda <- tibble(lambda=box$x, loglik=box$y) %>%
    top_n(1, loglik) %>% # get the largest log likelihood
    .[["lambda"]]
  return(lambda)
}
bc <- hipuf %>%
  # dplyr::select(RECID, E00200) %>%
  gather(vname, value, -RECID) %>%
  mutate(value=ifelse(value <=0, 1, value)) %>%
  group_by(vname) %>%
  summarise(n=n(),
            nz=sum(value==0),
            nnz=n - nz,
            mean=mean(value), 
            sd=sd(value), 
            kurt=e1071::kurtosis(value, type=2),
            skew=e1071::skewness(value, type=2),
            lambda=getlam(value))
bc
