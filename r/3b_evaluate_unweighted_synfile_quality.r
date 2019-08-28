
#****************************************************************************************************
#                Libraries ####
#****************************************************************************************************
library("magrittr")
library("plyr") # needed for ldply; must be loaded BEFORE dplyr
library("tidyverse")
options(tibble.print_max = 60, tibble.print_min = 60) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats

library("readxl") # readxl, for .xls and .xlsx files.

library("scales")

library("btools") # should not need this. If you do, it's at https://github.com/donboyd5/btools

library("knitr")

library("janitor") # for adorn_totals


#****************************************************************************************************
#                Includes ####
#****************************************************************************************************
source("./r/includes/globals_system_specific_boyd.r") # use a different version of this file if changing systems
source("./r/includes/globals_other.r")

source("./r/includes/functions_general.r")
source("./r/includes/functions_synthesis.r")


#****************************************************************************************************
#                Globals ####
#****************************************************************************************************



#****************************************************************************************************
#                Get puf.names, and puf or other files ####
#****************************************************************************************************
puf.vnames <- get_puf_vnames()
# puf.names

# sfname <- "synthpop4_all.rds"
# sfname <- "synthpop5_all.rds"
# sfname <- "synthpop6_all.rds"
sfname <- "synthpop8.rds"

# Max's files often start with synpuf not synthpop
# sfname <- "synpuf8.rds"

output <- readRDS(paste0(globals$tc.dir, "pufsyn_stack/", sfname))

names(output)
idvars <- output$idvars
xvars <- output$xvars
synvars.vs <- output$synvars.vs
visit.sequence <- output$synx$visit.sequence
# vist.sequence <- c(idvars, xvars, synvars.vs)

stack <- output$stack
glimpse(stack)
count(stack, ftype)
stack %>%
  group_by(ftype, m) %>%
  summarise(seq.min=min(pufseqn), seq.max=max(pufseqn), wt.min=min(wt), wt.max=max(wt), wage.min=min(e00200), wage.max=max(e00200))


#******************************************************************************************************************
#  utility analysis - logistic regression ####
#******************************************************************************************************************
#.. TODO: reconcile this stack with that above - I moved this from elsewhere ----
# sfname <- "synthpop5_all.rds"
# sfname <- "synthpop6_all.rds"

synprep <- readRDS(paste0(globals$tc.dir, str_remove(sfname, ".rds"), "_rwprep.rds"))

names(synprep)
# merge and then split
tcvars <- c("c00100", "taxbc") # taxcalc vars
mrgdf <- left_join(synprep$tc.base, synprep$tc.output %>% dplyr::select(RECID, tcvars))

stack <- mrgdf %>%
  # filter(MARS==2, c00100>=0, c00100<=25e3, m==1) %>%
  mutate(syn.ind=ifelse(ftype=="puf", 0, ifelse(ftype=="syn", 1, -9)))
count(stack, syn.ind, m)

glimpse(stack)
count(stack, ftype)
names(stack) %>% sort

dropvars <-  c("syn.ind", "RECID", "pufseqn", "ftype", "m", "MARS", "msname", "wt", "xfpt", "xfst", "xocah",
               "e00200p", "e00200s", "e00900p", "e00900s", "e00800", "e01500", "e02100p", "e02100s",
               "divratio", "e01500_minus_e01700")
(rhs.vars <- setdiff(names(stack), dropvars) %>% sort)

(frm <- as.formula(paste("syn.ind ~ ", paste(rhs.vars, collapse= "+"))))

mod1 <- glm(frm, data = stack %>% filter(syn.ind==0 | (syn.ind==1 & m==1)) %>% mutate(XTOT=as.factor(XTOT)), family = "binomial")
summary(mod1)
utility(mod1)
names(mod1)
mod1$coefficients
length(mod1$control)

cbind(coef(summary(mod1))[,2], coef(summary(mod1))[,4])
coefs <- coef(summary(mod1))[,2]
pvals <- coef(summary(mod1))[,4]
tibble(vname=names(coefs), coef=coefs, pval=pvals) %>%
  filter(!row_number()==1) %>%
  left_join(get_pufbase_sums() %>% rename(puf.billions=value)) %>%
  arrange(pval) %>%
  head(25) %>%
  kable(digits=c(0, 6, 6, 1, 0))


utility(mod1)


#****************************************************************************************************
#                Correlations ####
#****************************************************************************************************
corvars <- names(stack)[(names(stack) %>% str_sub(., 1, 1) %in% c("e", "p")) | names(stack)=="wt"]
cor1 <- stack %>%
  dplyr::select(ftype, corvars) %>%
  group_by(ftype) %>%
  do(cordf(.[, -1])) %>%
  do(trimcor(.)) %>%
  separate(combo, c("vname1", "vname2"), sep="-", remove=TRUE) %>%
  ungroup %>%
  left_join(puf.vnames %>% dplyr::select(vname1=vname, vdesc1=vdesc)) %>%
  left_join(puf.vnames %>% dplyr::select(vname2=vname, vdesc2=vdesc))
ht(cor1)

cor.comp <- cor1 %>%
  spread(ftype, value) %>%
  mutate(diff=syn - puf,
         adiff=abs(diff))

vsize <- 30
cor.comp %>%
  arrange(-adiff) %>%
  mutate(vdesc1=str_sub(vdesc1, 1, vsize), vdesc2=str_sub(vdesc2, 1, vsize)) %>%
  head(25) %>%
  kable(digits=3)

txt <- "wages"
txt <- "pension"
txt <- "dividend"
txt <- "contrib"
txt <- "medical"
cor.comp %>%
  filter(str_detect(vdesc1, coll(txt, ignore_case=TRUE)) |
           str_detect(vdesc2, coll(txt, ignore_case=TRUE))) %>%
  mutate(vdesc1=str_sub(vdesc1, 1, vsize), vdesc2=str_sub(vdesc2, 1, vsize)) %>%
  arrange(-adiff) %>%
  head(25) %>%
  kable(digits=3)

# we have not captured negative corrs between 
# e01000 and e02000, e26270
# e01000 seems to be the problem

#****************************************************************************************************
#                Measures of fits of individual variables ####
#****************************************************************************************************
comp <- stack %>%
  dplyr::select(ftype, c(idvars, synvars.vs)) %>%
  gather(vname, value, -ftype)
# glimpse(comp)
# count(comp, vname)

# don't use Anderson-Darling test - it takes forever

comp.ks <- comp %>%
  group_by(vname) %>%
  summarise(ks.p=ks.p(ftype, value), ks.D=ks.D(ftype, value)) %>%
  left_join(puf.vnames %>% dplyr::select(vname, vdesc)) %>%
  arrange(ks.p) # higher ks.p is better; ks.D lower is better

comp.ks %>% kable(digits=4) # dividend-related variables are bad, others are much better
comp.ks %>% select(-ks.p) %>% arrange(desc(ks.D)) %>% top_n(25) %>% kable(digits=4) # dividend-related variables are bad, others are much better

# https://www.researchgate.net/publication/269355453_The_multisample_Cucconi_test

# create bins for each variable ----
# names(stack)
dropvars <- c("RECID", "pufseqn", "wt", "syn.ind", "MARS", "XTOT", "divratio", "ftype", "m", "msname", "e01500_minus_e01700",
              "e00200p", "e00200s", "e00900p", "e00900s", "e02100p", "e02100s")
vars <- setdiff(names(stack), dropvars)
vars
nbins <- 1000
bins <- ldply(vars, getbins, stack, nbins, .progress="text")
glimpse(bins)

binstats <- bins %>%
  group_by(vname) %>%
  summarise(r2=rsq(puf_pct, syn_pct),
            sse=sse(puf_pct, syn_pct)) %>%
  left_join(puf.vnames %>% dplyr::select(vname, vdesc))

binstats %>%
  arrange(desc(sse)) %>%
  kable(digits=6)
# we seem to be giving e00650 in syn to bins that don't have it in puf, at least for low amounts

# compare the shares of a variable that fall in different ntiles  
bins %>%
  filter(vname=="e00650") %>%
  filter(puf_pct < 100) %>%
  ggplot(aes(puf_pct, syn_pct)) +
  geom_point(colour="blue") +
  scale_x_continuous(breaks=seq(0, 100, 5), limits=c(0, 65)) +
  scale_y_continuous(breaks=seq(0, 100, 5), limits=c(0, 65)) +
  geom_abline(slope=1, intercept = 0) +
  theme_bw() +
  ggtitle("puf vs. syn relative frequencies, using 1,000 puf-based bins",
          subtitle="e00650: Qualified dividends")

ulim <- 1
bins %>%
  filter(vname=="e00650") %>%
  filter(puf_pct <= ulim) %>%
  ggplot(aes(puf_pct, syn_pct)) +
  geom_point(colour="blue") +
  scale_x_continuous(breaks=seq(0, 100, .1), limits=c(0, ulim)) +
  scale_y_continuous(breaks=seq(0, 100, .1), limits=c(0, ulim)) +
  geom_abline(slope=1, intercept = 0) +
  theme_bw() +
  ggtitle("puf vs. syn relative frequencies where puf_pct <= 1, using 1,000 puf-based bins",
          subtitle="e00650: Qualified dividends")

ulim <- 1
bins %>%
  filter(vname=="e00650") %>%
  filter(puf_pct <= ulim) %>%
  mutate(ratio=syn_pct / puf_pct) %>%
  ggplot(aes(bin.num, ratio)) +
  geom_point(colour="blue") +
  geom_hline(yintercept = 1) +
  scale_x_continuous(breaks=seq(0, 1000, 20)) +
  scale_y_continuous(breaks=seq(0, 100, .05)) +
  theme_bw() +
  ggtitle("Ratio of syn to puf relative frequencies where puf_pct <= 1, using 1,000 puf-based bins",
          subtitle="e00650: Qualified dividends")

bins %>%
  filter(vname=="e00650") %>%
  filter(puf_pct > 50) %>%
  kable(digits=2, format.args=list(big.mark = ','))

bins %>%
  filter(vname=="e00650") %>%
  filter(puf_share ==0 ) %>%
  mutate_at(vars(puf_share, syn_share), funs(. * 100)) %>%
  kable(digits=2, format.args=list(big.mark = ','))

bins %>%
  filter(vname=="e00650") %>%
  mutate_at(vars(puf_share, syn_share), funs(. * 100)) %>%
  mutate(adiff=abs(syn_share - puf_share)) %>%
  arrange(-adiff) %>%
  head(25) %>%
  kable(digits=2, format.args=list(big.mark = ','))

bins %>%
  filter(vname=="e00650") %>%
  filter(puf_share==0 ) %>%
  summarise(syn=sum(syn), syn_share=sum(syn_share))



# combine the summary measures and bring in measure of total value
varstats <- left_join(comp.ks %>% dplyr::select(vname, ks.D),
                      binstats %>% dplyr::select(vname, sse)) %>%
  left_join(get_pufbase_sums() %>% rename(puf.billions=value)) %>%
  dplyr::select(vname, sse, ks.D, puf.billions, vdesc) %>%
  arrange(desc(sse))

varstats %>% 
  head(25) %>%
  kable(digits=c(0, 4, 3, 1, 0))


#****************************************************************************************************
#                Quantiles of individual variables ####
#****************************************************************************************************
var <- "e00650"
var <- "e00600"
var <- "divratio"
var <- "e26270"
stack %>%
  group_by(ftype) %>%
  do(qtiledf(.[[var]], probs=c(0, .01, .05, seq(.1, .8, .1), .85, .9, .95, .99, 1)))


#****************************************************************************************************
#                CDF plots of individual variables ####
#****************************************************************************************************
names(stack)
cdfplot.unwtd("wt")
cdfplot.unwtd("e00200")
cdfplot.unwtd("p04470")
cdfplot.unwtd("e04600")
cdfplot.unwtd("e01500")
cdfplot.unwtd("e01500_minus_e01700")
cdfplot.unwtd("e01700")
cdfplot.unwtd("e02400")
cdfplot.unwtd("e01000") # Net cap gain loss, needs work! MAYBE THE PROBLEM IS NEGATIVE VALUES??
cdfplot.unwtd("divratio")
cdfplot.unwtd("e00650") # bad on high end -- Qualified Dividends
# cdfplot.unwtd("e00600_minus_e00650")
cdfplot.unwtd("e00600")

var <- "e01400"
cdfplot.unwtd(var)

# variables that need work:
# e00900
# e01500
# e00600


#******************************************************************************************************************
#  regression analysis of puf vs syn ####
#******************************************************************************************************************
library("broom")
glimpse(stack)

psum <- get_pufbase_sums()
psum
psum %>% arrange(vname)

frm <- formula(e00200 ~ e00300 + e00600 + e00700 + e01500 + e02000 + e02100 + e02300 + e02400 + e02500)

frm <- formula(e00200 ~ e17500 + e18400 + e18500 + e19200 + e04600 + e19700)


# m1p <- lm(frm, data=stack %>% filter(ftype=="puf"))
# m12 <- lm(frm, data=stack %>% filter(ftype=="syn"))

mods <- stack %>%
  group_by(ftype) %>%
  do(lm(frm, data=.) %>% tidy)

mods %>%
  arrange(term, ftype)

stats <- stack %>%
  group_by(ftype) %>%
  do(lm(frm, data=.) %>% glance)

stats

# now do it by sample
mods <- stack %>%
  group_by(ftype, m) %>%
  do(lm(frm, data=.) %>% tidy)

mods %>%
  arrange(term, ftype, m)

