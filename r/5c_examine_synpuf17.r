# synpuf17 seems to be turning out really badly. try to figure out why.


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

# devtools::install_github("donboyd5/btools")
library("btools") # library that I created (install from github)

library("ipoptr")
library("nloptr")

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
source("./r/includes/functions_target_setup_and_analysis.r")
source("./r/includes/functions_ipopt.r")

# functions specific to the weighting from scratch approach:
source("./r/includes/functions_weight_from_scratch.r")


#****************************************************************************************************
#                Get puf variable names ####
#****************************************************************************************************
puf.vnames <- get_puf_vnames()


#****************************************************************************************************
#                debug ####
#****************************************************************************************************
#.. first, track the trail to make sure data are ok ----
sfname <- "synpuf17"

syn.from.max <- read_csv(paste0(globals$synd, "syntheses/", sfname, ".csv"), 
                         col_types = cols(.default= col_double()), 
                         n_max=-1)
glimpse(syn.from.max)

# I expect wages of 6581.0
sum(syn.from.max$E00200 * syn.from.max$S006 / 100) / 1e6 / 5
# good 6580952

#.. now look at the stacked file I created in 3a_run_tax_plans, before reweighting ----
df <- readRDS(paste0(globals$tc.dir, sfname, "_rwprep.rds"))
names(df)

tcvars <- c("c00100", "taxbc")
mrgdf <- left_join(df$tc.base, df$tc.output %>% dplyr::select(RECID, tcvars))
glimpse(mrgdf)
count(mrgdf, ftype)

# do I get the same wage value using wt.syn?
mrgdf %>%
  filter(ftype=="syn") %>%
  summarise(wages=sum(e00200 * wt.syn / 5) / 1e6) # good 6580952

# Can I reproduce the prior agi table, shown below? (based on puf and syn only)

#   |agirange        |     puf|     syn|    wtfs| syn_diff| wtfs_diff| syn_pdiff| wtfs_pdiff|vname  |vdesc                              |
#   |:---------------|-------:|-------:|-------:|--------:|---------:|---------:|----------:|:------|:----------------------------------|
#   |[-Inf,0)        |   -59.0|  -119.5|   -59.0|    -60.6|       0.0|     102.7|        0.0|c00100 |Adjusted gross income (calculated) |
#   |[0,2.5e+04)     |   717.6|   976.7|   655.9|    259.1|     -61.7|      36.1|       -8.6|c00100 |Adjusted gross income (calculated) |
#   |[2.5e+04,5e+04) | 1,242.7| 1,004.1| 1,115.5|   -238.5|    -127.1|     -19.2|      -10.2|c00100 |Adjusted gross income (calculated) |
#   |[5e+04,7.5e+04) | 1,172.7| 1,297.5| 1,082.0|    124.9|     -90.7|      10.6|       -7.7|c00100 |Adjusted gross income (calculated) |
#   |[7.5e+04,1e+05) | 1,037.7| 1,104.5|   990.8|     66.8|     -46.9|       6.4|       -4.5|c00100 |Adjusted gross income (calculated) |
#   |[1e+05,2e+05)   | 1,984.8| 1,916.2| 1,946.7|    -68.6|     -38.2|      -3.5|       -1.9|c00100 |Adjusted gross income (calculated) |
#   |[2e+05,5e+05)   | 1,101.8| 1,179.0| 1,077.7|     77.2|     -24.0|       7.0|       -2.2|c00100 |Adjusted gross income (calculated) |
#   |[5e+05,1e+06)   |   417.1|   545.8|   393.1|    128.7|     -24.0|      30.9|       -5.8|c00100 |Adjusted gross income (calculated) |
#   |[1e+06,1e+07)   |   654.3|   763.6|   642.2|    109.3|     -12.1|      16.7|       -1.8|c00100 |Adjusted gross income (calculated) |
#   |[1e+07, Inf)    |   226.7|   298.6|   228.5|     71.9|       1.8|      31.7|        0.8|c00100 |Adjusted gross income (calculated) |
#   |Total           | 8,496.4| 8,966.4| 8,073.4|    470.1|    -423.0|       5.5|       -5.0|-      |-                                  |

# summaries by income range
agiranges <- c(-Inf, 0, 25e3, 50e3, 75e3, 100e3, 200e3, 500e3, 1e6, 10e6, Inf)
vlist <- c("c00100", "e00200", "e00300", "e00600", "e01700", "p23250", "taxbc")

dfsums <- mrgdf %>%
  mutate(wt=ifelse(ftype=="puf", wt.puf, wt.syn / 5)) %>%
  mutate(agirange=cut(c00100, agiranges, right=FALSE),
         wtone=1e9) %>%
  dplyr::select(ftype, agirange, wt, wtone, vlist) %>%
  gather(vname, value, -ftype, -agirange, -wt) %>%
  group_by(ftype, agirange, vname) %>%
  summarise(n=n(), wtsum.m=sum(wt) / 1e6, valsum.b=sum(wt * value) / 1e9) %>%
  left_join(puf.vnames %>% dplyr::select(vname, vdesc))
dfsums

dfsums %>%
  filter(vname=="c00100") %>%
  dplyr::select(ftype, agirange, stat="valsum.b", vname, vdesc) %>%
  spread(ftype, stat) %>%
  janitor::adorn_totals(where="row") %>%
  kable(digits=c(0, rep(1, 5), rep(1, 2)), format.args=list(big.mark = ','))
# yes -- good

# now, wages


# and how about taxbc for mgroup 2
#   |mgroup |agirange        |    puf|    syn|   wtfs| syn_diff| wtfs_diff| syn_pdiff| wtfs_pdiff|stat     |vname |vdesc                          |
#   |:------|:---------------|------:|------:|------:|--------:|---------:|---------:|----------:|:--------|:-----|:------------------------------|
#   |2      |[-Inf,0)        |   0.00|   0.00|   0.00|     0.00|      0.00|       NaN|        NaN|valsum.b |taxbc |Tax before credit (calculated) |
#   |2      |[0,2.5e+04)     |   0.20|   0.17|   0.22|    -0.03|      0.02|     -13.8|        8.0|valsum.b |taxbc |Tax before credit (calculated) |
#   |2      |[2.5e+04,5e+04) |  12.95|  16.49|  12.90|     3.55|     -0.04|      27.4|       -0.3|valsum.b |taxbc |Tax before credit (calculated) |
#   |2      |[5e+04,7.5e+04) |  43.43|  46.00|  42.57|     2.57|     -0.86|       5.9|       -2.0|valsum.b |taxbc |Tax before credit (calculated) |
#   |2      |[7.5e+04,1e+05) |  65.78|  69.34|  64.07|     3.56|     -1.71|       5.4|       -2.6|valsum.b |taxbc |Tax before credit (calculated) |
#   |2      |[1e+05,2e+05)   | 201.37| 183.85| 196.82|   -17.53|     -4.55|      -8.7|       -2.3|valsum.b |taxbc |Tax before credit (calculated) |
#   |2      |[2e+05,5e+05)   | 172.39| 176.95| 164.71|     4.56|     -7.68|       2.6|       -4.5|valsum.b |taxbc |Tax before credit (calculated) |
#   |2      |[5e+05,1e+06)   |  91.10| 115.11|  83.97|    24.01|     -7.13|      26.4|       -7.8|valsum.b |taxbc |Tax before credit (calculated) |
#   |2      |[1e+06,1e+07)   | 159.42| 179.69| 154.52|    20.27|     -4.91|      12.7|       -3.1|valsum.b |taxbc |Tax before credit (calculated) |
#   |2      |[1e+07, Inf)    |  48.27|  60.46|  48.30|    12.18|      0.03|      25.2|        0.1|valsum.b |taxbc |Tax before credit (calculated) |
#   |Total  |-               | 794.92| 848.06| 768.09|    53.14|    -26.83|       6.7|       -3.4|valsum.b |-     |-                              |

mrgdf %>%
  mutate(wt=ifelse(ftype=="puf", wt.puf, wt.syn / 5)) %>%
  mutate(agirange=cut(c00100, agiranges, right=FALSE),
         wtone=1e9) %>%
  dplyr::select(ftype, MARS, agirange, wt, wtone, vlist) %>%
  gather(vname, value, -ftype, -agirange, -wt, -MARS) %>%
  group_by(ftype, MARS, agirange, vname) %>%
  summarise(n=n(), wtsum.m=sum(wt) / 1e6, valsum.b=sum(wt * value) / 1e9) %>%
  filter(vname=="taxbc", MARS==2) %>%
  dplyr::select(ftype, agirange, stat="valsum.b", vname) %>%
  spread(ftype, stat) %>%
  janitor::adorn_totals(where="row") %>%
  kable(digits=c(0, rep(1, 5), rep(1, 2)), format.args=list(big.mark = ','))

#   |MARS  |agirange        |vname |   puf|   syn|
#   |:-----|:---------------|:-----|-----:|-----:|
#   |2     |[-Inf,0)        |taxbc |   0.0|   0.0|
#   |2     |[0,2.5e+04)     |taxbc |   0.2|   0.2|
#   |2     |[2.5e+04,5e+04) |taxbc |  12.9|  16.5|
#   |2     |[5e+04,7.5e+04) |taxbc |  43.4|  46.0|
#   |2     |[7.5e+04,1e+05) |taxbc |  65.8|  69.3|
#   |2     |[1e+05,2e+05)   |taxbc | 201.4| 183.8|
#   |2     |[2e+05,5e+05)   |taxbc | 172.4| 176.9|
#   |2     |[5e+05,1e+06)   |taxbc |  91.1| 115.1|
#   |2     |[1e+06,1e+07)   |taxbc | 159.4| 179.7|
#   |2     |[1e+07, Inf)    |taxbc |  48.3|  60.5|
#   |Total |-               |-     | 794.9| 848.1|


# yes, so it does not appear that anything went wrong prior to weighting

# now let's pick a group where 17 did poorly and 8 did well and see if we can reproduce the objective function


#****************************************************************************************************
#                weight one group ####
#****************************************************************************************************
# make a bigger stacked file with puf, synpuf17, and synthpop8
tcvars <- c("c00100", "taxbc")

df2 <- readRDS(paste0(globals$tc.dir, "synpuf17", "_rwprep.rds"))
mrgdf.17 <- left_join(df2$tc.base, df2$tc.output %>% dplyr::select(RECID, tcvars))

df3 <- readRDS(paste0(globals$tc.dir, "synthpop8", "_rwprep.rds"))
mrgdf.8 <- left_join(df3$tc.base, df3$tc.output %>% dplyr::select(RECID, tcvars))

setdiff(names(mrgdf.17), names(mrgdf.8)) # "wt.puf" "e00600_minus_e00650" "wt.syn"  on 17 but not 8
setdiff(names(mrgdf.8), names(mrgdf.17)) # wt" "divratio" "e09600"   "msname"  on 8 but not 17

# Let’s use group 67 -- mgroup 1, imin 1.6m, imax 2m, obj.17 28363.2, obj.8 3.4 puf 525 records
# also Group 177 mgroup 2 imin 3.6m imax 4m obj 17 16,361.9, obj 8 3.0
grp <- expression(MARS==1 & c00100 > 1.6e6 & c00100 <= 2e6)
igrp <- 67

grp <- expression(MARS==2 & c00100 > 3.6e6 & c00100 <= 4e6)
igrp <- 177

sub17 <- mrgdf.17 %>%
  filter(eval(grp)) %>%
  mutate(wt=ifelse(ftype=="puf", wt.puf, wt.syn)) %>%
  dplyr::select(-e00600_minus_e00650, -wt.puf, -wt.syn)
names(sub17) %>% sort
count(sub17, ftype)

sub8 <- mrgdf.8 %>%
  filter(eval(grp)) %>%
  dplyr::select(-divratio, -e09600, -msname)
count(sub8, ftype)

setdiff(names(sub17), names(sub8))
setdiff(names(sub8), names(sub17))

check17 <- readRDS(paste0(globals$tc.dir, "weight_pieces/optim_group_", igrp, ".rds"))
names(check17)

check8 <- readRDS(paste0("D:/tcdir/weight_pieces/synthpop8_500iter/optim_group_", igrp, ".rds"))
names(check8)

nrow(check17$puf); nrow(check17$syn) # good
nrow(check8$puf); nrow(check8$syn) # good

check8$result$objective
check17$result$objective

# ok, can we reproduce things? ----

recipe <- get_recipe_long(get_weighting_recipe("recipe5")) %>%
  filter(vname %in% names(sub17)) %>%
  dplyr::select(vname, vname, fn)
recipe

# run the functions below and then run this
reclist17 <- getrec(puf=sub17 %>% filter(ftype=="puf"),
                    syn=sub17 %>% filter(ftype=="syn"),
                    recipe, puf.vnames)
names(reclist17)
reclist17$recipe.use # 56 elements (not 57??)
reclist17$recipe.flagged

reclist8 <- getrec(puf=sub8 %>% filter(ftype=="puf"),
                    syn=sub8 %>% filter(ftype=="syn"),
                    recipe, puf.vnames)
names(reclist8)
reclist8$recipe.use # 56 elements (not 57??)
reclist8$recipe.flagged # 96

inputs17 <- getinplist(syn=sub17 %>% filter(ftype=="syn"),
                       recipe.use=reclist17$recipe.use)
names(inputs17)

inputs8 <- getinplist(syn=sub8 %>% filter(ftype=="syn"),
                       recipe.use=reclist8$recipe.use)
names(inputs8)

# run # 17
puf <- sub17 %>% filter(ftype=="puf")
syn <- sub17 %>% filter(ftype=="syn")
inputs <- inputs17

# run # 8
puf <- sub8 %>% filter(ftype=="puf")
syn <- sub8 %>% filter(ftype=="syn")
inputs <- inputs8



# bounds on the weights
xlb <- rep(1, nrow(syn))
xub <- rep(1.5*max(puf$wt), nrow(syn)) # FIX THIS djb

# starting point:
x0 <- (xlb + xub) / 2
x0 <- x0 * sum(puf$wt / sum(x0))

opts <- list("algorithm"="NLOPT_LD_MMA",
             "xtol_rel"=1.0e-8,
             "maxeval"=500)

a <- proc.time()
result <- nloptr(x0, 
                 eval_f=eval_f_wtfs,
                 eval_grad_f = eval_grad_f_wtfs,
                 lb = xlb, ub = xub,
                 opts = opts, inputs=inputs)
b <- proc.time()
b - a

# result17 <- result
# or
# result8 <- result


#****************************************************************************************************
#                look at results ####
#****************************************************************************************************
# Let’s use group 67 -- mgroup 1, imin 1.6m, imax 2m, obj.17 28363.2, obj.8 3.4 puf 525 records

# objective fun:
#  synpuf17 obj=28362  num.evals=1
#  synthpop8 obj=3.41
t8 <- tibble(wt=result8$solution, ftype="syn8")
t17 <- tibble(wt=result17$solution, ftype="syn17")

wts <- bind_rows(t8, t17)
wts %>%
  group_by(ftype) %>%
  do(qtiledf(.$wt, probs=c(0, .10, .25, .5, .75, .8, .85, .9, .92, .91, .93, .94, .95, .99, .999, 1)))

wcuts <- c(0, 1, 1.25, 1.5, 1.75, 2, 3, 4, 5, 10, 15, 20)
wcuts <- c(1, 1.0001, 1.1, seq(1.2, 2, .1), 2.5, 3, 4, 5, 10, 15, 20, 30)
wts %>%
  #mutate(wbreaks=cut(wt, 40, right=FALSE)) %>%
  mutate(wbreaks=cut(wt, wcuts, right=FALSE)) %>%
  group_by(ftype, wbreaks) %>%
  summarise(n=n()) %>%
  spread(ftype, n)

wts %>%
  ggplot(aes(wt, y=..density..)) +
  geom_histogram(fill="blue", binwidth = .05) +
  facet_wrap(~ftype, ncol=1)
  

# ok this seems to proveout what we expected

puf8 <- sub8 %>% filter(ftype=="puf")
syn8 <- sub8 %>% filter(ftype=="syn")

summary(syn8)

puf17 <- sub17 %>% filter(ftype=="puf")
syn17 <- sub17 %>% filter(ftype=="syn")

summary(puf8)
summary(puf17)

# compare inputs -- ok, puf looks good
bind_rows(puf8 %>% mutate(file="syn8"),
          puf17 %>% mutate(file="syn17")) %>%
  dplyr::select(-pufseqn) %>%
  gather(variable, value, -ftype, -file, -wt) %>%
  group_by(ftype, file, ftype, variable) %>%
  summarise(wsum6=sum(wt * value) / 1e6) %>%
  spread(file, wsum6) %>%
  mutate(diff=syn17 - syn8,
         pdiff=diff / syn8 * 100) %>%
  kable(digits=2)

# syn weighted sums are very different
bind_rows(syn8 %>% mutate(file="syn8"),
          syn17 %>% mutate(file="syn17")) %>%
  dplyr::select(-pufseqn) %>%
  gather(variable, value, -ftype, -file, -wt) %>%
  group_by(ftype, file, ftype, variable) %>%
  summarise(wsum6=sum(wt * value) / 1e6) %>%
  spread(file, wsum6) %>%
  mutate(diff=syn17 - syn8,
         pdiff=diff / syn8 * 100) %>%
  kable(digits=2)

# ok, now look at correlations with agi
gstack <- bind_rows(puf8 %>% mutate(ftype="puf"),
                    syn8 %>% mutate(ftype="syn8"),
                    syn17 %>% mutate(ftype="syn17"))

gstack %>%
  dplyr::select(-pufseqn) %>%
  gather(variable, value, -ftype, -wt) %>%
  group_by(ftype, ftype, variable) %>%
  summarise(wsum6=sum(wt * value) / 1e6) %>%
  spread(ftype, wsum6) %>%
  mutate(diff17=syn17 - puf,
         diff8=syn8 - puf,
         pdiff17=diff17 / puf * 100,
         pdiff8=diff8 / puf * 100,
         apworst=pmax(abs(pdiff17), abs(pdiff8))) %>%
  arrange(-apworst) %>%
  kable(digits=1)

corvars <- names(gstack)[(names(gstack) %>% str_sub(., 1, 1) %in% c("c", "e", "p", "t"))]
corvars <- corvars[-c(1, 4, 58:63)]
corvars

cor1 <- gstack %>%
  dplyr::select(ftype, corvars) %>%
  setNames(str_to_lower(names(.))) %>%
  group_by(ftype) %>%
  do(cordf(.[, -1])) %>%
  do(trimcor(.)) %>%
  separate(combo, c("vname1", "vname2"), sep="-", remove=TRUE) %>%
  ungroup %>%
  left_join(puf.vnames %>% dplyr::select(vname1=vname, vdesc1=vdesc)) %>%
  left_join(puf.vnames %>% dplyr::select(vname2=vname, vdesc2=vdesc))
ht(cor1)
count(cor1, ftype)

cor1 %>% 
  filter(vname1=="e00200" | vname2=="e00200") %>%
  arrange(vname2, ftype)

cor1 %>% 
  filter(vname1=="e00200" | vname2=="e00200") %>%
  arrange(desc(abs(value)))

cor.comp <- cor1 %>%
  mutate(ftype=factor(ftype, levels=c("puf", "syn8", "syn17"))) %>%
  spread(ftype, value) %>%
  filter(!is.na(puf), !is.na(syn8), !is.na(syn17)) %>% # ????
  mutate(diff8=syn8 - puf,
         diff17=syn17 - puf,
         diff.worst=pmax(abs(diff8), abs(diff17)),
         worst=ifelse(abs(diff8) > abs(diff17), "syn8", "syn17"))

cor.comp %>% filter(vname1=="c00100") %>% arrange(-diff.worst)

cor.comp %>%
  group_by(worst) %>%
  summarise(n=n(), adiff=sum(abs(diff.worst)))

count(cor.comp, worst) %>%
  janitor::adorn_totals() %>%
  mutate(pct=n / n[worst=="Total"] * 100) %>%
  kable(digits=1)

vsize <- 30
cor.comp %>%
  arrange(desc(diff.worst)) %>%
  mutate(vdesc1=str_sub(vdesc1, 1, vsize), vdesc2=str_sub(vdesc2, 1, vsize)) %>%
  head(25) %>%
  rename(cor.puf=puf, cor17=syn17, cor8=syn8) %>%
  kable(digits=3)
  


#****************************************************************************************************
#                density plots ####
#****************************************************************************************************

kdplot2 <- function(var){
  sq10 <- c(0, 1e3, 10e3, 25e3, 50e3, 100e3, 250e3, 500e3, 750e3, 1e6,
            1.5e6, 2e6, 3e6, 4e6, 5e6, 10e6, 25e6, 50e6, 100e6)
  xlabs <- scales::comma(sq10 / 1e3)
  
  xscale.l10 <- scale_x_log10(name=paste0(var, " in $ thousands, log 10 scale"), breaks=sq10, labels=xlabs)
  
  vdesc <- puf.vnames$vdesc[puf.vnames$vname==var]
  
  gtitle <- paste0("Kernel density plot for: ", vdesc, " (", var, ")")
  
  p <- gstack %>%
    dplyr::select(ftype, value=var) %>%
    filter(value > 1000) %>%
    mutate(ftype=factor(ftype, levels=c("puf", "syn8", "syn17"))) %>%
    ggplot(aes(value, colour=ftype)) +
    geom_density(size=1.5) +
    scale_colour_manual(values=c("blue", "red", "darkgreen")) +
    # scale_x_continuous(name="value of variable") +
    xscale.l10 +
    theme_bw() +
    ggtitle(gtitle, subtitle="Only includes values above $1,000") +
    theme(plot.title=element_text(size=12)) +
    theme(axis.text.x=element_text(angle=45, size=10, hjust=1, colour="black"))
  return(p)
}

count(gstack, ftype)

# e00200, e01500, e01700, e02400 4 largest variables
kdplot2("c00100")
kdplot2("e00200")
kdplot2("e00300")
kdplot2("e00600")
kdplot2("taxbc")

#****************************************************************************************************
#                functions below ####
#****************************************************************************************************
getrec <- function(puf, syn, recipe, puf.vnames){
    # adjust a previously created recipe
    tscale <- 1
    
    # use get rather than do.call to find a function as do.call does not seem to work in parallel
    recipe$target <- NA_real_
    for(i in 1:nrow(recipe)){
      recipe$target[i] <- get(recipe$fn[i])(df=puf, var=recipe$vname[i], puf$wt)
    }
    
    recipe <- recipe %>%
      mutate(scale=ifelse(target!=0, abs(target / tscale), 1/ tscale),
             obj.element=paste0(vname, "_", fn)) %>%
      dplyr::select(obj.element, vname, fn, scale, target) %>%
      arrange(vname, fn)
    
    #..weed out unnecessary elements of the recipe ----
    # if the target is 0 for negative values we can drop the neg versions AS LONG AS WE HAVE NO SYNTH NEG VERSIONS
    # if the val.pos and val.sum versions are identical then we can drop the val.sum version
    # can drop the "neg" and "pos" versions
    recipe.flagged <- recipe 
    recipe.flagged$syn.unwtd <- NA_real_
    for(i in 1:nrow(recipe.flagged)){
      recipe.flagged$syn.unwtd[i] <- get(recipe.flagged$fn[i])(df=syn, var=recipe.flagged$vname[i], weight=rep(1, nrow(syn)))
    }
    recipe.flagged <- recipe.flagged %>%
      group_by(vname) %>%
      mutate(flag.dropneg=ifelse(str_detect(fn, "neg") & target==0 & syn.unwtd==0, 1, 0),
             flag.dropdupsum=ifelse(target==target[match("val.sum", fn)] & (fn=="val.pos"), 1, 0),
             flag.dropdupn=ifelse(target==target[match("n.sum", fn)] & (fn=="n.pos"), 1, 0)) %>%
      mutate_at(vars(starts_with("flag")), funs(naz)) %>%
      ungroup %>%
      arrange(vname, fn)
    
    # remove recipe elements where the target is zero
    recipe.use <- recipe.flagged %>%
      filter(!(flag.dropneg | flag.dropdupsum | flag.dropdupn)) %>%
      filter(target!=0) %>%
      dplyr::select(obj.element, vname, fn, scale, target)
    
    # finally, add priority weights
    recipe.use <- recipe.use %>%
      mutate(priority.weight=case_when(vname %in% c("wt", "c00100", "e00200", "taxbc") ~ 100,
                                       fn %in% c("n.sum", "val.sum") ~ 100,
                                       TRUE ~ 1))  %>%
      left_join(puf.vnames %>% dplyr::select(vname, vdesc))
    return(list(recipe.use=recipe.use, recipe.flagged=recipe.flagged))
  }
  
  getinplist <- function(syn, recipe.use){
    inputs <- list()
    inputs$recipe <- recipe.use
    # inputs$synsub <- syn[, unique(inputs$recipe$vname)] %>% mutate(wt=1)
    inputs$synsub <- syn %>% dplyr::select(unique(inputs$recipe$vname)) %>% mutate(wt=1)
    synlong <- inputs$synsub %>%
      dplyr::mutate(wtnum=row_number()) %>%
      gather(vname, value, -wtnum)
    
    # create a data frame with one row for each weight and obj.element combination
    coeffs <- expand.grid(wtnum=1:nrow(inputs$synsub), 
                          obj.element=inputs$recipe$obj.element, stringsAsFactors = FALSE) %>%
      ungroup %>%
      left_join(inputs$recipe %>% dplyr::select(obj.element, vname, fn, scale, priority.weight, target)) %>%
      left_join(synlong) %>%
      mutate(coeff=case_when(fn=="val.sum" ~ value,
                             fn=="val.pos" ~ value*(value>0),
                             fn=="val.neg" ~ value*(value<0),
                             fn=="n.sum" ~ 1,
                             fn=="n.pos" ~ 1*(value>0),
                             fn=="n.neg" ~ 1*(value<0),
                             TRUE  ~ 0)) %>%
      dplyr::select(obj.element, vname, fn, wtnum, scale, priority.weight, value, coeff, target)
    
    inputs$coeffs <- coeffs
    return(inputs)
  }
  







