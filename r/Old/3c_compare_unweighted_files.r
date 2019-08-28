# Next up, compare:
#   synpuf17.csv -- Max's latest RF file on 100% of sample, with seeds MARS, DSI, XTOT, S006, F6251, MIDR, FDED and 5x records
#   synpuf20.csv -- Max --  200 trees and 5x records
#   synthpop10.csv -- my latest CART file, 9 seeds, 5x records, wages density, cart.minbucket 10 for MARS2


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

library("caret") # for confustionMatrix


#****************************************************************************************************
#                Includes ####
#****************************************************************************************************
source("./r/includes/globals_system_specific_boyd.r") # use a different version of this file if changing systems
source("./r/includes/globals_other.r")

source("./r/includes/functions_general.r")
source("./r/includes/functions_synthesis.r")


#****************************************************************************************************
#                Get summary data ####
#****************************************************************************************************
puf.vnames <- get_puf_vnames()

psum <- get_pufbase_sums()
psum
psum %>% arrange(vname)


#****************************************************************************************************
#                Get comparison files ####
#****************************************************************************************************
# get pufcart
# synname <- "synthpop7"
synname <- "synthpop10"
# dir <- globals$synd
dir <- globals$tc.dir
pufcart <- read_csv(paste0(dir, synname, "_stack.csv"), 
         col_types = cols(ftype=col_character(),
                          msname=col_character(),
                          .default= col_double()),
         n_max=-1)
glimpse(pufcart)
names(pufcart) %>% sort

puf <- pufcart %>% filter(ftype=="puf")
cart <- pufcart %>% filter(ftype=="syn") %>% mutate(ftype="cart")

# get the RF files
rf <- read_csv(paste0(globals$synd, "syntheses/", "synpuf20", ".csv"), 
                         col_types = cols(.default= col_double()), 
                         n_max=-1) %>% 
  mutate(ftype="rf") %>%
  setNames(change_case(names(.)))
glimpse(rf)


#****************************************************************************************************
#                Compare available variables ####
#****************************************************************************************************
not.rf <- setdiff(names(cart), names(rf))
not.rf %>% sort

not.cart <- setdiff(names(rf), names(cart))
not.cart %>% sort

psum %>%
  filter(vname %in% not.rf)

psum %>%
  filter(vname %in% not.cart)

psum %>%
  filter(vname %in% names(rf)) %>%
  arrange(vname) %>% tail(50)


#****************************************************************************************************
#                Create stack file ####
#****************************************************************************************************
# get minimum set of names
usenames <- intersect(names(cart), names(rf))
usenames %>% sort

stack <- bind_rows(puf %>% dplyr::select(usenames),
                   cart %>% dplyr::select(usenames),
                   rf %>% dplyr::select(usenames)) %>%
  mutate(ftype=factor(ftype, levels=c("puf", "cart", "rf"))) # so that ftype sorts as we want
count(stack, ftype)
names(stack)


#****************************************************************************************************
#                Summary statistics ####
#****************************************************************************************************
epvars <- names(stack)[(names(stack) %>% str_sub(., 1, 1) %in% c("e", "p"))]
epvars <- setdiff(epvars, "e01500_minus_e01700")
epvars

a <- proc.time()
long <- stack %>%
  dplyr::select(ftype, epvars) %>%
  # pivot_longer(cols=-ftype, names_to="vname", values_to="value") # 75 secs
  gather(vname, value, -ftype) # 3 secs
b <- proc.time()
b - a
glimpse(long)

a <- proc.time()
stats <- long %>%
  group_by(ftype, vname) %>%
  summarise(n=n(),
            n.NA=sum(is.na(value)),
            mean=mean(value, na.rm=TRUE),
            median=median(value, na.rm=TRUE), 
            sd=sd(value, na.rm=TRUE),
            kurtosis=e1071::kurtosis(value, type=2, na.rm=TRUE),
            skewness=e1071::skewness(value, type=2, na.rm=TRUE))
b <- proc.time()
b - a
stats %>% ht

wstats <- stats %>%
  gather(stat, value, -ftype, -vname) %>%
  spread(ftype, value) %>%
  mutate(stat=factor(stat, levels=c("n", "n.NA", "mean", "median", "sd", "kurtosis", "skewness"))) %>%
  arrange(vname, stat) %>%
  dplyr::select(vname, stat, puf, cart, rf) %>%
  mutate(cart.diff=cart - puf,
         rf.diff=rf - puf,
         worst=case_when(abs(cart.diff) > abs(rf.diff) ~ "cart",
                         abs(rf.diff) > abs(cart.diff) ~ "rf",
                         TRUE ~ "tie"))
ht(wstats)

wstats %>%
  filter(!stat %in% c("n", "n.NA")) %>%
  group_by(stat, worst) %>%
  summarise(n=n()) %>%
  spread(worst, n) %>%
  mutate_at(vars(-stat), ~naz(.)) %>%
  rename(cart.nworst=cart, rf.nworst=rf)


#****************************************************************************************************
#                Correlations ####
#****************************************************************************************************
corvars <- names(stack)[(names(stack) %>% str_sub(., 1, 1) %in% c("e", "p"))]

a <- proc.time()
cor1 <- stack %>%
  dplyr::select(ftype, corvars) %>%
  setNames(str_to_lower(names(.))) %>%
  group_by(ftype) %>%
  do(cordf(.[, -1])) %>%
  do(trimcor(.)) %>%
  separate(combo, c("vname1", "vname2"), sep="-", remove=TRUE) %>%
  ungroup %>%
  left_join(puf.vnames %>% dplyr::select(vname1=vname, vdesc1=vdesc)) %>%
  left_join(puf.vnames %>% dplyr::select(vname2=vname, vdesc2=vdesc))
b <- proc.time()
b - a

ht(cor1)
count(cor1, ftype)

cor1 %>% 
  filter(vname1=="e00200" | vname2=="e00200") %>%
  arrange(vname2, ftype)

cor1 %>% 
  filter(vname1=="e00200" | vname2=="e00200") %>%
  arrange(desc(abs(value)))

cor.comp <- cor1 %>%
  mutate(ftype=factor(ftype, levels=c("puf", "cart", "rf"))) %>%
  spread(ftype, value) %>%
  mutate(diff.cart=cart - puf,
         diff.rf=rf - puf,
         diff.worst=pmax(abs(diff.cart), abs(diff.rf)),
         worst=ifelse(abs(diff.cart) > abs(diff.rf), "cart", "rf"))

count(cor.comp, worst) %>%
  janitor::adorn_totals() %>%
  mutate(pct=n / n[worst=="Total"] * 100) %>%
  kable(digits=1)

vsize <- 30
cor.comp %>%
  arrange(desc(diff.worst)) %>%
  mutate(vdesc1=str_sub(vdesc1, 1, vsize), vdesc2=str_sub(vdesc2, 1, vsize)) %>%
  head(25) %>%
  rename(cor.puf=puf, cor.cart=cart, cor.rf=rf) %>%
  kable(digits=3)

cor.comp %>%
  filter(vname1=="e00200" | vname2=="e00200") %>%
  arrange(desc(diff.worst)) %>%
  mutate(vdesc1=str_sub(vdesc1, 1, vsize), vdesc2=str_sub(vdesc2, 1, vsize)) %>%
  head(25) %>%
  rename(cor.puf=puf, cor.cart=cart, cor.rf=rf) %>%
  kable(digits=3)

comp2 <- cor1 %>%
  mutate(ftype=factor(ftype, levels=c("puf", "cart", "rf"))) %>%
  spread(ftype, value) %>%
  gather(syn, value, cart, rf) %>%
  mutate(diff=value - puf,
         adiff=abs(diff))

probs <- c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1)
comp2 %>%
  group_by(syn) %>%
  do(qtiledf(.$diff, probs)) %>%
  kable(digits=3)

comp2 %>%
  group_by(syn) %>%
  do(qtiledf(.$adiff, probs)) %>%
  gather(ptile, value, -syn, -n, -n.notNA) %>%
  mutate(ptile=factor(ptile, levels=unique(ptile))) %>%
  spread(syn, value) %>%
  dplyr::select(ptile, cart, rf) %>%
  kable(digits=3)

# comp2 %>%
#   ggplot(aes(diff, colour=syn)) +
#   geom_density(size=1.5) +
#   scale_colour_manual(values=c("blue", "red")) +
#   scale_x_continuous(name="Syn correlation minus puf correlation",
#                      breaks=seq(-.4, .4, .05),
#                      limits=c(-.2, .2)) +
#   geom_vline(xintercept = 0, linetype="dashed", size=1) +
#   theme_bw() +
#   ggtitle("Kernel density: Difference between correlation of variable pairs in each syn file\nand correlation in puf file", 
#           subtitle="1540 pairs of variables") +
#   theme(plot.title=element_text(size=12))


#****************************************************************************************************
#                Bivariate comparisons ####
#****************************************************************************************************



#****************************************************************************************************
#                fits of individual variables ####
#****************************************************************************************************
glimpse(stack)

# variable by marital status e00200 e01500
var <- "e00200"
var <- "e01500"
var <- "e19200"
vdesc <- puf.vnames$vdesc[puf.vnames$vname==var]
stack %>%
  dplyr::select(ftype, MARS, value=var) %>%
  filter((var=="e00200") | (value > 0)) %>%
  mutate(ftype=factor(ftype, levels=c("puf", "cart", "rf")),
         marsgroup=case_when(MARS==1 ~ "single",
                             MARS==2 ~ "married",
                             MARS %in% 3:4 ~ "other"),
         marsgroup=factor(marsgroup, levels=c("single", "married", "other")),
         lvalue=log(value + 1)) %>%
  ggplot(aes(ftype, lvalue, colour=ftype)) +
  scale_colour_manual(values=c("blue", "red", "darkgreen")) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Log 10 of value by marital status and file type",
          subtitle = vdesc) +
  facet_wrap(~marsgroup, nrow=1, scales="free_y")

# dev.off()

var <- "e00200"
var <- "e01500"
var <- "e19200"
vdesc <- puf.vnames$vdesc[puf.vnames$vname==var]
stack %>%
  dplyr::select(ftype, XTOT, value=var) %>%
  filter((var=="e00200") | (value > 0)) %>%
  mutate(ftype=factor(ftype, levels=c("puf", "cart", "rf")),
         lvalue=log(value + 1)) %>%
  ggplot(aes(ftype, lvalue, colour=ftype)) +
  scale_colour_manual(values=c("blue", "red", "darkgreen")) +
  geom_boxplot() +
  theme_bw() +
  ggtitle("Log 10 of value by XTOT and file type",
          subtitle = vdesc) +
  facet_wrap(~XTOT, nrow=1, scales="free_y")


#****************************************************************************************************
#                binning ####
#****************************************************************************************************
ncuts <- function(v, nbins=100) {
  cuts <- quantile(v, probs = seq(0, 1, length = nbins), na.rm = TRUE, type = 2)
  # replace endpoints
  cuts[1] <- -Inf
  cuts[length(cuts)] <- Inf
  # cuts <- c(-Inf, as.vector(cuts), Inf)
  return(cuts)
}

ncutsdf <- function(v){
  tibble(bincut=ncuts(v), cutnum=1:length(bincut))
}


# get bin cutpoints based upon puf
glimpse(stack)
epvars <- names(stack)[(names(stack) %>% str_sub(., 1, 1) %in% c("e", "p"))]
epvars <- setdiff(epvars, "e01500_minus_e01700")
epvars
binvars <- epvars[1:5]
binvars <- epvars

bins.all <- stack %>%
  filter(ftype=="puf") %>%
  dplyr::select(binvars) %>%
  # sample_frac(size=.1) %>%
  gather(vname, value) %>%
  group_by(vname) %>%
  do(ncutsdf(.$value))
length(unique(bins.all$vname))

# collapse the bins, as there may be multiple bins for values of zero (for example)
bins <- bins.all %>%
  group_by(vname, bincut) %>%
  summarise(cutnum.min=min(cutnum), cutnum.max=max(cutnum)) %>%
  group_by(vname) %>%
  mutate(cutnum=row_number())


getbin <- function(df){
  brks <- bins$bincut[bins$vname==df$vname[1]]
  df$valgroup <- cut(df$value, brks)
  df$binnum <- as.integer(df$valgroup)
  df$valgroup <- as.character(df$valgroup) # factors don't work once we combine with other variables with different levels
  return(df)
}

# df <- stack %>%
#   sample_n(100) %>%
#   dplyr::select(ftype, binvars) %>%
#   gather(vname, value, -ftype) %>%
#   group_by(ftype, vname) 
# 
# df <- df %>% filter(ftype=="cart", vname=="e00200")

# determine the bin for each observation
a <- proc.time()
bindat <- stack %>%
  # sample_n(100) %>%
  dplyr::select(ftype, binvars) %>%
  gather(vname, value, -ftype) %>%
  group_by(ftype, vname) %>%
  do(getbin(.)) %>%
  ungroup %>%
  arrange(ftype, vname, binnum)
b <- proc.time()
b - a # 2-3 minutes

# collapse by ftype, vname, and bin
binned <- bindat %>%
  group_by(ftype, vname) %>%
  mutate(ntot=n()) %>%
  group_by(ftype, vname, binnum, valgroup) %>%
  summarise(n=n(), ntot=first(ntot)) %>%
  mutate(pct=n / ntot * 100)

# get ssd by ftype, vname
bincomp <- binned %>%
  group_by(vname, binnum, valgroup) %>%
  mutate(diff=pct - pct[ftype=="puf"],
         diffsq=diff^2) %>%
  group_by(ftype, vname) %>%
  summarise(ssd=sum(diffsq)) %>%
  filter(ftype!="puf")

# get file-total summaries, across all variables
bincomp %>%
  group_by(vname) %>%
  mutate(worst=case_when(ftype=="cart" & ssd > ssd[ftype=="rf"] ~ 1,
                         ftype=="rf" & ssd > ssd[ftype=="cart"] ~ 1,
                         TRUE ~ 0)) %>%
  group_by(ftype) %>%
  summarise(n.worst=sum(worst), ssd.mdn=median(ssd), ssd.mean=mean(ssd)) %>%
  kable(digits=2)

bincomp %>%
  spread(ftype, ssd) %>%
  mutate(abs_diff=abs(rf - cart),
         worst=ifelse(rf > cart, "rf", "cart")) %>%
  rename(cart_ssd=cart, rf_ssd=rf) %>%
  adorn_totals() %>%
  left_join(psum %>% dplyr::select(vname, vdesc, puf.billions=value)) %>%
  # arrange(-abs_diff) %>%
  arrange(-puf.billions) %>%
  kable(digits=1, format.args=list(big.mark = ','))

# why do wages do so poorly??
f <- function(x, y) (x - y)^2
binned %>%
  filter(vname=="e00200") %>%
  dplyr::select(ftype, vname, binnum, valgroup, pct) %>%
  spread(ftype, pct) %>%
  dplyr::select(vname, binnum, valgroup, puf, cart, rf) %>%
  mutate_at(vars(cart, rf), list(ssd=~f(., puf))) %>%
  kable(digits=2)

# wages: cart 203.1 of 203.1 ssd is in 1st 2 brackets, 0-$150, 0 is elsewhere
#        rf of 38.9 of 81.9 is in 1st 2, 43 is elsewhere
# so if we can live with cart's pushing returns into >0, <=$150, cart is better

# var <- "e00200"
# binplot <- function(var, llim=0, ulim=1, vnames=puf.vnames){
#   vdesc <- vnames$vdesc[match(var, vnames$vname)]
#   
#   gtitle <- paste0("puf vs. syn relative frequencies for ", var, ": ", vdesc)
#   gsub <- paste0("using 1,000 puf-based bins; plot includes bins where puf_pct ranged from ", llim, " to ", ulim)
#   # ylab <- paste0("Cumulative proportion of the sum of unweighted ", var)
#   
#   nvals <- 10
#   p <- bins %>%
#     filter(vname==var) %>%
#     # filter(bin!="[0,1e-09)") %>%
#     filter(puf_pct >= llim, puf_pct <= ulim) %>%
#     gather(syntype, syn_pct, cart_pct, rf_pct) %>%
#     ggplot(aes(puf_pct, syn_pct, colour=syntype)) +
#     geom_point() +
#     scale_x_continuous(breaks=seq(0, 100, (ulim - llim)/nvals), limits=c(llim, ulim)) +
#     scale_y_continuous(breaks=seq(0, 100, (ulim - llim)/nvals), limits=c(llim, ulim)) +
#     geom_abline(slope=1, intercept = 0) +
#     theme_bw() +
#     ggtitle(gtitle, subtitle=gsub)
#   return(p)
# }
# 
# # compare the shares of a variable that fall in different ntiles
# binplot("e00200")
# binplot("e00200", llim=.038, ulim=.23)
# 
# binplot("e00300")


# ksf <- function(vname, ftype.in){
#   # return the Kolmogorov-Smirnov test p.value
#   # null hypothesis that x and y were drawn from the same continuous distribution is performed.
#   # reject null when D is large
#   # the larger the test statistic D is, the greater the distance between the cdfs of the distributions for
#   # x and y (puf and syn)
#   # we want D to be small, and p to be high
#   
#   x <- stack %>% filter(ftype=="puf") %>% .[[vname]]
#   y <- stack %>% filter(ftype==ftype.in) %>% .[[vname]]
#   D <- ks.test(x, y)$statistic
#   return(D)
# }


# ksvars <- names(stack)[(names(stack) %>% str_sub(., 1, 1) %in% c("e", "p"))]
# 
# D.cart <- laply(ksvars, ksf, ftype.in="cart")
# D.rf <- laply(ksvars, ksf, ftype.in="rf")

#****************************************************************************************************
#                selected kernel density plots ####
#****************************************************************************************************
glimpse(stack)

kdplot <- function(var){
  sq10 <- c(0, 1e3, 10e3, 25e3, 50e3, 100e3, 250e3, 500e3, 750e3, 1e6,
            1.5e6, 2e6, 3e6, 4e6, 5e6, 10e6, 25e6, 50e6, 100e6)
  xlabs <- scales::comma(sq10 / 1e3)
  
  xscale.l10 <- scale_x_log10(name=paste0(var, " in $ thousands, log 10 scale"), breaks=sq10, labels=xlabs)
  
  vdesc <- puf.vnames$vdesc[puf.vnames$vname==var]
  
  gtitle <- paste0("Kernel density plot for: ", vdesc, " (", var, ")")
  
  p <- stack %>%
    dplyr::select(ftype, value=var) %>%
    filter(value > 1000) %>%
    mutate(ftype=factor(ftype, levels=c("puf", "cart", "rf"))) %>%
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

# 4 largest income items e00200, e01500, e01700, e02400
kdplot("e00200")
kdplot("e01500")
kdplot("e01700")
kdplot("e02400")

# 3 large deduction items e19200 e18400 e19800
kdplot("e19200")
kdplot("e18400")
kdplot("e19800")



f <- function(vname){
  p <- kdplot(vname)
  ggsave(paste0("./results/kdplot_", vname, ".png"), width=8, height=8, dpi=300)
}

f("e00200")
f("e01500")
f("e01700")
f("e02400")


#****************************************************************************************************
#                selected CDF plots ####
#****************************************************************************************************
cdfplot2 <- function(var, stackdf=stack, vnames=puf.vnames){
  # print(var)
  vdesc <- vnames$vdesc[match(var, vnames$vname)]
  df <- stackdf %>%
    dplyr::select(ftype, value=var) %>%
    group_by(ftype) %>%
    arrange(value) %>%
    mutate(cum.pct=cumsum(value) / sum(value)) %>%
    ungroup
  
  if(nrow(df) < 100){
    p <- paste0(var, " has only ", nrow(df), " rows.")
    return(p)
  }
  
  # find a good minimum x-axis value to start the plot on -- based on a desired cumulative percentage
  cum.pct.threshold <- .01
  iminval <- min(which(df$cum.pct[df$ftype=="puf"] > cum.pct.threshold))
  minval <- df$value[df$ftype=="puf"][iminval]
  # minval
  
  capt <- "- x-axis is log10 scale\n- For display purposes x-axis is truncated at left to start at puf's cum.pct=1%"
  gtitle <- paste0("Cumulative distribution of unweighted ", var, ": ", vdesc)
  gsub <- "Aggregate records excluded"
  ylab <- paste0("Cumulative proportion of the sum of unweighted ", var)
  
  # define x scale break points and associated labels
  sq10 <- c(0, 1e3, 10e3, 25e3, 50e3, 100e3, 250e3, 500e3, 750e3, 1e6,
            1.5e6, 2e6, 3e6, 4e6, 5e6, 10e6, 25e6, 50e6, 100e6)
  xlabs <- scales::comma(sq10 / 1e3)
  xscale.l10 <- scale_x_log10(name=paste0(var, " in $ thousands"), breaks=sq10, labels=xlabs)
  
  p <- df %>%
    filter(value > minval) %>%
    mutate(ftype=factor(ftype, levels=c("puf", "cart", "rf"))) %>%
    ggplot(aes(value, cum.pct, colour=ftype)) + 
    geom_line(size=1.5) +
    scale_colour_manual(values=c("blue", "red", "darkgreen")) +
    theme_bw() +
    ggtitle(gtitle, subtitle=gsub) +  labs(caption=capt) +
    scale_y_continuous(name=ylab, breaks=c(seq(0, .9, .05), seq(.92, 1, .02))) +
    xscale.l10 +
    theme(axis.text.x=element_text(angle=45, size=10, hjust=1, colour="black")) +
    theme(plot.caption = element_text(hjust=0, size=rel(.8)))
  return(p)
}


f <- function(vname){
  p <- cdfplot2(vname)
  ggsave(paste0("./results/cdfplot_", vname, ".png"), width=8, height=8, dpi=300)
}

f("e00200")
f("e01500")
f("e01700")
f("e02400")



#******************************************************************************************************************
#  NOT USEFUL:utility analysis - logistic regression ####
#******************************************************************************************************************

glimpse(stack)
count(stack, ftype)
names(stack) %>% sort

stack.cart <- stack %>% 
  filter(ftype %in% c("puf", "cart")) %>%
  mutate(syn.ind=ifelse(ftype=="cart", 1, 0))
count(stack.cart, ftype, syn.ind)
glimpse(stack.cart)
#  %>% mutate(XTOT=as.factor(XTOT))

stack.rf <- stack %>% 
  filter(ftype %in% c("puf", "rf")) %>%
  #bind_rows(stack %>% filter(ftype=="rf")) %>%
  mutate(syn.ind=ifelse(ftype=="rf", 1, 0))
count(stack.rf, ftype, syn.ind)
glimpse(stack.rf)


factor.vars <- c("DSI", "MARS", "MIDR", "n24", "XTOT", "EIC", "FDED", "f2441", "f6251")
dropvars <- c("ftype", factor.vars)

(rhs.vars <- setdiff(names(stack), dropvars) %>% sort)
(frm <- as.formula(paste("syn.ind ~ ", paste(rhs.vars, collapse= "+"))))

(rhs.vars <- c(setdiff(names(stack), dropvars),
               "factor(XTOT)", "factor(MIDR)", "factor(MARS)", "factor(FDED)") %>% sort)
(frm <- as.formula(paste("syn.ind ~ ", paste(rhs.vars, collapse= "+"))))


mod1 <- glm(frm, data = stack.cart, family = "binomial")
summary(mod1)
p <- predict(mod1, type="response")
p_class=ifelse(p > .5, 1, 0) %>% as.numeric
cart.cm <- confusionMatrix(factor(p_class), factor(stack.cart$syn.ind))
cart.cm


mod2 <- glm(frm, data = stack.rf, family = "binomial")
summary(mod2)
p <- predict(mod2, type="response")
p_class=ifelse(p > .66667, 1, 0) %>% as.numeric
rf.cm <- confusionMatrix(factor(p_class), factor(stack.rf$syn.ind))
rf.cm

summary(mod1)
summary(mod2)

utility(mod1) # need to take # obs into account!
utility(mod2)


# now, include every variable, and every variable crossed with wages

base <- c(setdiff(names(stack), dropvars)) %>% sort
wcross <- paste0(setdiff(base, "e00200"), ":", "e00200")
factvars <- c("factor(XTOT)", "factor(MIDR)", "factor(MARS)", "factor(FDED)")
(rhs.vars <- c(base, wcross, factvars))

(frm.cross <- as.formula(paste("syn.ind ~ ", paste(rhs.vars, collapse= "+"))))


mod1c <- glm(frm.cross, data = stack.cart, family = "binomial")
summary(mod1c)
p <- predict(mod1c, type="response")
p_class=ifelse(p > .5, 1, 0) %>% as.numeric
cart.cm.c <- confusionMatrix(factor(p_class), factor(stack.cart$syn.ind))
cart.cm.c


mod2c <- glm(frm, data = stack.rf, family = "binomial")
summary(mod2c)
p <- predict(mod2c, type="response")
p_class=ifelse(p > .66667, 1, 0) %>% as.numeric
rf.cm.c <- confusionMatrix(factor(p_class), factor(stack.rf$syn.ind))
rf.cm.c

summary(mod1)
summary(mod2)

utility(mod1) # need to take # obs into account!
utility(mod2)



ufrm.cross <- formula(syn.ind ~ v1 + v2 + v3 + v4 + v5 + 
                        v1:v2 + v1:v3 + v1:v4 + v1:v5 +
                        v2:v3 + v2:v4 + v2:v5 +
                        v3:v4 + v3:v5 +
                        v4:v5)




# names(mod1)
# mod1$coefficients
# length(mod1$control)
# 
# cbind(coef(summary(mod1))[,2], coef(summary(mod1))[,4])
# coefs <- coef(summary(mod1))[,2]
# pvals <- coef(summary(mod1))[,4]
# tibble(vname=names(coefs), coef=coefs, pval=pvals) %>%
#   filter(!row_number()==1) %>%
#   left_join(get_pufbase_sums() %>% rename(puf.billions=value)) %>%
#   arrange(pval) %>%
#   head(25) %>%
#   kable(digits=c(0, 6, 6, 1, 0))
# 
# 
# utility(mod1)


#****************************************************************************************************
#                quantiles ####
#****************************************************************************************************
qvars <- names(stack)[(names(stack) %>% str_sub(., 1, 1) %in% c("e", "p"))]

tmp <- stack %>%
  dplyr::select(ftype, qvars[1:5]) %>%
  gather(vname, value, -ftype) %>%
  filter(value>=1000) %>%
  group_by(vname, ftype) %>%
  do(qtiledf(.$value))

tmp %>%
  kable(digits=1, format.args=list(big.mark = ','))




#****************************************************************************************************
#                NOT USED regression ####
#****************************************************************************************************
library("broom")
glimpse(stack)


frm <- formula(e00200 ~ e00300 + e00600 + e00700 + e01500 + e02000 + e02100 +
                 e02300 + e02400 + e00900 + e18500 + e19800 + factor(XTOT) + factor(DSI))

#frm <- formula(e00200 ~ e17500 + e18400 + e18500 + e19200 + e04600 + e19700)

mods <- stack %>%
  group_by(ftype) %>%
  do(lm(frm, data=.) %>% tidy)

mods %>%
  dplyr::select(ftype, term, estimate) %>%
  spread(ftype, estimate) %>%
  dplyr::select(term, puf, cart, rf) %>%
  mutate(pdiff.cart=cart / puf * 100 - 100,
         pdiff.rf=rf / puf * 100 - 100,
         worst=ifelse(abs(cart - puf) > abs(rf - puf), "cart", "rf")) %>%
  kable(digits=1)


mods %>%
  arrange(term, ftype)

count(stack, ftype, DSI)


# mods %>%
#   dplyr::select(ftype, term, estimate, p.value) %>%
#   gather(stat, value, estimate, p.value) %>%
#   unite(fstat, ftype, stat) %>%
#   spread(fstat, value) %>%
#   dplyr::select(term, puf_estimate, contains("estimate"), puf_p.value, contains("p.value")) %>%
#   kable(digits=3)

# stats <- stack %>%
#   group_by(ftype) %>%
#   do(lm(frm, data=.) %>% glance)
# 
# stats


                 
                 