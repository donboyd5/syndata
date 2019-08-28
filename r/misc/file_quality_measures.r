

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

library("mvtnorm")

library("synthpop")

library("Matrix")

library("knitr")

library("broom")

#****************************************************************************************************
#                functions ####
#****************************************************************************************************

mindist <- function(adf, bdf){
  # compute distances from file a to nearest rec in file b
  # NOTE first variable must be ftype, all others must be numeric
  ftypea <- adf$ftype[1]
  ftypeb <- bdf$ftype[1]
  
  df.scale <- bind_rows(adf, bdf) %>%
    group_by(ftype) %>%
    do(as.data.frame(scale(.[, -1]))) %>%
    mutate(rn=row_number(), lab=paste0(ftype, "_", rn)) %>%
    ungroup %>%
    dplyr::select(ftype, rn, lab, everything())
  
  # now that we are prepped for distances, do the calculations
  m <- as.matrix(df.scale[, -c(1:3)])
  rownames(m) <- df.scale$lab
  
  # get the distance matrix
  dm <- as.matrix(dist(m))
  
  dmdf <- as_tibble(dm, rownames="lab") %>%
    gather(brec, dist, -lab) %>%
    separate(lab, c("afile", "rn1"), sep="_") %>%
    separate(brec, c("bfile", "rn2"), sep="_") %>%
    mutate(rn1=as.numeric(rn1), rn2=as.numeric(rn2)) %>%
    filter(afile==ftypea, bfile!=ftypea)
  
  # finally, retrieve the minimum distances
  dmin <- dmdf %>%
    group_by(rn1, bfile) %>%
    arrange(dist) %>%
    filter(row_number()<=1) %>%
    ungroup %>%
    arrange(afile, rn1, dist)
  
  return(list(df.scale=df.scale, dmin=dmin))
}


rnoise <- function(x, frac=.1){
  # add noise to a vector, with mean 0, sd a fraction of the variable's mean
  # set.seed(1234) # allow noise to vary from variable to variable
  err <- rnorm(length(x), mean=0, sd=sd(x) * frac)
  return(x + err)
}


#****************************************************************************************************
#                globals ####
#****************************************************************************************************

source("D:/Dropbox/RPrograms PC/OSPC/syndata4/r/includes/functions_general.r")

# define file types we may want to include
ft1 <- c("a", "b.lownoise", "b.syn", "b.rf", "b.midnoise", "b.hinoise") # all of them
ft2 <- c("a", "b.lownoise", "b.syn", "b.rf", "b.hinoise")
ft3 <- c("a", "b.lownoise", "b.syn", "b.hinoise")

ft <- ft1

seedval <- 10
nrows <- 300

#****************************************************************************************************
#                create data ####
#****************************************************************************************************

# define characteristics
# correlation matrix
sigma <- matrix(0, nrow=4, ncol=4)
sigma[upper.tri(sigma)] <- c(-.5, .4, .6, .3, -.2, -.5)
sigma <- sigma + t(sigma)
diag(sigma) <- 1
sigma

# make positive definite
sigma.pd <- nearPD(sigma, corr=TRUE)$mat %>% as.matrix # to be safe
sigma
sigma.pd
(sigma.pd - sigma) %>% round(3)

mvals <- c(1, 2, 3, 4) # means of variables

# create the a data file
set.seed(seedval)
# create a data frame of 4 multivariate normal random variables with means and covariance as defined above
a <- rmvnorm(n=nrows, mean=mvals, sigma=sigma.pd) %>%
  as_tibble(.name_repair=function(x) paste0("v", 1:4)) %>%
  mutate(ftype="a") %>%
  dplyr::select(ftype, everything())
a
cor(a %>% dplyr::select(-ftype))

# construct a 5th variable with a linear relationship to the other 4, plus noise
set.seed(seedval)
err.rnorm <- rnorm(n=nrows, mean=0, sd=1)
coefs.true <- c(10, 3, -2, 2, -.5)

sdfrac <- .05
quantile(err.rnorm*sdfrac)

a2 <- a %>%
  mutate(err=err.rnorm * sdfrac, 
         lhs.true=coefs.true[1] + (coefs.true[2]*v1) + (coefs.true[3]*v2) + (coefs.true[4]*v3) + (coefs.true[5]*v4),
         lhs.obs=lhs.true + err)
a2 %>% ht
a2 %>% ggplot(aes(lhs.true, lhs.obs)) + geom_point() + geom_abline(slope=1, intercept=0)
cor(a2 %>% dplyr::select(-ftype))
colMeans(a2 %>% dplyr::select(-ftype))
summary(a2)

mod0 <- lm(lhs.true ~ v1 + v2 + v3 + v4, data=a2)
summary(mod0)
mod1 <- lm(lhs.obs ~ v1 + v2 + v3 + v4, data=a2)
summary(mod1)

#.. create files -- a is our base, b files have noise or are synthesized ----
adf <- a2 %>%
  dplyr::select(ftype, v1, v2, v3, v4, v5=lhs.obs)
adf

# synthesized
synx <- syn(adf, method="cart", seed=seedval)
b.cart <- synx$syn %>% mutate(ftype="b.cart")
# b.cart

synx <- syn(adf, method="rf", seed=seedval)
b.rf <- synx$syn %>% mutate(ftype="b.rf")
# b.rf


#.. different noise choices ----
# b files will be adf plus random noise around each variable
# vnoise <- c(.2, .7, 1)
# vnoise <- c(.2, .8, 1.2)
# vnoise <- c(.2, .9, 1.4)
# vnoise <- c(.3, 1, 3)

# Before:
# funs(name = f(.)
#      
# After: 
# list(name = ~f(.))


vnoise <- c(.25, .4, 1)

set.seed(seedval)
b.lownoise <- adf %>%
  mutate(ftype="b.lownoise") %>%
  mutate_at(vars(-ftype), list(~rnoise(., vnoise[1]))) # new dplyr syntax
# b.lownoise

set.seed(seedval)
b.midnoise <- adf %>%
  mutate(ftype="b.midnoise") %>%
  mutate_at(vars(-ftype), funs(rnoise(., vnoise[2])))
# b.midnoise

set.seed(seedval)
b.hinoise <- adf %>%
  mutate(ftype="b.hinoise") %>%
  mutate_at(vars(-ftype), funs(rnoise(., vnoise[3])))
# b.hinoise


stack <- bind_rows(adf, b.lownoise, b.midnoise, b.cart, b.rf, b.hinoise)


# examine the files ----

#.. minimum distances ----
# each synfile vs. true
dlist.low <- mindist(b.lownoise, adf)
dlist.mid <- mindist(b.midnoise, adf)
dlist.high <- mindist(b.hinoise, adf)
dlist.cart <- mindist(b.cart, adf)
dlist.rf <- mindist(b.rf, adf)

# names(dlist.low)
# names(dlist.low$dmin)

dmin <- bind_rows(dlist.low$dmin,
                  dlist.mid$dmin,
                  dlist.high$dmin,
                  dlist.cart$dmin,
                  dlist.rf$dmin)
# glimpse(dmin)
# count(dmin, afile, bfile)

dmin %>% 
  group_by(afile) %>%
  summarise(dist.mdn=median(dist), dist.sd=sd(dist)) %>%
  arrange(dist.mdn)

dmin %>%
  mutate(bfile=factor(afile, 
                      levels=c("b.cart", "b.rf", "b.lownoise", "b.midnoise", "b.hinoise"))) %>%
  ggplot(aes(bfile, dist, colour=bfile)) +
  geom_boxplot(size=.7) +
  theme_bw() +
  scale_y_continuous(name="Euclidean distance", breaks=seq(0, 10, .5)) +
  scale_x_discrete(name=NULL) +
  ggtitle("Distances from records in synthetic file to nearest record in true file")


#.. correlations ----
# glimpse(stack)
cor1 <- stack %>%
  group_by(ftype) %>%
  do(cordf(.[, -1])) %>%
  do(trimcor(.)) %>%
  separate(combo, c("vname1", "vname2"), sep="-", remove=TRUE) %>%
  ungroup
ht(cor1)
count(cor1, ftype)

cor1 %>%
  group_by(vname1, vname2) %>%
  mutate(diff=value - value[ftype=="a"],
         adiff=abs(diff)) %>%
  ungroup %>%
  arrange(vname1, vname2, adiff) %>%
  filter(ftype!="a") %>%
  dplyr::select(ftype, vname1, vname2, adiff) %>%
  group_by(vname1, vname2) %>%
  filter(row_number()==1) %>%
  group_by(ftype) %>%
  summarise(n.best=n())

cor1 %>%
  group_by(vname1, vname2) %>%
  mutate(diff=value - value[ftype=="a"],
         adiff=abs(diff)) %>%
  ungroup %>%
  filter(ftype!="a") %>%
  unite(pair, vname1, vname2) %>%
  dplyr::select(pair, ftype, adiff) %>%
  spread(pair, adiff) %>%
  mutate(adiff.median=apply(.[, -1], 1, median)) %>%
  arrange(adiff.median) %>%
  kable(digits=3)



# create a correlation measure - median of absolute differences
cordiffs <- cor1 %>%
  group_by(vname1, vname2) %>%
  mutate(diff=value - value[ftype=="a"],
         adiff=abs(diff)) %>%
  ungroup %>%
  arrange(vname1, vname2, adiff) %>%
  filter(ftype!="a") %>%
  group_by(ftype) %>%
  summarise(n=n(), adiff.mean=mean(adiff), adiff.mdn=median(adiff))

cordiffs %>%
  kable(digits=3)

#..bin distributions ----
# v <- 1:100
# nbins <- 10
ncuts <- function(v, nbins=10) {
  cuts <- quantile(v, probs = seq(0, 1, length = nbins - 1), na.rm = TRUE, type = 2)
  # replace endpoints
  cuts[1] <- -Inf
  cuts[length(cuts)] <- Inf
  # cuts <- c(-Inf, as.vector(cuts), Inf)
  return(cuts)
}

ncutsdf <- function(v){
  tibble(bincut=ncuts(v), cutnum=1:length(bincut))
}

bins <- stack %>%
  filter(ftype=="a") %>%
  gather(vname, value, -ftype) %>%
  group_by(vname) %>%
  do(ncutsdf(.$value))

f <- function(df){
  df$valgroup <- cut(df$value, breaks=bins$bincut[bins$vname==df$vname])
  df$binnum <- as.integer(df$valgroup)
  return(df)
}

bindat <- stack %>%
  gather(vname, value, -ftype) %>%
  group_by(ftype, vname) %>%
  do(f(.))

binned <- bindat %>%
  group_by(ftype, vname) %>%
  mutate(ntot=n()) %>%
  group_by(ftype, vname, binnum, valgroup) %>%
  summarise(n=n(), ntot=first(ntot)) %>%
  mutate(pct=n / ntot * 100)

bincomp <- binned %>%
  group_by(vname, binnum, valgroup) %>%
  mutate(diff=pct - pct[ftype=="a"],
         diffsq=diff^2) %>%
  group_by(ftype, vname) %>%
  summarise(ssd=sum(diffsq)) %>%
  filter(ftype!="a")

bin1 <- bincomp %>%
  mutate(vname=paste0("ssd.", vname)) %>%
  spread(vname, ssd) %>%
  ungroup %>%
  mutate(ssd.median=apply(.[, -1], 1, median),
         ssd.mean=(ssd.v1 + ssd.v2 + ssd.v3 + ssd.v4 + ssd.v5) / 5) %>%
  arrange(ssd.median)

bin1 %>% dplyr::select(-ssd.mean) %>% kable(digits=1)

#.. regression ----
frm <- formula(v5 ~ v1 + v2 + v3 + v4)
terms <- stack %>%
  group_by(ftype) %>%
  do(lm(frm, data=.) %>% tidy)
# terms %>%
#   arrange(term, ftype)

# rmse2 <- function(mod) sqrt(mean(mod$residuals^2))
# rmse2(mod1)
# summary(mod1)
# summary(mod1)$sigma

# create dataframe with estimated coefs
coefdf <- terms %>%
  dplyr::select(ftype, term, estimate) %>%
  mutate(term=ifelse(term=="(Intercept)", "intercept", paste0(term, "_coef"))) %>%
  spread(term, estimate)
# coefdf

stats <- stack %>%
  group_by(ftype) %>%
  do(lm(frm, data=.) %>% glance)

stats %>%
  dplyr::select(ftype, sigma, r.squared) %>%
  arrange(ftype!="a", -r.squared) %>%
  # left_join(coefdf) %>%
  kable(digits=2)


#.. propensity scores ----
ftypes <- c("b.cart", "b.rf", "b.lownoise", "b.midnoise", "b.hinoise")

# ufrm <- formula(syn.ind ~ v1 + v2 + v3 + v4 + v5)

ufrm <- formula(syn.ind ~ v1 + v2 + v3 + v4 + v5)

ufrm.cross <- formula(syn.ind ~ v1 + v2 + v3 + v4 + v5 + 
                  v1:v2 + v1:v3 + v1:v4 + v1:v5 +
                  v2:v3 + v2:v4 + v2:v5 +
                  v3:v4 + v3:v5 +
                  v4:v5)

ufrm.squared <- formula(syn.ind ~ v1 + v2 + v3 + v4 + v5 +
                          I(v1^2) + I(v2^2) + I(v3^2) + I(v4^2) + I(v5^2))
ufrm.sqcube <- formula(syn.ind ~ v1 + v2 + v3 + v4 + v5 +
                         I(v1^2) + I(v2^2) + I(v3^2) + I(v4^2) + I(v5^2) +
                         I(v1^3) + I(v2^3) + I(v3^3) + I(v4^3) + I(v5^32))

ufrm.xsq <- formula(syn.ind ~ v1 + v2 + v3 + v4 + v5 +
                  v1:v2 + v1:v3 + v1:v4 + v1:v5 +
                          v2:v3 + v2:v4 + v2:v5 +
                                  v3:v4 + v3:v5 +
                                          v4:v5 +
                  I(v1^2) + I(v2^2) + I(v3^2) + I(v4^2) + I(v5^2))

mod1 <- glm(ufrm, data = stack %>% filter(ftype %in% c("a", ftypes[1])) %>% mutate(syn.ind=ifelse(ftype=="a", 0, 1)),
            family = "binomial")
# summary(mod1)
utility(mod1)

uf <- function(ftype.in, frm){
  frmname <- deparse(substitute(frm))
  mod <- glm(frm, 
             data = stack %>% 
               filter(ftype %in% c("a", ftype.in)) %>% 
               mutate(syn.ind=ifelse(ftype=="a", 0, 1)),
             family = "binomial")
  df <- tibble(ftype=ftype.in, frmname=frmname, utility=utility(mod))
  return(df)
}

udf <- ldply(ftypes, uf, ufrm)
udf.cross <- ldply(ftypes, uf, ufrm.cross)
udf.squared <- ldply(ftypes, uf, ufrm.squared)
udf.xsq <- ldply(ftypes, uf, ufrm.xsq)
udf.sqcube <- ldply(ftypes, uf, ufrm.sqcube)

udf1 <- bind_rows(udf, udf.cross) %>%
  spread(frmname, utility)

udf1 %>% arrange(-ufrm.cross) %>% kable(digits=3)




#.. show all measures ----
d1 <- dmin %>% 
  group_by(ftype=afile) %>%
  summarise(dist.mdn=median(dist), dist.sd=sd(dist))

reg1 <- stats %>%
  dplyr::select(ftype, sigma, r.squared) %>%
  arrange(ftype!="a", -r.squared) %>%
  dplyr::select(ftype, r.squared, sigma) %>%
  left_join(coefdf) 

bin1 <- bincomp %>%
  mutate(vname=paste0("ssd.", vname)) %>%
  spread(vname, ssd) %>%
  mutate(ssd.mean=(ssd.v1 + ssd.v2 + ssd.v3 + ssd.v4 + ssd.v5) / 5) %>%
  arrange(ssd.mean)

d1 %>% dplyr::select(-dist.sd) %>%
  left_join(cordiffs %>% dplyr::select(ftype, cor.adiff.mdn=adiff.mdn)) %>%
  left_join(bin1 %>% dplyr::select(ftype, ssd.median)) %>%
  left_join(reg1 %>% dplyr::select(ftype, sigma)) %>%
  left_join(udf1 %>% dplyr::select(ftype, ufrm.cross)) %>%
  arrange(dist.mdn) %>%
  kable(digits=3)


#.. kd plots ----
stack %>%
  gather(vname, value, -ftype) %>%
  mutate(ftype=factor(ftype, levels=c("a", "b.lownoise", "b.cart", "b.rf", "b.midnoise", "b.hinoise"))) %>%
  ggplot(aes(value, colour=ftype)) +
  geom_density(size=1) +
  theme_bw() +
  ggtitle("Kernel densities") +
  theme(plot.title=element_text(size=12)) +
  facet_wrap(~vname, ncol = 2, scales="free")

# theme(axis.text.x=element_text(angle=45, size=10, hjust=1, colour="black"))



# ---- BREAK ----






#****************************************************************************************************
#                test data ####
#****************************************************************************************************

#mu <- rep(0,4)

sigma <- matrix(.7, nrow=4, ncol=4) + diag(4)*.3
sigma

mvals <- c(1, 2, 3, 4)
nrows <- 100

set.seed(1234)
a <- rmvnorm(n=nrows, mean=mvals, sigma=sigma)
colnames(a) <- paste0("c", 1:4)
a
colMeans(a)

# create a synthetic a
adf <- data.frame(a)
vs <- names(adf)
methods <- rep("cart", length(vs))
synx <- syn(adf, method=methods, visit.sequence=vs, seed=1234)
a.syn <- synx$syn
a.syn

# add random noise to a
nvals <- nrow(a) * ncol(a)
set.seed(1234)
noise <- matrix(rnorm(nvals, mean=0, sd=.4), ncol=ncol(a))
noise
noise.adj <- noise %*% diag(mvals)
colMeans(noise.adj)
apply(noise.adj, 2, sd)

b <- a + noise.adj
b - a
b

a.syn <- a + 2*noise.adj

# get mean, sd
# apply(a, 2, mean)
colMeans(a)
colMeans(b)
colMeans(a.syn)

# sd
apply(a, 2, sd)
apply(b, 2, sd)
apply(a.syn, 2, sd)


# get a and b correlations
cor(a)
cor(b)
cor(a.syn)
cdiffs <- cor(b) - cor(a)
cdiffs

# now compute distances across files - make df's and scale
# scale
df <- bind_rows(as.data.frame(a) %>% mutate(file="a"),
                as.data.frame(b) %>% mutate(file="b"),
                as.data.frame(a.syn) %>% mutate(file="a.syn"))
df

df2 <- as.data.frame(scale(df[, -5])) %>%
  mutate(file=df$file)
df2

df2 <- df %>%
  group_by(file) %>%
  do(as.data.frame(scale(.[, -5]))) %>%
  mutate(rn=row_number(), lab=paste0(file, "_", rn)) %>%
  ungroup %>%
  dplyr::select(file, rn, lab, everything())
ht(df2)

m <- as.matrix(df2[, -c(1:3)])
rownames(m) <- df2$lab # paste0("r", 1:rows)
# colnames(m) <- paste0("c", 1:cols)
ht(m)
d <- dist(m)
dm <- as.matrix(d)
ht(dm)

dmdf <- as_tibble(dm, rownames="lab") %>%
  gather(altrec, dist, -lab) %>%
  separate(lab, c("file1", "rn1"), sep="_") %>%
  separate(altrec, c("file2", "rn2"), sep="_") %>%
  mutate(rn1=as.numeric(rn1), rn2=as.numeric(rn2)) %>%
  filter(file1=="a", file2!="a")
ht(dmdf)
count(dmdf, file1, file2)

d.min <- dmdf %>%
  group_by(rn1, file2) %>%
  arrange(dist) %>%
  filter(row_number()<=1) %>%
  ungroup %>%
  arrange(file1, rn1, dist)

d.comp <- d.min %>%
  dplyr::select(-rn2) %>%
  spread(file2, dist) %>%
  mutate(worst=ifelse(a.syn > b, "a.syn", "b"))
count(d.comp, worst)


d.min %>%
  ggplot(aes(file2, dist, colour=file2)) +
  geom_boxplot()


df %>%
  ggplot(aes(c2, colour=file)) +
  geom_density(size=1.5)

# now make a function of it ----
f <- function(nrows=100, nvars=4, sd.noise=.4, pcor=.7, seed=1234){
  sigma <- matrix(pcor, nrow=nvars, ncol=nvars)
  diag(sigma) <- 1
  mvals <- c(1, 2, 3, 4)
  
  set.seed(seed)
  a <- rmvnorm(n=nrows, mean=mvals, sigma=sigma)
  colnames(a) <- paste0("c", 1:nvars)
  
  set.seed(seed)
  noise <- matrix(rnorm(nrows * nvars, mean=0, sd=sd.noise), ncol=nvars)
  noise.adj <- noise %*% diag(mvals)
  b <- a + noise.adj
  
  # create a synthetic a
  adf <- data.frame(a)
  vs <- names(adf)
  methods <- rep("cart", length(vs))
  synx <- syn(adf, method=methods, visit.sequence=vs, seed=seed)
  a.syn <- synx$syn
  
  # now compute distances across files - make dfs and scale
  # scale
  df <- bind_rows(as.data.frame(a) %>% mutate(file="a"),
                  as.data.frame(b) %>% mutate(file="b"),
                  as.data.frame(a.syn) %>% mutate(file="a.syn"))
  df2 <- df %>%
    group_by(file) %>%
    do(as.data.frame(scale(.[, -ncol(.)]))) %>%
    mutate(rn=row_number(), lab=paste0(file, "_", rn)) %>%
    ungroup %>%
    dplyr::select(file, rn, lab, everything())
  
  # now that we are prepped for distances, do the calculations
  m <- as.matrix(df2[, -c(1:3)])
  rownames(m) <- df2$lab
  d <- dist(m)
  dm <- as.matrix(d)
  
  dmdf <- as_tibble(dm, rownames="lab") %>%
    gather(altrec, dist, -lab) %>%
    separate(lab, c("file1", "rn1"), sep="_") %>%
    separate(altrec, c("file2", "rn2"), sep="_") %>%
    mutate(rn1=as.numeric(rn1), rn2=as.numeric(rn2)) %>%
    filter(file1=="a", file2!="a")
  
  # finally, retrieve the minimum distances
  dmin <- dmdf %>%
    group_by(rn1, file2) %>%
    arrange(dist) %>%
    filter(row_number()<=1) %>%
    ungroup %>%
    arrange(file1, rn1, dist)
  
  return(list(df=df, dmin=dmin))
}

dlist <- f(nrows=100, nvars=4, sd.noise=.4, pcor=.7, seed=1234)

f2 <- function(seed.in, sd.noise=.4){
  f(seed=seed.in)$dmin %>% mutate(seed=seed.in)
}

f3 <- function(seed.in, sd.noise=.4){
  f(seed=seed.in)$df %>% mutate(seed=seed.in)
}
f2(1)

dmin <- ldply(20:25, f2, sd.noise=.01)
ht(dmin)

d.comp <- dmin %>%
  dplyr::select(-rn2) %>%
  spread(file2, dist) %>%
  mutate(worst=ifelse(a.syn > b, "a.syn", "b"))
count(d.comp, worst)
count(d.comp, seed, worst) %>% 
  spread(worst, n) 


d.min %>%
  ggplot(aes(file2, dist, colour=file2)) +
  geom_boxplot()

df.all <- ldply(20:25, f3, sd.noise=.01)
ht(df.all)

# df.all %>%
#   group_by(seed, file) %>%
  
