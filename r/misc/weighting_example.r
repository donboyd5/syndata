


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

library("nloptr")

library("Matrix")

library("knitr")

library("broom")

library("SimMultiCorrData")


#****************************************************************************************************
#                Includes ####
#****************************************************************************************************
#****************************************************************************************************
source("./r/includes/globals_system_specific_boyd.r") # use a different version of this file if changing systems
source("./r/includes/globals_other.r")

source("./r/includes/functions_general.r")
source("./r/includes/functions_target_setup_and_analysis.r")
source("./r/includes/functions_ipopt.r")

# functions specific to the weighting from scratch approach:
source("./r/includes/functions_weight_from_scratch.r")



#****************************************************************************************************
#                functions ####
#****************************************************************************************************
# define characteristics

# correlation matrix
msize <- 3
cormat <- matrix(0, nrow=msize, ncol=msize)
cormat[upper.tri(sigma)] <- c(-.5, .4, .6)
cormat <- cormat + t(cormat)
diag(cormat) <- 1
cormat

# make positive definite
cormat.pd <- nearPD(cormat, corr=TRUE)$mat %>% as.matrix # to be safe
cormat
cormat.pd
(cormat.pd - cormat) %>% round(3)


covmat <- function(cormat, sdvec){
  b <- sdvec %*% t(sdvec)
  covmat <- b * cormat
  dimnames(covmat) <- dimnames(cormat)
  return(covmat)
}


# now make 2 covariance matrices
gvals1 <- c(100, 50, 20) # means of variables for grp1
sd1 <- .05 * gvals1
gvals2 <- c(150, 65, 25) # means of variables for grp2
sd2 <- .05 * gvals2

covmat1 <- covmat(cormat.pd, sd1)
covmat2 <- covmat(cormat.pd, sd2)


nrows <- 5000
seedval <- 1234

# create the a data file
set.seed(seedval)
mat1 <- rmvnorm(n=nrows / 2, mean=gvals1, sigma=covmat1)
mat2 <- rmvnorm(n=nrows / 2, mean=gvals2, sigma=covmat2)
mat <- rbind(mat1, mat2)
mat

# create a data frame of 3 multivariate normal random variables with means and covariance as defined above
a <- mat %>%
  as_tibble(.name_repair=function(x) paste0("v", 1:3)) %>%
  mutate(ftype="a", 
         grp=factor(c(rep(1, nrows/2), rep(2, nrows/2))),
         wt=c(seq(10, 40, length.out = nrows/2), seq(20, 80, length.out = nrows/2)),
         v4=v1 + v2 - v3 / 2) %>%
  dplyr::select(ftype, grp, wt, everything())

ht(a)



#****************************************************************************************************
#                describe the data ####
#****************************************************************************************************
cor(a %>% dplyr::select(-ftype, -grp))

# histograms
a %>%
  gather(variable, value, -ftype, -grp) %>%
  ggplot(aes(x=value,  y = ..density.., fill=grp)) +
  geom_histogram(bins=20) +
  facet_wrap(~variable, scales = "free")

# get weighted sums as targets
targets.grp <- a %>%
  group_by(grp) %>%
  mutate(wt.sum=1) %>%
  summarise_at(vars(wt.sum, v1, v2, v3, v4), ~sum(. * wt))
targets.grp

targets <- a %>%
  mutate(wt.sum=1) %>%
  summarise_at(vars(wt.sum, v1, v2, v3, v4), ~sum(. * wt))
targets

obj <- function(w, inputs){
  obj <- (sum(w) - inputs$targets$wt.sum)^2 +
    (sum(w * inputs$data$v1) - inputs$targets$v1)^2 +
    (sum(w * inputs$data$v2) - inputs$targets$v2)^2 +
    (sum(w * inputs$data$v3) - inputs$targets$v3)^2 +
    (sum(w * inputs$data$v4) - inputs$targets$v4)^2
  return(obj)
}


grad <- function(w, inputs){
  # grad={2 * coeff * priority.weight * (calc - target)} / {scale^2})
  # g[1] <- 2*v1[1]*(calc1 - t1) + 2*v2[2]*(calc2 - t2)
  # ...
  # g[40]
  g <- rep(0, inputs$ngrads)
  diff <- rep(0, 5)
  diff[1] <- sum(w) - inputs$targets$wt.sum
  diff[2] <- sum(w * inputs$data$v1) - inputs$targets$v1
  diff[3] <- sum(w * inputs$data$v2) - inputs$targets$v2
  diff[4] <- sum(w * inputs$data$v3) - inputs$targets$v3
  diff[5] <- sum(w * inputs$data$v4) - inputs$targets$v4
  
  for(i in 1:inputs$ngrads){
    g[i] <- 2 * diff[1] +
      2 * inputs$data$v1[i] * diff[2] +
      2 * inputs$data$v2[i] * diff[3] +
      2 * inputs$data$v3[i] * diff[4] +
      2 * inputs$data$v4[i] * diff[5]
  }
  return(g)
}

inputs <- list()
inputs$ngrads <- nrows
inputs$targets <- targets
inputs$data <- a

xlb <- rep(0.5*min(a$wt), nrows)
xub <- rep(2*max(a$wt), nrows)
x0 <- (xlb + xub) / 2

# grad(x0, inputs)

# set.seed(seedval)
# x0 <- runif(nrows, min = 0, max = 2*max(a$wt))

# t1 <- mma(x0, 
#           fn=obj,
#           gr=grad,
#           lower=xlb, upper=xub,
#           nl.info = FALSE, 
#           inputs=inputs)

# obj(t1$par, inputs)
# obj(a$wt, inputs)

opts <- list("algorithm"="NLOPT_LD_MMA",
             "xtol_rel"=1.0e-8,
             "maxeval"=500)

opts <- list("algorithm"="NLOPT_LD_LBFGS",
             "xtol_rel"=1.0e-8,
             "maxeval"=500)

t.a <- proc.time()
t3 <- nloptr(x0, 
             eval_f=obj,
             eval_grad_f = grad,
             lb = xlb, ub = xub,
             opts = opts, inputs=inputs)
t.b <- proc.time()
t.b - t.a

names(t3)
t3$termination_conditions
t3$message
t3$iterations
t3$objective
t3$status
t3$solution %>% ht
# t3$options

wt.new <- t3$solution


a2 <- a %>% 
  mutate(wt2=wt.new)
ht(a2)

a2 %>%
  gather(wtype, weight, wt, wt2) %>%
  group_by(ftype, wtype) %>%
  mutate(wt.sum=1) %>%
  summarise_at(vars(wt.sum, v1, v2, v3, v4), ~sum(. * weight))

a2 %>%
  gather(wtype, weight, wt, wt2) %>%
  group_by(ftype, grp, wtype) %>%
  mutate(wt.sum=1) %>%
  summarise_at(vars(wt.sum, v1, v2, v3, v4), ~sum(. * weight))

a2 %>%
  ggplot(aes(wt, wt2, colour=grp)) + geom_point()

inputs$targets

# tout <- t1
# names(tout)
# tout$par
# tout$value
# tout$iter
# tout$convergence
# tout$message
# 
# wt.new <- t1$par


#****************************************************************************************************
#                synthesize 2 new data files ####
#****************************************************************************************************





