


library("magrittr")
library("plyr") # needed for ldply; must be loaded BEFORE dplyr
library("tidyverse")
options(tibble.print_max = 60, tibble.print_min = 60) # if more than 60 rows, print 60 - enough for states

source("./r/includes/globals_system_specific_boyd.r") # use a different version of this file if changing systems
source("./r/includes/globals_other.r")

source("./r/includes/functions_general.r")
library("syntpop")

psums <- get_pufbase_sums()


df <- get_puf.base()
glimpse(df)
names(df) %>% sort

set.seed(3)
df2 <- df %>%
  sample_n(20) %>%
  dplyr::select(pufseqn, e00100, e00200, e00300) %>%
  mutate(row=row_number()) %>%
  dplyr::select(row, e00100, e00200, e00300)
df2
head(df2)

visit.sequence <- names(df2) # setdiff(names(df2), "row")
methods <- c("", "", "cart", "cart")
methods <- c("", "", "norm", "cart")
methods <- c("", "", "cart", "norm")
cbind(visit.sequence, methods)

test <- syn(df2, visit.sequence=visit.sequence, method=methods, seed=1234, m=0, k=nrow(df2), proper=FALSE) # m: 1=32, 2=64, 3=92
test$predictor.matrix
pm <- test$predictor.matrix
pm[, "row"] <- 0
pm

synx <- syn(df2, visit.sequence=visit.sequence, method=methods, seed=1, m=2, k=nrow(df2), proper=FALSE, predictor.matrix = pm) # m: 1=32, 2=64, 3=92
# names(synx) 
# synx$method
# synx$predictor.matrix
syndf <- bind_rows(synx$syn) %>% mutate(m=rep(1:synx$m, each=synx$k))

left_join(syndf %>% dplyr::rename(e00200.syn=e00200, e00300.syn=e00300),
          df2) %>% 
  dplyr::select(row, m, e00100, e00200, e00300, e00200.syn, e00300.syn) %>%
  arrange(row, m)


# no proper
synx <- syn(df2, visit.sequence=visit.sequence, method=methods, seed=1234, m=2, k=nrow(df2), proper=FALSE, predictor.matrix = pm) # m: 1=32, 2=64, 3=92
syndf <- bind_rows(synx$syn) %>% mutate(m=rep(1:synx$m, each=synx$k))
noprop <- left_join(syndf %>% dplyr::rename(e00200.syn=e00200, e00300.syn=e00300),
               df2) %>% 
  dplyr::select(row, m, e00100, e00200, e00200.syn, e00300, e00300.syn) %>%
  mutate(type="noprop") %>%
  arrange(row, m)
# proper
synx <- syn(df2, visit.sequence=visit.sequence, method=methods, seed=1234, m=2, k=nrow(df2), proper=TRUE, predictor.matrix = pm) # m: 1=32, 2=64, 3=92
syndf <- bind_rows(synx$syn) %>% mutate(m=rep(1:synx$m, each=synx$k))
prop <- left_join(syndf %>% dplyr::rename(e00200.syn=e00200, e00300.syn=e00300),
                    df2) %>% 
  dplyr::select(row, m, e00100, e00200, e00200.syn, e00300, e00300.syn) %>%
  mutate(type="prop") %>%
  arrange(row, m)
pnp <- bind_rows(noprop, prop)
write_csv(pnp, "./pnp.csv")



# why is computed 2013 agi different from agi on file?
# get a puf that is in same order as synthpop4
puf1 <- get_puf.base()
glimpse(puf1)

puf2 <- puf1 %>% select(pufseqn, wt.puf=wt, e00100)

spuf2 <- stack %>% filter(ftype=="puf") %>% select(pufseqn, wt.spuf=wt, c00100)

tmp <- left_join(puf2, spuf2)
ht(tmp)

f <- function(var, wt){sum(var * wt) / 1e9}
tmp %>%
  mutate(diff=c00100 - e00100,
         pdiff=diff / e00100 * 100,
         apdiff=abs(pdiff)) %>%
  summarise(e00100.s=f(e00100, wt.puf),
            c00100.s=f(c00100, wt.puf),
            diff.s=f(diff, wt.puf),
            dgt0=sum(diff>0),
            dlt0=sum(diff<0),
            deq0=sum(diff==0),
            pdiff=median(pdiff, na.rm=TRUE))



puf3 <- puf2 %>%
  select(seqn, RECID, wt, MARS, XTOT, e00100, e00200)
ht(puf3)
sum(puf3$wt)
sum(puf3$wt * puf3$e00200)

puf3 %>% filter(round(wt, 2)==1425.93, e00200==5870) %>% kable(digits=2)


s4 <- readRDS(paste0(globals$tc.dir, "synthpop4", "_reweighted_stackedfiles.rds"))
ps4 <- s4 %>% 
  filter(ftype=="puf") %>%
  mutate(seqn=row_number()) %>%
  select(seqn, RECID, wt, MARS, XTOT, c00100, e00200)
ht(ps4)
sum(ps4$wt)
sum(ps4$wt * ps4$e00200)


sum(puf3$wt * puf3$e00200) - sum(ps4$wt * ps4$e00200)


sfname <- "synthpop4"
sfname <- "synpuf8"

s8 <- readRDS(paste0(globals$tc.dir, "synpuf8", "_reweighted_stackedfiles.rds"))
s4 <- readRDS(paste0(globals$tc.dir, "synthpop4", "_reweighted_stackedfiles.rds"))

sfname <- "synthpop4"
sfname <- "synpuf8"
stack <- readRDS(paste0(globals$tc.dir, sfname, "_reweighted_stackedfiles.rds"))

puf.vnames %>%
  filter(vname %in% setdiff(names(s8), names(s4)))

puf.vnames %>%
  filter(vname %in% setdiff(names(s4), names(s8)))


glimpse(s4)
count(s4, ftype)

p <- get_puf.base()
glimpse(p)

p2 <- p %>% select(RECID, wt, e00100, e00200)
ht(p2)

ps4 <- s4 %>% 
  filter(ftype=="puf") %>%
  select(RECID, wt, MARS, XTOT, c00100, e00200)
ht(ps4)
sum(ps4$wt)


sfname <- "synthpop4_all"
tmp <- readRDS(paste0(globals$tc.dir, sfname, "_rwprep.rds"))
names(tmp)
glimpse(tmp$tc.base)

df <- tmp$tc.base %>% filter(ftype=="puf")



# 400 secs =6.7 mins
result$objective
result$solution[1:10]


a <- proc.time()
t2 <- mma(x0, fn=eval_f_wtfs, gr = eval_grad_f_wtfs,
          lower=xlb, upper=xub,
          nl.info = FALSE, inputs=inputs)
b <- proc.time()
b - a
names(t2)
# 181 secs=3 mins
t2$value
t2$par[1:10]


library("BB")
a <- proc.time()
t3 <- BBoptim(par=x0, fn=eval_f_wtfs, gr = eval_grad_f_wtfs,
          lower=xlb, upper=xub,inputs=inputs)
b <- proc.time()
b - a
names(t3)
# 76 secs
t3$value
t3$par[1:10]


# NLOPT_LD_MMA good
# NLOPT_LD_CCSAQ NOT AVAIL
# NLOPT_LD_SLSQP NO NO NO
# NLOPT_LD_LBFGS_NOCEDAL NOT AVAIL
# NLOPT_LD_LBFGS ok maybe, but NOT very good
# NLOPT_LD_VAR1 not very good
# NLOPT_LD_VAR2 not very good
# NLOPT_LD_TNEWTON not very good
# NLOPT_LD_TNEWTON_RESTART
# NLOPT_LD_TNEWTON_PRECOND
# NLOPT_LD_TNEWTON_PRECOND_RESTART NO

# NLOPT_LN_NEWUOA_BOUND NO
# NLOPT_LN_COBYLA NO
# NLOPT_LN_NELDERMEAD NO

# NLOPT_GD_STOGO NO NO NO
# NLOPT_GD_STOGO_RAND NO NO NO
# NLOPT_GN_ISRES

# use a local solver with these:
# NLOPT_LD_AUGLAG
# NLOPT_LN_AUGLAG_EQ
# NLOPT_LD_AUGLAG_EQ

a <- proc.time()
res0 <- nloptr( x0=x0,
                eval_f=eval_f_wtfs,
                # eval_grad_f=eval_grad_f_wtfs,
                lb = xlb,
                ub = xub,
                opts = list("algorithm"="NLOPT_GN_ISRES",
                            "print_level"=1,
                            "maxeval"=10,
                            "xtol_abs"=1e-10),
                inputs=inputs)
b <- proc.time()
b - a
# names(res0)
res0$message
res0$objective
res0$solution[1:10]
# res0$iterations


local_opts <- list(algorithm = "NLOPT_LD_MMA", xtol_rel= 1.0e-7 )
a <- proc.time()
res1 <- nloptr( x0=x0,
                eval_f=eval_f_wtfs,
                eval_grad_f=eval_grad_f_wtfs,
                lb = xlb,
                ub = xub,
                opts = list("algorithm"="NLOPT_LD_AUGLAG",
                            local_opts = local_opts,
                            "print_level"=1,
                            "maxeval"=100,
                            "xtol_abs"=1e-10),
                inputs=inputs)
b <- proc.time()
b - a
res1$message
res1$objective
res1$solution[1:10]


nloptr.print.options()$algorithm

# nlminb not very good
a <- proc.time()
port <- nlminb(start=x0, objective=eval_f_wtfs, gradient = eval_grad_f_wtfs,
               inputs=inputs,
               scale = 1, control = list(iter.max=100, trace=0), lower = xlb, upper = xub)
b <- proc.time()
b - a
# names(port)
port$message
port$objective
port$par[1:10]



a <- proc.time()
t2a <- mma(t3$par, fn=eval_f_wtfs, gr = eval_grad_f_wtfs,
           lower=xlb, upper=xub,
           nl.info = FALSE, inputs=inputs)
b <- proc.time()
b - a
names(t2a)
t2a$value


# compare the improvement in convergence when bounds are specified
BBoptim(par=p0, fn=rosbkext, lower=0) 

# identical to spg() with defaults
BBoptim(par=p0, fn=rosbkext, method=3, control=list(M=10, trace=TRUE)) 


# check out binning ----

# tmp <- bins %>% filter(vname=="e00650")
vname.in <- "e00650"
df <- stack %>%
  dplyr::select(ftype, value=vname.in)

# construct puf-based bins
bins <- df %>%
  filter(ftype=="puf") %>%
  mutate(ntile=ntile(value, nbins)) %>%
  group_by(ntile) %>%
  summarise(pufn=n(), valmin=min(value), valmax=max(value)) %>%
  group_by(valmin) %>%
  summarise(pufn=sum(pufn)) %>%
  ungroup

brks <- c(-Inf, unique(bins$valmin2), max(bins$valmax), Inf) %>% unique # unique, in case maxval is in the 2nd term

# replace first min with -Inf, and add Inf
# explicitly add - Inf, 0, small positive, and Inf, get unique values, and sort
brks <- c(-Inf, 0, 1e-9, Inf, bins$valmin) %>% unique %>% sort
brks %>% ht

levels(binned$bin) %>% ht
sum(binned$syn)

# use these puf-based breaks to bin both the puf and the syn data and then get relative frequencies
binned <- df %>%
  mutate(bin=cut(value, brks, include.lowest=TRUE, right=FALSE), bin.num=as.numeric(bin)) %>%
  group_by(ftype, bin.num, bin) %>%
  summarise(n=n()) %>%
  spread(ftype, n) %>%
  mutate_at(vars(puf, syn), funs(naz)) %>%
  ungroup %>%
  mutate_at(vars(puf, syn), funs(pct=. / sum(.) * 100)) %>%
  mutate(vname=vname.in) %>%
  dplyr::select(vname, everything()) %>%
  arrange(vname, bin.num)


tmp2 <- getbins("e00650", stack, 1000)

t3 <- tmp2 %>%
  group_by(vname) %>%
  mutate(error=puf_share*100 - syn_share*100,
         e2=error^2)
  summarise(r2=rsq(puf_share, syn_share),
            sse=sse(puf_share*100, syn_share*100)) %>%
  left_join(puf.vnames %>% dplyr::select(vname, vdesc))
t3

