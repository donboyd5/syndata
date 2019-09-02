# !diagnostics off  # to get rid of errant tibble warnings about unknown column

# https://github.com/tidyverse/tibble/issues/450
# or set bad variable to NULL tbl %>% has_name("c")

#****************************************************************************************************
#                Libraries ####
#****************************************************************************************************
library("magrittr")
library("plyr") # needed for ldply; must be loaded BEFORE dplyr
library("tidyverse")
options(tibble.print_max = 60, tibble.print_min = 60) # if more than 60 rows, print 60 - enough for states
# ggplot2 tibble tidyr readr purrr dplyr stringr forcats

# library("readxl") # readxl, for .xls and .xlsx files
# library("haven") # haven, for SPSS, SAS and Stata files
# library("vctrs")
# library("precis")

# library("grDevices")
library("knitr")

# library("nloptr")
library("survey")
library("ipoptr")

# devtools::install_github("donboyd5/btools")
library("btools") # library that I created (install from github)


#****************************************************************************************************
#                Includes ####
#****************************************************************************************************
source("./r/includes/globals_system_specific_boyd.r") # use a different version of this file if changing systems
source("./r/includes/globals_other.r")

source("./r/includes/functions_general.r")

source("./r/includes/functions_ipopt.r")

# functions specific to the weighting from scratch approach:
source("./r/includes/functions_weight_from_scratch.r")


#****************************************************************************************************
#                globals ####
#****************************************************************************************************
gd <- "c:/Users/donbo/Google Drive/"

pufd <- paste0(gd, "synpuf/puf/")
synd <- paste0(gd, "synpuf/syntheses/")

puf_fn <- "puf2011.csv"
synrf_fn <- "synpuf20.csv"
synrfnd_fn <- "synpuf20_no_disclosures.csv"


#****************************************************************************************************
#                ONETIME PER RUN - set PATH ####
#****************************************************************************************************
# note - might just be able to add Anaconda3/Scripts to the system path
# Temporarily change the path environment variable so tc will work from R
Sys.getenv("PATH")
# # Caution: these seem to be added to PATH when Anaconda powershell is open
p1 <- "C:/ProgramData/Anaconda3"
p2 <- paste0(p1, "/", "bin")
p3 <- paste0(p1, "/", "condabin")
p4 <- paste0(p1, "/", "Scripts")
p5 <- paste0(p1, "/", "Library/bin")
p6 <- paste0(p1, "/", "Library/usr/bin")
p7 <- paste0(p1, "/", "Library/mingw-w64/bin")
Sys.setenv(
  PATH = paste(p1, p2, p3, p4, p5, p6, p7,
               Sys.getenv("PATH"),
               sep = ";"
  )
)
shell("echo %PATH% ", intern= TRUE)


#****************************************************************************************************
#                functions - general ####
#****************************************************************************************************
ns <- function(df) {names(df) %>% sort } # names sort 

combine <- function(prefix, suffix) {as.vector(t(outer(prefix, suffix, paste, sep="_")))}

flag_dups <- function(df) {duplicated(df) | duplicated(df, fromLast=TRUE)} # this gets all dups

showvars <- function(sort="vname", usevars=syn_info$vname){
  if(sort=="vname") temp <- syn_info %>% arrange(vname) else
    if(sort=="sum") temp <- syn_info %>% arrange(-abs(sum))
    
    temp %>% 
      filter(vname %in% usevars) %>%
      kable(digits=0, format.args = list(big.mark=","))
}


#****************************************************************************************************
#                functions for constraints ####
#****************************************************************************************************

nnz <- function(x, wt0) {(x != 0) * 1}
nz <- function(x, wt0) {(x == 0) * 1}
npos <- function(x, wt0) {(x > 0) * 1}
nneg <- function(x, wt0) {(x < 0) * 1}

sumnz <- function(x, wt0) {(x != 0) * x}
sumpos <- function(x, wt0) {(x > 0) * x}
sumneg <- function(x, wt0) {(x < 0) * x}


#****************************************************************************************************
#                functions for optimization ####
#****************************************************************************************************

calc_constraints <- function(weights, data, constraint_vars){
  # weights: numeric vector
  # data: df with 1 column per constraint variable (plus other info); 1 row per person
  colSums(weights * select(data, constraint_vars))
}

pdiff <- function(weights, data, constraints) {
  # percent difference between calculated constraints and targets
  calculated_constraints <- calc_constraints(weights, data, names(constraints))
  (calculated_constraints - constraints) / constraints * 100
}


objfn <- function(weights, data, constraints, p=2){
  # this objfn scales the data and constraint values used in the obj calculation
  calculated_constraints <- calc_constraints(weights, data, names(constraints))
  # scale the constraints and the calculations by the constraint values
  diff_to_p <- (calculated_constraints / constraints - constraints / constraints)^p
  obj <- sum(diff_to_p)
  return(obj)
}


define_jac_g_structure_dense <- function(n_constraints, n_variables){
  # list with n_constraints elements
  # each is 1:n_variables
  lapply(1:n_constraints, function(n_constraints, n_variables) 1:n_variables, n_variables)
} 
# eval_jac_g_structure_dense <- define_jac_g_structure_dense(n_constraints=2, n_variables=4)


eval_jac_g_dense <- function(x, inputs){
  # the Jacobian is the matrix of first partial derivatives of constraints (these derivatives may be constants)
  # this function evaluates the Jacobian at point x
  
  # This is the dense version that returns a vector with one element for EVERY item in the Jacobian (including zeroes)
  # Thus the vector has n_constraints * n_variables elements
  # first, all of the elements for constraint 1, then all for constraint 2, etc...
  
  # because constraints in this problem are linear, the derivatives are all constants
  
  # ipoptr requires that ALL functions receive the same arguments, so the inputs list is passed to ALL functions
  
  return(inputs$cc_dense)
}

eval_g_dense <- function(x, inputs){
  unname(calc_constraints(inputs$wt * x, inputs$data, inputs$constraint_vars))
}

# eval_g_dense(x0, inputs)
# eval_jac_g_dense(x0, inputs)
# weights <- data$wt0

# now create the dense hessian
# eval_h_dense <- function(x, obj_factor, hessian_lambda, inputs){
#   # The Hessian matrix has many zero elements and so we set it up as a sparse matrix
#   # We only keep the (potentially) non-zero values that run along the diagonal.
#   
#   # http://www.derivative-calculator.net/
#   # w{x^p + x^(-p) - 2}                                 objective function
#   # w{px^(p-1) - px^(-p-1)}                             first deriv
#   # p*w*x^(-p-2)*((p-1)*x^(2*p)+p+1)                    second deriv
#   
#   # make it easier to read:
#   p <- inputs$p
#   w <- inputs$wt
#   
#   hess <- obj_factor * 
#     { p*w*x^(-p-2) * ((p-1)*x^(2*p)+p+1) }
#   
#   return(hess)
# }

# get_tols <- function(tolerances, constraints, default=0){
#   # default should be in 0:Inf
#   # tolerances is data frame with cols variable and tol
#   # constraints is 1-row data frame with cols for constraints
#   # row has constraint values
#   # return: named vector with a tolerance for each constraint
#   tols <- rep(default, ncol(constraints))
#   names(tols) <- names(constraints)
#   
#   var_indexes <- match(names(constraints), tolerances$variable)
#   tol_subset <- tolerances[var_indexes, ]
#   tols[tol_subset$variable] <- tol_subset$tol
#   # tols
#   return(tols)
# }
# get_tols(tolerances, constraints, .01)


calibrate_reweight <- function(weights, data, constraints){
  # id = ~0 means no clusters
  sdo <- svydesign(id = ~0, weights = weights, data = data)
  
  pop_totals <- constraints %>% 
    gather(vname, value) %>%
    deframe
  
  eps <- abs(pop_totals * .005)
  
  # pop_totals <- pop_totals[-1] # djb skipping this step
  
  frm_string <- paste0("~ ", paste(names(pop_totals), collapse=" + "), " - 1")
  frm <- as.formula(frm_string)
  
  # note that logit calibration requires finite bounds
  fn <- c("linear","raking","logit")
  calibrated <- calibrate(sdo, 
                          formula=frm, 
                          population=pop_totals, 
                          calfun=fn[1], 
                          # eta=weights, 
                          # bounds=c(-Inf,Inf),
                          # bounds=c(0.1,Inf),
                          bounds=c(0.00001, 1000),
                          bounds.const=FALSE, # default FALSE
                          # trim=c(0, Inf),
                          epsilon=eps, # default 1e-7
                          sparse=TRUE, 
                          force=FALSE)
  
  return(calibrated)
}


uni_opt <- function(method, wts0, data, constraints, objfn=NULL, ...){
  # call any of the methods and get a uniform return value
  
  methods <- c("simann", "mma", "calibrate", "ipopt")
  if(!method %in% methods){
    err_msg1 <- paste0("STOP: Method ", method, " not supported. Method must be one of:")
    err_msg2 <- paste(methods, collapse=", ")
    print(err_msg1)
    print(err_msg2)
    return(NULL)
  }
  
  args <- list(...) # only gets the dots values; match.call() gets all names, but not the values
  if(hasArg(maxiter)) maxiter <- args$maxiter else maxiter <- 1000 # always define even if not needed
  
  if(method=="simann"){
    set.seed(seed)
    a <- proc.time()
    sa_result <- sim_ann(objfn, wts0, data, constraints, 
                         p=2, niter = maxiter, step = 0.1, 
                         phase_iter=500, max_recshare=.3, max_wtshare=1)
    b <- proc.time()
    
    result <- list()
    result$method <- method
    result$obj <- sa_result$best_obj
    result$weights <- sa_result$best_weights
    result$iter <- sa_result$iterations
    result$elapsed_seconds <- unname((b - a)[3])
    result$opt_object <- sa_result
    
  } else if(method=="mma"){
    a <- proc.time()
    mma_result <- mma_nloptr(objfn, wts0, data, constraints, niter=maxiter)
    b <- proc.time()
    
    result <- list()
    result$method <- method
    result$obj <- mma_result$value
    result$weights <- mma_result$par
    result$iter <- mma_result$iter
    result$elapsed_seconds <- unname((b - a)[3])
    result$opt_object <- mma_result
    
  } else if(method=="calibrate"){
    a <- proc.time()
    calib_result <- calibrate_reweight(wts0, data, constraints)
    b <- proc.time()
    
    result <- list()
    result$method <- method
    result$obj <- objfn(weights(calib_result), data, constraints)
    result$weights <- unname(weights(calib_result))
    result$iter <- 0
    result$elapsed_seconds <- unname((b - a)[3])
    result$opt_object <- NULL
    result$notes <- "(1) obj is NOT from the calibrate function, it is the general objfn; (2) result object not returned - too big"
  }
  return(result)
}


getwts <- function(df, all_constraints, good_constraints){
  ftype.in <- first(df$ftype)
  ugroup.in <- first(df$ugroup)
  
  convars <- good_constraints %>%
    filter(ftype==ftype.in,
           ugroup==ugroup.in) %>%
    .[["good_constraint"]]
  
  constraints <- all_constraints %>%
    ungroup %>%
    filter(ftype=="puf", # IMPORTANT!
           ugroup==ugroup.in) %>%
    select(convars)
  
  calib_result <- calibrate_reweight(df$wt0, df[, convars], constraints)
  
  fname <- paste0("calib", ugroup.in, ".rds")
  saveRDS(calib_result, paste0("d:/temp/", fname))
  df$wt1 <- unname(weights(calib_result))
  return(df)
}



#****************************************************************************************************
#                Ipopt functions ####
#****************************************************************************************************
ipopt_reweight <- function(weights, data, constraint_df, ugroup.in){
  # arguments for all functions passed to ipoptr must be x, inputs
  inputs <- list()
  inputs$p <- 2
  inputs$wt <- weights
  inputs$data <- data
  inputs$constraint_vars <- names(constraint_df %>% select(-variable))
  inputs$cc_dense <- c(as.matrix(data[, inputs$constraint_vars] * weights)) # flattens the cc matrix
  
  xlb <- rep(0, nrow(data))
  xub <- rep(1000, nrow(data))
  x0 <- rep(1, nrow(data))
  
  constraints <- constraint_df %>% filter(variable=="constraint") %>% select(-variable) # do I need this??
  
  clb <- constraint_df %>% filter(variable=="clb") %>% select(-variable) %>% t %>% as.vector
  cub <- constraint_df %>% filter(variable=="cub") %>% select(-variable) %>% t %>% as.vector
  
  eval_jac_g_structure_dense <- define_jac_g_structure_dense(
    n_constraints=ncol(data[, names(constraints)]), 
    n_variables=nrow(data[, names(constraints)]))
  
  eval_h_structure <- lapply(1:length(inputs$wt), function(x) x) # diagonal elements of our Hessian
  
  # ma86 was faster in one test I did
  opts <- list("print_level" = 0,
               "file_print_level" = 5, # integer
               "linear_solver" = "ma86", # mumps pardiso ma27 ma57 ma77 ma86 ma97
               "max_iter"=200,
               "nlp_scaling_max_gradient" = 1e-3, # 10 improved things -- default 100
               "obj_scaling_factor" = 1, # default 1
               # "derivative_test"="first-order",
               # "derivative_test_print_all"="yes",
               "output_file" = "syntarget.out")
  
  result <- ipoptr(x0 = weights,
                   lb = xlb,
                   ub = xub,
                   eval_f = eval_f_xtop, # arguments: x, inputs
                   eval_grad_f = eval_grad_f_xtop,
                   eval_g = eval_g_dense, # constraints LHS - a vector of values
                   eval_jac_g = eval_jac_g_dense,
                   eval_jac_g_structure = eval_jac_g_structure_dense,
                   eval_h = eval_h_xtop, # the hessian is essential for this problem
                   eval_h_structure = eval_h_structure,
                   constraint_lb = clb,
                   constraint_ub = cub,
                   opts = opts,
                   inputs = inputs)
  
  return(result)
}


getwts_ipopt <- function(df, all_constraints, bounds){
  # get weights using ipoptr
  ftype.in <- first(df$ftype)
  ugroup.in <- first(df$ugroup)
  
  # create a 3-record data frame from bounds, of constraints, clb, and cub
  # so that all are in exactly the same order
  constraint_df <- bounds %>%
    filter(ftype==ftype.in, ugroup==ugroup.in)%>%
    select(good_constraint, constraint=target, clb, cub) %>%
    gather(variable, value, -good_constraint) %>%
    spread(good_constraint, value)
  convars <- names(constraint_df %>% select(-variable))
  
  print(paste0("Starting group: ", ugroup.in))
  ipopt_result <- ipopt_reweight(df$wt0, df[, convars], constraint_df, ugroup.in)
  
  if(ipopt_result$status != 0){
    print(ipopt_result$message)
    print("\n")
    fname <- paste0("ipopt_bad_", ugroup.in, ".rds")
    saveRDS(ipopt_result, paste0("d:/temp/", fname))
    
    tmp <- read_file("syntarget.out")
    fname <- paste0("d:/temp/", "ipopt_bad_", ugroup.in, ".out")
    write_file(tmp, fname)
  }
  
  df$wt1 <- ipopt_result$solution * df$wt0
  return(df)
}





#****************************************************************************************************
#                get the puf, syn, and syn-nodups file ####
#****************************************************************************************************
syn_info <- readRDS("./data/syn_info.rds")


showvars()
showvars("sum")

syn <- read_csv(paste0(synd, synrf_fn))
glimpse(syn)
ns(syn) # we have RECID and S006 on this file

syn_nd <- read_csv(paste0(synd, synrfnd_fn))
glimpse(syn_nd)
ns(syn_nd)

puf <- read_csv(paste0(pufd, puf_fn))
ns(puf)
puf2 <- puf[, names(syn)] %>% filter(MARS!=0) # this is the one we'll use

identical(names(syn), names(puf2))
identical(names(syn_nd), names(puf2))


#****************************************************************************************************
#                define disclosure-review variables, save with descriptions and puf sums ####
#****************************************************************************************************
# df1 <- read_csv(paste0(synd, synrfnd_fn), n_max=0)
# glimpse(df1)
# ns(df1)
nodup_vars <- setdiff(names(syn_nd), c("RECID", "S006"))
nodup_vars # 64 variables

puf_vnames <- get_puf_vnames() %>% select(vname, vdesc)
tc_vnames <- tribble(
  ~vname, ~vdesc,
  "c17000", "Sch A: Medical expenses deducted (calculated)",
  "c18300", "Sch A: State and local taxes deducted (calculated)",
  "c21060", "Itemized deductions before phase-out (zero for non-itemizers) (calculated)",
  "standard", "standard Standard deduction (zero for itemizers) (calculated)",
  "c04800", "Regular taxable income (calculated)",
  "taxbc", "regular tax before credits (calculated)",
  "c09600", "Alternative Minimum Tax (AMT) liability (calculated)",
  "c05800", "taxbc plus AMT liability (calculated)"
  )
var_vnames <- bind_rows(puf_vnames, tc_vnames)


f <- function(x, wt) sum(x * wt / 100)
syn_info <- puf2 %>%
  select(nodup_vars, S006) %>%
  summarise_at(vars(-S006), list(~f(., S006))) %>%
  gather(vname, sum) %>%
  left_join(puf_vnames %>% 
              mutate(vname=str_to_upper(vname)) %>% 
              group_by(vname) %>%
              filter(row_number()==1) %>%
              ungroup %>%
              select(vname, vdesc))
saveRDS(syn_info, "./data/syn_info.rds")





#****************************************************************************************************
#                ONETIME: stack the 3 files, prepare for tax calculator, and run  ####
#****************************************************************************************************
stack <- bind_rows(puf2 %>% mutate(ftype="puf"),
                   syn %>% mutate(ftype="syn"),
                   syn_nd %>% mutate(ftype="syn_nd")) %>%
  setNames(change_case(names(.))) %>% # Tax-Calculator expects mostly lower-case names
  # do(impose_variable_rules(.)) %>% # not needed for synpuf5 and later
  do(prime_spouse_splits(.)) %>%
  mutate(RECID_original=RECID,
         RECID=row_number()) # RECID is the one way we can be SURE we match the output
write_csv(stack, paste0(globals$tc.dir, "pufsyn.csv"))
count(stack, ftype, MARS)


set.seed(1234)
stack_samp <- stack %>%
  group_by(ftype) %>%
  sample_n(1000)
write_csv(stack_samp, paste0(globals$tc.dir, "pufsyn_samp.csv"))


# pick one of these
#    taxplan.fn -- character variable holding the name (without directory) of the taxplan json file e.g., "rate_cut.json"

tc.fn <- "pufsyn.csv"
# tc.fn <- "pufsyn_samp.csv"

cmd <- tc.wincmd(tc.fn, globals$tc.dir, globals$tc.cli, taxyear=2013, taxplan.fn=NULL, taxplans.dir=globals$taxplans.dir)
cmd    # a good idea to look at the command

a <- proc.time()
system(cmd) # CAUTION: this will overwrite any existing output file that had same input filename!
# consider system2
b <- proc.time()
b - a  # it can easily take 5-10 minutes depending on the size of the input file


#****************************************************************************************************
#                ONETIME: get the results, combine as needed  ####
#****************************************************************************************************
# c00100 AGI
# c17000 Sch A: Medical expenses deducted
# c18300 Sch A: State and local taxes deducted
# c21060 Itemized deductions before phase-out (zero for non-itemizers)
# standard Standard deduction (zero for itemizers)
# c04800 Regular taxable income
# taxbc regular tax before credits
# c09600 Alternative Minimum Tax (AMT) liability
# c05800 taxbc plus AMT liability

stack <- read_csv(paste0(globals$tc.dir, "pufsyn.csv"))

tc.outfn <- paste0(tools::file_path_sans_ext(tc.fn),
                   "-", 13, "-#-#-#.csv")
tc.output <- read_csv(paste0(globals$tc.dir, tc.outfn),
                      col_types = cols(.default= col_double()),
                      n_max=-1)
glimpse(tc.output)
ns(tc.output)

tcvars <- c("RECID", "c00100", "taxbc", "c09600", "c05800") # taxcalc vars we want to get
outfile <- "stack_for_weighting.rds"
stack2 <- stack %>%
  left_join(tc.output %>% select(tcvars))
saveRDS(stack2, paste0(globals$tc.dir, outfile))


#****************************************************************************************************
#                START: weight the files ####
#****************************************************************************************************
ffw <- readRDS(paste0(globals$tc.dir, "stack_for_weighting.rds")) # file for weighting
glimpse(ffw)
ffw %>%
  group_by(ftype) %>%
  do(qtiledf(.$s006))

# create a starting-point weight, wt0, that has the right (i.e., puf) total
wt_totals <- ffw %>% 
  group_by(ftype) %>%
  summarise(wt=sum(s006) / 100) %>%
  mutate(ratio=wt[ftype=="puf"] / wt)
wt_totals

ffw <- ffw %>%
  left_join(wt_totals %>% select(ftype, ratio)) %>%
  mutate(wt0=s006 / 100 * ratio,
         n=1 / wt0, pop=1)

ffw %>%
  group_by(ftype) %>%
  summarise(wt_m=sum(wt0) / 1e6,
            c00100_b=sum(c00100 * wt0) / 1e9,
            taxbc_b=sum(taxbc * wt0) / 1e9) %>%
  as.data.frame


#****************************************************************************************************
#               automate ####
#****************************************************************************************************
#.. 1. Define subsets that have approximately 500-1000 records ----
puf <- ffw %>%
  filter(ftype=="puf")
count(puf, MARS)

ntarget <- 2000
nmin <- 1500

ngroups <- puf %>%
  mutate(agige0=ifelse(c00100 >= 0, 1, 0)) %>%
  select(agige0, MARS, c00100) %>%
  group_by(agige0, MARS) %>%
  arrange(c00100) %>%
  mutate(ngroups=n() / ntarget,
         grp=ntile(n=ngroups)) %>%
  group_by(agige0, MARS, grp) %>%
  summarise(n=n(), 
            agimin=min(c00100), 
            agimax=max(c00100)) %>%
  mutate(grp=ifelse(n < nmin, 
                    ifelse(grp==1, grp + 1, grp - 1), # can cause grp numbers to skip
                    grp)) %>%
  group_by(agige0, MARS, grp) %>%
  summarise(n_final=sum(n),
            agimin=min(agimin),
            agimax=max(agimax)) %>%
  group_by(agige0, MARS) %>%
  mutate(grp=row_number(), ngrps=max(grp)) %>% # ensure that grp numbers start at 1 and rise sequentially
  ungroup
ngroups %>% print(n = Inf)
quantile(ngroups$n_final)

group_agibreaks <- ngroups %>%
  group_by(agige0, MARS) %>%
  mutate(agimax_prior=lag(agimax),
         agimin_next=lead(agimin),
         agilow_ge=case_when(agige0==0 & grp==1 ~ -Inf,
                             agige0==1 & grp==1 ~ 0,
                             TRUE ~  (agimin + agimax_prior) / 2),
         agihigh_lt=case_when(agige0==0 & grp==max(grp) ~ 0,
                              agige0==1 & grp==max(grp) ~ Inf,
                              TRUE ~ (agimax + agimin_next) / 2),
         range=agihigh_lt - agilow_ge) %>%
  select(-agimax_prior, -agimin_next) %>%
  ungroup %>%
  arrange(MARS, agilow_ge) %>%
  mutate(ugroup=row_number()) %>% # universal group (group within the entire universe)
  select(agige0, MARS, grp, ngrps, ugroup, everything())
  
# checks:
group_agibreaks %>% 
  print(n = Inf) %>% 
  kable(digits=0, format.args = list(big.mark=","))

# look at first and last rec in each agege0, MARS grouping
group_agibreaks %>% 
  filter(grp==1) %>% 
  kable(digits=0, format.args = list(big.mark=","))

group_agibreaks %>% 
  filter(grp==ngrps) %>% 
  kable(digits=0, format.args = list(big.mark=","))

group_agibreaks %>%
  group_by(MARS) %>%
  summarise(ngrps=n(), nsum=sum(n_final), nmin=min(n_final), nmax=max(n_final)) %>%
  mutate(totgrps=sum(ngrps), totn=sum(nsum))

# now loop through the data and put ugroup on each record
count(ffw, ftype)

getgroup <- function(agivec, MARSval){
  # MARSval <- 1
  # agivec <- c(-Inf, -1000, 0, 1000, 20e3, Inf)
  
  gindex <- function(agival, agilow_ge_vec) {
    ifelse(agival < Inf, 
           min(which(agival < agilow_ge_vec)) - 1,
           length(agilow_ge_vec) - 1)
  }
  
  breaks <- group_agibreaks %>%
    filter(MARS==MARSval) %>%
    arrange(agilow_ge) %>%
    select(MARS, ugroup, agilow_ge, agihigh_lt)
  
  agilow_ge_vec <- c(breaks %>% .[["agilow_ge"]], Inf)
  
  indexes <- gindex(1.5e6, agilow_ge_vec)
  indexes <- sapply(agivec, gindex, agilow_ge_vec)
  ugroup <- breaks$ugroup[indexes]
  return(ugroup)
}

a <- proc.time()
ffw2 <- ffw %>%
  group_by(MARS) %>%
  mutate(ugroup=getgroup(c00100, first(MARS)))
b <- proc.time()
b - a # 18 secs

#..2. Now that groups are defined, get constraints for each group ----
# get size-ordered list of vnames
(size_vars <- syn_info %>%
    mutate(vname=change_case(vname)) %>%
    arrange(-abs(sum)) %>% .[["vname"]])
showvars("sum")

#.. define variable groupings ----
pcatvars <- c("XTOT", "DSI", "EIC", "FDED", "MIDR", "n24", "f6251", "f2441") # don't include MARS as it is used for grouping
(continuous_vars <- setdiff(size_vars, c(pcatvars, "MARS")))
tcvars <- c("c00100", "taxbc", "c09600", "c05800")
# priority levels
p1 <- continuous_vars[1:10]
p2 <- continuous_vars[11:20]
p3 <- continuous_vars[21:30]
p4 <- continuous_vars[31:40]
p5 <- continuous_vars[41:50]
p6 <- continuous_vars[51:length(continuous_vars)]
cbasevars <- c(tcvars, p1, p2, p3, p4, p5, p6)
#.. end define variable groupings ----

# get constraint coefficients
ccoef <- ffw2 %>%
  ungroup %>%
  select(ftype, RECID, MARS, ugroup, wt0, cbasevars) %>%
  mutate_at(vars(cbasevars), list(npos = ~npos(., wt0),
                                  nz = ~nneg(., wt0),
                                  sumpos = ~sumpos(., wt0),
                                  sumneg = ~sumneg(., wt0)))
glimpse(ccoef) # 228 constraints plus ftype and ugroup

all_constraint_vars <- ccoef %>%
  select(contains("_n"), contains("_sum")) %>%
  names(.) # length 228

a <- proc.time()
all_constraint_vals <- ccoef %>%
  # filter(ugroup %in% 1:2) %>%
  group_by(ftype, ugroup) %>%  
  do(calc_constraints(.$wt0, ., all_constraint_vars) %>% 
       enframe %>%
       spread(name, value))
b <- proc.time()
b - a # about a min
glimpse(all_constraint_vals) # 228 constraints plus ftype ugroup; 246 obs (3 ftypes x 82 groups)


# fdiff <- function(x, puf) x - puf
# fpdiff <- function(x, puf) fdiff(x, puf) / puf * 100
# pdiffs <- all_constraint_vals %>%
#   gather(variable, value, -ftype, -ugroup) %>%
#   spread(ftype, value) %>%
#   mutate_at(vars(syn, syn_nd),
#             list(diff = ~fdiff(., puf),
#                  pdiff= ~fpdiff(., puf)))
# glimpse(pdiffs) # 82 groups
# length(unique(pdiffs$variable)) # length 228

# now find good constraints for each group, within each non-puf ftype:
# a) drop constraints where 

# b) drop duplicate constraints that have identical concoefs (keep the first)
getgoodcols <- function(data){
  # get names of non-duplicated constraint coefficients
  dupcols <- which(duplicated(as.matrix(data), MARGIN = 2))
  df <- tibble(good_constraint=setdiff(names(data), names(dupcols)))
  return(df)
}

a <- proc.time()
good_con <- ccoef %>%
  filter(ftype != "puf") %>%
  # filter(ugroup %in% 1:2) %>%
  group_by(ftype, ugroup) %>%
  do(getgoodcols(.[, all_constraint_vars]))
b <- proc.time()
b - a # only 13 secs
glimpse(good_con)
length(unique(good_con$good_constraint)) # only 132 good constraints
count(good_con, ftype, ugroup)

# get sums within groups of constraint coefficients to help with setting tolerances
# str_extract(all_constraint_vars, "[^_]+")

ccsums <- all_constraint_vals %>%
  gather(constraint_var, file_value, -ftype, -ugroup) %>%
  group_by(ugroup, constraint_var) %>%
  mutate(target=file_value[ftype=="puf"],
         diff=file_value - target,
         pdiff=diff / target * 100) %>%
  left_join(group_agibreaks %>% select(ugroup, agilow_ge, agihigh_lt)) %>%
  select(ugroup, constraint_var, agilow_ge, agihigh_lt, ftype, target, file_value, diff, pdiff) %>%
  mutate(vdesc=var_vnames$vdesc[match(str_extract(constraint_var, "[^_]+"), var_vnames$vname)]) %>%
  ungroup
glimpse(ccsums)
ccsums %>% filter(ugroup==3, ftype=="syn_nd")
# ,nnz_puf=sum(file_value[ftype=="puf"]!=0
# ccsums <- pdiffs %>%
#   group_by(ugroup, variable) %>%
#   summarise_at(vars(puf, syn, syn_nd, syn_diff, syn_nd_diff), ~sum(.))


# combine them
targets <- ccsums %>% 
  rename(good_constraint=constraint_var) %>%
  right_join(good_con, by=c("ftype", "ugroup", "good_constraint"))
  # left_join(good_pdiffs) %>%
  # filter(abs(pdiff) < 100) %>%
targets
targets %>% filter(ftype=="syn_nd", ugroup==3)
targets %>% filter(ftype=="syn_nd", target==0, file_value!=0) # 825 in the entire file
targets %>% filter(ftype=="syn_nd", target==!0, file_value==0) # none in the file


#.3. Prepare the constraint bounds ----
# look at priority groupings
showvars(usevars=str_to_upper(p1))
showvars(usevars=str_to_upper(p2))
showvars(usevars=str_to_upper(p3))
showvars(usevars=str_to_upper(p4))
showvars(usevars=str_to_upper(p5))
showvars(usevars=str_to_upper(p6))

# e09800 e58990 e03400 e07240 p08000 e07600 e24518
# create priority groupings and assign tolerances

tolerances <- tibble(vname=unique(str_extract(all_constraint_vars, "[^_]+"))) %>%
  left_join(syn_info %>% mutate(vname=change_case(vname))) %>%
  mutate(sum=case_when(vname=="c00100" ~ 1e20,
                       vname=="taxbc"  ~ 1e19,
                       TRUE ~ sum)) %>%
  arrange(-abs(sum)) %>%
  mutate(rank=row_number(),
         # establish default tolerances
         tol=case_when(rank %in% 1:10 ~ .001,
                       rank %in% 11:22 ~ .05,
                       rank %in% 23:40 ~ .25,
                       TRUE ~ Inf),
         # zero_frac is the max fraction of the value to keep if target is zero
         zero_frac=case_when(rank %in% 1:10 ~ .2,
                       rank %in% 11:22 ~ .4,
                       rank %in% 23:40 ~ .6,
                       TRUE ~ Inf))
tolerances
vars <- c("e58990", "e03400", "p08000", "e07240", "e07600", "e24518")
tolerances %>% filter(vname %in% vars)

# tolerances <- tibble(variable=all_constraint_vars) %>%
#   mutate(tol=case_when(variable %in% grp1 ~ .001,
#                        variable %in% grp2 ~ .05,
#                        variable %in% grp3 ~ .25,
#                        variable %in% grp4 ~ .5,
#                        variable %in% grp5 ~ .9,
#                        TRUE ~ 1.5),
#          # zero_frac is the max fraction of the value to keep if target is zero
#          zero_frac=case_when(variable %in% grp1 ~ .3,
#                              variable %in% grp2 ~ .4,
#                              variable %in% grp3 ~ .5,
#                              variable %in% grp4 ~ .6,
#                              variable %in% grp5 ~ .7,
#                              TRUE ~ .8),
#          zero_frac=1)
# tolerances

bounds <- targets %>%
  mutate(vname=str_extract(good_constraint, "[^_]+")) %>%
  left_join(tolerances %>% select(vname, rank, tol, zero_frac)) %>%
  mutate(clb=ifelse(target==0, -zero_frac * abs(diff), target - tol * abs(target)),
         cub=ifelse(target==0, +zero_frac * abs(diff), target + tol * abs(target))) %>%
  # ensure that logical inconsistencies cannot occur
  # .. counts cannot be negative, and sumpos cannot be negative
  mutate(clb=ifelse(str_detect(good_constraint, "_n") & (clb < 0), 0, clb),
         cub=ifelse(str_detect(good_constraint, "_n") & (cub < 0), 0, cub),
         clb=ifelse(str_detect(good_constraint, "_sumpos") & (clb < 0), 0, clb),
         cub=ifelse(str_detect(good_constraint, "_sumneg") & (cub > 0), 0, cub)) %>%
  select(-vdesc, everything(), vdesc)
bounds
bounds %>% filter(good_constraint=="taxbc_npos", ftype=="syn_nd", ugroup==7)
bounds %>% filter(ftype=="syn_nd", ugroup==7) %>% print(n = Inf)


#.4. Run the optimization ----
a <- proc.time()
opt <- ccoef %>%
  filter(ftype=="syn_nd") %>%
  # filter(ugroup %in% 1:5) %>%
  group_by(ftype, ugroup) %>%
  do(getwts_ipopt(., all_constraint_vals, bounds)) %>%
  ungroup
b <- proc.time()
b - a

saveRDS(opt, "d:/temp/opt.rds")

opt <- readRDS("d:/temp/opt.rds")
# glimpse(opt)
probs <- c(0, .01, .05, .1, .25, .5, .75, .9, .95, .99, 1)
quantile(opt$wt0, probs) %>% round(2)
quantile(opt$wt1, probs) %>% round(2)
quantile(opt$wt1 / opt$wt0, probs) %>% round(2)

tol <- .01
check <-  opt %>%
  # filter(ugroup==1) %>%
  select(ftype, MARS, ugroup, wt0, wt1, contains("_n"), contains("_sum")) %>%
  gather(variable, value, -ftype, -MARS, -ugroup, -wt0, -wt1) %>%
  group_by(ftype, MARS, ugroup, variable) %>%
  summarise(sum0=sum(value * wt0), sum1=sum(value * wt1)) %>%
  ungroup %>%
  left_join(bounds %>%
              filter(ftype=="syn_nd") %>% 
              select(ugroup, variable=good_constraint, target, clb, cub, vdesc) %>% 
              mutate(range=cub - target)) %>%
  mutate(cnum=row_number(),
         pdiff0=sum0 / target * 100 - 100,
         pdiff1=sum1 / target * 100 - 100,
         diff0=sum0 - target,
         diff1=sum1 - target,
         viol=!(sum1>=(clb-abs(clb*tol)) & sum1<=(cub+abs(cub*tol)))) %>%
  select(ftype, MARS, ugroup, variable, sum0, target, clb, sum1, cub, range, diff1, pdiff1, vdesc)
check

# tol <- .001
# check <- check %>%
#   mutate(viol=!(sum1>=(clb-abs(clb*tol)) & sum1<=(cub+abs(cub*tol)))) %>%
#   select(ftype, MARS, ugroup, variable, sum0, target, clb, sum1, cub, good, viol, diff1, pdiff1, vdesc)
# glimpse(check)
# glimpse(bounds)

check %>%
  filter(is.infinite(pdiff1) | (abs(pdiff1)>50)) %>%
  group_by(variable, vdesc) %>%
  summarise(n=n()) %>%
  arrange(-n)

check %>%
  filter(str_detect(variable, "_")) %>%
  # filter(ugroup==1) %>%
  filter(target!=0) %>%
  # pick one of the sorts
  arrange(-abs(diff1)) %>%
  # arrange(-abs(pdiff1)) %>%
  head(300) %>%
  mutate(vdesc=str_sub(vdesc, end=25)) %>%
  kable(digits=c(rep(0, 11), 1, 0), format.args = list(big.mark=","))

check %>%
  filter(str_detect(variable, "_")) %>%
  group_by(ftype, variable, vdesc) %>%
  summarise_at(vars(sum0, target, clb, sum1, cub), ~sum(.)) %>%
  mutate(diff0=sum0 - target,
         diff1=sum1 - target,
         pdiff0=sum0 / target * 100 - 100,
         pdiff1=sum1 / target * 100 - 100) %>%
  # pick one of the sorts
  arrange(-abs(diff1)) %>%
  # arrange(-abs(pdiff1)) %>%
  head(300) %>%
  mutate(vdesc=str_sub(vdesc, end=25)) %>%
  select(-diff0, -pdiff0) %>%
  kable(digits=c(rep(0, 9), 1), format.args = list(big.mark=","))


# find constraints where every puf is zero so target is zero but we have nonzero syn records -- 


check %>%
  filter(str_detect(variable, "_")) %>%
  filter(ugroup==7) %>%
  mutate(cnum=row_number(),
         pdiff0=sum0 / target * 100 - 100,
         pdiff1=sum1 / target * 100 - 100,
         diff0=sum0 - target,
         diff1=sum1 - target) %>%
  select(-vdesc, cnum, everything(), vdesc) %>%
  arrange(-abs(diff1)) %>%
  head(100) %>%
  mutate(vdesc=str_sub(vdesc, end=25)) %>%
  kable(digits=c(rep(0, 9), 1, 1, 0, 0, 0), format.args = list(big.mark=","))


#****************************************************************************************************
#               TEST:  parallelization ####
#****************************************************************************************************
# devtools::install_github("tidyverse/multidplyr")
# library("multidplyr")
# cluster <- new_cluster(4)
# 
# ccoef2 <- ccoef %>% 
#   filter(ftype != "puf") %>%
#   group_by(ftype, ugroup) %>%
#   partition(cluster)
# 
# good_con2 <- ccoef2 %>%
#   select(ftype, ugroup, all_constraint_vars) %>%
#   do(getgoodcols(.)) %>%
#   collect()


# library("doParallel")
# cl <- makeCluster(6)
# registerDoParallel(cl)
# 
# # define recipe if new one desired
# recipe <- get_recipe_long(get_weighting_recipe("recipe5")) %>%
#   filter(vname %in% names(puf)) %>%
#   dplyr::select(vname, vname, fn)
# 
# packages <- c("magrittr", "tidyverse", "dplyr", "nloptr")
# # CAUTION:  DO NOT PASS large items as function arguments
# # instead export them
# xport <- c("globals", "idfile", "recipe", "puf.vnames", "mrgdf",
#            "n.neg", "n.pos", "n.sum", "val.neg", "val.pos", "val.sum",
#            "naz", "eval_f_wtfs", "eval_grad_f_wtfs") 
# popts <- list(.packages=packages, .export=xport)
# popts
# djb ----
