
#****************************************************************************************************
#                ####
#****************************************************************************************************
synpart <- function(puf.subset.rules, visit.sequence, methods, smoothing.rules=NULL, save.name=NULL, test=TRUE){
  # function to synthesize a part of the file, based on a logical variable passed to it 
  
  # puf extract
  pufx <- puf.base %>%
    filter(eval(puf.subset.rules)) %>%
    mutate(divratio=e00650 / e00600,
           divratio=ifelse(is.na(divratio), 0, divratio),
           penratio=e01700 / e01500,
           penratio=ifelse(is.na(penratio), 0, penratio))
  
  # construct variables of interest
  pufx.base <- pufx %>%
    dplyr::select(visit.sequence) # note that RECID is dropped; if we want it in dataset, we should remove it from predictors
  
  seed <- 1234
  if(test) {
    a <- proc.time()
    synx <- syn(pufx.base, visit.sequence=visit.sequence, method=methods, smoothing=smoothing.rules, seed=seed, m=0)
    b <- proc.time()
    print(b - a)      
  } else {
    a <- proc.time()
    synx <- syn(pufx.base, visit.sequence=visit.sequence, method=methods, smoothing=smoothing.rules, seed=seed, m=1)
    b <- proc.time()
    print(b - a)
  }
  return(synx)
}
