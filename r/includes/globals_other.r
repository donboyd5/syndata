

if(!exists("globals")) globals <- list()

# these globals are not system-specific

# variables that I believe Tax-Calculator expects to be upper case
# note that FLPDYR and RECID are not synthesized
globals$upper_case_vars <- c("RECID", "MARS", "XTOT", "DSI", "EIC", "FLPDYR", "MIDR", "FDED") 

# agibrks <- c(-Inf, 0, 25e3, 50e3, 75e3, 100e3, 200e3, 500e3, Inf)



# synvars.cats <- c("dsi", "eic", "f2441", "f6251", "fded", "flpdyr", "mars", "midr", "n24", 
#                  "recid", "s006", "xtot") # generally (but not entirely) categorical-type variables