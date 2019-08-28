
# make a data frame with wtvar (variable to be targeted) and subgroup, of the form
# wtvar  subgroup                        
# s006   c00100 <= 0                     
# c00100 c00100 <= 0                     
# e00200 c00100 <= 0                     
# s006   c00100 > 0 & c00100 <= 25e3     
# c00100 c00100 > 0 & c00100 <= 25e3   

iranges <- c(
  "c00100 <= 0", 
  "c00100 > 0 & c00100 <= 25e3",
  "c00100 > 25e3 & c00100 <= 50e3",
  "c00100 > 50e3 & c00100 <= 75e3",
  "c00100 > 75e3 & c00100 <= 100e3",
  "c00100 > 100e3 & c00100 <= 200e3",
  "c00100 > 200e3 & c00100 <= 500e3",
  "c00100 > 500e3 & c00100 <= 1e12")
iranges
paste0("(", iranges, ")")

# vars to weight
vtw <- c("s006", "c00100", "e00200")
vtw <- c("s006", "c00100", "e00200", "e00200p")
vtw <- c("s006", "c00100", "e00200", "e00200p", "e00600", "e00650")
vtw <- c("s006", "c00100", "e00200", "e00200p", "e00600", "e00650", "p23250")
vtw <- c("s006", "c00100", "e00200", "e00200p", "e00600", "e00650", "p23250", "e26270")

# create target rules
target.rules.s <- expand.grid(wtvar=vtw, subgroup=iranges, stringsAsFactors = FALSE)
# target.rules



