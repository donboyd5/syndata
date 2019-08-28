
# http://blog.aicry.com/multidplyr-dplyr-meets-parallel-processing/index.html
# http://blog.aicry.com/r-parallel-computing-in-5-minutes/

# devtools::install_github("hadley/multidplyr")


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

library("multidplyr")


#****************************************************************************************************
#                Test ####
#****************************************************************************************************

library(nycflights13)
glimpse(flights)
count(flights, year)
count(flights, month)
summary(flights$arr_delay)

f2x <- function(x) 2 * x
f3x <- function(x) 3 * x
fxsq <- function(x) x^2


recipe <- tibble(fn=sample(c("f2x", "fxsq"), 12, replace=TRUE), month=1:12)
recipe

# mutate(syn.unwtd=do.call(fn, list(syn, vname, rep(1, nrow(syn))))) 
get_precs <- function(df){
  require(dplyr)
  require(magrittr)
  df <- df %>%
    left_join(recipe) %>%
    rowwise() %>%
    mutate(fval=do.call(fn, list(x=arr_delay)))
  # df <- df %>%
  #   left_join(recipe) %>%
  #   mutate(fval=do.call(fn, list(x=arr_delay)))
  write_csv(df, "d:/df.csv")
  return(df)
}

cl <- create_cluster(3)
set_default_cluster(cl)
cluster_copy(cl, get_precs)
cluster_ls(cl)

cluster_library(cl, c("magrittr", "tidyverse"))

cl %>%
  cluster_library("magrittr") %>%
  cluster_eval(search())

# cluster_assign_value(cl, 'get_precs', get_precs)

# cluster_get(cl, c("get_precs", "fxsq", "f2x"))

# cluster_copy(cl, c(get_precs, fxsq, f2x))
# cluster_copy(cl, c("get_precs", "fxsq", "f2x"))

cluster_copy(cl, get_precs)
cluster_copy(cl, fxsq)
cluster_copy(cl, f2x)
cluster_copy(cl, recipe)

df2 <- flights %>% 
  filter(!is.na(arr_delay)) %>%
  dplyr::select(month, arr_delay) %>%
  partition(month) %>%
  do(get_precs(.))

names(df2)
df3 <- df2 %>% collect()

ht(df3)
count(df3, fn)
recipe


s <- partition(mtcars)
s %>% mutate(cyl2 = 2 * cyl)
s %>% filter(vs == 1)
s %>% summarise(n())
s %>% select(-cyl)

if (require("nycflights13")) {
  planes <- partition(flights, tailnum)
  summarise(planes, n())
  
  month <- partition(flights, month)
  month %>% group_by(day) %>% summarise(n())
}





