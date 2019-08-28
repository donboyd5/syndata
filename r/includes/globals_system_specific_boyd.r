
# change these variables when moving to a new system 

if(!exists("globals")) globals <- list()

globals$pufdir <- "D:/Dropbox/OSPC - Shared/IRS_pubuse_2011/" # location of puf2011.csv
globals$synd <- "D:/Google Drive/synpuf/"

# globals$tc.cli <- "C:/ProgramData/Anaconda3/Scripts/tc" # location of Tax-Calculator command-line interface
globals$tc.cli <- "C:/Users/donbo/Anaconda3/Scripts/tc" # location of Tax-Calculator command-line interface

# private directory for Tax-Calculator record-level output that we don't want moved from this machine
globals$tc.dir <- "D:/tcdir/"

globals$taxplans.dir <- "D:/Dropbox/RPrograms PC/OSPC/syndata4/tax_plans/"
