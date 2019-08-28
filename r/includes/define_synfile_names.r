
#.. define named vector of names of synthesized files ----
# for now, I include all of the synthesized files - we will only look at a subset
# for meanings of the files, see the Google sheet:
#    https://drive.google.com/open?id=1qTQJd2DGMm5zXnFxyP2Rw-8rszOkLNFrobg-NIikIsw

# files that Max created, generally with random forests
ttfnames <- c("puf_10p_sample_train.csv", "puf_10p_sample_test.csv")
names(ttfnames) <- c("train", "test")

synums <- 1:8
synfnames <- paste0("synpuf", synums, ".csv")
names(synfnames) <- paste0("synpuf", synums)

# files that Boyd created, generally with CART via R's synthpop
spnums <- 1:3
spnames <- paste0("synthpop", spnums, ".csv")
names(spnames) <- paste0("synthpop", spnums)

# put them together in one big vector
synfile.fnames <- paste0(synd, c(ttfnames, synfnames, spnames))
names(synfile.fnames) <- names(c(ttfnames, synfnames, spnames))
# fnames

# in some of the early synthesized files, we need to round MARS to an integer, so define those files for later use
bad.MARS.files <- paste0("synpuf", 1:4)

