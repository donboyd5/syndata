

glimpse(pufcart)

xvars <- c("pufseqn", "wt", "e01500_minus_e01700", "divratio", "m", "msname")

length(setdiff(names(pufcart), xvars))
# only rownum



# mark dup puf records
vset <- setdiff(names(pufcart), c(xvars, "rownum"))
tmp <- puf %>% select(vset)
puf_dups <- duplicated(tmp) | duplicated(tmp, fromLast=TRUE)

sum(puf_dups)  
# 19527 dups
ipuf_unique <- !puf_dups
sum(ipuf_unique)

upuf <- puf[ipuf_unique, ]
names(upuf)

# first, a simple test - stack the cart and upuf and look for dups
ctmp <- cart %>% select(vset)
icart_unique <- !(duplicated(ctmp) | duplicated(ctmp, fromLast=TRUE))
# icart_unique <- icart_dups
sum(!icart_unique)
ucart <- cart[icart_unique, ]
stack <- bind_rows(upuf, ucart)

stack2 <- stack %>% select(vset)

dups <- duplicated(stack2) | duplicated(stack2, fromLast=TRUE)
sum(dups)
names(stack2) %>% sort



