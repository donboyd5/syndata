
library(bestNormalize)
x <- rgamma(100, 1, 1)
plot(x)

d <- tibble(z=rnorm(100))
# tibble(x=rnorm(100)) %>%
d %>%
  ggplot(aes(z)) +
  stat_function(fun=dnorm,
                color="red",
                args=list(mean=mean(d$z), 
                          sd=sd(d$z))) +
  geom_histogram(aes(y=..density..), fill="blue")

# tibble(z=rnorm(100)) %>%
tibble(z=rgamma(100, 1, 1)) %>%
  ggplot(aes(z)) +
  geom_histogram(aes(y=..density..), fill="cyan") +
  stat_function(fun=dnorm,
                color="red", size=1.5)

# With Repeated CV
x <- rgamma(10e3, 1, 1)
x <- hipuf$E00200
x <- hipuf %>% filter(E00200!=0) %>% .[["E00200"]]

x
tibble(z=x) %>%
  ggplot(aes(z)) +
  geom_histogram(aes(y=..density..), fill="cyan") +
  stat_function(fun=dnorm,
                color="red", size=1.5)

BN_obj <- bestNormalize(x)
BN_obj
p <- predict(BN_obj) # get the normal version
x2 <- predict(BN_obj, newdata = p, inverse = TRUE) # get x values from quasi-normal p values
all.equal(x2, x)

tibble(z=p) %>%
  ggplot(aes(z)) +
  geom_histogram(aes(y=..density..), fill="cyan") +
  stat_function(fun=dnorm,
                color="red", size=1.5)

# prediction -- creating wages
p3 <- rnorm(length(p), mean=mean(p), sd=sd(p))
x3 <- predict(BN_obj, newdata = p3, inverse = TRUE)
x3

tibble(z=x3) %>%
  ggplot(aes(z)) +
  geom_histogram(aes(y=..density..), fill="cyan") +
  stat_function(fun=dnorm,
                color="red", size=1.5)

tibble(rn=1:length(x), x.true=x, x.syn=x3) %>%
  gather(variable, value, -rn) %>%
  ggplot(aes(value, colour=variable)) +
  geom_density()

# repeat without dropping zeros
wages <- hipuf$E00200
wages_bn <- bestNormalize(wages)
wages_bn
wages_norm <- predict(wages_bn) # get the normal version

# not very normal looking when we keep zeros
tibble(z=wages_norm) %>%
  ggplot(aes(z)) +
  geom_histogram(aes(y=..density..), fill="cyan") +
  stat_function(fun=dnorm,
                color="red", size=1.5)

# prediction -- creating wages
wages_norm_pred <- rnorm(length(wages_norm), mean=mean(wages_norm), sd=sd(wages_norm))
wages_pred <- predict(wages_bn, newdata = wages_norm_pred, inverse = TRUE)
wages_pred

# looks very different from true data, predicts negative wages
tibble(z=wages_pred) %>%
  ggplot(aes(z)) +
  geom_histogram(aes(y=..density..), fill="cyan") +
  stat_function(fun=dnorm,
                color="red", size=1.5)

tibble(rn=1:length(wages), x.true=wages, x.syn=wages_pred) %>%
  gather(variable, value, -rn) %>%
  ggplot(aes(value, colour=variable)) +
  geom_density()


library(MultiVarMI) 
library(mvnfast)



## Not run:
# With leave-one-out CV
BN_obj <- bestNormalize(x, loo = TRUE)
BN_obj
p <- predict(BN_obj)
x2 <- predict(BN_obj, newdata = p, inverse = TRUE)
all.equal(x2, x)
## End(Not run)
# Without CV
BN_obj <- bestNormalize(x, allow_orderNorm = FALSE, out_of_sample = FALSE)
BN_obj
p <- predict(BN_obj)
x2 <- predict(BN_obj, newdata = p, inverse = TRUE)
all.equal(x2, x)