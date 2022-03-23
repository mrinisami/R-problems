#Problems taken from Scientific Programming and Simulation using R
#Ch 11, prob 2 and 3

f1 <- function(x) {2}
f2 <- function(x) {4 * (pi ^ 2) - (x - (2 * pi)) ^ 2}
f1squared <- function(x) {2 ^ 2}
f2squared <- function(x) {(4 * pi - x) ^ 2 * (x ^ 2)}

pol_int <- function(n, fn) {
  tot = 0
  for(i in 0:(n - 1)){
    el1 <- (2 * pi * i) / n
    item1 <- fn(el1)
    el2 <- (2 * pi * (i + 1)) / n
    item2 <- fn(el2)
    tot <- tot + ((pi * item1 * item2) / n)
  }
  return(tot)
}

#1st function

pol_int(10000, f1)
integrate(Vectorize(f1squared), lower = 0, upper = 2 * pi)[[1]] * 0.5

#No convergence since we have a scalar function

#2nd function

vec <- v()
for(i in seq(10, 10000, 500)){
  val <- pol_int(i, f2)
  vec <- append(vec, val)
}
pol_int(100000, f2)
integrate(f2squared, lower = 0, upper = 2 * pi)[[1]] * 0.5

show(vec)

#Convergence towards 2611.368 when N tends to infinity 

#7.3

rm(list = ls())

library(spuRs)

phi <- function(x) return(exp(-x^2/2)/sqrt(2*pi))
Phi <- function(z) {
  if (z < 0) {
    return(0.5 - integrate(phi, z, 0)[[1]])
  } else {
    return(0.5 + integrate(phi, 0, z)[[1]])
  }
}

find_z <- function(p, max.iter=1000){
  zup <- 5
  zdown <- 0
  p_up <- Phi(zup)
  while(p_up != p){
    dif <- zup - zdown
    p_up <- Phi(zup)
    if(p_up > p){zup <- zup - (dif / 2)}
    else if(p_up < p){
      zdown <- zup
      zup <-  zdown + (dif / 2)
    }
  }
  return(zup)
}

test <- c(0.5, 0.95, 0.975, 0.99)
ans <- sapply(test, find_z)
quant <- sapply(test, qnorm)

#Bissection method to retrieve the quantile function of a distribution

quant
ans
