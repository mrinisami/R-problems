#Problems taken from Scientific Programming and Simulation using R
#Ch 12, prob 6 and 7

library(numDeriv)

fx <- function(x, y){
  -(x ^ 2 + y ^ 2 - 2) * (x ^ 2 + y ^ 2 -1) * (x ^ 2 + y ^ 2) * (x ^ 2 + y ^ 2 + 1) * (x ^ 2 + y ^ 2 + 2) * (2 - sin(x ^ 2 - y ^ 2)
                                                                                                             * (cos(y - exp(y))))
}

grads <- function(x, y){
  fun <- expression((-(x^2 + y^2 - 2)*(x^2 + y^2 - 1)*(x^2 + y^2)*(x^2 + y^2 + 1)*(x^2 + y^2 + 2)*(2 - sin(x^2 - y^2)*cos(y - exp(y)))))
  dx <- D(fun, 'x')
  dy <- D(fun, 'y')
  xi <- eval(dx, list(x, y))
  yi <- eval(dy, list(x, y))
  return(c(xi, yi))
}

steep <- function(iter){
  draw <- runif(2, -1.5, 1.5)
  x <- draw[1]
  y <- draw[2]
  gradient <- grads(x, y)
  for (i in 1:iter){
    find_d <- function(d){
      x_d <- x + d * gradient[1]
      y_d <- y + d * gradient[2]
      return(fx(x_d, y_d))
    }
    opt_d <- optim(0, find_d, method="BFGS", control=list(fnscale=-1))
    x <- x + opt_d[[1]] * gradient[1]
    y <- y + opt_d[[1]] * gradient[2]
    gradient <- grads(x, y)
    if (is.na(gradient)) {
      draw <- runif(2, -1.5, 1.5)
      x <- draw[1]
      y <- draw[2]
    }
  }
  return(c(x, y))
}

ans <- vector()
for (i in 1:50){
  data <- steep(200)
  x <- data[1]
  y <- data[2]
  max_found <- fx(x, y)
  ans <- append(ans, max_found)
}

ans

#Global maxima is 7.6135. You can show ans to figure out the local maxima.  


#Prob 7

library(spuRs)
direc <- .libPaths()[1]
data <- read.csv(file = paste0(direc,"/spuRs/resources/data/trees.csv"))

split.screen(c(2, 1))
trees <- data.frame(data)
sites <- unique(data[,1])               
loc <- data[2]

richards <- function(t, theta){
  theta[1] * (1 - exp(-theta[2] * t)) ^ theta[3]
}

loss.L2 <- function(theta, age, vol){
  sum((vol - richards(age, theta)) ^ 2)
}

theta.L2 <- optim(theta0, loss.L2, age=data$Age, vol=data$Vol)

screen(1)
plot(data$Age, data$Vol, type = "p", xlab = "Age", ylab = "Volume", text(trees$Age, data$Vol, substr(data$ID, 1, 1)))

sites[1]

theta0 <- c(1000, 0.1, 3)

xs <- vector()
ys <- vector()
site <- vector()

#plot1

for (el in sites){
  tree <- trees[data$ID == el, 2:3]
  theta.L2 <- optim(theta0, loss.L2, age=tree$Age, vol=tree$Vol)
  param <- theta.L2[[1]]
  a <- param[1]
  b <- param[2]
  xs <- append(xs, a)
  ys <- append(ys, b)
  site <- append(site, substr(el, 1, 1))
  colours <- c('red', 'blue', 'yellow', 'pink', 'grey')
  lines(tree$Age, richards(tree$Age, theta.L2$par), col=colours[strtoi(substr(el, 1, 1))])
}
legend('topleft', legend = c('Site 1', 'Site 2', 'Site 3', 'Site 4', 'Site 5'), text.col = colours)

#plot2

min_x <- min(xs)
max_x <- max(xs)
min_y <- min(ys)
max_y <- max(ys)

screen(2)
plot(c(min_x, max_x), c(min_y, max_y), xlab='Max size', ylab='Growth rate')

for (i in 1:length(xs)) {
  #points(xs[i], ys[i])
  text(xs[i], ys[i], site[i])
}

#Je n'ai pas mis les points pour le 2e graph car c'est encore plus difficile de lire les sites, puisque les chiffres sont centres sur les points,
#je me suis dis que c'etait mieux comme ca
#On dirait qu'il existe une relation entre le site et les parametres. Les arbres provenant des sites 1-2 semblent avoir un plus petit growth rate
#et une plus grande size.