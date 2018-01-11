library(MASS) # Generation of multidimensional normal distribution
library(rgl) 
norm_vec <- function(x) sqrt(sum(x^2))

estimate_mu <- function(objects) {
  rows <- dim(objects)[1]
  cols <- dim(objects)[2]
  mu <- matrix(NA, 1, cols)
  for (col in 1:cols) {
    mu[1, col] = mean(objects[ ,col])
  }
  return(mu)
}

# Covariation matrix of normal distribution
estimate_cov_matrix <- function(objects, mu) {
  rows <- dim(objects)[1]
  sigma <- rep(0, ncol(objects))
  for (i in 1:rows) {
    sigma <- sigma + (objects[i, ] - mu)^2
  }
  sigma <- sigma / rows
  
  return(sigma)
}

NG <- function(x, mu, sigma)
{
  n <- length(mu)
  ro <- 0
  for(d in 1:n)
  {
    ro <- ro + ((x[d] - mu[d])^2)/sigma[d]
  }
  
  N <- ((2*pi)^(-n/2))*(prod((sigma^(-0.5))))
  return(N*exp((-0.5)*ro))
}
p <- function(x, w, mus, sigmas)
{
  n <- length(w)
  res <- 0
  for(i in 1:n)
  {
    res <- res + w[i]*NG(x, mus[[i]], sigmas[[i]])
  }
  return (res)
}

EM <- function(X, R, m0, delta)
{
  m <- nrow(X)
  k <- 1#количество компонент
  w <- c(1)#веса компонент
  mus <- list()
  mus[[1]] <- estimate_mu(X)
  sigmas <- list()
  sigmas[[1]] <- estimate_cov_matrix(X, mus[[1]])
  
  while(TRUE)
  {
    if(k == 1){
      a<-c()
    }
    ps <- c()
    for(i in 1:m)
    {
      ps <- c(ps, p(X[i, ], w, mus, sigmas))
    }
    maxp <- max(ps) / R
    
    toDel <- c()
    for(i in 1:m)
    {
      if(ps[i] < maxp) toDel <- c(toDel, i)
    }
    U <- X
    for(i in 1:length(toDel))
    {
      U <- U[-toDel[length(toDel) - i + 1], ]
    }
    
    if(length(U) < 3) break
    if(nrow(U) < m0) break
    k <- k + 1
    
    w <- c(w, nrow(U)/m)
    for(j in 1:k-1)
    {
      w[j] <- w[j]*(1 - w[k])
    }
    mus[[k]] <- estimate_mu(U)
    sigmas[[k]] <- estimate_cov_matrix(U, mus[[k]]) / w[k]
    
    g <- matrix(rep(0, m*k), nrow = m, ncol = k, byrow = TRUE)
    while(TRUE)
    {
      #E step
      g0 <- g
      for(i in 1:m)
      {
        tmpsum1 <- 0
        for(j in 1:k)
        {
          tmpsum1 <- tmpsum1 + w[j]*NG(X[i, ], mus[[j]], sigmas[[j]])
        }
        for(j in 1:k)
        {
          g[i, j] = w[j]*NG(X[i, ], mus[[j]], sigmas[[j]]) / tmpsum1
        }
      }
      #M step
      for(j in 1:k)
      {
        w[j] <- sum(g[1:m, j]) / m
        tmpsum2 <- rep(0, 2)
        for(i in 1:m)
        {
          tmpsum2 <- X[i, ] * g[i, j] + tmpsum2
        }
        mus[[j]] <- tmpsum2 / (m*w[j])
        tmpsum3 <- 0
        for(i in 1:m)
        {
          tmpsum3 <- ((X[i, ] - mus[[j]])^2) * g[i, j] + tmpsum3
        }
        sigmas[[j]] <- tmpsum3 / (m*w[j])
      }
      
      tmp <- c()
      for(i in 1:m)
      {
        tmp <- c(tmp, sqrt(sum((g[i] - g0[i])^2)))
      }
      
      crntdelta <- max(tmp)
      if(crntdelta < delta) break
    }
  }
  
  return(list(k, mus, w, sigmas))
}

# Count of objects in each class
objects_count <- 100 

# Generation of test data
Sigma11 <- matrix(c(0.5, 0, 0, 0.5), 2, 2)
Sigma12 <- matrix(c(0.5, 0, 0, 0.5), 2, 2)
Sigma13 <- matrix(c(1, 0, 0, 1), 2, 2)
Sigma21 <- matrix(c(1, 0, 0, 1), 2, 2)
Sigma22 <- matrix(c(1, 0, 0, 1), 2, 2)
Sigma31 <- matrix(c(0.7, 0, 0, 1), 2, 2)
Mu11 <- c(-1, 0)
Mu12 <- c(1, 7)
Mu13 <- c(5, 2)
Mu21 <- c(15, 0)
Mu22 <- c(17, 6)
Mu31 <- c(10, 12)
xy1 <- mvrnorm(n = objects_count, Mu11, Sigma11)
xy2 <- mvrnorm(n = objects_count, Mu21, Sigma21)
xy3 <- mvrnorm(n = objects_count, Mu31, Sigma31)
xy1 <- rbind(xy1, mvrnorm(n = objects_count, Mu12, Sigma12))
xy1 <- rbind(xy1, mvrnorm(n = objects_count, Mu13, Sigma13))
xy2 <- rbind(xy2, mvrnorm(n = objects_count, Mu22, Sigma22))

m1 <- nrow(xy1)
m2 <- nrow(xy2)
m3 <- nrow(xy3)
# Assembling 2 classes in one sample xl
xl <- rbind(cbind(xy1, 1), cbind(xy2, 2), cbind(xy3, 3))

# Drawing the training sample
colors <- c("blue2", "green3", "red")
plot(xl[ , 1], xl[ , 2], pch = 21, bg = colors[xl[ ,3]], asp = 1, xlab = "x", ylab = "y")

m0 <- 50 #минимальная длина выборки, по которой можно восстановить плотность
R <- 2#максимальный допустимый разброс правдоподобия объектов
delta <- 0.3#критерий остонова

res1 <- EM(xy1, R, m0, delta - 0.2)

k1 <- res1[[1]]
mus1 <- res1[[2]]
w1 <- res1[[3]]
sigmas1 <- res1[[4]]

res2 <- EM(xy2, R, m0, delta)

k2 <- res2[[1]]
mus2 <- res2[[2]]
w2 <- res2[[3]]
sigmas2 <- res2[[4]]

res3 <- EM(xy3, R, m0, delta)

k3 <- res3[[1]]
mus3 <- res3[[2]]
w3 <- res3[[3]]
sigmas3 <- res3[[4]]

netCol <- 300
netRow <- 180
x <- -5 + (1:netCol) * 0.1
y <- -3 + (1:netRow) * 0.1
z <- c()
colors <- matrix(NA, netCol, netRow)
for(i in 1:netRow)
{
  for(j in 1:netCol)
  {
    value1 <- p(c(-5 + j*0.1, -3 + i*0.1), w1, mus1, sigmas1) * 10
    value2 <- p(c(-5 + j*0.1, -3 + i*0.1), w2, mus2, sigmas2) * 10
    value3 <- p(c(-5 + j*0.1, -3 + i*0.1), w3, mus3, sigmas3) * 10
    z <- c(z, max(c(value1, value2, value3)))
    if(value1 > value2 && value1 > value3)
    {
      colors[j, i] = "blue"
    }
    else if(value2 > value1 && value2 > value3)
    {
      colors[j, i] = "green"
    }
    else
    {
      colors[j, i] = "red"
    }
  }
}

z <- matrix(z, netRow, netCol)

open3d()
surface3d(x, y, z, col=colors, back = "lines")

#plot3d(y, x, z, col = rainbow(1000), back = "lines")
