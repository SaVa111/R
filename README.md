# Сеть радиальных базисных функций

Радиально базисными функцами называется класс функций являющихся Гауссианами с диагональной матрицей ковариации и имеет вид: 
![RBF](https://github.com/SaVa111/R/blob/master/Images/image.png)
![RBF](https://github.com/SaVa111/R/blob/master/Images/RBFnet.png)
```
EM <- function(X, R, m0, delta)
{
  m <- nrow(X)
  k <- 1 #count of components
  w <- c(1) #weights of components
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
```
![RBF](https://github.com/SaVa111/R/blob/master/Images/RBF2d.png)
![RBF](https://github.com/SaVa111/R/blob/master/Images/RBF3d.png)
![RBF](https://github.com/SaVa111/R/blob/master/Images/RBF2d2.png)
![RBF](https://github.com/SaVa111/R/blob/master/Images/RBF3d2.png)
