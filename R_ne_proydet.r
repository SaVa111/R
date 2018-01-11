library(MASS)

estimateMu <- function(objects)
{
  rows <- dim(objects)[1]
  cols <- dim(objects)[2]
  
  mu <- matrix(NA, 1, cols)
  for (col in 1:cols)
  {
    mu[1, col] = mean(objects[,col])
  }
  
  return(mu)
}

estimateCovarianceMatrix <- function(objects, mu)
{
  rows <- dim(objects)[1]
  cols <- dim(objects)[2]
  sigma <- matrix(0, cols, cols)
  
  for (i in 1:rows)
  {
    sigma <- sigma + (t(objects[i,] - mu) %*% (objects[i,] - mu)) / (rows - 1)
  }
  
  return (sigma)
}
getPlugInDiskriminantCoeffs <- function(mu1, sigma1, mu2,sigma2)
{
  invSigma1 <- solve(sigma1)
  26
  invSigma2 <- solve(sigma2)
  
  f <- log(abs(det(sigma1))) - log(abs(det(sigma2))) +
    mu1 %*% invSigma1 %*% t(mu1) - mu2 %*% invSigma2 %*%
    t(mu2);
  
  alpha <- invSigma1 - invSigma2
  a <- alpha[1, 1]
  b <- 2 * alpha[1, 2]
  c <- alpha[2, 2]
  
  beta <- invSigma1 %*% t(mu1) - invSigma2 %*% t(mu2)
  d <- -2 * beta[1, 1]
  e <- -2 * beta[2, 1]
  return (c("x^2" = a, "xy" = b, "y^2" = c, "x" = d, "y" = e, "1" = f))
}

ObjectsCount <- 600

Sigma1 <- matrix(c(0.1, 0, 0, 0.1), 2, 2)
Sigma2 <- matrix(c(1, 0, 0, 9), 2, 2)
Mu1 <- c(1, 0)
Mu2 <- c(15, 0)
xy1 <- mvrnorm(n=100, Mu1, Sigma1)
xy2 <- mvrnorm(n=ObjectsCount, Mu2, Sigma2)

xl <- rbind(cbind(xy1, 1), cbind(xy2, 2))

colors <- c("blue2", "green3")
plot(xl[,1], xl[,2], pch = 21, bg = colors[xl[,3]], asp = 1, xlab = "x", ylab = "y")

objectsOfFirstClass <- xl[xl[,3] == 1, 1:2]
objectsOfSecondClass <- xl[xl[,3] == 2, 1:2]

mu1 <- estimateMu(objectsOfFirstClass)
mu2 <- estimateMu(objectsOfSecondClass)
sigma1 <- estimateCovarianceMatrix(objectsOfFirstClass, mu1)
sigma2 <- estimateCovarianceMatrix(objectsOfSecondClass, mu2)
coeffs <- getPlugInDiskriminantCoeffs(mu1, sigma1, mu2, sigma2)


x <- y <- seq(-10, 20, len=100)
z <- outer(x, y, function(x, y) coeffs["x^2"]*x^2 +
             coeffs["xy"]*x*y
           + coeffs["y^2"]*y^2 + coeffs["x"]*x
           + coeffs["y"]*y + coeffs["1"])
contour(x, y, z, levels=0, drawlabels=FALSE, lwd = 3, col = "red", add = TRUE)