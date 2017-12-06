euclideanDistance <- function(u, v)
{
  sqrt(sum((u - v)^2))
}
LOO <- function(matrica, k)
{
    error <- 0
    for(i in 1:nrow(matrica))
    {
        if(KNN(matrica[i,3:4], matrica[-i,], k) != matrica[i, 5])
        {
            error <- error + 1
        }
    }
    return (error / nrow(matrica))
}

epanenchenkov_core <- function(r) {
  if (abs(r) > 1)
    return (0.)
  return (0.75 * (1 - r^2)) 
}
gaus_core <- function(r) {
  if (abs(r) > 1)
    return (0.)
  return (2 * pi)^(-0.5) * exp^(-0.5 * r^2)
}
rectangle_core <- function(r) {
  if (abs(r) > 1)
    return (0.)
  return (0.5)
}
triangle_core <- function(r) {
  if (abs(r) > 1)
    return (0.)
  return (1 - abs(r))
}
quartic_core <- function(r) {
  if (abs(r) > 1)
    return (0.)
  return (15/16 * (1 - r^2)^2)
}

ConstParzen <- function(point, matrica, h)
{
  values <- rep(0, 3)
  for(i in 1:nrow(matrica))
  {
    distance <- euclideanDistance(c(matrica[i,3], matrica[i,4]), point)
    
    value <- rectangle_core(distance / h)
    values[matrica[i,5]] <- values[matrica[i,5]] + value
  }
  if (max(values) == 0)
    return(0)
  return(UnName(which.max(values)))
}

VairableParzen <- function(point, matrica, k)
{
  
  distances <- matrix(data = NA ,nrow = 1)
  for(i in 1:nrow(matrica))
  {
    distances[i] <- euclideanDistance(c(matrica[i,3], matrica[i, 4]), point)
  }
  namedDistances <- data.frame(matrica[,5], distances)
  sortedNamedDistances <- namedDistances[order(distances),]
  if(sortedNamedDistances[1,2] == 0.)
    return(sortedNamedDistances[1, 1])
  #print(sortedNamedDistances)
  values <- rep(0, 3)
  for(i in 1:k)
  {
    value <- rectangle_core(sortedNamedDistances[i,2] / sortedNamedDistances[k+1,2])
    values[sortedNamedDistances[i,1]] <- values[sortedNamedDistances[i,1]] + value
  }
  if (max(values) == 0)
    return(0)
  return(UnName(which.max(values)))
}

UnName <- function(id)
{
  if (id == 1)
    return ("setosa")
  if (id == 2)
    return ("versicolor")
  if (id == 3)
    return ("virginica")
  return("unknown")
}


KNN <- function(point, matrica, k, q = 1)
{
  distances <- matrix(data = NA ,nrow = 1)
  for(i in 1:nrow(matrica))
  {
    distances[i] <- euclideanDistance(c(matrica[i,3], matrica[i,4]), point)
  }
  namedDistances <- data.frame(matrica[,5], distances)
  sortedNamedDistances <- namedDistances[order(distances),]
  
  aaa <- rep(0, 3)
  aaa[sortedNamedDistances[1:k, 1]] <- 0
  for(i in 1:k)
  {
    aaa[sortedNamedDistances[i, 1]] <- aaa[sortedNamedDistances[i, 1]] + q^i
  }
  
  max <- 0
  id <- 0
  
  for (i in 1:length(aaa)) {
    if (aaa[i] > max) {
      max = aaa[i]
      id = i
    }
  }
  
  return(UnName(id))
}

CParzen_plot <- function(h)
{
  colors <- c("setosa" = "red", "versicolor" = "green3", "virginica" = "blue")
  plot(iris[, 3:4], pch = 21, bg = colors[iris$Species], col = colors[iris$Species])
  for(i in 0:70)
  {
    for(j in 0:25)
    {
      point <- c(1 + 0.1 * i, 0.1 * j)
      points(point[1], point[2], pch = 1, col = colors[ConstParzen(point, iris, h)], asp = 1)
    }
  }
}

VParzen_plot <- function(k)
{
  colors <- c("setosa" = "red", "versicolor" = "green3", "virginica" = "blue")
  plot(iris[, 3:4], pch = 21, bg = colors[iris$Species], col = colors[iris$Species])
  for(i in 0:70)
  {
    for(j in 0:25)
    {
      point <- c(1 + 0.1 * i, 0.1 * j)
      points(point[1], point[2], pch = 1, col = colors[VairableParzen(point, iris, k)], asp = 1)
    }
  }
}

KNN_plot <- function(k, q)
{
  colors <- c("setosa" = "red", "versicolor" = "green3", "virginica" = "blue")
  plot(iris[, 3:4], pch = 21, bg = colors[iris$Species], col = colors[iris$Species])
  for(i in 0:70)
  {
    for(j in 0:25)
    {
      point <- c(1 + 0.1 * i, 0.1 * j)
      points(point[1], point[2], pch = 1, col = colors[KNN(point, iris, k, q)], asp = 1)
    }
  }
}

LOO4K_plot <- function()
{
  plot(0, 0, type = "n", pch = 21, asp = 0, main = "LOO",xlim = c(0, nrow(iris)), ylim = c(0, 1), xlab = "k", ylab = "value")
  x <- c()
  y <- c()
  for(k in 1:nrow(iris))
  {
    x <- c(x, k)
    y <- c(y, LOO(iris, k))
    lines(x[k-1:k], y[k-1:k], pch = 8, bg = "black", col = "green3")
  }
}