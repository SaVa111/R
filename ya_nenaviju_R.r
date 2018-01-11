preventParalysis <- function(normXl){
  return (min(normXl) + (max(normXl) - min(normXl))  / 2)
}

prepareSample <- function(){
  xl <- iris[iris$Species == "setosa" | iris$Species == "versicolor",]
  x_sign1 <- xl[,3] - preventParalysis(xl[,3])
  x_sign2 <- xl[,4] - preventParalysis(xl[,4])
  class_sign <- xl[,5]
  y_sign <- c() 
  for(i in 1:nrow(xl)){
    if(xl[i,5] == "setosa"){
      y_sign[i] <- -1
    }else if(xl[i,5] == "versicolor"){
      y_sign[i] <- 1
    }
    
  }
  return(data.frame(x_sign1, x_sign2, y_sign,class_sign)) 
}

majorantFun <- function(currParam){
  return((1 - currParam)^2)
}

majorantToDiff <- function(currParam){
  return(2*currParam - 2)
}

stochGrad <- function(xl, eta, lambda){
  wghtVec <- c()
  wghtVecNew <- c()
  wghtVec <- runif(nrow(xl), -1/nrow(xl), 1/nrow(xl))
  Q <- 0
  eps <- 0.999
  

  for(i in 1:nrow(xl)){
    currSmplPoint <- c(xl[i,1],xl[i,2],rep(0,length(wghtVec) - 2))
    Q <- Q + majorantFun(crossprod(currSmplPoint,wghtVec) * xl[i,3])
    
  }
  print(Q)
  
  while(Q  >= eps){
    rndElXl <- sample(1:nrow(xl), 1, replace = TRUE)
    
    currRndXL <- c(xl[rndElXl, 1],xl[rndElXl, 2],rep(0,length(wghtVec) - 2))
    currError <- majorantFun(crossprod(currRndXL,wghtVec) * xl[rndElXl, 3])
    prodXiYi <- currRndXL  * xl[rndElXl,3]
    wghtVecNew <- wghtVec - as.vector(eta * (majorantToDiff(crossprod(currRndXL, wghtVec)))) * prodXiYi
    Q <- (1 - lambda)*Q + lambda*currError
    
    abline(wghtVecNew[3]/wghtVecNew[2],-wghtVecNew[1]/wghtVecNew[2], col = "black")
    
  }
  abline(wghtVecNew[3]/wghtVecNew[2],-wghtVecNew[1]/wghtVecNew[2], col = "red", lwd = 3)
  
}

drawGrad <- function(xl){
  colors <- c("setosa" =  "red", "versicolor"  =  "green")
  plot(xl[,1:2], pch = 21, xlab = 'x', ylab = 'y', col = colors[xl$class_sign], asp = 1, bg = colors[xl$class_sign])
  stochGrad(xl, 0.5, 0.5) 
}

xl_sample <- prepareSample()
drawGrad(xl_sample)