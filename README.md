# Метрические алгоритмы классификации
Метрические методы обучения -- методы, основанные на анализе сходства объектов.
Мерой близости называют функцию расстояния. Чем меньше расстояние между объектами, тем больше объекты похожи друг на друга.
Метрические алгоритмы классификации опираются на гипотезу компактности: схожим объектам соответствуют схожие ответы.
## KNN Метод К ближайших соседей

Метод прост, и основывается на гипотизе компактности. Оценивая К соседних элементов мы относим текущий элемент к тому классу - элементов которого больше среди его соседей.
Пример работы алгоритма для ирисов Фишера и 6 ближайших соседей.
![KNN](https://github.com/SaVa111/R/blob/master/Images/6NN.png)


### LOO
Метод скользящего контроля используется для подбора параметров, для необходимого алгоритма.
мы перебераем возможные варианты параметра, и для каждого элемента выборки проверяем ошибается ли алгоритм, если данный элемент не включать в выбоку. Проверяя возможные варианты параметра мы можем выбрать тот у которого наименьший коэффициент ошибки.
### LOO для KNN
Результат работы скользящего контроля для KNN.
![KNN LOO](https://github.com/SaVa111/R/blob/master/Images/LOOKNN.png)


Как видно из графика, алгоритм деёт наилучший результат при k=6.
### Взвешенный KNN
Основная проблема KNN в том, что наиближайший элемент и самый дальний от классифицируемого, абсолютно равнозначны. Хотя очевидно, что в большенстве случаев чем дальше объект находится тем меньше он похож на классифицируемый.

Эту проблему решает взвешенный KNN который помима самих соседей оценивает, так же и расстояние до объектов.
### Результат работы:
![KNN](https://github.com/SaVa111/R/blob/master/Images/W6NN.png)


Исходный код:

```
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
```

## Парзеновское окно

Возьмём весовую функцию - некоторую невозрастающую на [0, ∞] функцию ядра.
и получим классификатор который имеет вид: 
![LOO Parzen](https://github.com/SaVa111/R/blob/master/Images/image3.png)
основной смысл в том, что мы просматриваем сферическую окрестность(окно) и попадающин в нее элементы определяют к какому из классов относится элемент.
![CParzen](https://github.com/SaVa111/R/blob/master/Images/CParzen.png)

Исходный код:
```
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
```

### Сравним различные ядра.
#### Гауссово ядро.
![LOO Parzen](https://github.com/SaVa111/R/blob/master/Images/LOOgaus.png)
#### Ядро Епаненченкова.
![LOO Parzen](https://github.com/SaVa111/R/blob/master/Images/LOOepanenchenkov.png)
#### Прямоугольное ядро.
![LOO Parzen](https://github.com/SaVa111/R/blob/master/Images/LOOrectangle.png)
#### Треугольное ядро.
![LOO Parzen](https://github.com/SaVa111/R/blob/master/Images/LOOtriangle.png)

Наилучший результат даёт Гаусово ядро со значением ошибки 0.1, в то время как у остальных он равен 0.4.

## Парзеновское окно переменной ширины
Чтобы избежать случаев, когда у нас окно не попадает, ни на один элемент или наоборот охватывает слишком большую окрестность, введем зависимость размера окна от попадающих в него элементов. Будем выбирать такой размер окна, чтобы в него всегда попадало k соседей.
![VParzen](https://github.com/SaVa111/R/blob/master/Images/VParzen.png)

```
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
```
# Байесовские методы классификации

## Plug-in алгоритм
Подстановачный метод. Мы предпологаем, что каждый класс образован нормальным распределением. Тогда, можем оценить математическое ожидание и дисперсию для каждого класса и попытаться востановить первоначальное распределение подстановкой этих параметров в оптимальный байесавский классификатор. Таким образом мы получим кривую второго порядка разделяющую классы.

#### Пример работы:
![plug-in](https://github.com/SaVa111/R/blob/master/Images/plug-in.png)
![plug-in](https://github.com/SaVa111/R/blob/master/Images/plug-in2.png)
## Сеть радиальных базисных функций

Радиально базисными функцами называется класс функций являющихся Гауссианами с диагональной матрицей ковариации и имеет вид: 
![RBF](https://github.com/SaVa111/R/blob/master/Images/image.png)

Где N - это нормировочный множитель

Алгоритм классификации имеет вид суперпозиции, состоящей из трех уровней или слоёв. Первый слой представляет собой набор из радиально базисных функций для каждой компоненты, каждого класса. второй слой состоит из M сумматоров, вычисляющих принадлежность объекта x каждому из классов. И третий слой определяющий к какому из классов отнести элемент.

![RBF](https://github.com/SaVa111/R/blob/master/Images/RBFnet.png)

Компоненты, их матожидания, дисперсии и веса получаются с помощью ЕМ-алгоритма.
ЕМ-алгоритм основывается на вычислении неизвестных апостериорных вероятностей и исспользования их для востановления матожидания и дисперсии.

Алгоритм состоит из 2х шагов:
1. На Е шаге вычисляются приближенные апостериорные вероятности каждого элемента выборки.
2. На М шаге вычисляются матожидание и дисперсия для каждой из компонент.
 Ниже предаставлен код для ЕМ-алгоритма
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
## Пример работы:
![RBF](https://github.com/SaVa111/R/blob/master/Images/RBF2d.png)
## Трёхмерная модель для наглядного демонстрирования полученного распределения:
![RBF](https://github.com/SaVa111/R/blob/master/Images/RBF3d.png)
## Аналогично для трёх классов состаящих из различного числа компонент.
![RBF](https://github.com/SaVa111/R/blob/master/Images/RBF2d2.png)
## Результат работы.
![RBF](https://github.com/SaVa111/R/blob/master/Images/RBF3d2.png)

# Линейные методы классификации
## Метод стохастического градиента
Метод заключается в минимизации импирического риска, путем градиентного спуска.
Имеется некоторый вектор весов, который регулирует работу нашего классификатора. Будем изменять данный вектор двигаясь против градиента функции ошибки, тем самым её минимизируя.
![grad](https://github.com/SaVa111/R/blob/master/Images/grad.png)

Исходный код:
```
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
```

