#############################################
# Programa noharm4 parcialmente portado a R #
#############################################

rm(list = ls(all = TRUE))

miNoharm = function(N,M,NS,IN,EX,IV,PR,C,F=NULL,P=NULL,data,maxiteraciones = 100,epsilon = 10^-8){
  
  ## sección de control ## 
  if(missing(N)){
      stop("'Falta el número de items 'N'")
  }
  if(missing(M)){
    stop("'Falta el número de dimensiones 'M'")
  }
  if(missing(NS)){
    stop("'Falta el número de observaciones 'NS'")
  }
  if(missing(IN)){
    stop("'Falta el parámetro de tipo de (input) 'IN'. 0 si las entradas son las respuestas
         dicotómicas, 1 si se ingresa la matriz de momentos")
  }
  if(missing(EX)){
    stop("'Falta el parámetro de análisis exploratorio 'EX'")
  }
  if(missing(IV)){
    stop("Falta el parámetro 'IV' de inicialización de F y P")
  }
  if(missing(PR)){
    stop("Falta parámetro de Impresión 'PR'")
  }
  
  ## Construcción de F y de P en caso de análisis exploratorio ##
  F = matrizEscalon(N,M)  
  P = diag(1,M)
  PP <<-P
  
  # Matrix de primer momento
  B <<- (1/NS) * t(data) %*% data
  
  # Matriz de covarianzas
  itemCov = cov(data)
 
  # vector d (distancias)
  d <<- rep(0,N)

  for(i in 1:N){
    if(M > 1){
      d[i] <<- dj(F[i,],P)      
    }else{
      d[i] <<- dj(F[i],P)
    }
   
  }
 
  #Vector e
  e = ej(d)
  A.est = B
  
  # Calculo de b  
  b0 <<- diag(A.est)
  f0.inicial = qnorm(b0)  
  #f0.inicial = qnorm((b0-C)/(1-C)) * e  
  b1 <<- (1-C)*dnorm(f0.inicial/e)*(d/e)
  b2 <<- (1-C)*(f0.inicial/e)*dnorm(f0.inicial/e)*(d/e)^2 / sqrt(2)
  b3 <<- (1-C)*((f0.inicial/e)^2 -1)*dnorm(f0.inicial/e)*(d/e)^3 / sqrt(3)
   
  NN <<- N
  opt = optim(F,fn,NULL,method = "BFGS",
              control = list(maxit=maxiteraciones,abstol=epsilon))
  
  if(opt$convergence != 0){
    stop(cat("El proceso no converge con ", 100, " iteraciones a una tolerancia de "
             ,epsilon))
  }
  
  F = opt$par
  A.est = matrix(rep(0,N*N),ncol=N)
  for(i in 1:N){
    for(j in 1:i){
      A.est[j,i] = b0[i]*b0[j]+
        b1[i]*b1[j] * ((t(F[i,]) %*% P %*% F[j,]) / d[i]*d[j]) +
        b2[i]*b2[j] * ((t(F[i,]) %*% P %*% F[j,]) / d[i]*d[j])^2 +
        b3[i]*b3[j] * ((t(F[i,]) %*% P %*% F[j,]) / d[i]*d[j])^3
      A.est[i,j] = A.est[j,i]
    }
  }
  
  
  list("Primer.momento"=B,
       "Covarianzas" = itemCov,
       "initf_0" = f0.inicial,
       "Residuals" = B - A.est,
       "FinalF" = opt$par,
       "SumCuadRes" = opt$value,
       "Iter" = opt$counts
  )
}

fn =function(F){
  N = NN
  if(is.vector(F)){
    F = matrix(F,nrow=5)
  }
  A.est = matrix(rep(0,N*N),ncol=N)
  for(i in 1:N){
    for(j in 1:i){
      A.est[j,i] = b0[i]*b0[j]+
        b1[i]*b1[j] * ((t(F[i,]) %*% PP %*% F[j,]) / d[i]*d[j]) +
        b2[i]*b2[j] * ((t(F[i,]) %*% PP %*% F[j,]) / d[i]*d[j])^2 +
        b3[i]*b3[j] * ((t(F[i,]) %*% PP %*% F[j,]) / d[i]*d[j])^3
      A.est[i,j] = A.est[j,i]
    }
  }
  resta = (B - A.est)^2
  for(i in 1:N){
    resta[i,i] = 0
  }
  sum(resta)
}

## polinomios de hermite ##
h2 = function(x){
  (x^2-1) / sqrt(2) 
}

h3 = function(x){
  (x^2-3*x) / sqrt(6)
}

## distancia d_j ##

dj = function(f,p){
  (t(f) %*% p %*% f) ^(1/2)
}

## calculo ej ##

ej = function(d){
  sqrt((1+d)^2)
}

## construcción matriz escalon ##

matrizEscalon = function(n,p){
  if(p > 1){
    A = matrix(rep(1,n*p),ncol=p)
    for(i in 2:p){
      for(j in 1:i-1){
        A[j,i] = 0
      }
    }
  }else{
    A = rep(1,n)
  }
  A
}

## Gradiente Conjugado para minimizar producto##

conjGrad2 = function(A,b,x0,epsilon=10^-8){
  r0 = b - A %*% x0
  tol = norm(r0)
  i=0
  while(tol > epsilon){
    r0 = b - A %*% x0
    rho = (t(r0) %*% r0) / (t(r0) %*% A %*% r0)      
    x0 = x0 + as.double(rho) * r0
    tol = norm(r0)
    i=i+1
  }
  x0
}

## Datos ##
library( mirt)

#eliminar todas las variables
#rm(list = ls(all = TRUE))


# Base LSAT7 viene agrupada por frecuencias
# así que se replica cada individuo por su frecuencia
data(LSAT7)
conv = function(LSAT7){
  LsatMio=matrix(rep(0,1000),ncol=5,nrow=1000)
  LSAT72 = as.matrix(LSAT7)
  m = 1
  for(i in 1:32){
    inds = LSAT7[i,6]
    for(j in 1:inds){
      LsatMio[m,] = LSAT72[i,-6]
      m = m + 1
    }
  }
  LsatMio 
}

LsatMio = conv(LSAT7)

## Llamado ##

C = rep(0,5)
salida = miNoharm(N = 5,M = 3,1000,1,1,0,0,0,C = C,data = LsatMio,maxiteraciones = 1000)

## Matriz de momentos
salida$Primer.momento

## Matriz de covarianzas
salida$Covarianzas

## f_0 inicial
salida$initf_0

## Residuales
salida$Residuals

## Suma de cuadrados de los residuales
salida$SumCuadRes

## Número de iteraciones necesarias para la convergencia
salida$Iter

## Estimación F
salida$FinalF

