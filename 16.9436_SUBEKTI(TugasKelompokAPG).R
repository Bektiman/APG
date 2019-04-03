#Nama  : Subekti
#Kelas : 3KS1
#NIM   : 16.9436
### Library
library(car)
library(ellipse)
data <- 
### Fungsi vektor mean
mv <- function(df) {
  M <- data.matrix(df)
  k <- ncol(M)
  n <- nrow(M)
  elm <- NULL
  for (i in  1:k) {
    val <- M[,i]
    rata <- 0
    for (j in 1:n) {
      jum <- val[j]
      rata <- jum +rata
      
    }
    elm[i] <- rata/n
  }
  vektor<- (matrix(t(elm)))
  return(vektor)
}
### Fungsi matrix var covar
vcm <- function(df) {
  M <- data.matrix(df)
  k <- ncol(M)
  n <- nrow(M)
  vec <- t(mv(df))
  M_mean <- matrix(data=1, nrow=n) %*% vec
  Dif <-M-M_mean
  VC <- (n-1)^-1 * t(Dif) %*% Dif
  return(VC)
}
### Coba Run

#### Fungsi matrix Corelasi R
kor <- function(df){
  vcm <- vcm(df)
  mv <- mv(df)
  M <- data.matrix(df)
  elm <- NULL
  for (i in  1:ncol(M)) {
    val <- M[,i]
    x <- 0
    dif2 <- 0
    for (j in 1:nrow(M)) {
      x <- val[j]
      dif <- (x-mv[i])^2
      
      dif2 <- dif + dif2
    }
    elm[i] <- dif2/(nrow(M)-1)
  }
  ds <- diag(elm^(1/2))
  corr <- solve(ds)%*%vcm%*%solve(ds)
  return(corr)
}

### Fungsi garis kontour 2Var

### untuk membuat kontur, disarankan menggunakan matriks bivariate normal


# function ----------------------------------------------------------------

###  1. Ini fungsi untuk generate matriks dengan baris dan kolom serta rentang nilai 
##      yang kita tentukan sendiri, tapi gak tau distribusi nya gimana


get_bvn <- function(n,meanx,meany,korelasi) {
  korelasi<-as.numeric(korelasi)
  s12<-korelasi*sqrt(meanx*meany)
  mat_sigma<-matrix(c(meanx,s12,s12,meany),nrow=2)
  
  rerata<-c(meanx,meany)
  data_bvn<-rmvnorm(n=n,mean = rerata, sigma = mat_sigma)
  return(data_bvn)
}

### 6. Fungsi untuk menggambar scatter plot, serta elips dengan level 20%,50%,75%,90%,95% dan 99% 

get_contour<-function(dataa){
  #  library(ellipse)
  plot(dataa, main = "Scatter Plot dan Garis Konturnya")
  rata<-mv(dataa)
  cov<-vcm(dataa)
  
  lines( ellipse( cov, centre = rata,level=0.20) , col='green')
  lines( ellipse( cov, centre = rata,level=0.50) , col='green')
  lines( ellipse( cov, centre = rata,level=0.75) , col='blue')
  lines( ellipse( cov, centre = rata,level=0.90) , col='blue')
  lines( ellipse( cov, centre = rata,level=0.95) , col='red')
  lines( ellipse( cov, centre = rata,level=0.99) , col='red')
}
#sigma <- c(10,-25,-25,100)
#df <-  Z <- matrix(rnorm(500),2,250)
#a <- c(0.5,0.25,0.90,0.95,0.75)
#kontour(df,a,sigma)


### Infenrnsia Multivariate
### Fungsi untuk uji hipotesis untuk Sigma diketahui
onessampleknown <- function(df,u0,sigmas){
  sigmas <- matrix(sigmas,(length(sigmas)^(1/2)),(length(sigmas)^(1/2)))
  rata <- matrix(t(u0))
  A <- mv(df)-rata
  z <- nrow(df)*t(A)%*%solve(sigmas)%*%A
  x <- qchisq(0.95 , df = ncol(df))
  
  uji <- if (z > x) {
    
    uji <- "Tolak H0"
  }else {
    uji <- "Gagal Tolak H0"
  }
  
  p <- list("Nilai Z" = z,"Nilai X" =x, "keputusan"=uji)
  return(p)
  
}
### Fungsi untuk uji hipotesis untuk Sigma tidak diketahui
onesampleunknown <- function(df,u0){
  rata <- matrix(t(u0))
  vcm <-vcm(df)
  A <- mv(df)-rata
  Thot <- nrow(df)*t(A)%*%solve(vcm)%*%A
  v <- nrow(df)-1
  p <- ncol(df)
  
  Fkonv <- qf(0.95, df1 = p, df2 = (v - p +1))*(v*p/(v-p+1))
  
  uji <- if (Thot > Fkonv) {
    
    uji <- "Tolak H0"
  }else {
    uji <- "Gagal Tolak H0"
  }
  
  p <- list("Nilai Thot" = Thot,"Nilai Fkonv" =Fkonv, "keputusan"=uji)
  return(p)
  
}
### Fungsi uji rata rata 2 samople multivariate
twosampleindep  <- function(df1,df2) {
  rata1 <- mv(df1)
  rata2 <- mv(df2)
  vcm1 <-vcm(df1)
  vcm2 <-vcm(df2)
  A <- rata1-rata2
  W1 <- (nrow(df1)-1)*vcm1
  W2 <- (nrow(df2)-1)*vcm2
  m <- nrow(df1)*nrow(df2)
  n <- nrow(df1)+nrow(df2)
  p <- ncol(df1)
  Spl <- (1/(nrow(df1)+nrow(df2)-2))*(W1 + W2)
  Thot <- (m/n)*t(A)%*%solve(Spl)%*%A
  
  Fkonv <- qf(0.95, df1 = p, df2 = (n - p - 1))*((n-2)*p/(n-p-1))
  
  uji <- if (Thot > Fkonv) {
    
    uji <- "Tolak H0"
  }else {
    uji <- "Gagal Tolak H0"
  }
  
  p <- list("Nilai Thot" = Thot,"Nilai Fkonv" =Fkonv, "keputusan"=uji)
  return(p)
}
### Fungsi uji rata rata 2 sample independen dengan data frame tidak diketahui
twosampleindep2  <- function(u1,u2,s1,s2,tot1,tot2,col) {
  u1 <- matrix(t(u1))
  u2 <- matrix(t(u2))
  s1 <- matrix(s1,length(s1)^(1/2),length(s1)^(1/2))
  s2 <- matrix(s2,length(s2)^(1/2),length(s2)^(1/2))
  A <-  u1-u2
  W1 <- (tot1-1)*s1
  W2 <- (tot2-1)*s2
  m <- tot1*tot2
  n <- tot1+tot2
  p <- col
  Spl <- (1/(n-2))*(W1 + W2)
  Thot <- (m/n)*t(A)%*%solve(Spl)%*%A
  
  Fkonv <- qf(0.99, df1 = p, df2 = (n - p - 1))*((n-2)*p/(n-p-1))
  
  uji <- if (Thot > Fkonv) {
    
    uji <- "Tolak H0"
  }else {
    uji <- "Gagal Tolak H0"
  }
  
  p <- list("Nilai Thot" = Thot,"Nilai Fkonv" =Fkonv, "keputusan"=uji)
  return(p)
}
### uji 2 sampel berpasangan
paired2 <- function(df1,df2){
  d <- matrix(df1)-matrix(df2)
  mv <- mv(d)
  vcm <-vcm(d)
  Thot <- nrow(df1)*t(mv)%*%solve(vcm)%*%mv
  v <- nrow(df1)-1
  p <- length(mv)
  
  Fkonv <- qf(0.95, df1 = p, df2 = (v - p +1))*(v*p/(v-p+1))
  
  uji <- if (Thot > Fkonv) {
    
    uji <- "Tolak H0"
  }else {
    uji <- "Gagal Tolak H0"
  }
  
  p <- list("Nilai Thot" = Thot,"Nilai Fkonv" =Fkonv, "keputusan"=uji)
  return(p)
  
}
### Uji 2 sampel berpasangan dataframe tidak diketahui
paired1 <- function(ud,sd,n){
  mv <- matrix(ud)
  vcm <-matrix(sd,length(sd)^(1/2),length(sd)^(1/2))
  Thot <- n*t(mv)%*%solve(vcm)%*%mv
  v <- n-1
  p <- length(mv)
  
  Fkonv <- qf(0.95, df1 = p, df2 = (v - p +1))*(v*p/(v-p+1))
  
  uji <- if (Thot > Fkonv) {
    
    uji <- "Tolak H0"
  }else {
    uji <- "Gagal Tolak H0"
  }
  
  p <- list("Nilai Thot" = Thot,"Nilai Fkonv" =Fkonv, "keputusan"=uji)
  return(p)
  
}
### Fungsi Likelihood ratio
lr <- function(df,u0){
  M <- data.matrix(df)
  rata <- matrix(u0)
  sigma1 <-vcm(df)
  n <- nrow(df)
  M_mean <- matrix(data=1, nrow=n) %*% t(rata)
  Dif <-M-M_mean
  sigma0 <- (n-1)^-1 * t(Dif) %*% Dif
  v <- nrow(df)-1
  p <- ncol(df)
  Thot <- (v*det(sigma0)/det(sigma1))-v
  Fkonv <- qf(0.95, df1 = p, df2 = (v - p +1))*(v*p/(v-p+1))
  
  uji <- if (Thot > Fkonv) {
    
    uji <- "Tolak H0"
  }else {
    uji <- "Gagal Tolak H0"
  }
  
  p <- list("Nilai Thot" = Thot,"Nilai Fkonv" =Fkonv, "keputusan"=uji)
  return(p)
  
}


## Uji Individual variable ketika Kondisi H0 ditolak 
### CI T Hotteling 1 Sampel data tidak diketahui
RentangT1 <- function(u,s,n) {
  u <- matrix(t(u))
  s <- matrix(s, length(s)^(1/2),length(s)^(1/2))
  p <- length(u)
  v <- n-1
  Fkonv <- qf(0.95, df1 = p, df2 = (n - p))*((n-1)*p/(n-p))
  ds <- (matrix(diag(s)))/n
  akards <- ds^(1/2)
  rentang <- Fkonv^(1/2)*akards
  bawah <- u - rentang
  atas <- u + rentang
  p <- list("batas bawah" = bawah, "batas atas" = atas)
  return(p)
}
### CI T hotelling satu sampel 
RentangT2 <- function(df) {
  
  u <- mv(matrix(df))
  s <- vcm(matrix(df))
  p <- length(u)
  n <- nrow(df) 
  Fkonv <- qf(0.95, df1 = p, df2 = (n - p))*((n-1)*p/(n-p))
  ds <- (matrix(diag(s)))/n
  akards <- ds^(1/2)
  rentang <- Fkonv^(1/2)*akards
  bawah <- u - rentang
  atas <- u + rentang
  p <- list("batas bawah" = bawah, "batas atas" = atas)
  return(p)
} 




### Uji Dua Populasi Independen
### statitik uji pakai T hoteling dengan S pooled, matriks ragam gabungan 
duapopindep <- function(df1,df2,u1,u2){
  u1 <- matrix(u1)
  u2 <- matrix(u2)
  rata1 <- mv(df1)
  rata2 <- mv(df2)
  vcm1 <-vcm(df1)
  vcm2 <-vcm(df2)
  A <- rata1-rata2
  B <- u1-u2
  C <- A-B
  W1 <- (nrow(df1)-1)*vcm1
  W2 <- (nrow(df2)-1)*vcm2
  m <- nrow(df1)*nrow(df2)
  n <- nrow(df1)+nrow(df2)
  p <- ncol(df1)
  Spl <- ((W1+W2)/n-2)
  Thot <- t(C)%*%solve((n/m)*Spl)%*%C
  
  Fkonv <- qf(0.95, df1 = p, df2 = (n - p - 1))*((n-2)*p/(n-p-1))
  
  uji <- if (Thot > Fkonv) {
    
    uji <- "Tolak H0"
  }else {
    uji <- "Gagal Tolak H0"
  }
  
  p <- list("Nilai Thot" = Thot,"Nilai Fkonv" =Fkonv, "keputusan"=uji)
  return(p)
  
  
}
### CI  T hotteling dua populasi independen
RentangTduapop <- function(df1,df2) {
  rata1 <- mv(df1)
  rata2 <- mv(df2)
  vcm1 <-vcm(df1)
  vcm2 <-vcm(df2)
  A <- rata1-rata2
  W1 <- (nrow(df1)-1)*vcm1
  W2 <- (nrow(df2)-1)*vcm2
  m <- nrow(df1)*nrow(df2)
  n <- nrow(df1)+nrow(df2)
  p <- ncol(df1)
  Spl <- ((W1+W2)/n-2)
  Fkonv <- qf(0.95, df1 = p, df2 = (n - p - 1))*((n-2)*p/(n-p-1))
  ds <- (matrix(diag(Spl)))*(n/m)
  akards <- ds^(1/2)
  rentang <- Fkonv^(1/2)*akards
  bawah <- A - rentang
  atas <- A + rentang
  p <- list("batas bawah" = bawah, "batas atas" = atas)
  return(p)
  
}
#### CI Benferoni dua populasi independen 
benferoniduapop  <- function(df1,df2,alpha) {
  rata1 <- mv(df1)
  rata2 <- mv(df2)
  vcm1 <-vcm(df1)
  vcm2 <-vcm(df2)
  A <- rata1-rata2
  W1 <- (nrow(df1)-1)*vcm1
  W2 <- (nrow(df2)-1)*vcm2
  m <- nrow(df1)*nrow(df2)
  n <- nrow(df1)+nrow(df2)
  p <- ncol(df1)
  Spl <- ((W1+W2)/n-2)
  t <- qt(1-(alpha/(2*p)),df = (n-2))
  ds <- (matrix(diag(Spl)))*(n/m)
  akards <- ds^(1/2)
  rentang <- t*akards
  bawah <- A - rentang
  atas <- A + rentang
  p <- list("batas bawah" = bawah, "batas atas" = atas)
  return(p)
}

a <- matrix(c(1,2,3,4),2,2)
b <- matrix(c(2,3,4,6),2,2)
c <- a-b
c
