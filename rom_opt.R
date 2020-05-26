rm(list = ls())
#======
mu <- 5/6
lambda <- mu
a1 <- 1/6
a2 <- 5/6
a3 <- a1
T <- 1
p0 <- 24
q0 <- 16
r0 <- q0
#======

cur_u <- function(t1, x1, t2, x2){
  p1 <- x1[1]; q1 <- x1[2]; r1 <- x1[3];
  p2 <- x2[1]; q2 <- x2[2]; r2 <- x2[3];
  u1 <- ((p2-p1)/(t2-t1)+.25*(lambda-mu)*(q2+q1)*(r2+r1))/a1
  u2 <- (mu*(q2-q1)/(t2-t1)+.25*(1-lambda)*(p2+p1)*(r2+r1))/a2
  u3 <- (lambda*(r2-r1)/(t2-t1)+.25*(mu-1)*(p2+p1)*(q2+q1))/a3
  return (c(u1 = u1, u2 = u2, u3 = u3))
}

J <- function(k, eps) {
  ucur <-  cur_u(t[k], x[,k], t[k+1], x[,k+1])
  uprev <- cur_u(t[k-1], x[,k-1], t[k], x[,k])
  # u <- sapply(1:N, function(i) cur_u(t[i], x[,i], t[i+1], x[,i+1]))
  # u[,k-1] <- uprev
  # u[,k] <- ucur
  u <- cbind(ucur, uprev)
  return (sum(apply(u, 2, function(x) sum(abs(x)+eps*(x^2))*tau)))
}


#======

Js <- data.frame(J=NA, eps=NA, N=NA, h=NA)

for (eps in c(0.1, 0.05, 2e-10, 2e-18, 0)){
  N <- 5
  tau <-T/N
  t <- (0:N)*tau
  if (eps==0.1){
    x <- array(dim = c(3, N+1), dimnames = list(c('p', 'q', 'r'), t))
    x['p',] <- p0*(-t+1)
    x['q',] <- q0*(-t+1)
    x['r',] <- r0*(-t+1)
  }
  while (N*2 < 100){
    N <- N*2
    tau <- T/N
    t <- (0:N)*tau
    if (eps < .1 & N == 10)
      x <- x[, as.character(t)]
    newx <- array(dim = c(3, N+1), dimnames = list(c('p','q','r'), t))
    newx[, colnames(x)] <- x
    difnames <- setdiff(colnames(newx), colnames(x))
    for (each in difnames){
      k <- which(each==colnames(newx))
      newx[,k] <- newx[,k-1]+(newx[,k+1]-newx[,k-1])*(t[k]-t[k-1])/(t[k+1]-t[k-1])
    }
    x <- newx
    u <- array(dim = c(3, N+1), dimnames = list(1:3, 0:N))
    h <- diag(0.1, nrow = 3, ncol = 3)
    while (h[1, 1] > 0.00001){
      Jcur <- 1000
      Jprev <- 1001
      while (Jprev-Jcur>1e-10){
        Jprev <- Jcur
        for (k in 2:N){
          u[, 1] <- cur_u(t[1], x[, 1], t[2], x[, 2])
          for (j in 1:3){
            x[k, ] <- x[k, ]+h[,j]
            J.pos <- get.J(eps)
            
            x[k, ] <- x[k, ]-2*h[,j]
            J.neg <- get.J(eps)
            
            x[k, ] <- x[k, ]+h[,j]
            J.null <- get.J(eps)
            
            if(J.pos < J.null & J.pos <= J.neg){x[k, ] <- x[k, ]+h[j, ]}
            else if(J.neg < J.null & J.neg<=J.pos){x[k, ] <- x[k, ]-h[j,]}
            
            u[k, ] <- get.u(t[k], x[k, ], t[k+1], x[k+1, ])
          }
        }
        res <- sapply(1:N, function(k) cur_u(t[k], x[,k], t[k+1], x[,k+1]))
        Jcur <- sum(apply(res, 2, function(u) sum(abs(u))*tau))
        cat(Jcur, ' ', eps, ' ', h[1, 1], ' ', N, '\n')
      }
      Js <- rbind(Js, c(J=Jcur, eps=eps, N=N, h = h[1,1]))
      h <- h/2
    }
    write.csv(t(u), paste(c(N, eps, '.csv'), collapse = '_'), row.names = F, quote = F)
  }
}

#=====
write.table(t(x), 'rom_x3.csv', row.names = F, quote = F)
write.table(t(u), 'rom_u3.csv', row.names = F, quote = F)
write.table(Js, 'rom_J3.csv', row.names = F, quote = F)
#=====
x <- t(read.csv('rom_x3.csv', sep=' '))
u <- t(read.csv('rom_u3.csv', sep=' '))
Js <- read.csv('rom_J3.csv', sep=' ')


plot(x[1,], type='l', ylim = c(-20, 30), col='red')
abline(h=0)
lines(x[2,], col='blue')
lines(x[3,], col='green')
legend(80, 30, legend = c('p', 'q', 'r'), col=c('red', 'blue', 'green'), lwd = 1)

plot(u[1,], type='l', col='red', ylim = c(-300, 100))
lines(u[2,], col='green')
lines(u[3,], col='blue')
par(xpd=T)
legend(80, 30, legend = c('u1', 'u2', 'u3'), col=c('red', 'green', 'blue'), lwd = 1)

