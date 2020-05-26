tn <- 3
h <- 0.01
m <- 50
s0 <- 100
r <- 0.05
sigma <- 0.2
N <- tn/h

w <- replicate(m, cumsum(rnorm(N, sd = h^(1/2))))

y <- matrix(nrow = N, ncol = m)
y[1, ] <- rep(s0, m)
for (i in 1:(N-1))
  y[i+1, ] <- y[i, ]+h*y[i, ]*r+sigma*y[i, ]*(w[i+1, ]-w[i, ])

ymean <- apply(y[,1:15], 1, mean)

t <- (1:300)/100
s <- s0*exp((r-0.5*sigma^2)*t+sigma*w[,50])
