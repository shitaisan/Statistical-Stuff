library(gtools)

df <- read.csv("Insult1.csv", sep = ';', dec = ',')[,1:4]
fs <- read.csv("Insult1.csv", sep = ';', dec = ',')[,5:6]

# for (i in 1:ncol(df)){
#   while (length(boxplot(df[,i], plot=FALSE)$out)>0)
#     df <- df[-which(df[,i] %in% boxplot(df[,i], plot = FALSE)$out),]
# }
means <- apply(df, 2, mean)
cov_matr <- cov(df)
cor_matr <- cor(df)
cov_inv <- solve(cov_matr)
covcor_inv <- solve(cor_matr)
sds <- diag(nrow = 4, ncol = 4)*apply(df, 2, sd)


lm4_all <- lm(Velocity ~ Retrak+lagtime+AreaCurve, data=df)
errors <- lm4_all$residuals
cor(df[,-1], errors)            
resid_var <- mean(errors^2)

plot(x=df$Velocity, y=lm4_all$fitted.values)

mult_cor <- cor(df$Velocity, lm4_all$fitted.values)

n <- ncol(df)
privcor <- matrix(nrow=n, ncol=n)
for (i in 1:n){
  for (j in 1:n){
    privcor[i,j] <- -cov_inv[i,j]/sqrt(cov_inv[i,i]*cov_inv[j,j])
  }
}

# cov11_inv <- solve(cov_matr[1:3, 1:3])
# cov12 <- cov_matr[1:3, 4:6]
# cov21 <- cov_matr[4:6, 1:3]
# cov22_inv <- solve(cov_matr[4:6, 4:6])
# canon_cor <- sqrt(eigen(cov11_inv%*%cov12%*%cov22_inv%*%cov21, only.values = TRUE)[[1]])

# principal components
C <- eigen(cov_matr)$vectors
tdf <- t(data.matrix(df))
tdf <- tdf-means
princ_comp <- t(C)%*%tdf
plot(x = princ_comp[1,], y = princ_comp[2,], col = fs[,1])

sum(apply(princ_comp, 1, var)) - sum(apply(tdf, 1, var)) # sum of vars saved

# anova
f1 <- quantcut(princ_comp[1,], q = 4)
summary(aov(df[,1]~f1))

# manova for f1
gr1 <- df[fs[,1]==1,1:4]
gr2 <- df[fs[,1]==2,1:4]
gr3 <- df[fs[,1]==3,1:4]

n1 <- nrow(gr1)
n2 <- nrow(gr2)
n3 <- nrow(gr3)

n <- nrow(df)
p <- ncol(df)
J <- max(fs[,1])

mean1 <- apply(gr1, 2, mean)
mean2 <- apply(gr2, 2, mean)
mean3 <- apply(gr3, 2, mean)

globmean <- apply(df[,1:4], 2, mean)
L1 <- cov(gr1)
L2<- cov(gr2)
L3 <- cov(gr3)
L <- ((n1-1)*L1+(n2-1)*L2+(n3-1)*L3)/(nrow(df)-3)

# Hotelling-Lawley test with chisq approximation
Q <- n1*(mean1-globmean)%*%t(mean1-globmean)+n2*(mean2-globmean)%*%t(mean2-globmean)+n3*(mean3-globmean)%*%t(mean3-globmean)
sqT <- sum(diag(Q%*%solve(L)))
# n1*t(mean1-globmean)%*%solve(L)%*%(mean1-globmean)+n2*t(mean2-globmean)%*%solve(L)%*%(mean2-globmean)+n3*t(mean3-globmean)%*%solve(L)%*%(mean3-globmean)
1-pchisq(sqT, p*(J-1))
summary(manova(as.matrix(df)~as.factor(fs[,1])), test = "Hotelling-Lawley")

# Hotelling-Lawley test with fisher approximation
g1 <- p*(J-1)*(n-J-p)/(n-J-p-(J-2)*(p-1))
g2 <- n-J-p+1
const <- (n-J-p-1)/((n-J)*(J-1)*p)*g2/(g2-2)
sqT_alt <- const*sqT
1-pf(sqT_alt, g1, g2)

# lambda-wilks test
lambda <- prod(1/(1+eigen(Q%*%solve(L))$values[1:min(J-1, p)]))
1-pchisq((p/2+J/2+1-n)*log(lambda), p*(J-1))
summary(manova(as.matrix(df)~as.factor(fs[,1])), test = 'Wilks')

# manova 1-2 groups (Scheffe method)
sqT12 <- const*n1*n2/(n1+n2)*t(mean1-mean2)%*%solve(L)%*%(mean1-mean2)
1-pf(sqT12, g1, g2)

# manova 1-3 groups (Scheffe method)
sqT13 <- const*n1*n3/(n1+n3)*t(mean1-mean3)%*%solve(L)%*%(mean1-mean3)
1-pf(sqT13, g1, g2)

# manova 2-3 groups (Scheffe method)
sqT23 <- const*n2*n3/(n2+n3)*t(mean2-mean3)%*%solve(L)%*%(mean2-mean3)
1-pf(sqT12, g1, g2)
