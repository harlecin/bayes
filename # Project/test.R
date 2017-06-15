data = as.data.frame(Seatbelts)

n <- nrow(data)
x = data$PetrolPrice
X = matrix(c(rep(1, n), x), nrow=n)
y = data$kms
p <- ncol(X)



## sampling steps
draws <- 10000
burnin <- 1000

## initialize some space to hold results:
res <- matrix(as.numeric(NA), nrow=draws+burnin+1, ncol=p+1)
colnames(res) <- c(paste("beta", 0:(p-1), sep='_'), "sigma2")

## prior values for beta
B0 = diag(0.1^2, nrow = p, ncol = p)
b0 = c(5000, 90000)
B0inv <- solve(B0)

## prior value for sigma^2
c0 <- 1
C0 <- 1

## starting values:
res[1, ] = rep(1, p+1)

## precalculating values
B0_inv = solve(B0)
### (X'X+B0_inv)^(-1)
pre1 = solve(crossprod(X) + B0_inv)

### (X'y + B0_inv*b0)
pre2 = t(X) %*% y + B0_inv %*% b0

### cn = c0 + n/2 + p/2
pre3 = c0 + n/2 + p/2

## sampling steps
for (i in 2:(draws+burnin+1)) {
 if (i%%100 == 0) cat("Iteration", i, "done.\n")
 ## beta sampling
 Bn = res[i-1, p+1] * pre1
 bn = Bn %*% (pre2/res[i-1,p+1])
 ## block draw betas:
 res[i,1:p] <- mvtnorm::rmvnorm(1, bn, Bn)
 
 # draw sigma^2:
 Cn = C0 + .5*(crossprod(y - X %*% res[i, 1:p]) + t(res[i, 1:p] - b0) %*% B0_inv %*% (res[i, 1:p] - b0))
 res[i, p+1] <- 1/rgamma(1, shape = pre3, rate = Cn)
}

## throw away initial value and burn-in:
res <- res[-(1:(burnin+1)),]

dinvgamma <- function(x, a, b) {
 b^a/gamma(a) * x^(-a-1) * exp(-b/x)
}

## posterior mean:
print(colMeans(res))

##OLS estimate:
print(ols <- lm(y~X+0))


