# Gauss-Newton algorithm
# A simple implementation of the Gauss-Newton algorithm applied on the Chwirut1 dataset.

library(NISTnls) # library containing the Chwirut1 dataset

# Load the dataset
x = Chwirut1$x
y = Chwirut1$y

# Explore the dataset
h = hist(y) # compute the histogram
par(mfrow=c(1,2))
plot(x,y, xlab="metal object distance", ylab="ultrasonic response", main="Chwirut1 dataset")
hist(y, xlab="response classes", ylab="frequencies", main="Chwirut1 histogram")
lines(h$mids, h$counts, col="red")

# Gauss-Newton algorithm

# Useful variables and coefficients initialization
iterNum = 0 # iteration number
maxNumIter = 10 # maximum number of iterations
eps = 1e-10 # threshold
e = 1

b = c(1.00E-01, 1.00E-02, 1.00E-03)

# Model training
while ((iterNum < maxNumIter) && (e > eps))
{
  cat("iteration = ", iterNum <- iterNum + 1, "\n")

  Jb1 = -x * exp(-b[1]*x) / (b[2] + b[3]*x)
  Jb2 = -exp(-b[1]*x) / (b[2] + b[3]*x)^2
  Jb3 = -x * exp(-b[1]*x) / (b[2] + b[3]*x)^2
  
  J = matrix(c(Jb1,Jb2,Jb3), nrow=length(x), ncol=length(b))
  K = tcrossprod(solve(crossprod(J,J)),J)
  
  r = y - (exp(-b[1]*x) / (b[2] + b[3]*x))
  b = b + K %*% r
  
  e = sum((2 * t(J) * r)^2)
}

# Model evaluation
# Compute y-hat (estimate) and plot it
yHat = exp(-b[1]*x) / (b[2] + b[3]*x)

par(mfrow=c(1,2))
plot(x,y, xlim=range(x), ylim=range(y), xlab="metal object distance", ylab="ultrasonic response", main="Predictions")
# Plot in red the predicted model
lines(x[order(x)], yHat[order(x)], xlim=range(x), ylim=range(yHat), col="red")
# Residuals
plot(x, y - yHat, xlab="metal object distance", ylab="residuals", main="Residuals")