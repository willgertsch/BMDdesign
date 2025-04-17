# model and gradient functions

# logistic #################################################################
logistic.fun = function(x, theta) {
  eta = theta[1] + theta[2] * x
  return(1/(1+exp(-eta)))
}

logistic.grad= function(x, theta) {

  eta = theta[1] + theta[2] * x
  sigma = exp(eta)/(1 + exp(eta))^2
  grad = sigma * c(1, x)
  return(grad)
}

logistic.bmdgrad = function(r, theta) {

  beta0 = theta[1]
  beta1 = theta[2]

  g1 = (exp(beta0) + 1)*r*(exp(beta0)*(r+1)+r-1)/(beta1*(exp(beta0)*r+r-1)*(exp(beta0)*(r+1)+r))
  g2 = suppressWarnings(- log(- (exp(beta0)*(exp(beta0)*r+r-1))/(exp(beta0)*(r+1)+r))/beta1^2)
  return(c(g1, g2))

}

# quantal linear ###############################################################
qlinear.fun = function(x, theta) {
  g = 1/(1 + exp(-theta[1]))
  b = theta[2]
  return(g + (1 - g) * (1 - exp(-b * x)))
}

qlinear.grad = function(x, theta) {
  g = 1/(1 + exp(-theta[1]))
  b = theta[2]
  g1 = exp(-b*x)
  g2 = (1 - g) * x * g1
  return(c(g1, g2))
}

# using extra risk definition
qlinear.bmd = function(r, theta) {
  g = 1/(1 + exp(-theta[1]))
  b = theta[2]
  return(-log(1-r)/b)
}


qlinear.bmdgrad = function(r, theta) {
  g = 1/(1 + exp(-theta[1]))
  b = theta[2]
  g1 = 0
  g2 = log(1-r)/b^2
  return(c(g1, g2))
}

# 3 parameter log-logistic #####################################################
loglogistic.fun = function(x, theta) {
  g = 1/(1 + exp(-theta[1]))
  a = theta[2]
  b = theta[3]
  return(suppressWarnings(g + (1-g)/(1+exp(-a-b*log(x)))))
}

loglogistic.grad = function(x, theta) {
  g = 1/(1 + exp(-theta[1]))
  a = theta[2]
  b = theta[3]

  g1 = 1/(exp(a)*x^(b)+1)
  g2 = suppressWarnings((1-g)*exp(-a-b*log(x))/((exp(-a-b*log(x))+1)^2))
  g3 = suppressWarnings(g2 * log(x))
  return(c(g1, g2, g3))

}

loglogistic.bmdgrad = function(r, theta) {
  g = 1/(1 + exp(-theta[1]))
  a = theta[2]
  b = theta[3]
  g1 = 0
  g2 = -exp(-a/b)*(r/(1-r))^(1/b) / b
  g3 = exp(-a/b) * (r/(1-r))^(1/b) * (a-log(r/(1-r))) / (b^2)
  return(c(g1,g2,g3))
}

loglogistic.bmd = function(r, theta) {
  g = 1/(1 + exp(-theta[1]))
  a = theta[2]
  b = theta[3]
  bmd = (log(r/(1-r)) - a)/b
  return(bmd)
}

# Weibull ######################################################################
# using version found in BMDS
# P[dose] = g + (1 - g) * (1 - exp(-b * dose^a))
# g is reparameterized
# code assumes output from ToxicR fit $parameters
weibull.fun = function(x, theta) {
  g = 1/(1 + exp(-theta[1]))
  a = theta[2]
  b = theta[3]

  return(g + (1 - g) * (1 - exp(-b * x^a)))
}

weibull.grad = function(x, theta) {

  g = 1/(1 + exp(-theta[1]))
  a = theta[2]
  b = theta[3]

  g1 = exp(-b * x^a)
  g2 = -b * (g - 1) * x^a * log(x) * exp(-b * x^a)
  g3 = (g - 1) * x^a * (-exp(-b * x^a))
  return(c(g1, g2, g3))
}

weibull.bmdgrad = function(r, theta) {

  g = 1/(1 + exp(-theta[1]))
  a = theta[2]
  b = theta[3]

  g1 = 0
  g2 = suppressWarnings(- log(-log(1-r)/b)*(-log(1-r)/b)^(1/a) / a^2)
  g3 = suppressWarnings(- (-log(1-r)/b)^(1/a) / (a*b))
  return(c(g1, g2, g3))
}

# dichotomous hill #############################################################
hill.fun = function(x, theta) {
  g = 1 / (1 + exp(-theta[1]))
  n = 1 / (1 + exp(-theta[2]))
  a = theta[3]
  b = theta[4]
  return(suppressWarnings(g + (1 - g) * n * (1 / (1 + exp(-a - b * log(x))))))
}

hill.grad = function(x, theta) {
  g = 1 / (1 + exp(-theta[1]))
  n = 1 / (1 + exp(-theta[2]))
  a = theta[3]
  b = theta[4]
  s = suppressWarnings(1/(1+exp(-a-b*log(x))))
  g1 = 1 - n * s
  g2 = (1 - g) * s
  g3 = (1 - g) * n * s * (1-s)
  g4 = suppressWarnings(g3 * log(x))
  return(c(g1, g2, g3, g4))
}

# commented out version:
# close, but not exactly matching what the internal code does
# see https://github.com/NIEHS/ToxicR/blob/main/src/include/DichHillBMD_NC.h
# is their documentation incorrect?
# uncommented code is my derivation
# this gives the same result as their internal code when I test it
hill.bmd = function(r, theta) {
  g = 1 / (1 + exp(-theta[1]))
  n = 1 / (1 + exp(-theta[2]))
  a = theta[3]
  b = theta[4]
  #return(exp((-a - log(-(r - n + g*n)/(r)))/(b)))
  return(suppressWarnings(exp((-a-log(n/r - 1))/(b))))
}

hill.bmdgrad = function(r, theta) {
  g = 1 / (1 + exp(-theta[1]))
  n = 1 / (1 + exp(-theta[2]))
  a = theta[3]
  b = theta[4]
  bmd = suppressWarnings(exp((-a-log(n/r - 1))/(b)))
  g1 = 0
  g2 = -(1/b) * bmd * 1/(n-r)
  g3 = -(1/b) * bmd
  g4 = suppressWarnings(bmd * (a + log(n/r - 1))/(b^2))
  return(c(g1, g2,g3, g4))
}
