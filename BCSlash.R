#'@title Box-Cox Slash (BCSlash)
#'
#'@description The Box-Cox Slash distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Slash distribution.
#'
#'@usage dBCSlash(\eqn{x, \mu, \sigma, \nu, q})
#'
#'@usage pBCSlash(\eqn{x, \mu, \sigma, \nu, q})
#'
#'@usage qBCSlash(\eqn{p, \mu, \sigma, \nu, q})
#'
#'@usage rBCSlash(\eqn{n, \mu, \sigma, \nu, q})
#'
#'@param \eqn{\mu} a specific value of scale parameter (related to median)
#'
#'@param \eqn{\sigma} a specific value of relative dispersion parameter (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} a specific value of skew parameter (given by the power transformation to symmetry)
#'
#'@param \eqn{q} a specific value of kurtosis parameter
#'
#'@param n number of observations to be generate
#'
#'@param x a specific value of observed variable
#'
#'@param p a specific value of probability
#'
#'@details The probability density function of the  Box-Cox symmetric distribution is given by:
#'
#'@details \deqn{fx(x)=\frac{x^{\nu-1}}{\mu^{\nu}\sigma}\frac{r((\frac{1}{\sigma\nu}((\frac{x}{\mu})^\nu-1))^2)}{R(\frac{1}{\sigma|\nu|})}}, if \eqn{\nu != 0}
#'
#'@details \deqn{fx(x)=\frac{1}{x\sigma}r((\frac{1}{\sigma}log(\frac{x}{\mu}))^2)}, if \eqn{\nu = 0}
#'
#'@details for *x > 0*, where \eqn{R(w)=\int_{-\infty}^{w}r(u^2)du} for w that belongs to the real numbers
#'
#'@return dBCSlash() provides the probability density function, pBCSlash() provides the cumulative distribution function, qBCSlash() provides the quantiles, and rBCCSlash() generates random numbers.
#'
#'@section Warning: The BCCSlash distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, q>0, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the Slash distribution is:
#'
#'@note \deqn{r(u)=\frac{\psi(\frac{q+1}{2},\frac{u}{2})q2^{\frac{q}{2}-1}}{\sqrt{\pi}u^{\frac{q+1}{2}}}, q > 0}, where \eqn{\psi(a, x)=\int_{0}^{x}t^{a-1}e^{-t}dt} is the incomplete gamma function.
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@examples -
#'
#'rBCSlash(1,5,0.5,0.5,4)
#'
#'rBCSlash(100,5,0.5,0.5,4)
#'
#'@export rBCSlash

#random numbers generation

rBCSlash<- function(n,mu,sigma,nu,q){
  if(n<=0|n!=as.integer(n)) stop(paste("n must be positive integer"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  if(any(q<=0)) stop(paste("q must be positive"))
  if (nu!=0){
    z=rnorm(n,0,1)/((runif(n,0,1))^{1/q})
    x=mu*(1+sigma*nu*z)^(1/nu)
  }
  if (nu==0){
    z=rnorm(n,0,1)/((runif(n,0,1))^{1/q})
    x=mu*exp(sigma*z)
  }
  return(x)
}

#'@title Box-Cox Slash (BCSlash)
#'
#'@description The Box-Cox Slash distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Slash distribution.
#'
#'@usage dBCSlash(\eqn{x, \mu, \sigma, \nu, q})
#'
#'@usage pBCSlash(\eqn{x, \mu, \sigma, \nu, q})
#'
#'@usage qBCSlash(\eqn{p, \mu, \sigma, \nu, q})
#'
#'@usage rBCSlash(\eqn{n, \mu, \sigma, \nu, q})
#'
#'@param \eqn{\mu} a specific value of scale parameter (related to median)
#'
#'@param \eqn{\sigma} a specific value of relative dispersion parameter (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} a specific value of skew parameter (given by the power transformation to symmetry)
#'
#'@param \eqn{q} a specific value of kurtosis parameter
#'
#'@param n number of observations to be generate
#'
#'@param x a specific value of observed variable
#'
#'@param p a specific value of probability
#'
#'@details The probability density function of the  Box-Cox symmetric distribution is given by:
#'
#'@details \deqn{fx(x)=\frac{x^{\nu-1}}{\mu^{\nu}\sigma}\frac{r((\frac{1}{\sigma\nu}((\frac{x}{\mu})^\nu-1))^2)}{R(\frac{1}{\sigma|\nu|})}}, if \eqn{\nu != 0}
#'
#'@details \deqn{fx(x)=\frac{1}{x\sigma}r((\frac{1}{\sigma}log(\frac{x}{\mu}))^2)}, if \eqn{\nu = 0}
#'
#'@details for *x > 0*, where \eqn{R(w)=\int_{-\infty}^{w}r(u^2)du} for w that belongs to the real numbers
#'
#'@return dBCSlash() provides the probability density function, pBCSlash() provides the cumulative distribution function, qBCSlash() provides the quantiles, and rBCCSlash() generates random numbers.
#'
#'@section Warning: The BCCSlash distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, q>0, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the Slash distribution is:
#'
#'@note \deqn{r(u)=\frac{\psi(\frac{q+1}{2},\frac{u}{2})q2^{\frac{q}{2}-1}}{\sqrt{\pi}u^{\frac{q+1}{2}}}, q > 0}, where \eqn{\psi(a, x)=\int_{0}^{x}t^{a-1}e^{-t}dt} is the incomplete gamma function.
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@examples -
#'
#'dBCSlash(4,5,0.5,0.5,4)
#'
#'dBCSlash(4,5,0.5,0,4)
#'
#'@export dBCSlash

#Probability density function

dBCSlash<- function(x,mu,sigma,nu,q){
  n<-max(length(x),length(mu),length(sigma),length(nu))
  fx<-rep(NA,n)
  if(any(x<=0)) stop(paste("x must be positive"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  if(any(q<=0)) stop(paste("q must be positive"))
  for (i in 1:n){
    fx<-ifelse (nu>0 & (1/(sigma*nu))*(((x/mu)^nu-1))!=0,
                (x^{nu-1})/(sigma*mu^nu)*((gamma((q+1)/2)*(pgamma((((1/(sigma*nu))*(((x/mu)^nu-1)))^2)/2,(q+1)/2,scale=1)*q*2^{(q/2)-1})/(sqrt(pi)*(abs((1/(sigma*nu))*(((x/mu)^nu-1))))^{q+1}))/(as.numeric(integrate(function(u) q*(u^(q-1))*(pnorm(u/(sigma*abs(nu)),mean=0,sd=1)), lower = 0, upper = 1)$value))),
                ifelse (nu>0 & (1/(sigma*nu))*(((x/mu)^nu-1))==0,
                        (x^{nu-1})/(sigma*mu^nu)*((q/((q+1)*(sqrt(2*pi))))/(as.numeric(integrate(function(u) q*(u^(q-1))*(pnorm(u/(sigma*abs(nu)),mean=0,sd=1)), lower = 0, upper = 1)$value))),
                        ifelse (nu<0 & (1/(sigma*nu))*(((x/mu)^nu-1))!=0,
                                (x^{nu-1})/(sigma*mu^nu)*((gamma((q+1)/2)*(pgamma((((1/(sigma*nu))*(((x/mu)^nu-1)))^2)/2,(q+1)/2,scale=1)*q*2^{(q/2)-1})/(sqrt(pi)*(abs((1/(sigma*nu))*(((x/mu)^nu-1))))^{q+1}))/(as.numeric(integrate(function(u) q*(u^(q-1))*(pnorm(u/(sigma*abs(nu)),mean=0,sd=1)), lower = 0, upper = 1)$valeu))),
                                ifelse (nu<0 & (1/(sigma*nu))*(((x/mu)^nu-1))==0,
                                        (x^{nu-1})/(sigma*mu^nu)*((q/((q+1)*sqrt(2*pi)))/(as.numeric(integrate(function(u) q*(u^(q-1))*(pnorm(u/(sigma*abs(nu)),mean=0,sd=1)), lower = 0, upper = 1)$value))),
                                        ifelse (nu==0 & (1/sigma)*log(x/mu)!=0,
                                                (1/(sigma*x))*gamma((q+1)/2)*((pgamma(((1/sigma)*log(x/mu))^2/2,(q+1)/2,scale=1)*q*2^{(q/2)-1})/(sqrt(pi)*(abs((1/sigma)*log(x/mu)))^{q+1})),
                                                (1/(sigma*x))*(q/((q+1)*sqrt(2*pi))))))))
  }
  return(fx)
}
 
#'@title Box-Cox Slash (BCSlash)
#'
#'@description The Box-Cox Slash distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Slash distribution.
#'
#'@usage dBCSlash(\eqn{x, \mu, \sigma, \nu, q})
#'
#'@usage pBCSlash(\eqn{x, \mu, \sigma, \nu, q})
#'
#'@usage qBCSlash(\eqn{p, \mu, \sigma, \nu, q})
#'
#'@usage rBCSlash(\eqn{n, \mu, \sigma, \nu, q})
#'
#'@param \eqn{\mu} a specific value of scale parameter (related to median)
#'
#'@param \eqn{\sigma} a specific value of relative dispersion parameter (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} a specific value of skew parameter (given by the power transformation to symmetry)
#'
#'@param \eqn{q} a specific value of kurtosis parameter
#'
#'@param n number of observations to be generate
#'
#'@param x a specific value of observed variable
#'
#'@param p a specific value of probability
#'
#'@details The probability density function of the  Box-Cox symmetric distribution is given by:
#'
#'@details \deqn{fx(x)=\frac{x^{\nu-1}}{\mu^{\nu}\sigma}\frac{r((\frac{1}{\sigma\nu}((\frac{x}{\mu})^\nu-1))^2)}{R(\frac{1}{\sigma|\nu|})}}, if \eqn{\nu != 0}
#'
#'@details \deqn{fx(x)=\frac{1}{x\sigma}r((\frac{1}{\sigma}log(\frac{x}{\mu}))^2)}, if \eqn{\nu = 0}
#'
#'@details for *x > 0*, where \eqn{R(w)=\int_{-\infty}^{w}r(u^2)du} for w that belongs to the real numbers
#'
#'@return dBCSlash() provides the probability density function, pBCSlash() provides the cumulative distribution function, qBCSlash() provides the quantiles, and rBCCSlash() generates random numbers.
#'
#'@section Warning: The BCCSlash distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, q > 0, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the Slash distribution is:
#'
#'@note \deqn{r(u)=\frac{\psi(\frac{q+1}{2},\frac{u}{2})q2^{\frac{q}{2}-1}}{\sqrt{\pi}u^{\frac{q+1}{2}}}, q > 0}, where \eqn{\psi(a, x)=\int_{0}^{x}t^{a-1}e^{-t}dt} is the incomplete gamma function.
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@examples -
#'
#'pBCSlash(5,5,0.5,0.5,4)
#'
#'pBCSlash(5,5,0.5,0,4)
#'
#'@export pBCSlash

#Distribution function

pBCSlash<- function(x,mu,sigma,nu,q){
  n<-max(length(x),length(mu),length(sigma),length(nu), length(q))
  Fx<-rep(NA,n)
  if(any(x<=0)) stop(paste("x be positive"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  if(any(q<=0)) stop(paste("q must be positive"))
  for (i in 1:n){
    Fx<-ifelse(x>0 & mu>0 & sigma>0 & nu>0 & q>0,(((integrate((function(u) q*(u^(q-1))*(pnorm(((1/(sigma*nu))*(((x/mu)^nu)-1)),mean=0,sd=1))), lower = 0, upper = 1))$value)-(1-((integrate((function(u) q*(u^(q-1))*(pnorm(u/(sigma*abs(nu)),mean=0,sd=1))), lower = 0, upper = 1))$value)))/((integrate((function(u) q*(u^(q-1))*(pnorm(u/(sigma*abs(nu)),mean=0,sd=1))), lower = 0, upper = 1))$value),
               ifelse(x>0 & mu>0 & sigma>0 & nu<0 & q>0,((integrate((function(u) q*(u^(q-1))*(pnorm(((1/(sigma*nu))*(((x/mu)^nu)-1)),mean=0,sd=1))), lower = 0, upper = 1))$value)/((integrate((function(u) q*(u^(q-1))*(pnorm(u/(sigma*abs(nu)),mean=0,sd=1))), lower = 0, upper = 1))$value),
                      (integrate((function(u) q*(u^(q-1))*(pnorm(((1/sigma)*log(x/mu)),mean=0,sd=1))), lower = 0, upper = 1))$value))
  }
  return(Fx)
}

#'@title Box-Cox Slash (BCSlash)
#'
#'@description The Box-Cox Slash distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Slash distribution.
#'
#'@usage dBCSlash(\eqn{x, \mu, \sigma, \nu, q})
#'
#'@usage pBCSlash(\eqn{x, \mu, \sigma, \nu, q})
#'
#'@usage qBCSlash(\eqn{p, \mu, \sigma, \nu, q})
#'
#'@usage rBCSlash(\eqn{n, \mu, \sigma, \nu, q})
#'
#'@param \eqn{\mu} a specific value of scale parameter (related to median)
#'
#'@param \eqn{\sigma} a specific value of relative dispersion parameter (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} a specific value of skew parameter (given by the power transformation to symmetry)
#'
#'@param \eqn{q} a specific value of kurtosis parameter
#'
#'@param n number of observations to be generate
#'
#'@param x a specific value of observed variable
#'
#'@param p a specific value of probability
#'
#'@details The probability density function of the  Box-Cox symmetric distribution is given by:
#'
#'@details \deqn{fx(x)=\frac{x^{\nu-1}}{\mu^{\nu}\sigma}\frac{r((\frac{1}{\sigma\nu}((\frac{x}{\mu})^\nu-1))^2)}{R(\frac{1}{\sigma|\nu|})}}, if \eqn{\nu != 0}
#'
#'@details \deqn{fx(x)=\frac{1}{x\sigma}r((\frac{1}{\sigma}log(\frac{x}{\mu}))^2)}, if \eqn{\nu = 0}
#'
#'@details for *x > 0*, where \eqn{R(w)=\int_{-\infty}^{w}r(u^2)du} for w that belongs to the real numbers
#'
#'@return dBCSlash() provides the probability density function, pBCSlash() provides the cumulative distribution function, qBCSlash() provides the quantiles, and rBCCSlash() generates random numbers.
#'
#'@section Warning: The BCCSlash distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, q > 0, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the Slash distribution is:
#'
#'@note \deqn{r(u)=\frac{\psi(\frac{q+1}{2},\frac{u}{2})q2^{\frac{q}{2}-1}}{\sqrt{\pi}u^{\frac{q+1}{2}}}, q > 0}, where \eqn{\psi(a, x)=\int_{0}^{x}t^{a-1}e^{-t}dt} is the incomplete gamma function.
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@examples -
#'
#'qBCSlash(0.5,5,0.5,0.5,4)
#'
#'qBCSlash(0.5,5,0.5,0,4)
#'
#'@export qBCSlash

#Quantile

qBCSlash<- function(p,mu,sigma,nu,q){
  n<-max(length(p),length(mu),length(sigma),length(nu),length(q))
  xp<-rep(NA,n)
  if(any(p<0|p>1)) stop(paste("p must be between 0 and 1"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  if(any(q<=0)) stop(paste("q must be positive"))
  for (i in 1:n){
    xp<-ifelse(p>=0 & p<=1 & mu>0 & sigma>0 & nu>0 & q>0,mu*(1+sigma*nu*(qnorm(1-((1-p)*(pnorm(1/(sigma*abs(nu)),0,1)/((punif(1/(sigma*abs(nu)),0,1))^{1/q}))),0,1)/((qunif(1-((1-p)*(pnorm(1/(sigma*abs(nu)),0,1)/((punif(1/(sigma*abs(nu)),0,1))^{1/q}))))^{1/q}))))^(1/nu),
               ifelse(p>=0 & p<=1 & mu>0 & sigma>0 & nu<0 & q>0,mu*(1+sigma*nu*(qnorm(p*(pnorm(1/(sigma*abs(nu)),0,1)/((punif(1/(sigma*abs(nu)),0,1))^{1/q})))/((qunif(p*(pnorm(1/(sigma*abs(nu)),0,1)/((punif(1/(sigma*abs(nu)),0,1))^{1/q}))))^{1/q})))^(1/nu),
                      mu*exp(sigma*(qnorm(p)/((qunif(p))^{1/q})))))
  }
  return(xp)
}

