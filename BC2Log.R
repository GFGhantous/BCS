#'@title Box-Cox Type II Logistic (dBC2Log)
#'
#'@description The Box-Cox Type II Logistic distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Type II Logistic distribution.
#'
#'@usage dBC2Log(\eqn{x, \mu, \sigma, \nu})
#'
#'@usage pBC2Log(\eqn{q, \mu, \sigma, \nu})
#'
#'@usage qBC2Log(\eqn{p, \mu, \sigma, \nu})
#'
#'@usage rBC2Log(\eqn{n, \mu, \sigma, \nu})
#'
#'@param \eqn{\mu} vector of scale parameter values (related to median)
#'
#'@param \eqn{\sigma} vector of relative dispersion parameter values (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} vector of skew parameter values (given by the power transformation to symmetry)
#'
#'@param n number of observations to be generate
#'
#'@param q vector of quantiles
#'
#'@param x vector of values of the observed variable
#'
#'@param p vector of probabilities
#'
#'@details The probability density function of the  Box-Cox symmetric distribution is given by:
#'
#'@details \deqn{fx(x)=\frac{x^{\nu-1}}{\mu^{\nu}\sigma}\frac{r((\frac{1}{\sigma\nu}((\frac{x}{\mu})^\nu-1))^2)}{R(\frac{1}{\sigma|\nu|})}}, if \eqn{\nu != 0}
#'
#'@details \deqn{fx(x)=\frac{1}{x\sigma}r((\frac{1}{\sigma}log(\frac{x}{\mu}))^2)}, if \eqn{\nu = 0}
#'
#'@details for *x > 0*, where \eqn{R(w)=\int_{-\infty}^{w}r(u^2)du} for w that belongs to the real numbers
#'
#'@return dBC2Log() provides the probability density function, pBC2Log() provides the cumulative distribution function, qBC2Log() provides the quantiles, and rBC2Log() generates random numbers.
#'
#'@section Warning: The BC2Log distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the Type II Logistic distribution is:
#'
#'@note \deqn{r(u)=exp(-u^\frac{1}{2})(1+exp(-u^{\frac{1}{2}}))^{-2}}
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@examples- 
#'
#'rBC2Log(1,5,0.5,0.5)
#'
#'rBC2Log(100,5,0.5,0.5)
#'
#'@export rBC2Log

#random numbers generation

rBC2Log<- function(n,mu,sigma,nu){
  if(n<=0 | n!=as.integer(n)) stop(paste("n must be positive integer"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  if (nu!=0){
    z=rexp(n)*(1+rexp(n))^{-2}
    x=mu*(1+sigma*nu*z)^(1/nu)
  }
  if (nu==0){
    z=rexp(n)*(1+rexp(n))^{-2}
    x=mu*exp(sigma*z)
  }
  return(x)
}

#'@title Box-Cox Type II Logistic (dBC2Log)
#'
#'@description The Box-Cox Type II Logistic distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Type II Logistic distribution.
#'
#'@usage dBC2Log(\eqn{x, \mu, \sigma, \nu})
#'
#'@usage pBC2Log(\eqn{q, \mu, \sigma, \nu})
#'
#'@usage qBC2Log(\eqn{p, \mu, \sigma, \nu})
#'
#'@usage rBC2Log(\eqn{n, \mu, \sigma, \nu})
#'
#'@param \eqn{\mu} vector of scale parameter values (related to median)
#'
#'@param \eqn{\sigma} vector of relative dispersion parameter values (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} vector of skew parameter values (given by the power transformation to symmetry)
#'
#'@param n number of observations to be generate
#'
#'@param q vector of quantiles
#'
#'@param x vector of values of the observed variable
#'
#'@param p vector of probabilities
#'
#'@details The probability density function of the  Box-Cox symmetric distribution is given by:
#'
#'@details \deqn{fx(x)=\frac{x^{\nu-1}}{\mu^{\nu}\sigma}\frac{r((\frac{1}{\sigma\nu}((\frac{x}{\mu})^\nu-1))^2)}{R(\frac{1}{\sigma|\nu|})}}, if \eqn{\nu != 0}
#'
#'@details \deqn{fx(x)=\frac{1}{x\sigma}r((\frac{1}{\sigma}log(\frac{x}{\mu}))^2)}, if \eqn{\nu = 0}
#'
#'@details for *x > 0*, where \eqn{R(w)=\int_{-\infty}^{w}r(u^2)du} for w that belongs to the real numbers
#'
#'@return dBC2Log() provides the probability density function, pBC2Log() provides the cumulative distribution function, qBC2Log() provides the quantiles, and rBC2Log() generates random numbers.
#'
#'@section Warning: The BC2Log distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the Type II Logistic distribution is:
#'
#'@note \deqn{r(u)=exp(-u^\frac{1}{2})(1+exp(-u^{\frac{1}{2}}))^{-2}}
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@examples
#'
#'dBC2Log(5,5,0.5,0.5)
#'
#'dBC2Log(5,5,0.5,0)
#'
#'@export dBC2Log

#Probability density function

dBC2Log<- function(x,mu,sigma,nu){
  n<-max(length(x),length(mu),length(sigma),length(nu))
  fx<-rep(NA,n)
  if(any(x<=0)) stop(paste("x must be positive"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  for (i in 1:n){
    fx<-ifelse(x>0 & mu>0 & sigma>0 & nu>0, (((x^{nu-1})/(sigma*mu^nu))*exp(-(1/(sigma*nu))*(((x/mu)^nu-1)))*(1+exp(-(1/(sigma*nu))*(((x/mu)^nu-1))))^{-2})/(1-(1/(1+exp(1/(sigma*nu))))),
        ifelse(x>0 & mu>0 & sigma>0 & nu<0, (((x^{nu-1})/(sigma*mu^nu))*exp(-(1/(sigma*nu))*(((x/mu)^nu-1)))*(1+exp(-(1/(sigma*nu))*(((x/mu)^nu-1))))^{-2})/(1/(1+exp(1/(sigma*nu)))),
          (1/(sigma*x))*(exp(-((1/sigma)*log(x/mu)))*(1+exp(-((1/sigma)*log(x/mu))))^{-2})))
  }
  return(fx)
}

#'@title Box-Cox Type II Logistic (dBC2Log)
#'
#'@description The Box-Cox Type II Logistic distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Type II Logistic distribution.
#'
#'@usage dBC2Log(\eqn{x, \mu, \sigma, \nu})
#'
#'@usage pBC2Log(\eqn{q, \mu, \sigma, \nu})
#'
#'@usage qBC2Log(\eqn{p, \mu, \sigma, \nu})
#'
#'@usage rBC2Log(\eqn{n, \mu, \sigma, \nu})
#'
#'@param \eqn{\mu} vector of scale parameter values (related to median)
#'
#'@param \eqn{\sigma} vector of relative dispersion parameter values (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} vector of skew parameter values (given by the power transformation to symmetry)
#'
#'@param n number of observations to be generate
#'
#'@param q vector of quantiles
#'
#'@param x vector of values of the observed variable
#'
#'@param p vector of probabilities
#'
#'@details The probability density function of the  Box-Cox symmetric distribution is given by:
#'
#'@details \deqn{fx(x)=\frac{x^{\nu-1}}{\mu^{\nu}\sigma}\frac{r((\frac{1}{\sigma\nu}((\frac{x}{\mu})^\nu-1))^2)}{R(\frac{1}{\sigma|\nu|})}}, if \eqn{\nu != 0}
#'
#'@details \deqn{fx(x)=\frac{1}{x\sigma}r((\frac{1}{\sigma}log(\frac{x}{\mu}))^2)}, if \eqn{\nu = 0}
#'
#'@details for *x > 0*, where \eqn{R(w)=\int_{-\infty}^{w}r(u^2)du} for w that belongs to the real numbers
#'
#'@return dBC2Log() provides the probability density function, pBC2Log() provides the cumulative distribution function, qBC2Log() provides the quantiles, and rBC2Log() generates random numbers.
#'
#'@section Warning: The BC2Log distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the Type II Logistic distribution is:
#'
#'@note \deqn{r(u)=exp(-u^\frac{1}{2})(1+exp(-u^{\frac{1}{2}}))^{-2}}
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@examples
#'
#'pBC2Log(0.5,3,0.5,0.5)
#'
#'pBC2Log(0.5,3,0.5,0)
#'
#'@export pBC2Log

#Distribution function

pBC2Log<- function(q,mu,sigma,nu){
  n<-max(length(q),length(mu),length(sigma),length(nu))
  Fx<-rep(NA,n)
  if(any(q<=0)) stop(paste("q must be positive"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  for (i in 1:n){
    Fx<-ifelse(q>0 & mu>0 & sigma>0 & nu>0,((1/(1+exp(-((1/(sigma*nu))*(((q/mu)^nu-1))))))-(1-(1/(1+exp(-1/(sigma*nu))))))/(1-(1-(1/(1+exp(-1/(sigma*nu)))))),
               ifelse(q>0 & mu>0 & sigma>0 & nu<0, (1/(1+exp(-((1/(sigma*nu))*(((q/mu)^nu-1))))))/(1-(1/(1+exp(-1/(sigma*nu))))),
                      1/(1+exp(-((1/sigma)*log(q/mu))))))
  }
  return(Fx)
}

#'@title Box-Cox Type II Logistic (dBC2Log)
#'
#'@description The Box-Cox Type II Logistic distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Type II Logistic distribution.
#'
#'@usage dBC2Log(\eqn{x, \mu, \sigma, \nu})
#'
#'@usage pBC2Log(\eqn{q, \mu, \sigma, \nu})
#'
#'@usage qBC2Log(\eqn{p, \mu, \sigma, \nu})
#'
#'@usage rBC2Log(\eqn{n, \mu, \sigma, \nu})
#'
#'@param \eqn{\mu} vector of scale parameter values (related to median)
#'
#'@param \eqn{\sigma} vector of relative dispersion parameter values (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} vector of skew parameter values (given by the power transformation to symmetry)
#'
#'@param n number of observations to be generate
#'
#'@param q vector of quantiles
#'
#'@param x vector of values of the observed variable
#'
#'@param p vector of probabilities
#'
#'@details The probability density function of the  Box-Cox symmetric distribution is given by:
#'
#'@details \deqn{fx(x)=\frac{x^{\nu-1}}{\mu^{\nu}\sigma}\frac{r((\frac{1}{\sigma\nu}((\frac{x}{\mu})^\nu-1))^2)}{R(\frac{1}{\sigma|\nu|})}}, if \eqn{\nu != 0}
#'
#'@details \deqn{fx(x)=\frac{1}{x\sigma}r((\frac{1}{\sigma}log(\frac{x}{\mu}))^2)}, if \eqn{\nu = 0}
#'
#'@details for *x > 0*, where \eqn{R(w)=\int_{-\infty}^{w}r(u^2)du} for w that belongs to the real numbers
#'
#'@return dBC2Log() provides the probability density function, pBC2Log() provides the cumulative distribution function, qBC2Log() provides the quantiles, and rBC2Log() generates random numbers.
#'
#'@section Warning: The BC2Log distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the Type II Logistic distribution is:
#'
#'@note \deqn{r(u)=exp(-u^\frac{1}{2})(1+exp(-u^{\frac{1}{2}}))^{-2}}
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
#'qBC2Log(0.5,3,0.5,0.5)
#'
#'qBC2Log(0.5,3,0.5,0)
#'
#'@export qBC2Log

#Quantile

qBC2Log<- function(p,mu,sigma,nu){
  n<-max(length(p),length(mu),length(sigma),length(nu))
  xp<-rep(NA,n)
  if(any(p<0|p>1)) stop(paste("p must be between 0 and 1"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  for (i in 1:n){
    xp<-ifelse(p>=0 & p<=1 & mu>0 & sigma>0 & nu>0, mu*(1+sigma*nu*(log((1-((1-p)*(1/(1+exp(-1/(sigma*nu))))))/(1-(1-((1-p)*(1/(1+exp(-1/(sigma*nu))))))))))^(1/nu),
               ifelse(p>=0 & p<=1 & mu>0 & sigma>0 & nu<0, mu*(1+sigma*nu*(log((p*(1-(1/(1+exp(-1/(sigma*nu))))))/(1-(p*(1-(1/(1+exp(-1/(sigma*nu))))))))))^(1/nu),
                      mu*exp(sigma*(log(p/(1-p))))))
  }
  return(xp)
}
