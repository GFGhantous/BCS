#'@title Box-Cox Canonical Slash (BCCSlash)
#'
#'@description The Box-Cox Canonical Slash distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Canonical Slash distribution.
#'
#'@usage dBCCSlash(\eqn{x, \mu, \sigma, \nu})
#'
#'@usage pBCCSlash(\eqn{q, \mu, \sigma, \nu})
#'
#'@usage qBCCSlash(\eqn{p, \mu, \sigma, \nu})
#'
#'@usage rBCCSlash(\eqn{n, \mu, \sigma, \nu})
#'
#'@param \eqn{\mu} a specific value of scale parameter (related to median)
#'
#'@param \eqn{\sigma} a specific value of relative dispersion parameter (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} a specific value of skew parameter (given by the power transformation to symmetry)
#'
#'@param n number of observations to be generate
#'
#'@param q a specific value of quantil
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
#'@return dBCCSlash() provides the probability density function, pBCCSlash() provides the cumulative distribution function, qBCCSlash() provides the quantiles, and rBCCSlash() generates random numbers.
#'
#'@section Warning: The BCCSlash distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the canonical slash distribution is:
#'
#'@note \deqn{r(u)=\frac{1}{\sqrt{2\pi}u}(1-\exp(\frac{-u}{2}))}, for \eqn{u != 0} and \deqn{r(u)=\frac{1}{(2\sqrt{2\pi})}}, for \eqn{u = 0}.
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
#'rBCCSlash(1,5,0.5,0.5)
#'
#'rBCCSlash(100,5,0.5,0.5)
#'
#'@export rBCCSlash

#random numbers generation

rBCCSlash<- function(n,mu,sigma,nu){
  if(n<=0 | n!=as.integer(n)) stop(paste("n must be positive integer"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  if (nu!=0){
    q=1
    z=rnorm(n,0,1)/((runif(n,0,1))^{1/q})
    x=mu*(1+sigma*nu*z)^(1/nu)
  }
  if (nu==0){
    q=1
    z=rnorm(n,0,1)/((runif(n,0,1))^{1/q})
    x=mu*exp(sigma*z)
  }
  return(x)
}


#'@title Box-Cox Canonical Slash (BCCSlash)
#'
#'@description The Box-Cox Canonical Slash distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Canonical Slash distribution.
#'
#'@usage dBCCSlash(\eqn{x, \mu, \sigma, \nu})
#'
#'@usage pBCCSlash(\eqn{q, \mu, \sigma, \nu})
#'
#'@usage qBCCSlash(\eqn{p, \mu, \sigma, \nu})
#'
#'@usage rBCCSlash(\eqn{n, \mu, \sigma, \nu})
#'
#'@param \eqn{\mu} a specific value of scale parameter (related to median)
#'
#'@param \eqn{\sigma} a specific value of relative dispersion parameter (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} a specific value of skew parameter (given by the power transformation to symmetry)
#'
#'@param n number of observations to be generate
#'
#'@param q a specific value of quantil
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
#'@return dBCCSlash() provides the probability density function, pBCCSlash() provides the cumulative distribution function, qBCCSlash() provides the quantiles, and rBCCSlash() generates random numbers.
#'
#'@section Warning: The BCCSlash distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the canonical slash distribution is:
#'
#'@note \deqn{r(u)=\frac{1}{\sqrt{2\pi}u}(1-\exp(\frac{-u}{2}))}, for \eqn{u != 0} and \deqn{r(u)=\frac{1}{(2\sqrt{2\pi})}}, for \eqn{u = 0}.
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
#'dBCCSlash(4,5,0.5,0.5)
#'
#'dBCCSlash(4,5,0.5,0)
#'
#'@export dBCCSlash

#Probability density function

dBCCSlash<- function(x,mu,sigma,nu){
  n<-max(length(x),length(mu),length(sigma),length(nu))
  fx<-rep(NA,n)
  if(any(x<=0)) stop(paste("x must be positive"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  for (i in 1:n){
    fx<-ifelse (nu>0 & (1/(sigma*nu))*(((x/mu)^nu-1))!=0,
                (x^{nu-1})/(sigma*mu^nu)*(((1/(sqrt(2*pi)*((1/(sigma*nu))*(((x/mu)^nu-1)))^2))*(1-exp((-((1/(sigma*nu))*(((x/mu)^nu-1)))^2)/2)))/(as.numeric(integrate(function(u) 1*(u^(1-1))*(pnorm(u/(sigma*abs(nu)),mean=0,sd=1)), lower = 0, upper = 1)$value))),
                ifelse (nu>0 & (1/(sigma*nu))*(((x/mu)^nu-1))==0,
                        (x^{nu-1})/(sigma*mu^nu)*((1/(2*(sqrt(2*pi))))/(as.numeric(integrate(function(u) 1*(u^(1-1))*(pnorm(u/(sigma*abs(nu)),mean=0,sd=1)), lower = 0, upper = 1)$value))),
                        ifelse (nu<0 & (1/sigma)*(log(x/mu))!=0,
                                (x^{nu-1})/(sigma*mu^nu)*((1/(sqrt(2*pi)*((1/(sigma*nu))*(((x/mu)^nu-1)))^2))*(1-exp((-((1/(sigma*nu))*(((x/mu)^nu-1)))^2)/2))/(as.numeric(integrate(function(u) 1*(u^(1-1))*(pnorm(u/(sigma*abs(nu)),mean=0,sd=1)), lower = 0, upper = 1)$value))),
                                ifelse (nu<0 & (1/sigma)*(log(x/mu))==0,
                                        (x^{nu-1})/(sigma*mu^nu)*((1/(2*(sqrt(2*pi))))/(as.numeric(integrate(function(u) 1*(u^(1-1))*(pnorm(u/(sigma*abs(nu)),mean=0,sd=1)), lower = 0, upper = 1)$value))),
                                        ifelse (nu==0 & (1/sigma)*(log(x/mu))!=0,
                                                (1/(sigma*x))*(1/(sqrt(2*pi)*((1/sigma)*(log(x/mu)))^2))*(1-exp((-((1/sigma)*(log(x/mu)))^2)/2)),
                                                (1/(sigma*x))*(1/(2*(sqrt(2*pi)))))))))
  }
  return(fx)
}

#'@title Box-Cox Canonical Slash (BCCSlash)
#'
#'@description The Box-Cox Canonical Slash distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Canonical Slash distribution.
#'
#'@usage dBCCSlash(\eqn{x, \mu, \sigma, \nu})
#'
#'@usage pBCCSlash(\eqn{q, \mu, \sigma, \nu})
#'
#'@usage qBCCSlash(\eqn{p, \mu, \sigma, \nu})
#'
#'@usage rBCCSlash(\eqn{n, \mu, \sigma, \nu})
#'
#'@param \eqn{\mu} a specific value of scale parameter (related to median)
#'
#'@param \eqn{\sigma} a specific value of relative dispersion parameter (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} a specific value of skew parameter (given by the power transformation to symmetry)
#'
#'@param n number of observations to be generate
#'
#'@param q a specific value of quantil
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
#'@return dBCCSlash() provides the probability density function, pBCCSlash() provides the cumulative distribution function, qBCCSlash() provides the quantiles, and rBCCSlash() generates random numbers.
#'
#'@section Warning: The BCCSlash distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the canonical slash distribution is:
#'
#'@note \deqn{r(u)=\frac{1}{\sqrt{2\pi}u}(1-\exp(\frac{-u}{2}))}, for \eqn{u != 0} and \deqn{r(u)=\frac{1}{(2\sqrt{2\pi})}}, for \eqn{u = 0}.
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
#'pBCCSlash(5,5,0.5,0.5)
#'
#'pBCCSlash(5,5,0.5,0)
#'
#'@export pBCCSlash

#Distribution function

pBCCSlash<- function(q,mu,sigma,nu){
  n<-max(length(q),length(mu),length(sigma),length(nu))
  Fx<-rep(NA,n)
  if(any(q<=0)) stop(paste("q must be positive"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  for (i in 1:n){
    Fx<-ifelse(q>0 & mu>0 & sigma>0 & nu>0,((integrate((function(u) 1*(u^(1-1))*(pnorm(((1/(sigma*nu))*(((q/mu)^nu)-1)),mean=0,sd=1))), lower = 0, upper = 1)$value)-(1-((integrate((function(u) 1*(u^(1-1))*(pnorm(u/(sigma*abs(nu)),mean=0,sd=1))), lower = 0, upper = 1))$value)))/((integrate((function(u) 1*(u^(1-1))*(pnorm(u/(sigma*abs(nu)),mean=0,sd=1))), lower = 0, upper = 1))$value),
        ifelse(q>0 & mu>0 & sigma>0 & nu<0,((integrate((function(u) 1*(u^(1-1))*(pnorm(((1/(sigma*nu))*(((q/mu)^nu)-1)),mean=0,sd=1))), lower = 0, upper = 1))$value)/((integrate((function(u) 1*(u^(1-1))*(pnorm(u/(sigma*abs(nu)),mean=0,sd=1))), lower = 0, upper = 1))$value),
            (integrate(function(u) 1*(u^(1-1))*(pnorm(((1/sigma)*log(q/mu)),mean=0,sd=1)), lower = 0, upper = 1))$value))
  }
  return(Fx)
}

#'@title Box-Cox Canonical Slash (BCCSlash)
#'
#'@description The Box-Cox Canonical Slash distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Canonical Slash distribution.
#'
#'@usage dBCCSlash(\eqn{x, \mu, \sigma, \nu})
#'
#'@usage pBCCSlash(\eqn{q, \mu, \sigma, \nu})
#'
#'@usage qBCCSlash(\eqn{p, \mu, \sigma, \nu})
#'
#'@usage rBCCSlash(\eqn{n, \mu, \sigma, \nu})
#'
#'@param \eqn{\mu} a specific value of scale parameter (related to median)
#'
#'@param \eqn{\sigma} a specific value of relative dispersion parameter (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} a specific value of skew parameter (given by the power transformation to symmetry)
#'
#'@param n number of observations to be generate
#'
#'@param q a specific value of quantil
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
#'@return dBCCSlash() provides the probability density function, pBCCSlash() provides the cumulative distribution function, qBCCSlash() provides the quantiles, and rBCCSlash() generates random numbers.
#'
#'@section Warning: The BCCSlash distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the canonical slash distribution is:
#'
#'@note \deqn{r(u)=\frac{1}{\sqrt{2\pi}u}(1-\exp(\frac{-u}{2}))}, for \eqn{u != 0} and \deqn{r(u)=\frac{1}{(2\sqrt{2\pi})}}, for \eqn{u = 0}.
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
#'qBCCSlash(0.5,5,0.5,0.5)
#'
#'qBCCSlash(0.5,5,0.5,0)
#'
#'@export qBCCSlash

#Quantile

qBCCSlash<- function(p,mu,sigma,nu){
  n<-max(length(p),length(mu),length(sigma),length(nu))
  xp<-rep(NA,n)
  if(any(p<0|p>1)) stop(paste("p must be between 0 and 1"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  for (i in 1:n){
    xp<-ifelse(p>=0 & p<=1 & mu>0 & sigma>0 & nu>0,mu*(1+sigma*nu*(qnorm(1-((1-p)*(pnorm(1/(sigma*abs(nu)),0,1)/((punif(1/(sigma*abs(nu)),0,1))^{1/1}))),0,1)/((qunif(1-((1-p)*(pnorm(1/(sigma*abs(nu)),0,1)/((punif(1/(sigma*abs(nu)),0,1))^{1/1}))))^{1/1}))))^(1/nu),
               ifelse(p>=0 & p<=1 & mu>0 & sigma>0 & nu<0,mu*(1+sigma*nu*(qnorm(p*(pnorm(1/(sigma*abs(nu)),0,1)/((punif(1/(sigma*abs(nu)),0,1))^{1/1})))/((qunif(p*(pnorm(1/(sigma*abs(nu)),0,1)/((punif(1/(sigma*abs(nu)),0,1))^{1/1}))))^{1/1})))^(1/nu),
                      mu*exp(sigma*(qnorm(p)/((qunif(p))^{1/1})))))
  }
  return(xp)
}
