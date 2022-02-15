#'@title Box-Cox Power Exponential (BCPE)
#'
#'@description The Box-Cox Power Exponential distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Power Exponential distribution.
#'
#'@usage dBCPE(\eqn{x, \mu, \sigma, \nu, \tau})
#'
#'@usage pBCPE(\eqn{q, \mu, \sigma, \nu, \tau})
#'
#'@usage qBCPE(\eqn{p, \mu, \sigma, \nu, \tau})
#'
#'@usage rBCPE(\eqn{n, \mu, \sigma, \nu, \tau})
#'
#'@param \eqn{\mu} vector of scale parameter values (related to median)
#'
#'@param \eqn{\sigma} vector of relative dispersion parameter values (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} vector of skew parameter values (given by the power transformation to symmetry)
#'
#'@param \eqn{\tau} vector of kurtosis parameter values (given by the degrees of freedom)
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
#'@return dBCPE() provides the probability density function, pBCPE() provides the cumulative distribution function, qBCPE() provides the quantiles, and rBCPE() generates random numbers.
#'
#'@section Warning: The BCPE distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, \tau>0, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the BCPE distribution is:
#'
#'@note \deqn{r(u)=\frac{\tau\exp{\frac{\frac{-1}{2}u^{\frac{\tau}{2}}}{|p(\tau)|^\tau}}}{p(\tau)2^{1+\frac{1}{\tau}}\Gamma\left(\frac{1}{\tau}\right)}}, where \eqn{\tau > 0} and \eqn{p(\tau)^2=2^{\frac{-2}{\tau}}\Gamma\left(\frac{1}{\tau}\right)\left(\Gamma\left(\frac{3}{\tau}\right)\right)^{-1}}, when \eqn{\tau = 1} and \eqn{\tau = 2}, r(u) coincides with the density generating function of the double exponential and normal, respectively.
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references RIGBY, R. A.; STASINOPOULOS, D. M. Smooth centile curves for skew and kurtotic data modelled using the Box-Cox power exponential distribution. \eqn{Statistics in Medicine}, v. 23, n.19, p. 3053-3076, 2004.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@references VOUDOURIS, V. et al. Modelling skewness and kurtosis with the BCPE density in GAMLSS. \eqn{Journal of Applied Statistics}, v. 39, n.6, p. 1279-1293, 2012.
#'
#'@examples -
#'
#'rBCPE(1,5,0.5,0.5,4)
#'
#'rBCPE(100,5,0.5,0.5,4)
#'
#'@export rBCPE

#random numbers generation

rBCPE<- function(n,mu,sigma,nu,tau){
  if(n<=0 | n!=as.integer(n)) stop(paste("n must be positive integer"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  if(any(tau<=0)) stop(paste("tau must be positive"))
  library(gamlss)
  if (nu!=0){
    z=rPE(n,0,1,tau)
    x=mu*(1+sigma*nu*z)^(1/nu)
  }
  if (nu==0){
    z=rPE(n,0,1,tau)
    x=mu*exp(sigma*z)
  }
  return(x)
}

#'@title Box-Cox Power Exponential (BCPE)
#'
#'@description The Box-Cox Power Exponential distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Power Exponential distribution.
#'
#'@usage dBCPE(\eqn{x, \mu, \sigma, \nu, \tau})
#'
#'@usage pBCPE(\eqn{q, \mu, \sigma, \nu, \tau})
#'
#'@usage qBCPE(\eqn{p, \mu, \sigma, \nu, \tau})
#'
#'@usage rBCPE(\eqn{n, \mu, \sigma, \nu, \tau})
#'
#'@param \eqn{\mu} vector of scale parameter values (related to median)
#'
#'@param \eqn{\sigma} vector of relative dispersion parameter values (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} vector of skew parameter values (given by the power transformation to symmetry)
#'
#'@param \eqn{\tau} vector of kurtosis parameter values (given by the degrees of freedom)
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
#'@return dBCPE() provides the probability density function, pBCPE() provides the cumulative distribution function, qBCPE() provides the quantiles, and rBCPE() generates random numbers.
#'
#'@section Warning: The BCPE distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, \tau>0, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the BCPE distribution is:
#'
#'@note \deqn{r(u)=\frac{\tau\exp{\frac{\frac{-1}{2}u^{\frac{\tau}{2}}}{|p(\tau)|^\tau}}}{p(\tau)2^{1+\frac{1}{\tau}}\Gamma\left(\frac{1}{\tau}\right)}}, where \eqn{\tau > 0} and \eqn{p(\tau)^2=2^{\frac{-2}{\tau}}\Gamma\left(\frac{1}{\tau}\right)\left(\Gamma\left(\frac{3}{\tau}\right)\right)^{-1}}, when \eqn{\tau = 1} and \eqn{\tau = 2}, r(u) coincides with the density generating function of the double exponential and normal, respectively.
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references RIGBY, R. A.; STASINOPOULOS, D. M. Smooth centile curves for skew and kurtotic data modelled using the Box-Cox power exponential distribution. \eqn{Statistics in Medicine}, v. 23, n.19, p. 3053-3076, 2004.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@references VOUDOURIS, V. et al. Modelling skewness and kurtosis with the BCPE density in GAMLSS. \eqn{Journal of Applied Statistics}, v. 39, n.6, p. 1279-1293, 2012.
#'
#'@examples -
#'
#'dBCPE(4,5,0.5,0.5,4)
#'
#'dBCPE(4,5,0.5,0,4)
#'
#'@export dBCPE

#Probability density function


dBCPE<- function(x,mu,sigma,nu,tau){
  n<-max(length(x),length(mu),length(sigma),length(nu),
         length(tau))
  fx<-rep(NA,n)
  if(any(x<=0)) stop(paste("x must be positive"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  if(any(tau<=0)) stop(paste("tau must be positive"))
  library(gamlss)
  for (i in 1:n){
    fx<-ifelse(x>0 & mu>0 & sigma>0 & nu>0 & tau>0,((x^{nu-1}/(sigma*mu^nu))*dPE((1/(sigma*nu))*(((x/mu)^nu)-1),0,1,tau))/(1-pPE(-1/(sigma*nu),0,1,tau)),
               ifelse(x>0 & mu>0 & sigma>0 & nu<0 & tau>0,((x^{nu-1}/(sigma*mu^nu))*dPE((1/(sigma*nu))*(((x/mu)^nu)-1),0,1,tau))/(pPE(-1/(sigma*nu),0,1,tau)),
                      (1/(sigma*x))*dPE((1/sigma)*log(x/mu),0,1,tau)))
  }
  return(fx)
}

#'@title Box-Cox Power Exponential (BCPE)
#'
#'@description The Box-Cox Power Exponential distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Power Exponential distribution.
#'
#'@usage dBCPE(\eqn{x, \mu, \sigma, \nu, \tau})
#'
#'@usage pBCPE(\eqn{q, \mu, \sigma, \nu, \tau})
#'
#'@usage qBCPE(\eqn{p, \mu, \sigma, \nu, \tau})
#'
#'@usage rBCPE(\eqn{n, \mu, \sigma, \nu, \tau})
#'
#'@param \eqn{\mu} vector of scale parameter values (related to median)
#'
#'@param \eqn{\sigma} vector of relative dispersion parameter values (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} vector of skew parameter values (given by the power transformation to symmetry)
#'
#'@param \eqn{\tau} vector of kurtosis parameter values (given by the degrees of freedom)
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
#'@return dBCPE() provides the probability density function, pBCPE() provides the cumulative distribution function, qBCPE() provides the quantiles, and rBCPE() generates random numbers.
#'
#'@section Warning: The BCPE distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, \tau>0, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the BCPE distribution is:
#'
#'@note \deqn{r(u)=\frac{\tau\exp{\frac{\frac{-1}{2}u^{\frac{\tau}{2}}}{|p(\tau)|^\tau}}}{p(\tau)2^{1+\frac{1}{\tau}}\Gamma\left(\frac{1}{\tau}\right)}}, where \eqn{\tau > 0} and \eqn{p(\tau)^2=2^{\frac{-2}{\tau}}\Gamma\left(\frac{1}{\tau}\right)\left(\Gamma\left(\frac{3}{\tau}\right)\right)^{-1}}, when \eqn{\tau = 1} and \eqn{\tau = 2}, r(u) coincides with the density generating function of the double exponential and normal, respectively.
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references RIGBY, R. A.; STASINOPOULOS, D. M. Smooth centile curves for skew and kurtotic data modelled using the Box-Cox power exponential distribution. \eqn{Statistics in Medicine}, v. 23, n.19, p. 3053-3076, 2004.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@references VOUDOURIS, V. et al. Modelling skewness and kurtosis with the BCPE density in GAMLSS. \eqn{Journal of Applied Statistics}, v. 39, n.6, p. 1279-1293, 2012.
#'
#'@examples -
#'
#'pBCPE(5,5,0.5,0.5,4)
#'
#'pBCPE(5,5,0.5,0,4)
#'
#'@export pBCPE

#Distribution function

pBCPE<- function(q,mu,sigma,nu,tau){
  n<-max(length(q),length(mu),length(sigma),length(nu),
         length(tau))
  Fx<-rep(NA,n)
  if(any(q<=0)) stop(paste("q must be positive"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  if(any(tau<=0)) stop(paste("tau must be positive"))
  library(gamlss)
  for (i in 1:n){
    Fx<-ifelse(q>0 & mu>0 & sigma>0 & nu>0 & tau>0,(pPE((1/(sigma*nu))*(((q/mu)^nu)-1),0,1,tau)-pPE(-1/(sigma*nu),0,1,tau))/(1-pPE(-1/(sigma*nu),0,1,tau)),
        ifelse(q>0 & mu>0 & sigma>0 & nu<0 & tau>0, pPE((1/(sigma*nu))*(((q/mu)^nu)-1),0,1,tau)/pPE(-1/(sigma*nu),0,1,tau),
          pPE((1/sigma)*log(q/mu),0,1,tau)))
  }
  return(Fx)
}

#'@title Box-Cox Power Exponential (BCPE)
#'
#'@description The Box-Cox Power Exponential distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox Power Exponential distribution.
#'
#'@usage dBCPE(\eqn{x, \mu, \sigma, \nu, \tau})
#'
#'@usage pBCPE(\eqn{q, \mu, \sigma, \nu, \tau})
#'
#'@usage qBCPE(\eqn{p, \mu, \sigma, \nu, \tau})
#'
#'@usage rBCPE(\eqn{n, \mu, \sigma, \nu, \tau})
#'
#'@param \eqn{\mu} vector of scale parameter values (related to median)
#'
#'@param \eqn{\sigma} vector of relative dispersion parameter values (associated with the variation coefficient based in percentiles)
#'
#'@param \eqn{\nu} vector of skew parameter values (given by the power transformation to symmetry)
#'
#'@param \eqn{\tau} vector of kurtosis parameter values (given by the degrees of freedom)
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
#'@return dBCPE() provides the probability density function, pBCPE() provides the cumulative distribution function, qBCPE() provides the quantiles, and rBCPE() generates random numbers.
#'
#'@section Warning: The BCPE distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, \tau>0, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the BCPE distribution is:
#'
#'@note \deqn{r(u)=\frac{\tau\exp{\frac{\frac{-1}{2}u^{\frac{\tau}{2}}}{|p(\tau)|^\tau}}}{p(\tau)2^{1+\frac{1}{\tau}}\Gamma\left(\frac{1}{\tau}\right)}}, where \eqn{\tau > 0} and \eqn{p(\tau)^2=2^{\frac{-2}{\tau}}\Gamma\left(\frac{1}{\tau}\right)\left(\Gamma\left(\frac{3}{\tau}\right)\right)^{-1}}, when \eqn{\tau = 1} and \eqn{\tau = 2}, r(u) coincides with the density generating function of the double exponential and normal, respectively.
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references RIGBY, R. A.; STASINOPOULOS, D. M. Smooth centile curves for skew and kurtotic data modelled using the Box-Cox power exponential distribution. \eqn{Statistics in Medicine}, v. 23, n.19, p. 3053-3076, 2004.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@references VOUDOURIS, V. et al. Modelling skewness and kurtosis with the BCPE density in GAMLSS. \eqn{Journal of Applied Statistics}, v. 39, n.6, p. 1279-1293, 2012.
#'
#'@examples -
#'
#'qBCPE(0.5,5,0.5,0.5,4)
#'
#'qBCPE(0.5,5,0.5,0,4)
#'
#'@export qBCPE

#Quantile

qBCPE<- function(p,mu,sigma,nu,tau){
  n<-max(length(p),length(mu),length(sigma),length(nu),
         length(tau))
  xp<-rep(NA,n)
  if(any(p<0|p>1)) stop(paste("p must be between 0 and 1"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  if(any(tau<=0)) stop(paste("tau must be positive"))
  for (i in 1:n){
    xp<-ifelse(p>=0 & p<=1 & mu>0 & sigma>0 & nu>0 & tau>0,mu*((1+sigma*nu*(qPE(1-((1-p)*(pPE(1/(sigma*abs(nu)),0,1,tau))),0,1,tau)))^(1/nu)),
               ifelse(p>=0 & p<=1 & mu>0 & sigma>0 & nu<0, mu*(1+sigma*nu*(qPE(p*(pPE(1/(sigma*abs(nu)),0,1,tau)),0,1,tau)))^(1/nu),
                mu*(exp(sigma*(qPE(p,0,1,tau))))))
  }
  return(xp)
}
