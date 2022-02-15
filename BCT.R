#'@title Box-Cox t (BCT)
#'
#'@description The Box-Cox t distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox t distribution.
#'
#'@usage dBCT(\eqn{x, \mu, \sigma, \nu, \tau})
#'
#'@usage pBCT(\eqn{q, \mu, \sigma, \nu, \tau})
#'
#'@usage qBCT(\eqn{p, \mu, \sigma, \nu, \tau})
#'
#'@usage rBCT(\eqn{n, \mu, \sigma, \nu, \tau})
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
#'@return dBCT() provides the probability density function, pBCT() provides the cumulative distribution function, qBCT() provides the quantiles, and rBCT() generates random numbers.
#'
#'@section Warning: The BCT distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, \tau > 0, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the t-Student distribution is:
#'
#'@note \deqn{r(u)=\tau^\frac{\tau}{2}(B(\frac{1}{2}, \frac{\tau}{2}))^{-1}(\tau+u)^{-\frac{\tau+1}{2}}}, \eqn{\tau > 0}, where B(.,.) is the beta function.
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references FUMES-GHANTOUS, G.; FERRARI, S. L. P.; CORRENTE, J. E. Box-Cox t random intercept model for estimating usual nutrient intake distributions. \eqn{Statistical Methods and Applications}, v. 27, p. 715-734, 2018.
#'
#'@references RIGBY, R. A.; STASINOPOULOS, D. M. Using the Box-Cox t distribution in GAMLSS to model skewness and kurtosis. \eqn{Statistical Modelling}, v. 6, n.3, p. 209-229, 2006.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@examples -
#'
#'rBCT(1,1,0.5,0.5,4)
#'
#'rBCT(100,5,0.5,0.5,4)
#'
#'@export rBCT

#random numbers generation

rBCT<- function(n,mu,sigma,nu,tau){
  if(n<=0|n!=as.integer(n)) stop(paste("n must be positive integer"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  if(any(tau<=0)) stop(paste("tau must be positive"))
  if (nu!=0){
    z=rt(n,tau)
    x=mu*(1+sigma*nu*z)^(1/nu)
  }
  if (nu==0){
    z=rt(n,tau)
    x=mu*exp(sigma*z)
  }
  return(x)
}


#'@title Box-Cox t (BCT)
#'
#'@description The Box-Cox t distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox t distribution.
#'
#'@usage dBCT(\eqn{x, \mu, \sigma, \nu, \tau})
#'
#'@usage pBCT(\eqn{q, \mu, \sigma, \nu, \tau})
#'
#'@usage qBCT(\eqn{p, \mu, \sigma, \nu, \tau})
#'
#'@usage rBCT(\eqn{n, \mu, \sigma, \nu, \tau})
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
#'@return dBCT() provides the probability density function, pBCT() provides the cumulative distribution function, qBCT() provides the quantiles, and rBCT() generates random numbers.
#'
#'@section Warning: The BCT distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, \tau > 0, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the t-Student distribution is:
#'
#'@note \deqn{r(u)=\tau^\frac{\tau}{2}(B(\frac{1}{2}, \frac{\tau}{2}))^{-1}(\tau+u)^{-\frac{\tau+1}{2}}}, \eqn{\tau > 0}, where B(.,.) is the beta function.
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references FUMES-GHANTOUS, G.; FERRARI, S. L. P.; CORRENTE, J. E. Box-Cox t random intercept model for estimating usual nutrient intake distributions. \eqn{Statistical Methods and Applications}, v. 27, p. 715-734, 2018.
#'
#'@references RIGBY, R. A.; STASINOPOULOS, D. M. Using the Box-Cox t distribution in GAMLSS to model skewness and kurtosis. \eqn{Statistical Modelling}, v. 6, n.3, p. 209-229, 2006.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@examples -
#'
#'dBCT(5,5,0.5,0.5,4)
#'
#'dBCT(5,5,0.5,0,4)
#'
#'@export dBCT

#Probability density function

dBCT<- function(x,mu,sigma,nu,tau){
  n<-max(length(x),length(mu),length(sigma),length(nu),length(tau))
  fx<-rep(NA,n)
  if(any(x<=0)) stop(paste("x must be positive"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  if(any(tau<=0)) stop(paste("tau must be positive"))
  for (i in 1:n){
    fx<-ifelse(x>0 & mu>0 & sigma>0 & nu>0 & tau>0,((x^{nu-1}/(sigma*mu^nu))*dt((1/(sigma*nu))*(((x/mu)^nu)-1),tau))/(1-pt(-1/(sigma*nu),tau)),
        ifelse(x>0 & mu>0 & sigma>0 & nu<0 & tau>0,((x^{nu-1}/(sigma*mu^nu))*dt((1/(sigma*nu))*(((x/mu)^nu)-1),tau))/(pt(-1/(sigma*nu),tau)),
        (1/(sigma*x))*dt((1/sigma)*log(x/mu),tau)))
  }
  return(fx)
}

#'@title Box-Cox t (BCT)
#'
#'@description The Box-Cox t distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox t distribution.
#'
#'@usage dBCT(\eqn{x, \mu, \sigma, \nu, \tau})
#'
#'@usage pBCT(\eqn{q, \mu, \sigma, \nu, \tau})
#'
#'@usage qBCT(\eqn{p, \mu, \sigma, \nu, \tau})
#'
#'@usage rBCT(\eqn{n, \mu, \sigma, \nu, \tau})
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
#'@return dBCT() provides the probability density function, pBCT() provides the cumulative distribution function, qBCT() provides the quantiles, and rBCT() generates random numbers.
#'
#'@section Warning: The BCT distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, \tau > 0, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the t-Student distribution is:
#'
#'@note \deqn{r(u)=\tau^\frac{\tau}{2}(B(\frac{1}{2}, \frac{\tau}{2}))^{-1}(\tau+u)^{-\frac{\tau+1}{2}}}, \eqn{\tau > 0}, where B(.,.) is the beta function.
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references FUMES-GHANTOUS, G.; FERRARI, S. L. P.; CORRENTE, J. E. Box-Cox t random intercept model for estimating usual nutrient intake distributions. \eqn{Statistical Methods and Applications}, v. 27, p. 715-734, 2018.
#'
#'@references RIGBY, R. A.; STASINOPOULOS, D. M. Using the Box-Cox t distribution in GAMLSS to model skewness and kurtosis. \eqn{Statistical Modelling}, v. 6, n.3, p. 209-229, 2006.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@examples -
#'
#'pBCT(0.5,3,0.5,0.5,4)
#'
#'pBCT(0.5,3,0.5,0,4)
#'
#'@export pBCT

#Distribution function

pBCT<- function(x,mu,sigma,nu,tau){
  n<-max(length(x),length(mu),length(sigma),length(nu), length(tau))
  Fx<-rep(NA,n)
  if(any(x<=0)) stop(paste("x must be positive"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  if(any(tau<=0)) stop(paste("tau must be positive"))
  for (i in 1:n){
    Fx<-ifelse(x>0 & mu>0 & sigma>0 & nu>0 & tau>0,(pt(((1/(sigma*nu))*(((x/mu)^nu)-1)),tau)-pt(-1/(sigma*nu),tau))/(1-pt(-1/(sigma*nu),tau)),
        ifelse(x>0 & mu>0 & sigma>0 & nu<0 & tau>0, pt(((1/(sigma*nu))*(((x/mu)^nu)-1)),tau)/pt(-1/(sigma*nu),tau),
        pt(((1/sigma)*log(x/mu)),tau)))
  }
  return(Fx)
}

#'@title Box-Cox t (BCT)
#'
#'@description The Box-Cox t distribution belongs to Box-Cox symmetric class. The function returns random numbers generation from the density function, the probability density function, the cumulative distribution function and the quantile for the specific parameterization of the Box-Cox t distribution.
#'
#'@usage dBCT(\eqn{x, \mu, \sigma, \nu, \tau})
#'
#'@usage pBCT(\eqn{q, \mu, \sigma, \nu, \tau})
#'
#'@usage qBCT(\eqn{p, \mu, \sigma, \nu, \tau})
#'
#'@usage rBCT(\eqn{n, \mu, \sigma, \nu, \tau})
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
#'@details The probability density function of the Box-Cox symmetric distribution is given by:
#'
#'@details \deqn{fx(x)=\frac{x^{\nu-1}}{\mu^{\nu}\sigma}\frac{r((\frac{1}{\sigma\nu}((\frac{x}{\mu})^\nu-1))^2)}{R(\frac{1}{\sigma|\nu|})}}, if \eqn{\nu != 0}
#'
#'@details \deqn{fx(x)=\frac{1}{x\sigma}r((\frac{1}{\sigma}log(\frac{x}{\mu}))^2)}, if \eqn{\nu = 0}
#'
#'@details for *x > 0*, where \eqn{R(w)=\int_{-\infty}^{w}r(u^2)du} for w that belongs to the real numbers
#'
#'@return dBCT() provides the probability density function, pBCT() provides the cumulative distribution function, qBCT() provides the quantiles, and rBCT() generates random numbers.
#'
#'@section Warning: The BCT distribution is suitable for all combinations of the parameters within their ranges
#'
#'@section Warning: (i.e. \eqn{\mu > 0, \sigma > 0, -\infty < \nu < +\infty, \tau > 0, n > 1, x > 0 and  0 < p < 1})
#'
#'@note The density generating function for the t-Student distribution is:
#'
#'@note \deqn{r(u)=\tau^\frac{\tau}{2}(B(\frac{1}{2}, \frac{\tau}{2}))^{-1}(\tau+u)^{-\frac{\tau+1}{2}}}, \eqn{\tau > 0}, where B(.,.) is the beta function.
#'
#'@author(s) FUMES-GHANTOUS, G.; CORRENTE, J. E.; TATIS, A. F. G. S. L.; VASQUEZ, R. A. M.
#'
#'@references FERRARI, S. L. P.; FUMES, G. Box-Cox symmetric distributions and applications to nutritional data. \eqn{AStA-Advances in Statistical Analysis}, v. 101, p. 321-344, 2017.
#'
#'@references FUMES-GHANTOUS, G.; CORRENTE, J. E. Basic Functions for Computational Implementation of the Box-Cox Symmetric Class of Distributions. \eqn{OPEN JOURNAL OF STATISTICS}, v. 11, p. 1010-1016, 2021.
#'
#'@references FUMES-GHANTOUS, G.; FERRARI, S. L. P.; CORRENTE, J. E. Box-Cox t random intercept model for estimating usual nutrient intake distributions. \eqn{Statistical Methods and Applications}, v. 27, p. 715-734, 2018.
#'
#'@references RIGBY, R. A.; STASINOPOULOS, D. M. Using the Box-Cox t distribution in GAMLSS to model skewness and kurtosis. \eqn{Statistical Modelling}, v. 6, n.3, p. 209-229, 2006.
#'
#'@references VANEGAS, L. H.; PAULA, G. Log-symmetric distributions: Statistical properties and parameter estimation. \eqn{Brazilian Journal of Probability and Statistics}, v. 30, n. 2, p. 196-220, 2016.
#'
#'@examples -
#'
#'qBCT(0.5,5,0.5,0.5,4)
#'
#'qBCT(0.5,5,0.5,0,4)
#'
#'@export qBCT

#Quantile

qBCT<- function(p,mu,sigma,nu,tau){
  n<-max(length(p),length(mu),length(sigma),length(nu),length(tau))
  xp<-rep(NA,n)
  if(any(p<0|p>1)) stop(paste("p must be between 0 and 1"))
  if(any(mu<=0)) stop(paste("mu must be positive"))
  if(any(sigma<=0)) stop(paste("sigma must be positive"))
  if(any(tau<=0)) stop(paste("tau must be positive"))
  for (i in 1:n){
    xp<-ifelse(p>=0 & p<=1 & mu>0 & sigma>0 & nu>0 & tau>0,mu*(1+sigma*nu*(qt(1-((1-p)*(pt(1/(sigma*abs(nu)),tau))),tau)))^(1/nu),
               ifelse(p>=0 & p<=1 & mu>0 & sigma>0 & nu<0 & tau>0,mu*(1+sigma*nu*(qt(p*(pt(1/(sigma*abs(nu)),tau)),tau)))^(1/nu),
                      mu*exp(sigma*(qt(p,tau)))))
  }
  return(xp)
}
