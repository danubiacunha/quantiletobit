#################################################################
# "On a log-symmetric quantile tobit model applied to female    #   
#  labor supply data"                                           # 
# Author: Danúbia R. Cunha                                      #
# Date: 19-June-2021                                            #
#################################################################
# Fit for the log-symmetric quantile tobit regression model     #
#################################################################

#################################################################
## Required packages
#################################################################

require(rootSolve)
require(MASS)
require(maxLik)
require(VGAM) 
require(EnvStats)
require(normalp)
require(cubature)
require(miscTools)
require(ssym)
library(quantreg)
require(ctqr)

#################################################################
## Function to estimate
#################################################################

tobitlogsymqreg.fit = function(x,
                          w,
                          y,
                          q = 0.5,
                          xi,
                          status,
                          cens,
                          logcens=TRUE,
                          link.mu = "identity",
                          link.phi = "log",
                          family = "Normal",
                          method = "BFGS") 
  
{
  n           = length(y)
  p           = NCOL(x)
  l           = NCOL(w)
  
  if(link.mu=="identity"){
    l.mu <- function(y) y
    l.mu.i <- function(y) y
    l1.mu <- function(mu) matrix(1,length(mu),1)
    l2.mu <- function(mu) matrix(0,length(mu),1)	
    attr(l1.mu,"link") <- "identity"}else {if(link.mu=="log"){
      l.mu <- function(y) log(y)
      l.mu.i <- function(y) exp(y)
      l1.mu <- function(mu) mu
      l2.mu <- function(mu) -1/mu^2	
      attr(l1.mu,"link") <- "logarithmic"}
      if(link.mu=="exp"){
        l.mu <- function(y) exp(y)
        l.mu.i <- function(y) log(y)
        l1.mu <- function(mu) 1/mu
        l2.mu <- function(mu) mu	
        attr(l1.mu,"link") <- "exponential"}
      if(link.mu=="recip"){
        l.mu <- function(y) 1/y
        l.mu.i <- function(y) 1/y
        l1.mu <- function(mu) -mu^2
        l2.mu <- function(mu) 2/mu^3	
        attr(l1.mu,"link") <- "reciprocal"}		
    }
  
  if(link.phi=="log"){
    l.phi <- function(y) log(y)
    l.phi.i <- function(y) exp(y)
    l1.phi <- function(phi) matrix(1,length(phi),1)
    l2.phi <- function(phi) -1/phi^2	
    attr(l1.phi,"link") <- "logarithmic"} else{ l.phi <- function(y) y
    l.phi.i <- function(y) y
    l1.phi <- function(phi) 1/phi
    l2.phi <- function(phi) matrix(0,length(phi),1)	
    attr(l1.phi,"link") <- "identity"}
  
  
  ystar   = log(y)
  
  
  if(family=="Normal"){
    xi <- 0
    cdf <- function(z) pnorm(z)
    pdf <- function(z) dnorm(z)
    qf  <- function(z) qnorm(z)
    xix <- 1
  }
  if(family=="Student"){
    if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
    nu <- xi[1]
    cdf <- function(z) pt(z,nu)
    pdf <- function(z) dt(z,nu)
    qf  <- function(z) qt(z,nu)
    xix <- nu/(nu-2)
    if(nu<=2) xix <- as.null(xix)
  }
 
  if (family == "Powerexp") {
    if (xi[1] <= -1 | xi[1] >= 1) 
      stop("the extra parameter must be within the interval (-1, 1)!!", 
           call. = FALSE)
    kk <- xi[1]
    pp <- 2/(kk + 1)
    sigmap <- (1 + kk)^((kk + 1)/2)
    cdf <- function(z) pnormp(z, mu = 0, sigmap = sigmap, 
                              p = pp)
    pdf <- function(z) dnormp(z, mu = 0, sigmap = sigmap, 
                              p = pp)
    qf <- function(z) qnormp(z , mu=0, sigmap= sigmap, 
                             p= pp)
    dzg <- function(z) -(1/2) * exp(-(1/2) * z^(1/(1 + kk))) * (1/(1 + kk)) * z^(1/(1 + kk) - 1)
    xix <- 2^(1 + kk) * gamma(3 * (1 + kk)/2)/gamma((1 + 
                                                       kk)/2)
  }
  if (family == "Sinh-normal") {
    if (xi[1] <= 0) 
      stop("the extra parameter must be positive!!", call. = FALSE)
    alpha <- xi[1]
    cdf <- function(z) pnorm(2 * sinh(z)/alpha)
    pdf <- function(z) (2 * cosh(sqrt(z^2)) * exp(-2 * sinh(sqrt(z^2)) * 
                                                    sinh(sqrt(z^2))/alpha^2)/(sqrt(2 * pi) * alpha))
    qf <- function(z) asinh(alpha * qnorm(z)/2)
    dzg <- function(z) sinh(sqrt(z)) * exp(-(2/alpha^2) * sinh(sqrt(z))^2) * (alpha^2 - 4 * cosh(sqrt(z))^2)/(2 * alpha^2 * sqrt(z))
    fgf <- function(z) dshn(z) * z^2
    xix <- 2 * integrate(fgf, 0, 20)$value
  }
 
  res <- try(ssym.l( ystar ~ x| w, family = family, xi = xi,  link.mu = link.mu, link.phi =link.phi))

  res1 <- rq( ystar ~ x -1, tau=q, method= "br")
  initvaluescoef <- res1$coefficients
  names(initvaluescoef) <- colnames(x)
  
  start = list(beta = initvaluescoef, kappa = kappasinit)
  
  
  if (is.list(start)) 
    start = do.call("c", start)
  
  
  if(logcens==TRUE){cens1 <- cens} else{
    cens1 <- log(cens) 
  }
  
  
  fr    = function(vp)
  {                                         
    betab   = vp[1:p]
    etahat  = as.vector(x%*%betab)
    Xi      = l.mu.i(etahat)
    kappa  = vp[-(1:p)] 
    phihat  = as.vector(w%*%kappa)
    phi_es  = l.phi.i(phihat)
    zp      = qf(q)
    z_es      <- ((ystar - Xi + sqrt(phi_es)*zp))/sqrt(phi_es)
    z_cen      <- ((cens1 - Xi + sqrt(phi_es)*zp))/sqrt(phi_es)
    loglik <- sum( (1-status)*log(cdf(z_cen)) + status*(log(pdf(z_es))-(1/2)*log(phi_es) )   )  
    return(loglik)
    
  } 
  
  
  if(method == "Nelder-Mead"){
    opt = maxNM(fn = fr, start = start, control = list(iterlim = 10000))
  }
  if(method == "BFGS"){
    opt = maxBFGS(fn = fr, start = start, control = list(iterlim = 10000))
  }
  if(method == "SANN"){
    opt = maxSANN(fn = fr, start = start, control = list(iterlim = 10000))
  }
  if(method == "CG"){
    opt = maxCG(fn = fr, start = start, control = list(iterlim = 10000))
  }
  
  
  log.lik.est = opt$maximum
  estimates   = opt$estimate
  
  AIC   = - 2 * log.lik.est + 2 * (p+l)
  AICc  = AIC + (2 * (p+l) * ((p+l) + 1)) / (n - (p+l) - 1)
  BIC   = - 2 * log.lik.est + log(n) * (p+l)
  
  beta  = as.vector(estimates[1:p])
  eta   = as.vector(x%*%beta)
  Xi    = l.mu.i(eta)
  Q     = exp(Xi)
  
  kappa = as.vector(estimates[-(1:p)])
  kappahat = as.vector(w%*%kappa)
  phi    = l.phi.i(kappahat)
  
  fisher = -opt$hessian
  se     = sqrt(diag(solve(fisher)))
  
  confidence.int <- cbind(estimates - qnorm(0.975)*se, estimates + qnorm(0.975)*se)
  
  hess = as.matrix(opt$hessian)
  
  zstatbeta  = beta / se[1:p]
  zstatkappa   = kappa / se[-(1:p)] 
  pvalorbeta = 2 * pnorm(abs(zstatbeta), lower.tail = F)
  pvalorkappa  = 2 * pnorm(abs(zstatkappa), lower.tail = F)
  
  names(beta)  = colnames(x)
  names(kappa)   = colnames(w)
  
  sebeta = se[1:p]
  sekappa  = se[-(1:p)]
  matcofbeta  = cbind(beta,sebeta,zstatbeta,pvalorbeta)
  matcofkappa   = cbind(kappa,sekappa,zstatkappa,pvalorkappa)

  
  tb5 <- miscTools::coefTable(beta,sebeta, df=(n-(p+l)))
  tb6 <- miscTools::coefTable(kappa,sekappa, df=(n-(p+l)))
  
  
  zp          <- qf(q)
  zpadro      <- (((ystar - Xi + sqrt(phi)*zp))/sqrt(phi))
  znotpadro  <- ((ystar - Xi + sqrt(phi)*zp))
  zhat        <- cdf(((ystar - Xi + sqrt(phi)*zp))/sqrt(phi))
 
  
  Shat        <- 1-zhat
  GCSresidual <- -log(1-zhat)
  RQresidual  <- qnorm(zhat)
  Residuals   <- y -  Q
  
  z            <- ((ystar - Xi + sqrt(phi)*zp))/sqrt(phi)
  score.beta   <- colSums(x * v(z) * z / sqrt(phi))
  score.phi    <-  colSums(w * (v(z) * z * (z - zp) - 1))/ 2
  
  RMSE <- sqrt(mean((y-Q)^2))
  
  MAE <- (mean(abs(y-Q)))

  
  cat("\n")
  cat("--------------------------------------------------------------\n")
  cat("       Quantile log-symmetric tobit Model                     \n")
  cat("--------------------------------------------------------------\n")
  cat("--------------------------------------------------------------\n")
  cat("Log-symmetric case:", family, "\n")
  cat("Maximum Likelihood estimation \n")
  cat("Log-Likelihood:", log.lik.est, "\n")
  cat("AIC:", AIC, "BIC:", BIC, "\n")
  cat("Number of observations:", n, "(", n-sum(status), "censored and", sum(status), "observed", ")", "\n")
  cat("Quantile of interest:", q, "\n")
  cat("Extra parameter:", xi, "\n")
  cat("--------------------------------------------------------------\n")
  cat("Q - Coefficients:\n")
  printCoefmat(tb5, signif.stars = TRUE, signif.legend = TRUE, digits = 4)
  cat("--------------------------------------------------------------\n")
  cat("Phi - Coefficients :\n")
  printCoefmat(tb6, signif.stars = TRUE, signif.legend = TRUE, digits = 4)
  cat("--------------------------------------------------------------\n")
  
  rval = list(matcofbeta = matcofbeta,
              matcofkappa = matcofkappa,
              coefficients = list(beta = beta, kappa = kappa),
              se = se, 
              conf.int = confidence.int,
              pvalor = list(beta = pvalorbeta, kappa = pvalorkappa),
              xi = xi,
              Initial.values = start, 
              converged = opt$message, 
              information.criterions = list(aic = AIC,bic = BIC,aicc=AICc), 
              loglik = log.lik.est,
              n = n, p = p, l = l, 
              score = list(beta = score.beta, kappa = score.phi),
              Hessian = hess,
              RMSE     = RMSE,
              MAE      = MAE,
              Zpadronizado   = zpadro,
              Znotpadronizado  = znotpadro,
              Shat     = Shat,
              phihat  = phi,
              Residuals = Residuals,
              GCSresidual = as.vector(GCSresidual),
              RQresidual = as.vector(RQresidual),
              fitted.values.Q = structure(Q, .Names = names(y)))
  
  
  return(rval)
  
  
}




