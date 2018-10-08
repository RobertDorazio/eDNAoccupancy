#' Fit a Multi-scale Occupancy Model
#'
#' This function is called by \code{occModel} -- rather than directly by the user -- and fits a multi-scale occupancy model using Bayesian methods (specifically, a Markov chain Monte Carlo (MCMC) algorithm).
#'
#' 
#' @param niter no. iterations of MCMC algorithm
#' @param niterInterval no. iterations for reporting progress of MCMC algorithm
#' @param y    M x J matrix of numbers of detections per sample
#' @param K    M x J matrix of numbers of replicates per sample
#' @param X  matrix of regressors associated with model parameter beta
#' @param W  array of regressors associated with model parameter alpha
#' @param V  array of regressors associated with model parameter delta
#' @param siteEffectInW   logical indicator of whether model parameter alpha contains M elements (that is, one element per site)
#' @param colNamesOfX  vector of names of regressors in X
#' @param colNamesOfW  vector of names of regressors in W
#' @param colNamesOfV  vector of names of regressors in V
#'
#'
#' 
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @import stats
#' @keywords internal
## #' @export



fitOccModel = function(niter=100, niterInterval=10, y, K, X, W, V, siteEffectInW,
                           colNamesOfX, colNamesOfW, colNamesOfV) {

M = nrow(y)
J = ncol(y)
jind = !is.na(y)  # matrix of boolean indices for non-missing observations of y
Jvec = rowSums(jind)
negInf = -1E10  # arbitrarily large negative number used to avoid assignment of qnorm(0)=-Inf
posInf =  1E10  # arbitrarily large positive number used to avoid assignment of qnorm(1)=+Inf

    ## .... temporary assignments needed for error checking
    siteNames = dimnames(y)[[1]]
    dimnames(y) = NULL
    dimnames(K) = NULL
    


## .... define functions used in Metropolis-Hastings sampling

    logDensity.beta = function(beta, zvec, X, mu, sigma) {
        ##  unnormalized full conditional pdf for a vector of binary observations (zvec) with success probabilities (prob)
        ##  prob is computed using Gaussian cdf of X * beta
        logPrior = sum(dnorm(beta, mean=mu, sd=sigma, log=TRUE))
        prob = as.vector(pnorm(X %*% beta))
        ind = (zvec==1)  # index for successes
        logDensity = sum(log(prob[ind])) + sum(log(1-prob[!ind]))
        logDensity + logPrior
    }

    logDensity.delta = function(beta, yvec, Kvec, X, mu, sigma) {
        ##  unnormalized full conditional pdf for a vector of binomial observations (yvec) with success probabilities (prob)
        ##  prob is computed using Gaussian cdf of X * beta
        logPrior = sum(dnorm(beta, mean=mu, sd=sigma, log=TRUE))
        prob = as.vector(pnorm(X %*% beta))
        logDensity = sum(dbinom(yvec, size=Kvec, prob=prob, log=TRUE))
        logDensity + logPrior
    }
    
    ##  ... dMvn() is used internally in the bayesm R package.   It is much faster than dmvnorm() in the mvtnorm R package.
    dMvn <- function(X,mu,Sigma) {
        k <- ncol(X)
        rooti <- backsolve(chol(Sigma),diag(k))
        quads <- colSums((crossprod(rooti,(t(X)-mu)))^2)
        return(exp(-(k/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads))
    }



## Begin MCMC    

## ... assign prior parameter values
mu.beta = 0
sigma.beta = 1
mu.alpha = 0
sigma.alpha = 1
mu.delta = 0
sigma.delta = 1
    
a.psi = 1
b.psi = 1
a.theta = 1
b.theta = 1
a.p = 1
b.p = 1


## .... initialize Markov chain
## ...... initialize vector z; then initialize beta
ysum = rowSums(y, na.rm=TRUE)
z = as.integer(ysum>0)
XtransposeX = t(X) %*% X # .... check to ensure model is not overparameterized given data
if (any(eigen(XtransposeX, only.values=TRUE)$values<=1E-4)) {
    stop('Model of psi is overparameterized for your data')
}
fit.glm = glm(z~X-1, family=binomial(link='probit'))
beta = fit.glm$coefficients

## ...... initialize matrix a; then initialize alpha
a = matrix(as.integer(y>0), nrow=M)
alpha = rep(NA, dim(W)[3])
if (siteEffectInW & dim(W)[3]==M) {  # check for model of theta with site-level fixed parameters
        a.sum = rowSums(a, na.rm=TRUE)
        thetaVec = a.sum/Jvec
        alpha = qnorm(thetaVec)
        alpha[thetaVec==0] = negInf
        alpha[thetaVec==1] = posInf
}
else {
    ## ...... arrange elements of a and w vectors of occupied sites in column-major order; then initialize alpha

    zmat = matrix(rep(z,J), ncol=J)
    zvec = as.vector(zmat[jind])
    avec = as.vector(a[jind])
    wtemp = W[,,1]
    Wmat = matrix(as.vector(wtemp[jind]))
    if (dim(W)[3] > 1) {
        for (icov in 2:dim(W)[3]) {   
            wtemp = W[,,icov]
            Wmat = cbind(Wmat, as.vector(wtemp[jind]))
        }
    }

    zind = zvec==1
    avec = avec[zind]
    Wmat = matrix(Wmat[zind, ], ncol=dim(W)[3])                                    
    XtransposeX = t(Wmat) %*% Wmat # .... check to ensure model is not overparameterized given data
    if (any(eigen(XtransposeX, only.values=TRUE)$values<=1E-4)) {
        stop('Model of theta is overparameterized for your data')
    }
    fit.glm = glm(avec~Wmat-1, family=binomial(link='probit'))
    alpha = fit.glm$coefficients
}


    ## ......  arrange elements of y and v vectors of occupied samples in column-major order; then initialize delta
    yvec = as.vector(y[jind])
    Kvec = as.vector(K[jind])
    avec = as.vector(a[jind])
    vtemp = V[,,1]
    Vmat = matrix(as.vector(vtemp[jind]))
    if (dim(V)[3] > 1) {
        for (icov in 2:dim(V)[3]) {   
            vtemp = V[,,icov]
            Vmat = cbind(Vmat, as.vector(vtemp[jind]))
        }
    }
    
    aind = avec==1
    yvec = yvec[aind]
    Kvec = Kvec[aind]
    Vmat = matrix(Vmat[aind, ], ncol=dim(V)[3])
    XtransposeX = t(Vmat) %*% Vmat # .... check to ensure model is not overparameterized given data
    if (any(eigen(XtransposeX, only.values=TRUE)$values<=1E-4)) {
        stop('Model of p is overparameterized for your data')
    }
    fit.glm = glm(cbind(yvec,Kvec-yvec)~Vmat-1, family=binomial(link='probit'))
    delta = fit.glm$coefficients




## .... compute MCMC draws

continueGibbs = TRUE
draw = 0
ndraws = niter
CR = '\n'
cat('Begin MCMC sampling:', CR, CR)


beta.names = paste('beta', colNamesOfX, sep='.')
alpha.names = paste('alpha', colNamesOfW, sep='.')
delta.names = paste('delta', colNamesOfV, sep='.')

cat(c(beta.names, alpha.names, delta.names),  sep=',', file='mc.csv')
cat(CR, file='mc.csv', append=TRUE)

cat(as.character(1:M), sep=' ', file='mc.Z.txt')
cat(CR, file='mc.Z.txt', append=TRUE)

cat(as.character(1:sum(jind)), sep=' ', file='mc.ypred.txt')
    cat(CR, file='mc.ypred.txt', append=TRUE)
    
cat(as.character(1:sum(jind)), sep=' ', file='mc.ypredmean.txt')
cat(CR, file='mc.ypredmean.txt', append=TRUE)
    

start.time = Sys.time()

while(continueGibbs) {
  
  draw = draw + 1
  drawinterval = niterInterval
  if (draw == round(draw/drawinterval)*drawinterval)  {
      end.time = Sys.time()
      elapsed.time = difftime(end.time, start.time, units='mins')
      cat('..... drawing sample #', draw, ' after ', elapsed.time, ' minutes', CR)
  }


  ##  draw z
  psi = as.vector(pnorm(X %*% beta))
  theta = matrix(nrow=M, ncol=J)
  p = matrix(nrow=M, ncol=J)
  for (j in 1:J) {
    Wmat = matrix(W[,j,], ncol=dim(W)[3])
    theta[, j] = pnorm(Wmat %*% alpha)
    Vmat = matrix(V[,j,], ncol=dim(V)[3])
    p[, j] = pnorm(Vmat %*% delta)
    }
    asum = rowSums(a, na.rm=TRUE)
    for (i in 1:M) {
        if(asum[i]==0) {
            z.prob = ifelse(psi[i]==1, 1, psi[i] * prod(1-theta[i,jind[i,]]) / (psi[i] * prod(1-theta[i,jind[i,]])  + 1 - psi[i]))
            z[i] = rbinom(1, size=1, prob=z.prob)
        }
        else {
            z[i] = 1
        }
    }

    
    ## draw a
    probMat = z*theta*(1-p)^K / (z*theta*(1-p)^K + 1-z*theta)
    for (i in 1:M) {

        zeroInd = y[i,]==0 & jind[i, ]
        if (any(zeroInd)) {
            probVec = probMat[i, zeroInd]
            probVec[probVec>1] = 1   # corrects for possible round-off errors
            a[i, zeroInd] = rbinom(sum(zeroInd), size=1, prob=probVec)
        }
        
        positiveInd = y[i,]>0 & jind[i, ]
        if (any(positiveInd)) {
            a[i, positiveInd] = 1
        }
        
    }
    


    ## draw beta
    if (ncol(X)==1) {
        psi.scalar = rbeta(1, a.psi + sum(z), b.psi + M - sum(z))
        beta = qnorm(psi.scalar)
    }
    else {
  ## ... find mode and hessian of unnormalized conditional density function
    fit = glm(z~X-1, family=binomial(link='probit'))
    if (fit$converged & !fit$boundary) {
  
  ## ... draw candidate using multivariate normal distribution as proposal
  meanOfProposal = fit$coefficients
  SigmaOfProposal = vcov(fit)
  beta.cand = as.vector(rmvnorm(1, mean=meanOfProposal, sigma=SigmaOfProposal))

  logL.beta.cand = logDensity.beta(beta.cand, z, X, mu.beta, sigma.beta)
  logL.beta = logDensity.beta(beta, z, X, mu.beta, sigma.beta)
  ## logQ.beta.cand = dmvnorm(beta.cand, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
  ## logQ.beta = dmvnorm(beta, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
  logQ.beta.cand = log(dMvn(matrix(beta.cand,nrow=1), mu=meanOfProposal, Sigma=SigmaOfProposal))
  logQ.beta = log(dMvn(matrix(beta,nrow=1), mu=meanOfProposal, Sigma=SigmaOfProposal))
  logR = logL.beta.cand - logL.beta - logQ.beta.cand + logQ.beta
  if (logR >= 0 | runif(1,0,1) <= exp(logR)) {
      beta = beta.cand
  }
        }
  }

  
  
    ## draw alpha
    if (siteEffectInW & dim(W)[3]==M) { # check for model of theta with site-level fixed parameters
        a.sum = rowSums(a, na.rm=TRUE)
        thetaVec = rbeta(M, a.theta + z*a.sum, b.theta + z*(Jvec - a.sum))
        alpha = qnorm(thetaVec)
        alpha[thetaVec==0] = negInf
        alpha[thetaVec==1] = posInf
    }
    else {
        ##  ... extract y values and w vectors of occupied sites, and arrange them in column-major order
        
        zmat = matrix(rep(z,J), ncol=J)
        zvec = as.vector(zmat[jind])
        avec = as.vector(a[jind])
        wtemp = W[,,1]
        Wmat = matrix(as.vector(wtemp[jind]))
        if (dim(W)[3] > 1) {
            for (icov in 2:dim(W)[3]) {   
                wtemp = W[,,icov]
                Wmat = cbind(Wmat, as.vector(wtemp[jind]))
            }
        }

        zind = zvec==1
        avec = avec[zind]
        Wmat = matrix(Wmat[zind, ], ncol=dim(W)[3])
    
        ## ... find mode and hessian of unnormalized conditional density function
        fit = glm(avec~Wmat-1, family=binomial(link='probit'))
        if (fit$converged & !fit$boundary) {
  
            ## ... draw candidate using multivariate normal distribution as proposal
            meanOfProposal = fit$coefficients
            SigmaOfProposal = vcov(fit)
            if (length(meanOfProposal) != nrow(SigmaOfProposal))  browser()
            alpha.cand = as.vector(rmvnorm(1, mean=meanOfProposal, sigma=SigmaOfProposal))

            logL.alpha.cand = logDensity.beta(alpha.cand, avec, Wmat, mu.alpha, sigma.alpha)
            logL.alpha = logDensity.beta(alpha, avec, Wmat, mu.alpha, sigma.alpha)
            ## logQ.alpha.cand = dmvnorm(alpha.cand, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
            ## logQ.alpha = dmvnorm(alpha, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
            logQ.alpha.cand = log(dMvn(matrix(alpha.cand,nrow=1), mu=meanOfProposal, Sigma=SigmaOfProposal))
            logQ.alpha = log(dMvn(matrix(alpha,nrow=1), mu=meanOfProposal, Sigma=SigmaOfProposal))
            logR = logL.alpha.cand - logL.alpha - logQ.alpha.cand + logQ.alpha
            if (logR >= 0 | runif(1,0,1) <= exp(logR)) {
                alpha = alpha.cand
            }
        }
    }


  
  ## draw delta
    if (dim(V)[3]==1) {
        p.scalar = rbeta(1, a.p + sum(a*y, na.rm=TRUE), b.p + sum(a*(K-y), na.rm=TRUE) )
        delta = qnorm(p.scalar)
    }
    else { 
        ##  ... extract y values and v vectors of occupied samples, and arrange them in column-major order
        yvec = as.vector(y[jind])
        Kvec = as.vector(K[jind])
        avec = as.vector(a[jind])
        vtemp = V[,,1]
        Vmat = matrix(as.vector(vtemp[jind]))
        for (icov in 2:dim(V)[3]) {   
            vtemp = V[,,icov]
            Vmat = cbind(Vmat, as.vector(vtemp[jind]))
        }
    
        aind = avec==1
        yvec = yvec[aind]
        Kvec = Kvec[aind]
        Vmat = matrix(Vmat[aind, ], ncol=dim(V)[3])
    
        ## ... find mode and hessian of unnormalized conditional density function
        fit = glm(cbind(yvec,Kvec-yvec)~Vmat-1, family=binomial(link='probit'))
        if (fit$converged & !fit$boundary) {

            ## ... draw candidate using multivariate normal distribution as proposal
            meanOfProposal = fit$coefficients
            SigmaOfProposal = vcov(fit)
            delta.cand = as.vector(rmvnorm(1, mean=meanOfProposal, sigma=SigmaOfProposal))

            logL.delta.cand = logDensity.delta(delta.cand, yvec, Kvec, Vmat, mu.delta, sigma.delta)
            logL.delta = logDensity.delta(delta, yvec, Kvec, Vmat, mu.delta, sigma.delta)
            ## logQ.delta.cand = dmvnorm(delta.cand, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
            ## logQ.delta = dmvnorm(delta, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
            logQ.delta.cand = log(dMvn(matrix(delta.cand,nrow=1), mu=meanOfProposal, Sigma=SigmaOfProposal))
            logQ.delta = log(dMvn(matrix(delta,nrow=1), mu=meanOfProposal, Sigma=SigmaOfProposal))
            logR = logL.delta.cand - logL.delta - logQ.delta.cand + logQ.delta
            if (logR >= 0 | runif(1,0,1) <= exp(logR)) {
                delta = delta.cand
            }
        }
    }
        
  
    ##  Compute predictions of counts (ypred) conditional on current values of parameters
    ##  ... (predictions are used to computed model-selection criteria)
    apMat = matrix(nrow=M, ncol=J)
    if (dim(V)[3]==1) {
        apMat = a * pnorm(delta)
    }
    else {
        for (i in 1:M) {
            for (j in 1:J) {
                apMat[i,j] = ifelse(jind[i,j],  a[i,j] * pnorm(sum(V[i,j, ]*delta)), NA)
            }
        }
    }
    ypred = rbinom(sum(jind), size=as.vector(K[jind]), prob=as.vector(apMat[jind]))
    ypredmean = as.vector(apMat[jind])

  
  
  ## write draw to file
  cat(c(beta, alpha, delta), sep=',', file='mc.csv', append=TRUE)
  cat(CR, file='mc.csv', append=TRUE)
  cat(z, sep=' ', file='mc.Z.txt', append=TRUE)
  cat(CR, file='mc.Z.txt', append=TRUE)
  cat(ypred, sep=' ', file='mc.ypred.txt', append=TRUE)
  cat(CR, file='mc.ypred.txt', append=TRUE)
  cat(ypredmean, sep=' ', file='mc.ypredmean.txt', append=TRUE)
  cat(CR, file='mc.ypredmean.txt', append=TRUE)

  if (draw == ndraws) {
    cat('Completed ', ndraws, ' draws of MCMC algorithm', CR)
    continueGibbs = FALSE
  }

}  # end of while loop

### Return current state of Markov chain and data used to fit the model
    state = list(beta=beta, alpha=alpha, delta=delta, z=z, a=a)
    dimnames(y)[[1]] = siteNames
    dimnames(K)[[1]] = siteNames
    list(niterations=niter, state=state, y=y, K=K, X=X, W=W, V=V, siteEffectInW=siteEffectInW,
         colNamesOfX=colNamesOfX,
         colNamesOfW=colNamesOfW,
         colNamesOfV=colNamesOfV)
    
}  # end of fitting model
