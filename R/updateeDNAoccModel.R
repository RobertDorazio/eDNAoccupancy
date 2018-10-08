#' Update a Multi-scale Occupancy Model
#'
#' This function computes additional iterations of a Markov chain Monte Carlo algorithm that was used to fit a multi-scale occupancy model to a data set.
#'
#' 
#' @param fit  object of class occModel that contains data and previous state of the model's Markov chain
#' @param niter no. iterations of MCMC algorithm
#' @param niterInterval no. iterations for reporting progress of MCMC algorithm
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @import stats
#' @export
#'
#' 
#'
#' @return A list (object of class occModel) containing the following objects:
#' \describe{
#' \item{niterations:}{total number of iterations MCMC algorithm was run}
#' \item{state:}{the value of each model parameter after the last iteration of the MCMC algorithm}
#' \item{y:}{  M x J matrix of numbers of detections per sample}
#' \item{K:}{  M x J matrix of numbers of replicates per sample}
#' \item{X:}{ matrix of regressors associated with model parameter beta}
#' \item{W:}{ array of regressors associated with model parameter alpha}
#' \item{V:}{ array of regressors associated with model parameter delta}
#' \item{siteEffectInW:}{ logical indicator of whether model parameter alpha contains M elements (that is, one element per site)}
#' \item{colNamesOfX:}{ vector of names of regressors in X}
#' \item{colNamesOfW:}{ vector of names of regressors in W}
#' \item{colNamesOfV:}{ vector of names of regressors in V}
#' }
#' 
#' where M = no. sites and J = maximum no. samples per site
#'
#'
#'
#' @details  This function is used to update a multi-scale occupancy model fitted using \code{occModel}.  The model is fitted using Bayesian methods (specifically, a Markov chain Monte Carlo (MCMC) algorithm that is run for a finite number of iterations).  Output from the MCMC algorithm is stored in the file "mc.csv" and can be summarized using functions \code{plotTrace}, \code{plotACF}, and \code{posteriorSummary}.
#'
#' @seealso  \code{\link{occModel}}, \code{\link{scaleData}}
#'
#'
#'
#' @examples
#' data(gobyDetectionData)
#' detections = occData(gobyDetectionData, "site", "sample")
#' data(gobySurveyData)
#' gobySurveyData = scaleData(gobySurveyData)  # center and scale numeric covariates
#' 
#' fit1 = occModel(formulaSite          = ~ veg,
#'                 formulaSiteAndSample = ~ sal + twg,
#'                 formulaReplicate     = ~ sal + fish,
#'                 detectionMats        = detections,
#'                 siteData             = gobySurveyData,
#'                 niter                = 110,
#'                 niterInterval        = 10,
#'                 siteColName = 'site',
#'                 )
#'
#' # Update the Markov chain of the model specified in fit1
#' fit2 = updateOccModel(fit1, niter=50, niterInterval=10)
#' 




updateOccModel = function(fit, niter, niterInterval) {

    
    ## make sure fit is a occModel object
    if(!is.null(fit)){
        if(!inherits(fit, "occModel")){
            stop(paste(fit, "is not an occModel object"))
        }
    }

    y = fit$y
    K = fit$K
    X = fit$X
    W = fit$W
    V = fit$V
    siteEffectInW = fit$siteEffectInW
    colNamesOfX = fit$colNamesOfX
    colNamesOfW = fit$colNamesOfW
    colNamesOfV = fit$colNamesOfV

    ## make sure column names of file "mc.csv" correspond to those specified for fit
    beta.names = paste('beta', fit$colNamesOfX, sep='.')
    alpha.names = paste('alpha', fit$colNamesOfW, sep='.')
    delta.names = paste('delta', fit$colNamesOfV, sep='.')
    mc.names = c(beta.names, alpha.names, delta.names)
    ## .... remove parentheses from mc.names
    mc.names = gsub(pattern='(', replacement='.', mc.names, fixed=TRUE)
    mc.names = gsub(pattern=')', replacement='.', mc.names, fixed=TRUE)
    mcColumnNames = dimnames(read.csv('mc.csv'))[[2]]
    if (length(mc.names) != length(mcColumnNames)) {
        stop(paste("Column names in file 'mc.csv' do not match the model matrices of the occModel object"))
    }
    if (any(mc.names != mcColumnNames)) {
        stop(paste("Column names in file 'mc.csv' do not match the model matrices of the occModel object"))
    }
    


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

logDensity.beta = function(beta, z, X, mu.beta, sigma.beta) {
    logPrior = sum(dnorm(beta, mean=mu.beta, sd=sigma.beta, log=TRUE))
    psi = as.vector(pnorm(X %*% beta))
    ind = psi>0 & psi<1
    logDensity = sum(z[ind]*log(psi[ind]) + (1-z[ind])*log(1-psi[ind]))
    logDensity + logPrior
}

logDensity.alpha = function(alpha, avec, Wmat, mu.alpha, sigma.alpha) {
    logPrior = sum(dnorm(alpha, mean=mu.alpha, sd=sigma.alpha, log=TRUE))
    thetavec = as.vector(pnorm(Wmat %*% alpha))
    ind = thetavec>0 & thetavec<1
    logDensity = sum(avec[ind]*log(thetavec[ind]) + (1-avec[ind])*log(1-thetavec[ind]))
    logDensity + logPrior
}

logDensity.delta = function(delta, yvec, Kvec, Vmat, mu.delta, sigma.delta) {
    logPrior = sum(dnorm(delta, mean=mu.delta, sd=sigma.delta, log=TRUE))
    pvec = as.vector(pnorm(Vmat %*% delta))
    ind = pvec>0 & pvec<1
    logDensity = sum(dbinom(yvec[ind], size=Kvec[ind], prob=pvec[ind], log=TRUE))
    logDensity + logPrior
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
z = fit$state$z
a = fit$state$a
beta = fit$state$beta
alpha = fit$state$alpha
delta = fit$state$delta
    


## .... compute MCMC draws

continueGibbs = TRUE
draw = 0
ndraws = niter
CR = '\n'
cat('Begin MCMC sampling:', CR, CR)



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
        if (all.equal(probMat[i, zeroInd], rep(1, sum(zeroInd))) == TRUE) {
            a[i, zeroInd] <- 1   # this line corrects a bug in R
        }
        else {
            a[i, zeroInd] = rbinom(sum(zeroInd), size=1, prob=probMat[i, zeroInd])
        }
        positiveInd = y[i,]>0 & jind[i, ]
        a[i, positiveInd] = 1
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
  logQ.beta.cand = dmvnorm(beta.cand, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
  logQ.beta = dmvnorm(beta, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
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
  ind = z==1
  amat = a[ind, ]
  Warr = array(W[ind, ,], dim=c(sum(ind), J, dim(W)[3]))
  aind = jind[ind, ]
  avec = as.vector(amat[aind[,1], 1])
  Wmat = matrix(Warr[aind[,1] ,1, ], ncol=dim(W)[3])
  for (j in 2:J) {
    avec = c(avec, as.vector(amat[aind[,j], j]))
    Wmat = rbind(Wmat, matrix(Warr[aind[,j] ,j, ], ncol=dim(W)[3]))
  }

   ## ... find mode and hessian of unnormalized conditional density function
  fit = glm(avec~Wmat-1, family=binomial(link='probit'))
  if (fit$converged & !fit$boundary) {
  
  ## ... draw candidate using multivariate normal distribution as proposal
  meanOfProposal = fit$coefficients
  SigmaOfProposal = vcov(fit)
  alpha.cand = as.vector(rmvnorm(1, mean=meanOfProposal, sigma=SigmaOfProposal))

  logL.alpha.cand = logDensity.alpha(alpha.cand, avec, Wmat, mu.alpha, sigma.alpha)
  logL.alpha = logDensity.alpha(alpha, avec, Wmat, mu.alpha, sigma.alpha)
  logQ.alpha.cand = dmvnorm(alpha.cand, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
  logQ.alpha = dmvnorm(alpha, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
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
    aind = a==1 & jind
    yvec = as.vector(y[aind[,1], 1])
    Kvec = as.vector(K[aind[,1], 1])
    Vmat = matrix(V[aind[,1] ,1, ], ncol=dim(V)[3])
    for (j in 2:J) {
        yvec = c(yvec, as.vector(y[aind[,j], j]))
        Kvec = c(Kvec, as.vector(K[aind[,j], j]))
        Vmat = rbind(Vmat, matrix(V[aind[,j] ,j, ], ncol=dim(V)[3]))
        }
    
    ## ... find mode and hessian of unnormalized conditional density function
    fit = glm(cbind(yvec,Kvec-yvec)~Vmat-1, family=binomial(link='probit'))
    if (fit$converged & !fit$boundary) {

    ## ... draw candidate using multivariate normal distribution as proposal
    meanOfProposal = fit$coefficients
    SigmaOfProposal = vcov(fit)
    delta.cand = as.vector(rmvnorm(1, mean=meanOfProposal, sigma=SigmaOfProposal))

    logL.delta.cand = logDensity.delta(delta.cand, yvec, Kvec, Vmat, mu.delta, sigma.delta)
    logL.delta = logDensity.delta(delta, yvec, Kvec, Vmat, mu.delta, sigma.delta)
    logQ.delta.cand = dmvnorm(delta.cand, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
    logQ.delta = dmvnorm(delta, mean=meanOfProposal, sigma=SigmaOfProposal, log=TRUE)
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
    niterations = nrow(read.csv('mc.csv'))
    state = list(beta=beta, alpha=alpha, delta=delta, z=z, a=a)
    dimnames(y)[[1]] = siteNames
    dimnames(K)[[1]] = siteNames
    updateToFit = list(niterations=niterations, state=state, y=y, K=K, X=X, W=W, V=V, siteEffectInW=siteEffectInW,
         colNamesOfX=colNamesOfX,
         colNamesOfW=colNamesOfW,
         colNamesOfV=colNamesOfV)


    class(updateToFit) <- "occModel"
    return(updateToFit)
    
}  # end of fitting model

