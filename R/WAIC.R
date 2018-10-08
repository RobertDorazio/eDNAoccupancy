#' Widely Applicable Information Criterion (WAIC) of a Multi-scale Occupancy Model
#'
#' Computes the value of WAIC for a fitted multi-scale occupancy model.
#'
#' @param fit  object of class occModel that contains data and previous state of the model's Markov chain
#' @param burnin  initial no. iterations of Markov chain to be omitted from calculations
#' 
#' @importFrom utils read.csv
#' @export
#'
#'
#' @return value of WAIC and its goodness-of-fit and predictive-variance components.
#'
#' @details  This function computes the WAIC value used in model-selection once a multi-scale occupancy model has been fitted using  \code{occModel}.
#'
#'
#' @examples
#'
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
#'                 niter                = 1100,
#'                 niterInterval        = 100,
#'                 siteColName = 'site',
#'                 )
#' 
#' WAIC(fit1, burnin=100)




WAIC = function(fit, burnin = 1) {

    ## make sure fit is a occModel object
    if(!is.null(fit)){
        if(!inherits(fit, "occModel")){
            stop(paste(fit, "is not an occModel object"))
        }
    }
    
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
    
    


    ## Retrieve matrix of observed counts and vectorize non-missing values
    y = fit$y
    K = fit$K
    jind = !is.na(y)  # matrix of boolean indices for non-missing observations of y
    yvec = as.vector(y[jind])
    Kvec = as.vector(K[jind])


    ## Read in vector of predictions of counts for each draw of Markov chain
    outSave = -(1:burnin)
    mc.ymean = as.matrix(read.table("mc.ypredmean.txt", header=TRUE))
    mc.ymean = mc.ymean[ outSave, ]


    ## Estimate expectations for criterion
    nobs = sum(jind)
    logOfProb = matrix(nrow=nrow(mc.ymean), ncol=nobs)
    for (i in 1:nrow(logOfProb)) {
        logOfProb[i,] = dbinom(yvec, size=Kvec, prob=mc.ymean[i,], log=TRUE)
    }
    firstExpectation = colMeans(exp(logOfProb))
    secondExpectation = colMeans(logOfProb*logOfProb)
    thirdExpectation = colMeans(logOfProb)

    ## Compute WAIC
    lof = (-1/nobs) * sum( log(firstExpectation) )
    predVariance =  (1/nobs) * sum( secondExpectation - thirdExpectation^2 )

    list(criterion=lof+predVariance,  lackOfFit=lof,  predVariance=predVariance)
}
