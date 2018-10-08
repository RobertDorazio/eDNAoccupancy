#' Posterior Summary of a Multi-scale Occupancy Model's Sample-Occupancy Probabilities
#'
#' Estimates the posterior mean, median, and 95\% credible limits for a multi-scale occupancy model's sample-specific occurrence probabilities.
#'
#' 
#' @param fit  object of class occModel that contains data and previous state of the model's Markov chain
#' @inheritParams posteriorSummary
#' @importFrom utils read.csv
#' @export
#'
#' @return Computes estimates of summary statistics for a multi-scale occupancy model's sample-occupancy probabilities.  If \code{mcError}=TRUE, the Monte Carlo standard errors of these estimates are computed.   All posterior summaries are returned in a list.
#'
#' @details  This function estimates the posterior mean, median, and 95\% credible limits of the sample-specific occurrence probabilities of a multi-scale occupancy model.
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
#' theta = posteriorSummaryOfSampleOccupancy(fit1, burnin=100)
#' plot(gobySurveyData[, 'sal'], theta$median[,1])



posteriorSummaryOfSampleOccupancy <- function(
                                              fit,
                                              burnin = 1,
                                              mcError = FALSE
                                              ) {

    ## make sure fit is a occModel object
    if(!is.null(fit)){
        if(!inherits(fit, "occModel")){
            stop(paste(deparse(substitute(fit)), "is not an occModel object"))
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
    
    
        

    ## Compute theta vector for each draw of alpha in Markov chain
    W = fit$W
    colNamesOfW = fit$colNamesOfW
    alpha.names = make.names(paste('alpha', colNamesOfW, sep='.'))

    
    ## Read Markov chain from file
    mc = as.matrix(read.csv("mc.csv"))
    mc.alpha = as.matrix(mc[ , alpha.names])

    M =  dim(W)[1]
    J =  dim(W)[2]
    lenAlpha = dim(mc.alpha)[2]
    noMCMC =  dim(mc.alpha)[1]
    mc.theta = array(NA, dim = c( noMCMC, M, J))
    for(j in 1:J) {
        mc.theta[ , , j] = pnorm(t(W[ , j, ] %*% t(mc.alpha)))
    }
    
    


    ## ...compute posterior means and quantiles
    post.stats = array(NA, dim=c(M, J, 4))
    post.stats.MCSE = array(NA, dim=c(M, J, 4))
    for (j in 1:J) {
        postStats = EstimatePosteriorStats(mc.theta[, , j], burnin)
        post.stats[, j, ] = postStats$estimate
        post.stats.MCSE[, j, ] = postStats$MCerror
    }
    

    retVal = list(mean=post.stats[,,1], median=post.stats[,,2], lower=post.stats[,,3], upper=post.stats[,,4])
    if (mcError) {
        retVal = list(retVal, mean.MCSE=post.stats.MCSE[,,1], median.MCSE=post.stats.MCSE[,,2], lower.MCSE=post.stats.MCSE[,,3], upper.MCSE=post.stats.MCSE[,,4])
    }
    
    retVal
}
