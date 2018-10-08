#' Posterior Summary of a Multi-scale Occupancy Model's Detection Probabilities
#'
#' Estimates the posterior mean, median, and 95\% credible limits for a multi-scale occupancy model's sample-specific detection probabilities.
#'
#' 
#' @param fit  object of class occModel that contains data and previous state of the model's Markov chain
#' @inheritParams posteriorSummary
#' @importFrom utils read.csv
#' @export
#'
#' @return Computes estimates of summary statistics for a multi-scale occupancy model's sample-specific detection probabilities.  If \code{mcError}=TRUE, the Monte Carlo standard errors of these estimates are computed.   All posterior summaries are returned in a list.
#'
#' @details  This function estimates the posterior mean, median, and 95\% credible limits of the sample-specific detection probabilities of a multi-scale occupancy model.
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
#' p = posteriorSummaryOfSampleOccupancy(fit1, burnin=100)
#' plot(gobySurveyData[, 'sal'], p$median[,1])



posteriorSummaryOfDetection <- function(
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
    V = fit$V
    colNamesOfV = fit$colNamesOfV
    delta.names = make.names(paste('delta', colNamesOfV, sep='.'))

    
    ## Read Markov chain from file
    mc = as.matrix(read.csv("mc.csv"))
    mc.delta = as.matrix(mc[ , delta.names])

    M =  dim(V)[1]
    J =  dim(V)[2]
    lenDelta = dim(mc.delta)[2]
    noMCMC =  dim(mc.delta)[1]
    mc.p = array(NA, dim = c( noMCMC, M, J))
    for(j in 1:J) {
        mc.p[ , , j] = pnorm(t(V[ , j, ] %*% t(mc.delta)))
    }
    
    

    ## ...compute posterior means and quantiles
    post.stats = array(NA, dim=c(M, J, 4))
    post.stats.MCSE = array(NA, dim=c(M, J, 4))
    for (j in 1:J) {
        postStats = EstimatePosteriorStats(mc.p[, , j], burnin)
        post.stats[, j, ] = postStats$estimate
        post.stats.MCSE[, j, ] = postStats$MCerror
    }
    

    retVal = list(mean=post.stats[,,1], median=post.stats[,,2], lower=post.stats[,,3], upper=post.stats[,,4])
    if (mcError) {
        retVal = list(retVal, mean.MCSE=post.stats.MCSE[,,1], median.MCSE=post.stats.MCSE[,,2], lower.MCSE=post.stats.MCSE[,,3], upper.MCSE=post.stats.MCSE[,,4])
    }
    
    retVal
}
