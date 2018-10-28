#' Posterior Summary of a Multi-scale Occupancy Model's Site-Occupancy Probabilities
#'
#' Estimates the posterior mean, median, and 95\% credible limits for a multi-scale occupancy model's site-specific occurrence probabilities.
#'
#' 
#' @param fit  object of class occModel that contains data and previous state of the model's Markov chain
#' @inheritParams posteriorSummary
#' @importFrom utils read.csv
#' @export
#'
#' @return Computes estimates of summary statistics for a multi-scale occupancy model's site-occupancy probabilities.  If \code{mcError}=TRUE, the Monte Carlo standard errors of these estimates are computed.   All posterior summaries are returned in a list.
#'
#' @details  This function estimates the posterior mean, median, and 95\% credible limits of the site-specific occurrence probabilities of a multi-scale occupancy model.
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
#' psi = posteriorSummaryOfSiteOccupancy(fit1, burnin=100)
#' plot(gobySurveyData[, 'turb'], psi$median)




posteriorSummaryOfSiteOccupancy <- function(
                          fit,
                          burnin = 1,
                          mcError = FALSE
                          ) {


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
    
    

    ## Compute psi vector for each draw of beta in Markov chain
    X = fit$X
    colNamesOfX = fit$colNamesOfX
    beta.names = make.names(paste('beta', colNamesOfX, sep='.'))
    
    ## Read Markov chain from file
    mc = as.matrix(read.csv("mc.csv"))
    mc.beta = as.matrix(mc[ , beta.names])
    
    mc.psi = pnorm(t(X %*% t(mc.beta)))

    ## Estimate posterior means and quantiles
    postStats = EstimatePosteriorStats(mc.psi, burnin)
    post.names = dimnames(fit$y)[[1]]
    dimnames(postStats$estimate)[[1]] = post.names
    dimnames(postStats$MCerror)[[1]] = post.names
    
    

    retVal = list(mean=postStats$estimate[,1], median=postStats$estimate[,2], lower=postStats$estimate[,3], upper=postStats$estimate[,4])
    if (mcError) {
        retVal = list(mean=postStats$estimate[,1], median=postStats$estimate[,2], lower=postStats$estimate[,3], upper=postStats$estimate[,4], mean.MCSE=postStats$MCerror[,1], median.MCSE=postStats$MCerror[,2], lower.MCSE=postStats$MCerror[,3], upper.MCSE=postStats$MCerror[,4])
    }
    
    retVal
}
