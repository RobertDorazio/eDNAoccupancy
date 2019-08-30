#' Posterior Summary For Parameters of a Multi-scale Occupancy Model
#'
#' Estimates the posterior mean, median, and 95\% credible limits for each parameter of a multi-scale occupancy model.
#'
#' @param fit  object of class occModel that contains data and previous state of the model's Markov chain
#' @param burnin  initial no. iterations of Markov chain to be omitted from calculations
#' @param mcError logical switch to estimate Monte Carlo standard errors
#' @param outputSummary   logical switch to return values of posterior summary statistics
#' 
#' @importFrom utils read.csv
#' @export
#'
#'
#' @return Prints estimates of posterior summary statistics.  If \code{mcError}=TRUE, the Monte Carlo standard errors of these estimates are computed.   If \code{outputSummary}=TRUE, the posterior summaries are returned in a list.
#'
#' @details  This function estimates the posterior mean, median, and 95\% credible limits for each parameter of a multi-scale occupany model using a Markov chain created by \code{occModel}.
#'
#'
#' @examples
#' data(gobyDetectionData)
#' detections = occData(gobyDetectionData, "site", "sample")
#' fit = occModel(detectionMats=detections)
#' posteriorSummary(fit, burnin=100,  mcError = TRUE)




posteriorSummary <- function(
                          fit,
                          burnin = 1,
                          mcError = FALSE,
                          outputSummary = FALSE
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
    mcColumnNames = dimnames(read.csv('mc.csv'))[[2]]
    if (length(mc.names) != length(mcColumnNames)) {
        stop(paste("Column names in file 'mc.csv' do not match the model matrices of the occModel object"))
    }
    if (any(make.names(mc.names, unique=TRUE) != mcColumnNames)) {
        stop(paste("Column names in file 'mc.csv' do not match the model matrices of the occModel object"))
    }
    
    
    ## Read Markov chain from file
    mc = as.matrix(read.csv('mc.csv'))
    colnames(mc) = mc.names
    
    
    ## Estimate posterior means and quantiles
    postStats = EstimatePosteriorStats(mc, burnin)
    

    ## Print estimated summaries of posterior
    CR = '\n'
    cat ('Bayesian estimates of model parameters', CR)
    print(round(postStats$estimate, 3))

    if (mcError){
       cat (CR)
       cat('Monte Carlo SE of Bayesian estimates', CR)
       print(round(postStats$MCerror,4))
    }



    retVal = invisible()
    
    if (outputSummary) {
        if (mcError) {
            retVal = list(post.stats=postStats$estimate, post.stats.MCSE=postStats$MCerror)
        }
        else {
            retVal = list(post.stats=postStats$estimate)
        }
    }
    
    retVal
}
