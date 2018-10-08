#' Posterior Summary of the AUC of a Multi-scale Occupancy Model's Predictions
#'
#' Estimates the posterior mean, median, and 95\% credible limits for the AUC of a multi-scale occupancy model's predictions.
#'
#' 
#' @param fit  object of class occModel that contains data and previous state of the model's Markov chain
#' @inheritParams posteriorSummary
#' @importFrom utils read.csv read.table
#' @importFrom pROC roc
#' @export
#'
#'
#' @return Prints estimates of summary statistics for the posterior distribution of the AUC of a multi-scale occupancy model's predictions.  If \code{mcError}=TRUE, the Monte Carlo standard errors of these estimates are computed.   If \code{outputSummary}=TRUE, the posterior summaries are returned in a list.
#'
#' @details  This function estimates the posterior mean, median, and 95\% credible limits of the AUC (Area Under the Receiver Operating Characteristic curve) of a multi-scale occupancy model's predictions.  The AUC of a model is sometimes used in model selection with higher values of AUC favored over lower values.
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
#' posteriorSummaryOfAUC(fit1, burnin=100, mcError=TRUE)





posteriorSummaryOfAUC <- function(
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
    mc.beta = mc[ , beta.names]
    
    mc.psi = pnorm(t(X %*% t(mc.beta)))
    
    ## Read in Z vector for each draw in Markov chain
    mc.Z = as.matrix(read.table("mc.Z.txt", header=TRUE))

    
    ## Compute AUC for each draw in Markov chain
    mc.AUC = rep(NA, nrow(mc.Z))
    mc.sumOfZ = rowSums(mc.Z)
    for (i in 1:nrow(mc.Z)) {
        if (mc.sumOfZ[i]>0 & mc.sumOfZ[i]<ncol(mc.Z)) mc.AUC[i] = as.numeric(roc(mc.Z[i, ], mc.psi[i, ])$auc)
    }
    if (sum(is.na(mc.AUC))>0) {
        stop(paste('AUC cannot be computed for', sum(is.na(mc.AUC)), 'elements of Markov chain wherein site occupancy predictions were all zeros or all ones'))
    }
    

    

    ## ...compute posterior means and quantiles
    postStats = EstimatePosteriorStats(matrix(mc.AUC,ncol=1), burnin)
    post.stats = postStats$estimate
    post.stats.MCSE = postStats$MCerror

    ## Print estimated summaries of posterior
    CR = '\n'
    cat ('Bayesian estimates of AUC', CR)
    print(round(post.stats, 3))

    if (mcError){
      cat (CR)
      cat('Monte Carlo SE of Bayesian estimates', CR)
        print(round(post.stats.MCSE,4))
    }

    
    retVal = invisible()

    if (outputSummary) {
        if (mcError) {
            retVal = list(post.stats=post.stats, post.stats.MCSE=post.stats.MCSE)
        }
        else {
            retVal = list(post.stats=post.stats)
        }
    }
    
    retVal
}
