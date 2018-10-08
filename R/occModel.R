#' Specify and Fit a Multi-scale Occupancy Model
#'
#' This function specifies and fits a multi-scale occupancy model.
#' 
#' @param formulaSite model of occurrence at a site
#' @param formulaSiteAndSample model of occurrence in different samples of a site
#' @param formulaReplicate model of detection in replicates of a sample
#' @param detectionMats list containing
#' \describe{
#' \item{y:}{   M x J matrix of numbers of detections per sample}
#' \item{K:}{   M x J matrix of numbers of replicates per sample}
#' }
#' where M = no. sites and J = maximum no. samples per site
#' 
#' @param siteData data frame containing site-level covariates
#' @param siteAndSampleData data frame containing site- and sample-level covariates
#' @param niter no. iterations of MCMC algorithm
#' @param niterInterval no. iterations for reporting progress of MCMC algorithm
#' @param siteColName column of data frame siteData (or column of data frame siteAndSampleData) containing sites (i.e., primary sample units)
#' @param sampleColName column of data frame siteData (or column of data frame siteAndSampleData) containing sample numbers (i.e., integers 1, 2, ...) for each site
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
#' @details  This function is used to specify and fit a multi-scale occupancy model (with or without covariates).  The model is fitted using Bayesian methods (specifically, a Markov chain Monte Carlo (MCMC) algorithm that is run for a finite number of iterations).  Additional MCMC iterations may be computed using the function \code{updateOccModel}.   Output from the MCMC algorithm is stored in the file "mc.csv" and can be summarized using functions \code{plotTrace}, \code{plotACF}, and \code{posteriorSummary}.
#'
#' @seealso  \code{\link{updateOccModel}}, \code{\link{scaleData}}
#'
#'
#' @note  Before fitting the model, this function checks to ensure that the model specification is possible given the data files.  These checks include:
#' \itemize{
#' \item  Value of siteColName matches that in data frame siteData.
#' \item  Values of siteColName and sampleColName match those in data frame siteAndSampleData.
#' \item  Values in sampleColName contain consecutive positive integers (1, 2, ...) for each site.
#' \item  Only one data frame (siteData or siteAndSampleData) is specified for survey data.
#' }
#'
#' If any of these checks fail, the function returns an error message.
#'
#' 
#'
#' @examples
#' data(gobyDetectionData)
#' detections = occData(gobyDetectionData, "site", "sample")
#'
#' # Fit a model without covariates of occurrence or detection.
#' fit.simplest = occModel(detectionMats=detections)
#'
#'
#' # Fit a model which assumes that probability of occurrence in samples differs with site,
#' # but the differences are not a function of site-level covariates.
#' 
#' data(gobySurveyData)
#' fit = occModel(~1, ~factor(site)-1, ~1,
#'           detectionMats=detections,
#'           siteData=gobySurveyData,
#'           siteColName="site")
#'
#'
#' # Fit a model assuming occurrence and detection probabilities
#' # are functions of site-level covariates.
#' 
#' data(gobySurveyData)
#' gobySurveyData = scaleData(gobySurveyData)  # center and scale numeric covariates
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
#' 
#' # Update the Markov chain of the model specified in fit1
#' fit2 = updateOccModel(fit1, niter=50, niterInterval=10)




occModel <- function(
                     formulaSite = ~ 1,
                     formulaSiteAndSample = ~ 1,
                     formulaReplicate = ~ 1,
                     detectionMats,
                     siteData = NULL,
                     siteAndSampleData = NULL,
                     niter = 1100,
                     niterInterval = 100,
                     siteColName = 'site',
                     sampleColName = 'sample'
                     ){


    ## Make sure value of siteColName matches that in siteData data frame
    if (!is.null(siteData)) {
        if(!any(grepl(siteColName, names(siteData), fixed=TRUE))) {
            stop(paste("Column name of sites in data frame does not match", siteColName))
        }
    }
    
    ## Make sure values of siteColName and sampleColName match those in siteAndSampleData data frame
    if (!is.null(siteAndSampleData)) {
        if(!any(grepl(siteColName, names(siteAndSampleData), fixed=TRUE))) {
            stop(paste("Column name of sites in data frame does not match", siteColName))
        }
        if(!any(grepl(sampleColName, names(siteAndSampleData), fixed=TRUE))) {
            stop(paste("Column name of samples in data frame does not match", sampleColName))
        }
    }
    
        

    ## Make sure the values in sampleColName contain consecutive positive integers (1, 2, ...) for each site
    if (!is.null(siteAndSampleData)) {
        siteID = unique(siteAndSampleData[, siteColName])
        for (i in 1:length(siteID)) {
            ind = siteAndSampleData[, siteColName] == siteID[i]
            sampleNum = siteAndSampleData[ind, sampleColName]
            if (any(sort(sampleNum) != (1:length(sampleNum)) )) {
                errMsg = paste("Sample numbers of site", siteID[i],
                               "do not contain consecutive positive integers in data frame")
                stop(errMsg)
            }
        }
    }
    
    

    ## Make sure that siteData and siteAndSampleData are not both specified
    if( !is.null(siteData) & !is.null(siteAndSampleData)){
        errMsg = paste("Data frames of site-level and site-and-sample-level covariates cannot both be specified.\nPlease use only one data frame.")
        stop(errMsg)
    }

    



    
    ## Assign values of replicate-level data matrices
    y = detectionMats$y # observed data (0s and 1s)
    K = detectionMats$K # number of replicates per site-sample combination
    M = nrow(K) # number of sites
    J = ncol(K) # maximum number of samples per site

    ## Initialize value of flag for model of fixed site effects for theta parameter vector
    siteEffectInW = FALSE
    colNamesOfX = NULL
    colNamesOfW = NULL
    colNamesOfV = NULL


    
    ##  Assign values to the model matrices (X and W and V)
    
    ## ... Only 1s case: do not have Site or Site And Sample covariates
    if( is.null(siteAndSampleData) & is.null(siteData)) {
        
        X = matrix(1, nrow = M, ncol = 1)
        W = array( 1,  dim = c( M, J, 1))
        V = array( 1 , dim = c( M, J, 1))
        colNamesOfX = '(Intercept)'
        colNamesOfW = '(Intercept)'
        colNamesOfV = '(Intercept)'        
    }

    ## ... Only have Site covariates, but not Site And Sample covariates
    if( is.null(siteAndSampleData) & !is.null(siteData)) {
        
        mfSite <- model.frame(formula=formulaSite, data = siteData)
        X      <- model.matrix(attr(mfSite, "terms"), data = mfSite)
        colNamesOfX = dimnames(X)[[2]]
        

        mfSiteAndSample <- model.frame(formula=formulaSiteAndSample, data=siteData)
        Wmat     <- model.matrix(attr(mfSiteAndSample, "terms"), data = mfSiteAndSample)
        colNamesOfW = dimnames(Wmat)[[2]]

        wNames = attr(attr(mfSiteAndSample, "terms"), "term.labels")        
        if (length(wNames)>0) siteEffectInW = any(grepl(siteColName, wNames, fixed=TRUE))

        W = array( dim = c( M, J, ncol(Wmat)))
        for (j in 1:J) {
            W[,j,] = Wmat
        }

        
        mfReplicate <- model.frame(formula = formulaReplicate, data = siteData)
        Vmat        <- model.matrix(attr(mfReplicate, "terms"), data = mfReplicate)
        colNamesOfV = dimnames(Vmat)[[2]]

        V = array( dim = c( M, J, ncol(Vmat)))
        for (j in 1:J) {
            V[,j,] = Vmat
        }
    }

    ## ... Have covariates for both Site and Sample
    if( !is.null(siteAndSampleData) & is.null(siteData)) {

        siteID = dimnames(y)[[1]]  # need this to ensure that rows of y and covariate matrices are consistent

        ## Extract site-level covariates for model matrix X
        siteData <- siteAndSampleData[match(siteID, siteAndSampleData[, siteColName]),  ]

        mfSite <- model.frame(formula=formulaSite, data = siteData)
        X      <- model.matrix( attr( mfSite, "terms"), data = mfSite)
        colNamesOfX = dimnames(X)[[2]]


        ## Use site- and sample-level covariates for model matrices W and V

        mfSiteAndSample <- model.frame(formula =formulaSiteAndSample, data = siteAndSampleData)
        Wmat     <- model.matrix( attr(mfSiteAndSample, "terms"), data = mfSiteAndSample)
        colNamesOfW = dimnames(Wmat)[[2]]

        
        wNames = attr(attr(mfSiteAndSample, "terms"), "term.labels")        
        if (length(wNames)>0) siteEffectInW = any(grepl(siteColName, wNames, fixed=TRUE))

        siteIndex = match(siteAndSampleData[, siteColName], siteID)

        W = array( dim = c(M, J, ncol(Wmat)))
        for( i in 1:nrow(Wmat)){
            W[ siteIndex[i],  siteAndSampleData[i, sampleColName],  ] <- Wmat[i , ]
        }

        
        mfReplicate <- model.frame(formula = formulaReplicate, data = siteAndSampleData)
        Vmat     <- model.matrix( attr(mfReplicate, "terms"), data = mfReplicate)

        V = array( dim = c(M, J, ncol(Vmat)))
        colNamesOfV = dimnames(Vmat)[[2]]
        for( i in 1:nrow(Vmat)){
            V[ siteIndex[i],  siteAndSampleData[i, sampleColName], ] <- Vmat[i , ]
        }
    }



 
    ##  Fit occupancy model to data or update model if already fitted
    fit = fitOccModel(niter = niter, niterInterval = niterInterval,
                              y, K, X, W, V, siteEffectInW, colNamesOfX, colNamesOfW, colNamesOfV)
    class(fit) <- "occModel"
    return(fit)
}
