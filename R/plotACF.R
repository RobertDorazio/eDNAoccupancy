#' Plot ACF for One or More Parameters of a Multi-scale Occupancy Model
#'
#' Plots the ACF (autocorrelation function) of the Markov chain for each of 1-4 parameters of a multi-scale occupancy model.
#'
#' 
#' @param fit  object of class occModel that contains data and previous state of the model's Markov chain
#' @param paramName  vector of names of parameters to be plotted
#' @param burnin  initial no. iterations of Markov chain to be omitted from plot
#' 
#' @importFrom utils read.csv
#' @importFrom graphics plot par
#' @importFrom stats acf
#' @export
#'
#' 
#' @details  This function plots the ACF (autocorrelation function) of the Markov chain created by \code{occModel} for each of 1-4 parameters of a multi-scale occupancy model.
#'
#'
#' @examples
#' data(gobyDetectionData)
#' detections = occData(gobyDetectionData, "site", "sample")
#' fit = occModel(detectionMats=detections)
#' plotACF(fit, c('beta.(Intercept)', 'alpha.(Intercept)', 'delta.(Intercept)'), burnin=100)



plotACF = function(fit, paramName, burnin = 1) {
    
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
    outSave = -(1:burnin)
    mc = as.matrix(read.csv("mc.csv"))
    colnames(mc) = mc.names
    mc = mc[ outSave, ]
    
    ## ... make sure paramName is valid
    mc.params = dimnames(mc)[[2]]
    
    nparam = length(paramName)

    for (i in 1:nparam) {
        if ( !any(paramName[i]==mc.names) ) {
            errMsg = paste(paramName[i], 'is not a valid parameter name')
            stop(errMsg)
        }
    }

    if (nparam > 4)  stop('Number of parameters for plotting cannot exceed 4!')

    cexLab = 1.3

    if (nparam==1) {
        acf(mc[, paramName], main=paramName, cex.lab=cexLab)
    }
    else if (nparam==2) {
        par(mfrow=c(2,1))
        acf(mc[, paramName[1]], main=paramName[1], cex.lab=cexLab)
        acf(mc[, paramName[2]], main=paramName[2], cex.lab=cexLab)
    }
    else {
        par(mfrow=c(2,2))
        for (i in 1:nparam) {
            acf(mc[, paramName[i]], main=paramName[i], cex.lab=cexLab)
        }
    }
    
    invisible()
}
