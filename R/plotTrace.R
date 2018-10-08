#' Plot Trace of One or More Parameters of a Multi-scale Occupancy Model
#'
#' Plots the Markov chain for each of 1-4 parameters of a multi-scale occupancy model.
#'
#' 
#' @param fit  object of class occModel that contains data and previous state of the model's Markov chain
#' @param paramName  vector of names of parameters to be plotted
#' @param burnin  initial no. iterations of Markov chain to be omitted from plot
#' 
#' @importFrom utils read.csv
#' @importFrom graphics plot par
#' @export
#'
#' 
#' @details  This function plots the Markov chain created by \code{occModel} for each of 1-4 parameters of a multi-scale occupancy model.
#'
#'
#' @examples
#' data(gobyDetectionData)
#' detections = occData(gobyDetectionData, "site", "sample")
#' fit = occModel(detectionMats=detections)
#' plotTrace(fit, 'beta..Intercept.', burnin=100)
#' plotTrace(fit, c('beta..Intercept.', 'alpha..Intercept.', 'delta..Intercept.'), burnin=100)




plotTrace = function(fit, paramName, burnin = 1) {
    
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
    
    
    ## Read Markov chain from file
    outSave = -(1:burnin)
    mc = as.matrix(read.csv("mc.csv"))
    mc = mc[ outSave, ]
    mc.params = dimnames(mc)[[2]]
    
    ## ... make sure paramName is valid
    nparam = length(paramName)

    for (i in 1:nparam) {
        if ( !any(grepl(paramName[i], mc.params)) ) {
            errMsg = paste(paramName[i], 'is not a valid parameter name')
            stop(errMsg)
        }
    }

    if (nparam > 4)  stop('Number of parameters for plotting cannot exceed 4!')

    cexLab = 1.3

    if (nparam==1) {
        plot(1:nrow(mc), mc[, paramName], type='l',  xlab='Iteration', ylab=paramName, cex.lab=cexLab)
    }
    else if (nparam==2) {
        par(mfrow=c(2,1))
        plot(1:nrow(mc), mc[, paramName[1]], type='l',  xlab='Iteration', ylab=paramName[1], cex.lab=cexLab)
        plot(1:nrow(mc), mc[, paramName[2]], type='l',  xlab='Iteration', ylab=paramName[2], cex.lab=cexLab)
    }
    else {
        par(mfrow=c(2,2))
        for (i in 1:nparam) {
            plot(1:nrow(mc), mc[, paramName[i]], type='l',  xlab='Iteration', ylab=paramName[i], cex.lab=cexLab)
        }
    }
    
    invisible()
}
