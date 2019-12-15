#' Center and Scale Numeric Columns of a Data Frame
#'
#' Centers and scales all numeric- or integer-valued columns of a data frame object.
#' 
#' @param df  a data frame object
#' @param sampleColName column of data frame containing sample numbers (i.e., integers 1, 2, ...) for each site
#' @export
#'
#' @return The input data frame with all numeric- or integer-valued columns transformed to have mean zero and unit variance.
#'
#' @details   This function extracts all numeric- or integer-valued columns of a data frame object.  Each of these columns is centered by subtracting its mean and then scaled by dividing the result by the standard deviation of the original column values.  This linear transformation produces columns with mean zero and unit variance.  Both input and output data frame objects have the same column ordering.  sampleColName is only needed if data frame contains sample numbers.




scaleData = function(df, sampleColName = NULL) {

    df.names = colnames(df)

    ## Make sure value of sampleColName matches that in  file
    if (!is.null(sampleColName)) {
        if(!any(grepl(sampleColName, df.names, fixed=TRUE))) {
            stop(paste("Column name of samples in file does not match", sampleColName))
        }
    }

    
    ## Separate sample column from rest of data frame
    if (!is.null(sampleColName)) df.sample = data.frame(samp = df[, sampleColName])
    if (!is.null(sampleColName)) names(df.sample) = sampleColName
    if (!is.null(sampleColName)) df = df[, -match(sampleColName, df.names)]
        

    ## Split into numbers and integers vs non-numerics
    ind = unlist(lapply(df, class)) == "numeric" | unlist(lapply(df, class)) == "integer"
    xmat = as.matrix(df[ , ind])
    colnames(xmat) = names(df)[ind]

    ind =  unlist(lapply(df, class)) != "numeric" & unlist(lapply(df, class)) != "integer"
    notXmat =  data.frame(df[ , ind ])
    colnames(notXmat) = names(df)[ind]

    ## Center and scale columns of numbers or integers
    xMean = apply(xmat, 2, mean, na.rm=TRUE)
    xSD = apply(xmat, 2, sd, na.rm=TRUE)
    oneVec = matrix(1, nrow = nrow(xmat), ncol=1)
    xmat = xmat - kronecker(oneVec, matrix(xMean,nrow=1))
    xmat = xmat / kronecker(oneVec, matrix(xSD,nrow=1))

    ## Recombine numeric and non-numeric data, preserving original order of df's columns
    df.out = data.frame(notXmat, xmat)
    if (!is.null(sampleColName))  df.out = data.frame(df.out, df.sample)
    df.out = df.out[ , df.names]
    return(df.out)
}

