#' Compute Occupancy Data Matrices
#'
#' Computes occupancy data matrices from a data frame of tabular detection data.
#' 
#' @param d   data frame containing tabular detection data
#' @param siteColName  column of data frame d containing sites (i.e., primary sample units)
#' @param sampleColName  column of data frame d containing sample numbers (i.e., integers 1, 2, ...) for each site
#' @export
#'
#'
#' 
#' @return A list containing two matrices:
#' \describe{
#' \item{y:}{   M x J matrix of numbers of detections per sample}
#' \item{K:}{   M x J matrix of numbers of replicates per sample}
#' }
#' where M = no. sites and J = maximum no. samples per site
#'
#' 
#' @details This function takes a data frame of tabular detection data and returns a list containing a matrix of numbers of detections per sample and a matrix of numbers of replicates per sample.
#'
#' The number of samples per site may differ among sites.  In addition, the number of replicates per sample may differ among samples.  Such unbalanced designs will induce NAs in the detection matrix (y) and zeros in the replicate number matrix (K).
#'
#' 
#' @examples
#' data(gobyDetectionData)
#' 
#' occData(gobyDetectionData, "site", "sample")




occData = function(d, siteColName = 'site', sampleColName = 'sample')  {


    ## check and see if any columns are all NAs. This may
    if(any(apply(apply(d, 2, is.na), 2, sum) == nrow(d))){
        warning("Your data file appears to have a column of all NAs. You may want to double check your data. This may have indicated that your file editor (e.g., Excel introduced extra columns in your data.")
    }

    
    ## Make sure values of siteColName and sampleColName match those in  file
    if(!any(grepl(siteColName, names(d), fixed=TRUE))) {
        stop(paste("Column name of sites in detection file does not match", siteColName))
    }
    if(!any(grepl(sampleColName, names(d), fixed=TRUE))) {
        stop(paste("Column name of samples in detection file does not match", sampleColName))
    }


    ## Make sure the values in sampleColName contain consecutive positive integers (1, 2, ...) for each site
    siteID = unique(d[, siteColName])
    for (i in 1:length(siteID)) {
        ind = d[, siteColName] == siteID[i]
        sampleNum = d[ind, sampleColName]
        if (!identical(sort(sampleNum), 1:length(sampleNum))) {
            errMsg = paste("Sample numbers of site", siteID[i],
                           "do not contain consecutive positive integers in detection file")
            stop(errMsg)
        }
    }



    ## Assign values to data matrices y and K
    M = length(unique(d[, siteColName]))
    J = max(d[, sampleColName])
    y = matrix(nrow=M, ncol=J)
    K = matrix(nrow=M, ncol=J)

    siteID = unique(d[, siteColName])
    for (i in 1:length(siteID)) {
        ind = d[, siteColName] == siteID[i]
        sampleNum = d[ind, sampleColName]
        ymat = d[ind, -c(1:2)]
        y[i, sampleNum] = rowSums(ymat, na.rm=TRUE)
        K[i, sampleNum] = rowSums(!is.na(ymat))
    }
    indMiss = is.na(K)
    K[indMiss] = 0
    y[K==0] = NA
    dimnames(y)[[1]] = as.character(siteID)
    dimnames(K)[[1]] = as.character(siteID)

    list(y=y, K=K)
}


