#' Estimate posterior stats and Monte Carlo errors
#'
#'
#' @param mc      Markov chain
#' @param burnin  number of initial elements of Markov chain to discard
#'
#'  
## #' @importFrom mcse, mcse.q
#' @import stats mcmcse
#' @keywords internal
## #' @export


### Define function for estimating posterior stats and Monte Carlo errors

EstimatePosteriorStats = function(mc, burnin=1) {
    
    mc.names = colnames(mc)

## ...estimate posterior means and quantiles
    out = -(1:burnin)
    mc = matrix(mc[out, ], ncol=ncol(mc))

    prob.quantiles = c(.50, .025, .975)  # for credible limits
    prob.names = paste(as.character(100*prob.quantiles), '%', sep='')
    post.names = mc.names
    post.stats = matrix(nrow=ncol(mc), ncol=1+length(prob.quantiles))
    dimnames(post.stats) = list(post.names, c('Mean', prob.names))
    post.stats.MCSE = post.stats


    for (j in 1:ncol(mc)) {
        xvec = as.vector(mc[, j])
        mc.estimate = suppressWarnings(tryCatch( mcse(x=xvec, method='obm'), error=function(condition) TRUE ))
        if (!is.logical(mc.estimate)) {
            post.stats[j, 1] = mc.estimate$est
            post.stats.MCSE[j, 1] = mc.estimate$se
        }
        else {
            post.stats[j, 1] = mean(xvec)
            ##post.stats.MCSE[j, 1] = NULL
        }
        
        for (k in 1:length(prob.quantiles)) {
            mc.estimate = suppressWarnings(tryCatch( mcse.q(x=xvec, q=prob.quantiles[k], method='obm'), error=function(condition) TRUE ))
            if (!is.logical(mc.estimate)) {
                post.stats[j, 1+k] = mc.estimate$est
                post.stats.MCSE[j, 1+k] = mc.estimate$se
            }
            else {
                post.stats[j, 1+k] = quantile(xvec, probs=prob.quantiles[k])
                ##post.stats.MCSE[j, 1+k] = NULL
            }
        }
    }
    list(estimate=post.stats, MCerror=post.stats.MCSE)
}
