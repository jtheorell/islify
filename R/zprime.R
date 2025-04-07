#' Compute the Z'-factor quality score
#'
#' This function is derived in its entirety from the deprecated imageHTS 
#' package maintained by Gregoire Pau <gregoire.pau at embl.de>. 
#' All the documentation below as well as the code is identical to the
#' original text. Only the simplified example is new. So kudos to the original
#' authors! This re-publication is just there as no alternative has seemingly
#' been available in R for years and it is needed for the downstream analyses
#' with islify data. 
#' 
#' The Z'-factor is a popular metric measuring the separation of control
#' features in high-throughput screens. The original paper describing the
#' Z'-factor is Zhang, 1999, J Biomol Screen.
#' 
#' @param a A vector or matrix of control features. For example, positive
#' controls.
#' @param b A vector or matrix of different control values. For example negative
#' controls. 
#' @param method A character vector, indicating which method should be used to
#' compute the Z'-factor. Default is mahalanobis. See details.
#' @references J. H. Zhang, T. D. Chung, K. R. Oldenburg. A Simple
#' Statistical Parameter for Use in Evaluation and Validation of High
#' Throughput Screening Assays. J Biomol Screening, 1999.
#' @details
#' #' Several univariate Z'-factor scores exist. The original Z'-factor from
#' Zhang, 1999 is computed by Z' = 1 - 3*(sd(a)+sd(b))/abs(mean(a)-mean(b)). 
#' A more rigorous definition of the score, implemented by the method fixsd is 
#' given by Z' = 1 - 3*sqrt(var(a)+var(b))/abs(mean(a)-mean(b)), where the
#' pooled standard deviation is computed by the square root of the sum of the
#' control variances. A robust method, less sensitive to outliers, is computed
#' by the relation Z' = 1-3*(mad(a)+mad(b))/abs(median(a)-median(b)) where the
#' control dispersions are computed with the mad and the control locations
#' with the median.
#' 
#' A multivariate extension of the Z'-factor score can be designed by
#' linearly transforming the multivariate data to one dimension and computing
#' the standard (here, fixsd) Z'-factor. It can be shown that the linear
#' transform that maximizes the score is the LDA. Moreover, one can
#' demonstrate that the resulting Z'-factor score is equivalent of computing
#' Z' = 1 - 3/dMaha(mu_a, mu_b, Sigma_a + Sigma_b) where dMaha is the
#' Mahalanobis distance.
#' 
#' @examples 
#' 
#' #Generate some data
#' a <- rnorm(100, mean = 0, sd = 1)
#' b <- rnorm(100, mean = 7, sd = 1)
#' 
#' zprime(a, b)
#' @export zprime
zprime = function(a, b, method=c('mahalanobis', 'robust', 'fixsd', 'original')) {
    method = match.arg(method)
    
    if (method=='mahalanobis') {
        if (is.null(dim(a))) a = matrix(a, ncol = 1)
        if (is.null(dim(b))) b = matrix(b, ncol = 1)
        mua = apply(a, 2, mean)
        mub = apply(b, 2, mean)
        dm = try(mahalanobis(mua, mub, cov(a) + cov(b)))
        if (inherits(dm, what ='try-error')) NA
        else 1-3/sqrt(dm)
    }
    else {
        if (method=='robust') 1-3*(mad(a)+mad(b))/abs(median(a)-median(b))
        else if (method=='fixsd') 1-3*sqrt(var(a)+var(b))/abs(mean(a)-mean(b))
        else 1-3*(sd(a)+sd(b))/abs(mean(a)-mean(b))
    }
}