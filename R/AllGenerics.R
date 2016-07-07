#' summary.stego
#'
#' summarizes the results of a stego analysis
#'
#' @param object an object of class "stego"
#' @param ... further arguments passed to or from other methods.
#' @keywords keywords
#' @export
#' @return Summary description of stego S3 object
summary.stego <- function(object, ...){
    cat("Summary")
}
#' print.stego
#'
#' prints the results of a stego analysis
#'
#' @param object an object of class "stego"
#' @param ... further arguments passed to or from other methods.
#' @keywords keywords
#' @export
#' @return Summary description of stego S3 object
print.stego <- function(object, ...){
    cat(" stego Analysis\n\n")
    cat(ifelse(is.null(object$var_s_hap),"Unphased data\n","Phased data\n"))
    cat(paste(ncol(object$s_matrix_dip), "individuals\n\n"))
    cat("Top related Pairs:\n")
    print(object$summary[1:5])
}

#' plot.stego
#'
#' plots the results of a stego analysis
#'
#' @param object an object of class "stego"
#' @param ... further arguments passed to or from other methods.
#' @keywords keywords
#' @export
#' @return Summary description of stego S3 object
plot.stego <- function(object, ...){
    plotFromGSM(object, ...)
}