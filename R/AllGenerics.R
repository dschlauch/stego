#' summary.westgo
#'
#' summarizes the results of a WESTGO analysis
#'
#' @param object an object of class "westgo"
#' @param ... further arguments passed to or from other methods.
#' @keywords keywords
#' @export
#' @return Summary description of westgo S3 object
summary.westgo <- function(object, ...){
    cat("Summary")
}
#' print.westgo
#'
#' prints the results of a WESTGO analysis
#'
#' @param object an object of class "westgo"
#' @param ... further arguments passed to or from other methods.
#' @keywords keywords
#' @export
#' @return Summary description of westgo S3 object
print.westgo <- function(object, ...){
    cat(" WESTGO Analysis\n\n")
    cat(ifelse(is.null(object$var_s_hap),"Unphased data\n","Phased data\n"))
    cat(paste(ncol(object$s_matrix_dip), "individuals\n\n"))
    cat("Top related Pairs:\n")
    print(object$summary[1:5])
}