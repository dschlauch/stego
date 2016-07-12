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
    if(is.null(object$analysisType)||object$analysisType=="all.together"){
        cat(" stego Analysis\n\n")
        cat("Analysis Type: all.together\n")
        cat("Data Type:",ifelse(is.null(object$var_s_hap),"Unphased data\n","Phased data\n"))
        cat("Number of individuals:", ncol(object$s_matrix_dip), "\n\n")
        cat("Population Structure: p =", object$structurePValue, "\n")
        cat("Cryptic Relatedness: p =", object$crypticPValue, "\n\n")
        cat("Top related Pairs:\n")
        print(object$kinships[1:5])
    } else {
        cat(" stego Analysis\n\n")
        cat("Analysis Type:",object$analysisType,"\n")
        cat("Data Type:",ifelse(is.null(object$var_s_hap),"Unphased data\n","Phased data\n"))
        groups <- names(object)
        groups <- groups[which(groups!="analysisType")]
        numIndividuals <- lapply(groups, function(x){
            ncol(object[[x]]$s_matrix_dip)
        })
        cat("Total individuals:", sum(unlist(numIndividuals)),  "\n\n")
        cat("Groups:",paste(groups,collapse=","))
        
    }
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
    plotStegoHist(object, ...)
}