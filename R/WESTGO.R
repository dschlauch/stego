#' Calculate the similarity matrix
#'
#' 
#'
#' @param genotypes data object containing the phased or unphased genotypes by samples
#' @param cex character expansion for the text
#' @param mar margin paramaters; vector of length 4 (see \code{\link[graphics]{par}})
#'
#' @return Object containing similarity matrix and other results
#'
#' @examples
#' data(toyGenotypes)
#' 
#' res <- westgo(toyGenotypes)
#' plotFromGSM(res, plotname="All Samples")
#' 
#' labels <- paste("Group",c(LETTERS[rep(1:5,40)]))
#' res <- westgo(toyGenotypes, groups="each.separately", labels=labels)
#' plotFromGSM(res)
#' 
#' labels <- paste("Group",c(LETTERS[rep(1:5,20)],LETTERS[rep(6:10,20)]))
#' super <- c(rep("Super A",100), rep("Super B",100))
#' res <- westgo(toyGenotypes, groups="pairwise.within.superpop", labels=labels, super=super)
#' plotFromGSM(res)
#'
#' @export
westgo <- function(genotypes,
                    diploid=F,
                    groups="all.together",
                    labels=NA,
                    super=NA,
                    minVariants=5, 
                    ldPrune=NA,
                    computeFST=T,
                    outputDir='.'){

    if (groups=="all.together"){
        genotypes <- pruneGenotypes(genotypes, ldPrune)
        results <- calculateSMatrix(gt=genotypes, diploid=diploid, minVariants=5, saveResult="all_together_results.rds")
        results$analysisType <- groups
        return(results)
    }
    
    assert_that(length(labels)==ncol(genotypes))
    assert_that(length(unique(labels))>1)
    
    results <- list()
    if (groups=="each.separately"){
        for(subpop in unique(labels)){
            print(gc())
            print(subpop)
            
            filter <- labels%in%subpop
            if(diploid){
                filter <- rep(filter,each=2)
            }
            gt <- pruneGenotypes(genotypes[,filter,with=F], ldPrune)
            results[[subpop]] <- calculateSMatrix(gt=gt, diploid=diploid, minVariants=minVariants, saveResult="each_separately_results.rds")
        }
        names(results) <- unique(labels)
        results$analysisType <- groups
        return(results)
    }
    assert_that(length(super)==ncol(genotypes))
    if (groups=="pairwise.within.superpop"){
        for(continent in unique(super)){
            pairs <- combn(unique(labels[super==continent]),2)
            for(i in seq_len(ncol(pairs))){
                pair <- pairs[,i]
                print(gc())
                print(pair)
                filter <- labels%in%pair
                if(diploid){
                    filter <- rep(filter,each=2)
                }
                gt <- pruneGenotypes(genotypes[,labels%in%pair,with=F], ldPrune)
                
                
                results[[paste(pair,collapse="_")]] <- calculateSMatrix(gt=gt, diploid=diploid, minVariants=minVariants, outputDir=outputDir, saveResult="pairwise_within_superpop_results.rds")
                
#                 if(computeFST){
#                     fstData <- as.data.frame(t(as.matrix(as.data.frame(gt)[,c(T,F)] + as.data.frame(gt)[,c(F,T)])))
#                     fstData <- cbind(pop = pop[filterdip], fstData)
#                     results[[paste(pair,collapse="_")]]$FST <- genet.dist(fstData,diploid=F,method="Fst")
#                 } 
            }
        }
        results$analysisType <- groups
        return(results)
    }
    
    stop("Something went wrong, no output")    
}

#' Prune genotypes based on blocks of a specified length
#'
#' @param gt data object containing the phased or unphased genotypes by samples
#'
#' @return Object containing similarity matrix and other results
#'
#' @examples
#' data(toyGenotypes)
#' pruneGenotypes(toyGenotypes)
#'
#' @export
pruneGenotypes <-  function(gt, ldPrune=1){
    
    if(is.na(ldPrune)||ldPrune==1){
        return(gt)
    }
    names(gt) <- make.unique(names(gt))
    numSamples <- ncol(gt)
    numVariants <- nrow(gt)
    sumVariants <- rowSums(gt)
    
    
    # reverse so that MAF<.5
    gt[sumVariants>(numSamples/2),] <- 1-gt[sumVariants>(numSamples/2),]
    sumVariants <- rowSums(gt)
    
    # Intelligently LD prune
    numblocks <- numVariants/ldPrune +1
    blocks <- rep(1:numblocks, each=ldPrune)[1:numVariants]
    
    system.time(runningWhichMax <- running(sumVariants,width=ldPrune,fun=which.max, by=ldPrune))
    prunedIndices <- runningWhichMax + seq(0,ldPrune*(length(runningWhichMax)-1),ldPrune)
    system.time(gt <- gt[prunedIndices])
    gt
}

#' Calculate the similarity matrix
#'
#' 
#'
#' @param genotypes data object containing the phased or unphased genotypes by samples
#' @param cex character expansion for the text
#' @param mar margin paramaters; vector of length 4 (see \code{\link[graphics]{par}})
#'
#' @return Object containing similarity matrix and other results
#'
#' @examples
#' data(toyGenotypes)
#' calculateSMatrix(toyGenotypes)
#'
#' @export
calculateSMatrix <- function(gt, diploid=F, minVariants=5, scaleBySampleAF=F, outputDir=".", saveResult=NA, varcov=T){
    
    numSamples <- ncol(gt)
    numVariants <- nrow(gt)
    sumVariants <- rowSums(gt)
    
    #     # reverse so that MAF<.5
    invertMinorAllele <- sumVariants>(numSamples/2)
    gt[invertMinorAllele] <- 1-gt[invertMinorAllele]
    sumVariants <- rowSums(gt)
    
    # remove < n variants
    sumVariants <- rowSums(gt)
    gt <- gt[sumVariants>minVariants,,with=F]
    gt <- as.matrix(gt)
    
    print("Number of used variants")
    print(nrow(gt))
    numFilteredVariants <- nrow(gt)
    sumFilteredVariants <- rowSums(gt)
    varcovMat <- NULL
    if(varcov){
        #         varcovMat <- cov(t(scale(t(gt[,c(T,F)] + gt[,c(F,T)]))),use="pairwise.complete.obs")
        #         varcovMat <- cov(gt[,c(T,F)] + gt[,c(F,T)],use="pairwise.complete.obs")
        varcovMat <- cor(gt[,c(T,F)] + gt[,c(F,T)],use="pairwise.complete.obs")
    }
    totalPossiblePairs <- choose(numSamples,2)
    totalPairs <- choose(sumFilteredVariants,2)
    weights <- totalPossiblePairs/totalPairs
    p <- 1/weights
    
    # var_s_hap <- sum((1-p)/p)/(numFilteredVariants^2)
    # recalculated variance 3/23/16
    var_s_hap <- sum(weights-1)/(numFilteredVariants^2)
    
    print("variance of s (haploid)")
    print(var_s_hap)
    
    # Calculate expected values conditional on kinship
    pkweightsMean <- mean(((sumFilteredVariants-2)/numSamples)*weights)
    kinships <- seq(0,.25,.001)
    kinshipExpectation <- 1+kinships*(pkweightsMean-1)
    
    s_matrix_numerator <- t(gt*weights)%*%gt
    s_matrix_denominator <- numFilteredVariants
    s_matrix_hap <- s_matrix_numerator/s_matrix_denominator
    if (scaleBySampleAF){
        numAllelesPerSample <- colSums(gt)
        relativeAF <- sqrt(mean(numAllelesPerSample)/numAllelesPerSample)
        #         s_matrix_hap <- s_matrix_hap*tcrossprod(relativeAF)
    }
    colnames(s_matrix_hap) <- colnames(gt)
    rownames(s_matrix_hap) <- colnames(gt)
    
    print(mean(s_matrix_hap[row(s_matrix_hap)!=col(s_matrix_hap)]))
    print(median(s_matrix_hap[row(s_matrix_hap)!=col(s_matrix_hap)]))
    
    estimatedKinship <- (s_matrix_hap-1)/(pkweightsMean-1)
    
    # Collapse to diploid
    s_matrix_dip <- (s_matrix_hap[c(T,F),c(T,F)] + s_matrix_hap[c(F,T),c(T,F)] +s_matrix_hap[c(T,F),c(F,T)] + s_matrix_hap[c(F,T),c(F,T)])/4
    colnames(s_matrix_dip) <- colnames(gt)[c(T,F)]
    rownames(s_matrix_dip) <- colnames(gt)[c(T,F)]
    
    var_s_dip <- var_s_hap/4
    
    
    popResult <- list(s_matrix_dip=s_matrix_dip, s_matrix_hap=s_matrix_hap, pkweightsMean=pkweightsMean, var_s_dip=var_s_dip, var_s_hap=var_s_hap, varcovMat=varcovMat)
    if(!is.na(saveResult)){
        saveRDS(popResult, saveResult)
    }
    popResult
    
}

#' Plot the distribution of similarity statistics
#'
#' 
#'
#' @param genotypes data object containing the phased or unphased genotypes by samples
#' @param cex character expansion for the text
#' @param mar margin paramaters; vector of length 4 (see \code{\link[graphics]{par}})
#'
#' @return Object containing similarity matrix and other results
#'
#' @examples
#' data(toyGenotypes)
#' 
#' westgo(toyGenotypes)
#' labels <- c(rep("Group A",100), rep("Group B",100))
#' westgo(toyGenotypes, groups="each.separately", labels=labels)
#' 
#' labels <- paste("Group",c(LETTERS[rep(1:4,25)],LETTERS[rep(5:8,25)]))
#' super <- c(rep("Super A",100), rep("Super B",100))
#' res <- westgo(toyGenotypes, groups="pairwise.within.superpop", labels=labels, super=super)
#'
#' @export
plotFromGSM <- function(resObj, plotname="", alphaCutoff=.01){
    if(is.null(resObj$analysisType)||resObj$analysisType=="all.together"){
        gsm <- resObj$s_matrix_dip
        var_s <- resObj$var_s_dip
        pkweightsMean <- resObj$pkweightsMean        
    } else {
        numPlots <- length(resObj)-1
        columns <- floor(sqrt(numPlots))+1
        rows <- ceiling(numPlots/columns)
#         par(mfrow=c(rows,columns))
        plots <- lapply(seq_len(length(resObj)-1), function(i){
            plotFromGSM(resObj=resObj[[i]], plotname=names(resObj)[i])
        })
        do.call("grid.arrange", c(plots, ncol=columns))
    }

    print(mean(gsm[row(gsm)!=col(gsm)]))
    print(median(gsm[row(gsm)!=col(gsm)]))
    num_comparisons_dip <- choose(ncol(gsm),2)
    sample_IDs <- rownames(gsm)
    bonferroni_cutoff_dip <- qnorm((1-alphaCutoff/2)^(1/num_comparisons_dip), sd=sqrt(var_s)) + 1
    
    topValuesDip <- sort(gsm[row(gsm)>col(gsm)], decreasing=T)
    topValuesKinship <- (topValuesDip-1)/(pkweightsMean-1)
    
    # Display only those that are above the cutoff and among the top 5
    label_cutoff <- max(bonferroni_cutoff_dip, topValuesDip[1])
    
    pairs <- outer(sample_IDs, sample_IDs, paste)
    plotData <- data.frame(values=gsm[row(gsm)>col(gsm)], pairs=paste0("  ",pairs[row(pairs)>col(pairs)]))
    minDip <- min(plotData$values)
    maxDip <- max(plotData$values)#ifelse(max(plotData$values)>bonferroni_cutoff_dip,max(plotData$values),NA)
    xmin <- min(minDip, 1-(maxDip-1)*.5)
    xmax <- maxDip + (maxDip-minDip)*.4
    dipPlot <- ggplot(plotData, aes(values)) + 
        geom_histogram(color="blue",bins=40,fill=I("blue")) + 
        ggtitle(plotname)  + xlab("Similarity score") + 
        #         scale_x_continuous(expand=c(.4,0))+#
        xlim(xmin, xmax) +
        theme_classic() +
        theme(plot.title = element_text(rel(2)), axis.title.x = element_text(size = rel(1)), axis.title.y = element_blank(), axis.text.y=element_blank()) + 
        
        geom_vline(xintercept = bonferroni_cutoff_dip, color="red", linetype="longdash") + 
        geom_vline(data=subset(plotData, (values == maxDip & values>bonferroni_cutoff_dip)),aes(xintercept = values), color="blue", linetype="dotted") + 
        geom_vline(xintercept = median(gsm[row(gsm)!=col(gsm)]), color="black", linetype=1) + 
        
        geom_text(data=subset(plotData, values >= label_cutoff), aes(values,label=pairs), y=0, angle = 80, hjust=0, size=rel(3)) +
        geom_text(data=subset(plotData, (values == maxDip & values>bonferroni_cutoff_dip)), x=maxDip, y=Inf, label=paste0("hat(phi)==", round(topValuesKinship[1],3),"  "),parse = TRUE, color="blue", angle = 0, size = 6, vjust = 2, hjust = 0) +
        
        #         annotate("text", x=bonferroni_cutoff_dip, y=Inf, label=paste0("alpha==",format(alphaCutoff/num_comparisons_dip, digits=1)),parse = TRUE, color="red", angle = 0, size = 6, vjust = 2, hjust = 1) +
        #         annotate("text", x=bonferroni_cutoff_dip, y=Inf, label=paste0("alpha "),parse = TRUE, color="red", angle = 0, size = 6, vjust = 2, hjust = 1.5) +
        annotate("text", x=median(gsm[row(gsm)!=col(gsm)]), y=Inf, label=paste0("m=",round(median(gsm[row(gsm)!=col(gsm)]),3)," "), color="black", angle = 0, size=rel(3), vjust=1.5, hjust = 1) 
    
    
    dipPlot    
    
}

