#' Similarity Test for Genetic Outliers
#' 
#' This function runs the analysis.  A genotype matrix ([0,1] for phased, or [0,1,2] for unphased) is the only required argument.
#'
#' 
#'
#' @param genotypes data object containing the phased or unphased genotypes by samples
#' @param phased logical defining whether data exists as phased data, as opposed to unphased data
#' @param groups character specifying grouping of analysis.  Default is to run analysis all at once- one of "all.together", "each.separately" or "pairwise.within.superpop"
#' @param sampleNames character vector with unique identifiers for each sample
#' @param labels character covariates, such as population membership.  This is unused if groups is "all.together".
#' @param super character covariates, such as population membership.  This is used only if groups is "pairwise.within.superpop".
#' @param minVariants integer specifing a minimum number of occurrences of the minor allele for the variant to be included in analysis.  Default is 5, minimum allowed is 2.
#' @param blocksize integer specifying the number of consecutive rows in the data matrix to be considered LD blocks.  One variant will be chosen from each block in the analysis.  Default is NA (no LD pruning, equivalent to blocksize=1)
#' @param simFun function for similarity comparision, such as cor or cov.  Default is null.
#' @param saveDir file to save results output.  Default is no saving (saveDir=NA)
#' @param verbose logical indicating whether to output status updates during analysis run
#'
#' @return List with class "stego" containing
#' \item{summary}{Summary statistics, including p-values, FDR, kinship coefficient estimate between all pairs of individuals}
#' \item{s_matrix_dip}{A matrix of pairwise s statistics between all individuals}
#' \item{s_matrix_hap}{For phased data only, a matrix of pairwise s statistics between all haplotypes}
#' \item{var_s_dip}{numeric estimate of the variance of pairwise subject test statistics}
#' \item{var_s_hap}{numeric estimate of the variance of pairwise haploid test statistics. For phased data only}
#' \item{simMat}{if simFun is used, A similarity matrix between all individuals}
#' \item{analysisType}{character indicating what manner the subjects were grouped in the analysis}
#' \item{pkweightsMean}{numeric value for the whole dataset as a function of the observed allele frequencies}
#' 
#' @examples
#' data(toyGenotypes)
#' sampleNames <- paste("Sample",1:100)
#' 
#' res <- run_stego(toyGenotypes, sampleNames=sampleNames)
#' plot(res, plotname="All Samples")
#' 
#' labels <- paste("Group",c(LETTERS[rep(1:5,20)]))
#' res <- run_stegotoyGenotypes, groups="each.separately", labels=labels)
#' plot(res)
#' 
#' labels <- paste("Group",c(LETTERS[rep(1:5,10)],LETTERS[rep(6:10,10)]))
#' super <- c(rep("Super A",50), rep("Super B",50))
#' res <- run_stego(toyGenotypes, groups="pairwise.within.superpop", labels=labels, super=super)
#' plotFromGSM(res)
#'
#' @author Dan Schlauch \email{dschlauch@fas.harvard.edu}
#' @export
run_stego <- function(genotypes,
                    phased=T,
                    groups="all.together",
                    sampleNames=NULL,
                    labels=NA,
                    super=NA,
                    minVariants=5, 
                    blocksize=NA,
                    simFun=NULL,
                    saveDir=NA,
                    verbose=F,
                    cores=NULL){
    
    
    
    if(verbose){
        print("Preparing analysis...")
    }
    
    if(!is.matrix(genotypes)&&!is.data.frame(genotypes)){
        stop("genotypes must be matrix-like object")
    }
    if(phased==T){
        if(ncol(genotypes)%%2 == 1){
            stop("Odd number of columns for phased data.  Perhaps data is unphased.")
        }
        if(!is.null(sampleNames)&&length(sampleNames)!=ncol(genotypes)/2){
            stop("length of sampleNames must be 1/2 of number of columns of genotypes in phased data")
        }
        if(max(genotypes)!= 1||min(genotypes)!=0){
            stop("Non-binary values for phased data")
        }
    }
    if(phased!=T){
        if(max(genotypes)!= 2||min(genotypes)!=0){
            stop("Unphased data needs to be in [0,1,2]")
        }
        if(!is.null(sampleNames)&&length(sampleNames)!=ncol(genotypes)){
            stop("length of sampleNames must be equal number of columns of genotypes in unphased data")
        }
    }
    
    genotypes = tryCatch({
        as.data.table(genotypes)
    }, error = function(e) {
        stop("couldn't coerce genotypes to data.table")
    }, finally = {
    } )
    
    if (!is.null(sampleNames)){
        if(phased && length(sampleNames)==ncol(genotypes)/2){
            sampleNames <- make.unique(rep(sampleNames,each=2))
        }else{
            # Do nothing, sampleNames is already good
        }
    } else {
        if(is.null(colnames(genotypes))){
            if(phased){
                sampleNames <- make.unique(rep(paste("Sample",1:ncol(genotypes)),each=2))
            }else{
                sampleNames <- paste("Sample",1:ncol(genotypes))
            }
        } else {
            sampleNames <- colnames(genotypes)
        }
    }
    colnames(genotypes) <- sampleNames
    
    if (groups=="all.together"){
        genotypes <- pruneGenotypes(genotypes, blocksize)
        if(verbose){
            print(paste("Running STEGO on", nrow(genotypes), "genes and", length(sampleNames)," samples."))
        }
        
        results <- calculateSMatrix(gt=genotypes, phased=phased, minVariants=5, saveDir=saveDir, simFun=simFun, verbose=verbose)
        results$analysisType <- groups
        results$labels <- labels
        results$super <- super
        results$minVariants <- minVariants
        results$blocksize <- blocksize
        results$sampleNames <- sampleNames
        class(results) <- "stego"
        return(results)
    }
    
    numberOfSamples = ifelse(phased, ncol(genotypes)/2,ncol(genotypes))
    assertthat::assert_that(length(labels)==numberOfSamples)
    assertthat::assert_that(length(unique(labels))>1)
    
    results <- list()
    if (groups=="each.separately"){
        for(subpop in unique(labels)){
            if(verbose){
                print(paste("Computing on population:",subpop))
            }
            
            filter <- labels%in%subpop
            if(phased){
                filter <- rep(filter,each=2)
            }
            gt <- pruneGenotypes(genotypes[,filter,with=F], blocksize)
            results[[subpop]] <- calculateSMatrix(gt=gt, phased=phased, minVariants=minVariants, simFun=simFun, saveDir=saveDir, verbose=verbose)
        }
        names(results) <- unique(labels)
        results$analysisType <- groups
        results$labels <- labels
        results$super <- super
        results$minVariants <- minVariants
        results$blocksize <- blocksize
        results$sampleNames <- sampleNames
        class(results) <- "stego"
        return(results)
    }
    assertthat::assert_that(length(super)==numberOfSamples)
    assertthat::assert_that(length(unique(super))>1)
    
    if (groups=="pairwise.within.superpop"){
        
        
        pairs <- do.call(cbind, lapply(unique(super), function(x){
            combn(unique(labels[super==x]),2)
        }))
        
        
        for(i in seq_len(ncol(pairs))){
            pair <- pairs[,i]
            if(verbose){
                print(paste("Computing on populations:",paste(pair, collapse=",")))
            }
            
            filter <- labels%in%pair
            if(phased){
                filter <- rep(filter,each=2)
            }
            gt <- pruneGenotypes(genotypes[,filter,with=F], blocksize)
            
            if (!is.na(saveDir)){
                saveFile <- paste(pair, collapse="-")
            }
            results[[paste(pair,collapse="_")]] <- calculateSMatrix(gt=gt, phased=phased, minVariants=minVariants, outputDir=outputDir, simFun=simFun, saveDir=saveDir, saveFile=saveFile, verbose=verbose)
            
#                 if(computeFST){
#                     fstData <- as.data.frame(t(as.matrix(as.data.frame(gt)[,c(T,F)] + as.data.frame(gt)[,c(F,T)])))
#                     fstData <- cbind(pop = pop[filterdip], fstData)
#                     results[[paste(pair,collapse="_")]]$FST <- genet.dist(fstData,phased=F,method="Fst")
#                 } 
        }
        results$analysisType <- groups
        results$labels <- labels
        results$super <- super
        results$minVariants <- minVariants
        results$blocksize <- blocksize
        results$sampleNames <- sampleNames
        class(results) <- "stego"
        return(results)
    }
    
    stop("Something went wrong, no output")    
}


kinshipsFromSMatrix <-function(result){
    result$s_matrix_dip[row(result$s_matrix_dip)>=col(result$s_matrix_dip)] <- NA
    sStats <- melt(result$s_matrix_dip, na.rm=T)
    colnames(sStats) <- c("Individual 1","Individual 2", "similarity")
    sStats <- sStats[order(-sStats$similarity),]
    sStats$pValue <- 1-pnorm((sStats$similarity-1)/sd(sStats$similarity))
    sStats$bonferroniPValue <-p.adjust(sStats$pValue, method="bonferroni")
    sStats$FDR <-p.adjust(sStats$pValue, method="BH")
    sStats$kinshipEstimate <- (sStats$similarity-1)/(result$pkweightsMean-1)
    as.data.table(sStats)
}


#' Prune genotypes based on blocks of a specified length
#'
#' @param gt data object containing the genotypes by samples
#'
#' @return Object containing similarity matrix and other results
#'
#' @examples
#' data(toyGenotypes)
#' pruneGenotypes(toyGenotypes)
#'
#' @export
pruneGenotypes <-  function(gt, blocksize=1){
    
    if(is.na(blocksize)||blocksize==1){
        return(gt)
    }
    if(!is.numeric(blocksize)||blocksize<1||blocksize>nrow(gt)||blocksize%%1!=0){
        stop("Block size must be a positive integer less than the number of rows in the genotype data.")
    }
    names(gt) <- make.unique(names(gt))
    numSamples <- ncol(gt)
    numVariants <- nrow(gt)
    sumVariants <- rowSums(gt)
    
    
    # reverse so that MAF<.5
    gt[sumVariants>(numSamples/2),] <- 1-gt[sumVariants>(numSamples/2),]
    sumVariants <- rowSums(gt)
    
    # Intelligently LD prune
    numblocks <- numVariants/blocksize +1
    blocks <- rep(1:numblocks, each=blocksize)[1:numVariants]
    
    system.time(runningWhichMax <- running(sumVariants,width=blocksize,fun=which.max, by=blocksize))
    prunedIndices <- runningWhichMax + seq(0,blocksize*(length(runningWhichMax)-1),blocksize)
    system.time(gt <- gt[prunedIndices,,with=F])
    gt
}

#' Calculate the similarity matrix
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
calculateSMatrix <- function(gt, sampleNames=sampleNames, phased=T, minVariants=5, scaleBySampleAF=F, outputDir=".", saveDir=NA, saveFile=NA, simFun=NULL, verbose=F){
    
    require(data.table)
    numAlleles <- ifelse(phased, ncol(gt), 2*ncol(gt))
    numVariants <- nrow(gt)
    sumVariants <- rowSums(gt)
    
    gt <- as.matrix(gt)
    #     # reverse so that MAF<.5
    invertMinorAllele <- sumVariants>(numAlleles/2)
    if(phased){
        gt[invertMinorAllele,] <- 1-gt[invertMinorAllele,]
    } else {
        gt[invertMinorAllele,] <- 2-gt[invertMinorAllele,]
    }
    # remove < n variants
    sumVariants <- rowSums(gt)
    gt <- gt[sumVariants>minVariants,]
    
    numFilteredVariants <- nrow(gt)
    sumFilteredVariants <- rowSums(gt)
    simMat <- NULL
    if(!is.null(simFun)){
        if(verbose){
            print("Computing supplied similarity function matrix...")
        }
        if(phased){
            simMat <- simFun(gt[,c(T,F)] + gt[,c(F,T)],use="pairwise.complete.obs")
        } else {
            simMat <- simFun(gt, use="pairwise.complete.obs")
        }
    }
    totalPossiblePairs <- choose(numAlleles,2)
    totalPairs <- choose(sumFilteredVariants,2)
    weights <- totalPossiblePairs/totalPairs
    p <- 1/weights
    
    # var_s_hap <- sum((1-p)/p)/(numFilteredVariants^2)
    # recalculated variance 3/23/16
    var_s_hap <- sum(weights-1)/(numFilteredVariants^2)
    
    
    # Calculate expected values conditional on kinship
    pkweightsMean <- mean(((sumFilteredVariants-2)/numAlleles)*weights)
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
    
#     print(mean(s_matrix_hap[row(s_matrix_hap)!=col(s_matrix_hap)]))
#     print(median(s_matrix_hap[row(s_matrix_hap)!=col(s_matrix_hap)]))
    
    estimatedKinship <- (s_matrix_hap-1)/(pkweightsMean-1)
    popResult <- NULL
    # Collapse if phased, rename variables if unphased
    if(phased){
        s_matrix_dip <- (s_matrix_hap[c(T,F),c(T,F)] + s_matrix_hap[c(F,T),c(T,F)] +s_matrix_hap[c(T,F),c(F,T)] + s_matrix_hap[c(F,T),c(F,T)])/4
        colnames(s_matrix_dip) <- colnames(gt)[c(T,F)]
        rownames(s_matrix_dip) <- colnames(gt)[c(T,F)]
        
        var_s_dip <- var_s_hap/4
        popResult <- list(s_matrix_dip=s_matrix_dip, s_matrix_hap=s_matrix_hap, pkweightsMean=pkweightsMean, var_s_dip=var_s_dip, var_s_hap=var_s_hap, simMat=simMat)
        
    } else {
        s_matrix_dip <- s_matrix_hap/4
        colnames(s_matrix_dip) <- colnames(gt)
        rownames(s_matrix_dip) <- colnames(gt)
        
        var_s_dip <- var_s_hap
        rm(s_matrix_hap)
        rm(var_s_hap)
        popResult <- list(s_matrix_dip=s_matrix_dip, pkweightsMean=pkweightsMean, var_s_dip=var_s_dip, simMat=simMat)
        
    }
    popResult$kinships <- kinshipsFromSMatrix(popResult)
    s_vector <- popResult$s_matrix_dip[row(popResult$s_matrix_dip)>col(popResult$s_matrix_dip)]
    popResult$structurePValue <- ks.test((s_vector-1)/sd(s_vector), "pnorm", alternative = c("less"))$p.value
    popResult$crypticPValue <- popResult$kinships$bonferroniPValue[1]
    
    if(!is.na(saveDir)){
        saveRDS(popResult, file=paste0(saveDir,"/",saveFile,".rds"))
    }
    class(popResult) <- "stego"
    popResult
    
}

#' Plot the distribution of similarity statistics
#'
#' 
#'
#' @param genotypes data object containing the phased or unphased genotypes by samples
#' @param plotname title of the plot
#'
#' @return Object containing similarity matrix and other results
#'
#' @examples
#' data(toyGenotypes)
#' 
#' run_stegotoyGenotypes)
#' labels <- c(rep("Group A",100), rep("Group B",100))
#' run_stegotoyGenotypes, groups="each.separately", labels=labels)
#' 
#' labels <- paste("Group",c(LETTERS[rep(1:4,25)],LETTERS[rep(5:8,25)]))
#' super <- c(rep("Super A",100), rep("Super B",100))
#' res <- run_stegotoyGenotypes, groups="pairwise.within.superpop", labels=labels, super=super)
#'
#' @export
plotStegoHist <- function(resObj, plotname="", alphaCutoff=.01){
    if(is.null(resObj$analysisType)||resObj$analysisType=="all.together"){
        gsm <- resObj$s_matrix_dip
        var_s <- resObj$var_s_dip
        pkweightsMean <- resObj$pkweightsMean
        
        #         print(mean(gsm[row(gsm)!=col(gsm)]))
        #         print(median(gsm[row(gsm)!=col(gsm)]))
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
    } else {
        numPlots <- length(resObj)-1
        columns <- floor(sqrt(numPlots))+1
        rows <- ceiling(numPlots/columns)
        #         par(mfrow=c(rows,columns))
        plots <- lapply(seq_len(length(resObj)-1), function(i){
            plotStegoHist(resObj=resObj[[i]], plotname=names(resObj)[i])
        })
        do.call(gridExtra::grid.arrange, c(plots, ncol=columns))
    }
}

#' P-value for test of Population Structure
#' 
#' Get the population structure p-values from a stego object.  This is a test of the null hypotheses that the data is derived a homogeneous population.
#'
#' 
#'
#' @param stego a stego result object
#'
#' @return p-value for test or a list a p-values from each group
#'
#' @examples
#' data(toyGenotypes)
#' 
#' run_stegotoyGenotypes)
#' labels <- c(rep("Group A",100), rep("Group B",100))
#' run_stegotoyGenotypes, groups="each.separately", labels=labels)
#' 
#' labels <- paste("Group",c(LETTERS[rep(1:4,25)],LETTERS[rep(5:8,25)]))
#' super <- c(rep("Super A",100), rep("Super B",100))
#' res <- run_stegotoyGenotypes, groups="pairwise.within.superpop", labels=labels, super=super)
#' structurePValues(res)
#' @export
structurePValues <- function(stego){
    if(is.null(stego$analysisType)||stego$analysisType=="all.together"){
        return(stego$structurePValue)
    } else {
        groups <- names(stego)
        groups <- groups[which(groups!="analysisType")]
        pvals <- lapply(groups, function(x){
            stego[[x]]$structurePValue
        })
        names(pvals) <- groups
        unlist(pvals)
    }
}

#' P-value for test of Cryptic Relatedness
#' 
#' Get the cryptic relatedness p-values from a stego object.  This is a test of the null hypotheses that the data is derived a homogeneous population.
#'
#' 
#'
#' @param stego a stego result object
#'
#' @return p-value for test or a list a p-values from each group
#'
#' @examples
#' data(toyGenotypes)
#' 
#' run_stegotoyGenotypes)
#' labels <- c(rep("Group A",100), rep("Group B",100))
#' run_stegotoyGenotypes, groups="each.separately", labels=labels)
#' 
#' labels <- paste("Group",c(LETTERS[rep(1:4,25)],LETTERS[rep(5:8,25)]))
#' super <- c(rep("Super A",100), rep("Super B",100))
#' res <- run_stegotoyGenotypes, groups="pairwise.within.superpop", labels=labels, super=super)
#' crypticPValues(res)
#' @export
crypticPValues <- function(stego){
    if(is.null(stego$analysisType)||stego$analysisType=="all.together"){
        return(stego$crypticPValue)
    } else {
        groups <- names(stego)
        groups <- groups[which(groups!="analysisType")]
        pvals <- lapply(groups, function(x){
            stego[[x]]$crypticPValue
        })
        names(pvals) <- groups
        unlist(pvals)
    }
}

#' Plot the distribution of similarity statistics
#'
#' 
#'
#' @param genotypes data object containing the phased or unphased genotypes by samples
#' @param plotname title of the plot
#'
#' @return Object containing similarity matrix and other results
#'
#' @examples
#' data(toyGenotypes)
#' 
#' run_stegotoyGenotypes)
#' labels <- c(rep("Group A",100), rep("Group B",100))
#' run_stegotoyGenotypes, groups="each.separately", labels=labels)
#' 
#' labels <- paste("Group",c(LETTERS[rep(1:4,25)],LETTERS[rep(5:8,25)]))
#' super <- c(rep("Super A",100), rep("Super B",100))
#' res <- run_stegotoyGenotypes, groups="pairwise.within.superpop", labels=labels, super=super)
#'
#' @export
plotStegoHeatmap <- function(resObj, plotname="", ...){
    if(is.null(resObj$analysisType)||resObj$analysisType=="all.together"){
        # Needs sample names for colorbars
    } else {
        # Can't plot for each.separately or pairwise...
    }
}

