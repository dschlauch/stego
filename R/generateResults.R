generateToyResults <- function(){
    data(toyGenotypes)
    toyResults$toyAllTogether <-  stego(toyGenotypes)
    toyResults$toyUnphasedAllTogether <-   stego(toyGenotypes[,rep(c(T,F),100),,with=F]+toyGenotypes[,rep(c(F,T),100),with=F], groups="all.together", phased=F)
    
    labels <- paste("Group",c(LETTERS[rep(1:5,20)]))
    toyResults$toyEachSeparately <- stego(toyGenotypes, groups="each.separately", labels=labels)
    
    labels <- paste("Group",c(LETTERS[rep(1:5,10)],LETTERS[rep(6:10,10)]))
    super <- c(rep("Super A",50), rep("Super B",50))
    toyResults$toyPairwiseWithinSuperpop <- stego(toyGenotypes, groups="pairwise.within.superpop", labels=labels, super=super)

    save(toyResults, file="./data/toyResults.rda")
}