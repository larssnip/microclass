#' @name multinomTrain
#' @title Training multinomial model
#' 
#' @description Training the multinomial K-mer method on sequence data.
#' 
#' @param sequence Character vector of sequences.
#' @param taxon Character vector of taxon labels for each sequence.
#' @param K Word length (integer).
#' @param col.names Logical indicating if column names (K-mers) should be added to the trained model matrix.
#' @param n.pseudo Number of pseudo-counts to use (positive numerics, need not be integer). Special case -1
#' will only return word counts, not log-probabilities.
#' 
#' @details The training step of the multinomial method (Vinje et al, 2015) means counting K-mers
#' on all sequences and compute their multinomial probabilities for each taxon. 
#' \code{n.pseudo} pseudo-counts are added equally to all K-mers, before probabilities
#' are estimated. The optimal choice of \code{n.pseudo} will depend on \code{K} and the 
#' training data set.
#' 
#' Adding the actual K-mers as column names (\code{col.names = TRUE}) will slow down the
#' computations.
#' 
#' The relative taxon frequencies in the \code{taxon} input are also computed and
#' returned as an attribute to the probability matrix.
#' 
#' @return A matrix with the multinomial probabilities, one row for each
#' \code{taxon} and one column for each K-mer. The sum of each row is 1.0. No
#' probabilities are 0 if \code{n.pseudo} > 0.0.
#' 
#' The matrix has an attribute \code{attr("prior",)}, that contains the relative
#' taxon frequencies.
#' 
#' @author Kristian Hovde Liland and Lars Snipen.
#' 
#' @references Vinje, H, Liland, KH, Almøy, T, Snipen, L. (2015). Comparing K-mer based methods for
#' improved classification of 16S sequences. BMC Bioinformatics, 16:205.
#' 
#' @seealso \code{\link{KmerCount}}, \code{\link{multinomClassify}}.
#' 
#' @examples # See examples for multinomClassify
#' 
#' @export multinomTrain
#' 
multinomTrain <- function(sequence, taxon, K = 5, col.names = FALSE, n.pseudo = 1.0){
  tax.int <- taxon
  if(is.character(taxon)){
    tax.int  <- factor(taxon)
    tax.levels <- levels(tax.int)
    tax.int  <- as.integer(tax.int)
    prior <- as.numeric(table(taxon) / length(taxon))
  }
  classes.in <- lapply(1:max(tax.int), function(i)which(i == tax.int))
  multinom.prob <- multinomTrainCpp(charToInt(sequence), K, col.names, classes.in, -1)
  if(n.pseudo >= 0){
    multinom.prob <- CountsToMultinom(multinom.prob, n.pseudo)
  }
  if(is.character(taxon)){
    dimnames(multinom.prob) <- list(tax.levels, NULL) # Avoids copying
  }
  attr(multinom.prob, "prior") <- prior
  return(multinom.prob)
}



#' @name multinomClassify
#' @title Classifying with a Multinomial model
#' 
#' @description Classifying sequences by a trained Multinomial model.
#' 
#' @param sequence Character vector of sequences to classify.
#' @param multinom.prob A matrix of multinomial probabilities, see \code{\link{multinomTrain}}.
#' @param post.prob Logical indicating if posterior log-probabilities should be returned.
#' @param prior Logical indicating if classification should be done by flat priors (default)
#' or with empirical priors.
#' @param full.post.prob Logical indicating if full posterior probability matrix should be returned.
#' 
#' @details The classification step of the multinomial method (Vinje et al, 2015) means counting 
#' K-mers on all sequences, and computing the posterior probabilities for each
#' taxon given the trained model. The predicted taxon for each input sequence is
#' the one with the maximum posterior probability for that sequence.
#' 
#' By setting \code{post.prob = TRUE} you will get the log-probability of the
#' best and second best taxon for each sequence. This may be used for evaluating
#' the certainty in the classifications.
#' 
#' The classification is parallelized through RcppParallel
#' employing Intel TBB and TinyThread. By default all available
#' processing cores are used. This can be changed using the
#' function \code{\link{setParallel}}.
#' 
#' @return If \code{post.prob = FALSE} a character vector of predicted taxa is returned.
#' 
#' If \code{post.prob = TRUE} a \code{data.frame} with three columns is returned.
#' \itemize{
#'   \item taxon. The predicted taxa, one for each sequence in \code{sequence}.
#'   \item post_prob. The posterior log-probability of the assigned taxon.
#'   \item post_prob_2. The largest posterior log-probability of the other taxa.
#' }
#' 
#' @author Kristian Hovde Liland and Lars Snipen.
#' 
#' @references Vinje, H, Liland, KH, Almøy, T, Snipen, L. (2015). Comparing K-mer based methods for
#' improved classification of 16S sequences. BMC Bioinformatics, 16:205.
#' 
#' @seealso \code{\link{KmerCount}}, \code{\link{multinomTrain}}.
#' 
#' @examples 
#' data("small.16S")
#' seq <- small.16S$Sequence
#' tax <- sapply(strsplit(small.16S$Header,split=" "),function(x){x[2]})
#' \dontrun{
#' trn <- multinomTrain(seq,tax)
#' primer.515f <- "GTGYCAGCMGCCGCGGTAA"
#' primer.806rB <- "GGACTACNVGGGTWTCTAAT"
#' reads <- amplicon(seq, primer.515f, primer.806rB)
#' predicted <- multinomClassify(unlist(reads[nchar(reads)>0]),trn)
#' print(predicted)
#' }
#' 
#' @export multinomClassify
#' 
multinomClassify <- function(sequence, multinom.prob, post.prob = FALSE, prior = FALSE, full.post.prob = FALSE){
  int.list <- charToInt(sequence)
  if(prior){
    priors <- log2(attr(multinom.prob, "prior"))
  } else {
    priors <- rep(0, nrow(multinom.prob))
  }
  X <- multinomClassifyCpp(int.list, log(ncol(multinom.prob), 4), multinom.prob, priors, post.prob, full.post.prob)
  if(post.prob){
    return(data.frame(taxon = rownames(multinom.prob)[X$first_ind], post_prob = X$first,
                      post_prob_2 = X$second, stringsAsFactors = FALSE))
  } else if (full.post.prob){
    logprobs <- X$ProbMat
    const <- 700-apply(logprobs,1,max)
    adjusted <- logprobs+const # add a constant (to be canceled out) to log-posterior to avoid numerical issues.
    adjusted[adjusted < -700] = -700
    probs <- (exp(adjusted))/(rowSums(exp(adjusted)))
    probs[probs < 1e-10] = 0
    probs[probs > 1-1e-10] = 1
    probsmat <- Matrix::Matrix(probs,sparse=TRUE)
    return(probsmat)
  } else {
    return(rownames(multinom.prob)[X$first_ind])
  }
}


#' Set number of parallel threads
#' 
#' @description Simple function to set the number of threads to use in parallel
#' computations. The default equals all available logical cores. An integer is
#' interpreted as the number of threads. A numeric < 1 is interpreted as a proportion
#' of the avialable logical cores.
#'
#' @param C a scalar indicating the number of threads, default = NULL (#available logical cores)
#'
#' @return NULL, returned silently.
#' @export
#'
#' @examples
#' \dontrun{
#' setParallel() # Use all available logical cores.
#' }
#' 
#' @importFrom RcppParallel RcppParallelLibs setThreadOptions defaultNumThreads
setParallel <- function(C = NULL){
  if(is.null(C)){ # Reset to default
    setThreadOptions(numThreads = defaultNumThreads())
  } else {
    if(C < 1) # Fraction of available cores
      setThreadOptions(numThreads = max(1, round(defaultNumThreads()*C )))
    if(C >= 1)
      setThreadOptions(numThreads = C)
  }
  invisible(NULL)
}

# Unexported functions for transforming sparse count matrices to
# multinomial classification matrices
CountsToRDP <- function(X, sizes){
  log2(X + rep(attr(X, 'p') + 0.5, each = nrow(X)) / (sum(sizes) + 1)) - log2(sizes + 1)
}
CountsToMultinom <- function(X, n.pseudo){
  log2(X + n.pseudo / ncol(X)) - log2(rowSums(X) + n.pseudo)
}
CountsToMultinomMulti <- function(X, n.pseudo, rowsums){
  log2(X + n.pseudo / ncol(X)) - log2(rowsums + n.pseudo)
}

