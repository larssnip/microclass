#' @name clr
#' @title Centred log-ratio transform
#' 
#' @description Transforms readcount-matrix with Aitchisons transform.
#' 
#' @usage clr(readcount.mat, n.pseudo = 1)
#' 
#' @param readcount.mat matrix with readcount data.
#' @param n.pseudo number of pseudo-readcounts to add.
#' 
#' @details This is a standard implementation of the Aitchisons centered log-ratio
#' transform (Aitchison 1986) for compositional data. Readcount data can be seen as compositional data
#' since the total number of readcounts in a sample does not carry any information, but is simply
#' an effect of sequencing depth. Thus, the information in the data lies in the relative 
#' values, not the absolute.
#' 
#' The \code{readcount.mat} must have the samples in the rows and the taxa in the
#' columns. Transpose if necessary.
#' 
#' By transforming such data with this function, you get data who are better suited for a
#' number of downstream analyses, e.g. typically analyses making use of sum-of-squares type of
#' statistics, like PCA, PLS, ANOVA or clustering with euclidean distances.
#' 
#' @return A matrix of same size as the input, but with transformed readcounts.
#' 
#' @author Lars Snipen.
#' 
#' @references Aitchison J. The Statistical Analysis of Compositional Data. London, UK: Chapman & Hall; 1986.
#' 
#' @examples
#' 
#' @export clr
#' 
clr <- function(readcount.mat, n.pseudo = 1){
  clr.mat <- readcount.mat + n.pseudo
  library.size <- rowSums(clr.mat)
  lc <- log(clr.mat)
  clr.mat <- lc - rowMeans(lc)
  attr(clr.mat, "library.size") <- library.size
  return(clr.mat)
}

#' @name rarefy
#' @title Rarefying readcounts
#' 
#' @description Down-sampling readcounts to equal depths.
#' 
#' @usage rarefy(readcount.mat, depth = NULL)
#' 
#' @param readcount.mat matrix with readcount data.
#' @param depth fixed number of reads in each sample after.
#' 
#' @details To rarefy readcounts means down-sampling the data in all samples to 
#' a fixed target read depth. This is sometimes done prior to some diversity analyses,
#' the idea being that some taxa are absent due to low depth, and this will bias
#' diversity estimates. After the use \code{rarefy} all samples have equal 
#' (low) depth.
#' 
#' If \code{depth = NULL} the smallest readcount among the samples is used as the
#' target \code{depth}, but the user may supply any other positive integer value.
#' Note: If \code{depth} is larger than the minimum readcount, some samples
#' will be up-sampled. This means they will get more reads for the taxa where they
#' already have reads, but still no reads for the absent taxa. 
#' 
#' The \code{readcount.mat} must have the samples in the rows and the taxa in the
#' columns. Transpose if necessary.
#' 
#' This function minimizes the variance in the down-sampling as follows:
#' For sample \code{i}, we first compute the expected readcount given target 
#' \code{depth}. Let
#' \code{Ecount = depth * read.counts[i,]/sum(read.counts[i,])}. From this we get
#' the base count \code{base.count = floor(Ecount)}. The remainder is 
#' \code{remainder = Ecount - base.count}. The remaining reads 
#' \code{R = depth - sum(base.count)} are finally distributed over 
#' the taxa using the \code{remainder} as the probabilities. Thus, only the 
#' latter (few) reads will vary randomly between repeated use of this function
#' on the same data.
#' 
#' @return A matrix of same size as the input, but with down-sampled readcounts.
#' This matrix has an attribute \code{"original.depth"} with the original depths
#' for each sample.
#' 
#' @author Lars Snipen.
#' 
#' @examples
#' 
#' @export rarefy
#' 
rarefy <- function(read.counts, depth = NULL){
  N <- nrow(read.counts)
  P <- ncol(read.counts)
  samp.depths <- rowSums(read.counts)
  if(is.null(depth)){
    depth <- min(samp.depths)
  }
  X <- matrix(0, nrow = N, ncol = P)
  taxa <- colnames(read.counts)
  samples <- rownames(read.counts)
  colnames(X) <- taxa
  rownames(X) <- samples
  cfac <- factor(1:P)
  for(ss in 1:nrow(read.counts)){
    Ecounts <- depth * read.counts[ss,] / samp.depths[ss]
    counts <- floor(Ecounts)
    remainder <- Ecounts - counts
    if(max(remainder) == 0){
      rest <- rep(0, P)
    } else {
      rest <- table(sample(cfac, size = depth - sum(counts), replace = T, prob = remainder))
    }
    X[ss,] <- counts + rest
  }
  attr(X, "original.depth") <- samp.depths
  return(X)
}







# # The matrix read.counts must have
# # - samples in the rows
# # - clusters in the columns
# quantileNormalize <- function(read.counts){
#   X.rank <- apply(read.counts, 2, rank, ties.method = "min")
#   X.sorted <- apply(read.counts, 2, sort)
#   x.mean <- rowMeans(X.sorted)
# 
#   i2m <- function(index, mean){
#     return(mean[index])
#   }
# 
#   X.final <- apply(X.rank, 2, function(idx){x.mean[idx]})
#   rownames(X.final) <- rownames(read.counts)
#   return(X.final)
# }
# 
# # The matrix read.counts must have
# # - samples in the rows
# # - clusters in the columns
# rarE <- function(read.counts, depth = NULL){
#   if(is.null(depth)){
#     depth <- min(rowSums(read.counts))
#   }
#   X <- round(depth * microCLR(read.counts, C = 0))
#   return(X)
# }
