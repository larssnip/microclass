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
#' The \code{readcount.mat} must have the samples in the rows and the taxa (or genes) in the
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
# 
# # The matrix read.counts must have
# # - samples in the rows
# # - clusters in the columns
# rarify <- function(read.counts, depth = NULL){
#   N <- nrow(read.counts)
#   P <- ncol(read.counts)
#   samp.depths <- floor(rowSums(read.counts)) - 1
#   if(is.null(depth)){
#     depth <- min(samp.depths)
#   }
#   X <- matrix(NA, nrow = N, ncol = P)
#   clst <- colnames(read.counts)
#   samp <- rownames(read.counts)
#   colnames(X) <- clst
#   rownames(X) <- samp
#   for(ss in 1:nrow(read.counts)){
#     if(samp.depths[ss] >= depth){
#       svec <- unlist(lapply(1:P, function(i){rep(clst[i], read.counts[ss,i])}))
#       reads <- sample(svec, size = depth, replace = F)
#       X[ss,] <- table(factor(reads, levels = clst))
#     }
#   }
#   return(X)
# }
# 
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
