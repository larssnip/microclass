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
#' transform (Aitchison 1986) for compositional data. Readcount data can be seen as
#' compositional data since the total number of readcounts in a sample does not 
#' carry any information about the biology, but is simply an effect of sequencing
#' depth. Thus, the information in the data lies in the relative values, not the
#' absolute. By transforming such data with this function, you get data who are 
#' better suited for a number of downstream analyses, e.g. typically analyses 
#' making use of sum-of-squares type of statistics, like PCA, PLS, ANOVA or 
#' clustering with euclidean distances.
#' 
#' The \code{readcount.mat} must have the samples in the rows and the taxa in the
#' columns. Transpose if necessary.
#' 
#' The transform does not accept zeros in any cell of the \code{readcount.mat}. To
#' cope with this you add pseudo-counts. Bu default 1 additional read is assigned to 
#' all cells in \code{readcount.mat}. You may change this value, and it need not
#' be an integer. The rationale behind this is that we a priori assume a uniform
#' distribution of the taxa, and the more pseudo-counts you add, the more weight 
#' you give to this prior.
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
  lc <- log2(clr.mat)
  clr.mat <- lc - rowMeans(lc)
  attr(clr.mat, "library.size") <- library.size
  return(clr.mat)
}
