#' @name KmerCount
#' @aliases KmerCount
#' @title K-mer counting
#' 
#' @description Counting overlapping words of length K in DNA/RNA sequences.
#' 
#' @param sequences Vector of sequences (text).
#' @param K Word length (integer).
#' @param col.names Logical indicating if the words should be added as columns names.
#' @param codon Logical indicating if K-mers should start at the beginning of a codon (default FALSE, i.e. all K-mers).
#' 
#' @details For each input sequence, the frequency of every word of length \code{K} is counted. 
#' Counting is done with overlap. The counting itself is done by a C++ function.
#' 
#' With \code{col.names = TRUE} the K-mers are added as column names, but this makes the
#' computations slower.
#' 
#' Counting of K-mers can also be restricted to the reading frame starting at the first nucleotide
#' by setting \code{codon = TRUE}, i.e. leading to a skip of 3 instead of 1 when moving along the sequence.
#' 
#' @return A matrix with one row for each sequence in \code{sequences} and one column for 
#' each possible word of length\code{K}.
#' 
#' @author Kristian Hovde Liland and Lars Snipen.
#' 
#' @seealso \code{\link{multinomTrain}}, \code{\link{multinomClassify}}.
#' 
#' @examples 
#' KmerCount("ATGCCTGAACTGACCTGC", K = 2)
#' 
#' @export KmerCount
#' 
KmerCount <- function(sequences, K = 1, col.names = FALSE, codon = FALSE){
  int.list <- charToInt(sequences)
  X <- Kmer_count(int.list, K, col.names, codon)
  rownames(X) <- names(sequences)
  return(X)
}

# #' @rdname KmerCount
# #' @export
# KmerCountParallel <- function( sequences, K = 1, col.names=FALSE ){
#   int.list <- charToInt( sequences )
#   X <- Kmer_parallel( int.list, K, col.names )
#   rownames( X ) <- names( sequences )
#   return(X)
#}

