#' @name KmerCount
#' @aliases KmerCount
#' @title K-mer counting
#' 
#' @description Counting overlapping words of length K in DNA/RNA sequences.
#' 
#' @param sequences Vector of sequences (text).
#' @param K Word length (integer).
#' @param col.names Logical indicating if the words should be added as columns names.
#' @param type Character selecting type of K-mer: "nucleotide" (default), "nucleotide 3" (skip 3), 
#' "amino acid", "codon".
#' 
#' @details For each input sequence, the frequency of every word of length \code{K} is counted. 
#' Counting is done with overlap. The counting itself is done by a C++ function.
#' 
#' With \code{col.names = TRUE} the K-mers are added as column names, but this makes the
#' computations slower.
#' 
#' Counting of nucleotide K-mers can be restricted to the reading frame starting at the first nucleotide
#' by setting \code{type = "nucleotide 3"}, i.e. leading to a skip of 3 instead of 1 when moving 
#' along the sequence.
#' 
#' @return A matrix with one row for each sequence in \code{sequences} and one column for 
#' each possible word of length\code{K}.
#' 
#' @author Kristian Hovde Liland and Lars Snipen.
#' 
#' @seealso \code{\link{multinomTrain}}, \code{\link{multinomClassify}}.
#' 
#' @examples 
#' # Nucleotides
#' KmerCount("ATGCCTGAACTGACCTGC", K = 2)
#' KmerCount("ATGCCTGAACTGACCTGC", K = 3, col.names = TRUE)
#' 
#' # Amino Acids
#' KmerCount("ATGCCTGAACTGACCTGC", K = 1, col.names = TRUE, type = "amino acid")
#' 
#' # Codons
#' length(KmerCount("ATGCCTGAACTGACCTGC", K = 3, col.names = TRUE, type = "codon"))
#' # 64^3 = 262144 for word length 3
#' 
#' @export KmerCount
#' 
KmerCount <- function(sequences, K = 1, col.names = FALSE, type = "nucleotide"){
  if(type == "nucleotide"){
    int.list <- charToInt(sequences)
    if(K > 15)
      stop('Maximum value for K = 15')
    X <- Kmer_count(int.list, K, col.names, FALSE)
  }
  if(type == "nucleotide 3"){
    int.list <- charToInt(sequences)
    if(K > 15)
      stop('Maximum value for K = 15')
    X <- Kmer_count(int.list, K, col.names, TRUE)
  }
  if(type == "amino acid"){
    int.list <- charToIntAminoAcid(translate(sequences))
    if(K > 6)
      stop('Maximum value for K = 6')
    X <- Kmer_count_amino_acid(int.list, K, col.names)
  }
  if(type == "codon"){
    int.list <- charToIntCodon(translate(sequences))
    if(K > 5)
      stop('Maximum value for K = 5')
    X <- Kmer_count_codon(int.list, K, col.names)
  }
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

