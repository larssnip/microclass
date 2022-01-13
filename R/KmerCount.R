#' @name KmerCount
#' @aliases KmerCount
#' @title K-mer counting
#' 
#' @description Counting overlapping words of length K in DNA/RNA sequences.
#' 
#' @param sequences Vector of sequences (text).
#' @param K Word length (integer >= 0, see Details).
#' @param col.names Logical indicating if the words should be added as columns names.
#' @param type Character selecting type of K-mer: "nucleotide" (default), 
#' "amino acid", "codon", "nucleotide x" (see Details).
#' 
#' @details For each input sequence, the frequency of every word of length \code{K} is counted. 
#' Counting is done with overlap by default. The counting itself is done by a C++ function.
#' \code{K=0} is interpreted as AT/GC.
#' 
#' With \code{col.names = TRUE} the K-mers are added as column names, but this makes the
#' computations slower.
#' 
#' Counting of nucleotide K-mers can be restricted to a particular reading frame 
#' by setting \code{type = "nucleotide 1"}, \code{type = "nucleotide 2"} or \code{type = "nucleotide 3"}, 
#' i.e. leading to a skip of 3 instead of 1 when moving along the sequence.
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
#' # AT/GC
#' KmerCount("ATGCCTGAACTGACCTGC", K = 0, col.names = TRUE)
#' 
#' # In-frame counting
#' KmerCount("ATGCCTGAACTGACCTGC", K = 2, col.names = TRUE, type = "nucleotide 2") # 2-mer, 2nd reading frame
#' KmerCount("ATGCCTGAACTGACCTGC", K = 0, col.names = TRUE, type = "nucleotide 3") # AT/GC, 3rd reading frame
#' 
#' @export KmerCount
#' 
KmerCount <- function(sequences, K = 1, col.names = FALSE, type = "nucleotide"){
  if(type == "nucleotide"){
    if(K > 15)
      stop('Maximum value for K = 15')
    if(K==0){
        int.list <- charToIntGC(sequences)
        X <- Kmer_countGC(int.list, 1, col.names, FALSE)
    } else {
      int.list <- charToInt(sequences)
      X <- Kmer_count(int.list, K, col.names, FALSE)
    }
  }
  if(type %in% c("nucleotide 1","nucleotide 2","nucleotide 3")){
    if(K > 15)
      stop('Maximum value for K = 15')
    frame <- as.integer(strsplit(type, " ")[[1]][2])
    if(K==0){
      int.list <- charToIntGC(sequences)
    } else {
      int.list <- charToInt(sequences)
    }
    if(frame > 1)
      int.list <- lapply(int.list, function(x)x[-(1:(frame-1))])
    if(K==0){
      X <- Kmer_countGC3(int.list, 1, col.names, TRUE)
    } else {
      X <- Kmer_count3(int.list, K, col.names, TRUE)
    }
  }
  if(type == "amino acid"){
    int.list <- charToIntAminoAcid(translate(sequences))
    if(K > 6)
      stop('Maximum value for K = 6')
    X <- Kmer_count_amino_acid(int.list, K, col.names)
  }
  if(type == "codon"){
    int.list <- charToIntCodon(translate(sequences, codon=TRUE))
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

