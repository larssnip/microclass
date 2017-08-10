#' @name blastClassify16S
#' @title Classifying using BLAST
#' 
#' @description A 16S based classification based on BLAST.
#' 
#' @param sequence Character vector of 16S sequences to classify.
#' @param bdb Name of BLAST data base, see \code{\link{blastDbase16S}}.
#' 
#' @details A vector of 16S sequences (DNA) are classified by first using BLAST \code{blastn} against
#' a database of 16S DNA sequences, and then classify according to the nearest-neighbour principle.
#' The nearest neighbour of a query sequence is the hit with the largest bitscore. The BLAST+
#' software  \url{https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download}
#' must be installed on the system. Type \code{system("blastn -help")} in the Console window,
#' and a sensible Help-text should appear.
#' 
#' The database must contain 16S sequences where the Header starts with a token specifying the taxon. 
#' More specifically, the tokens must look like:
#' 
#' <taxon>_1
#' 
#' <taxon>_2
#' 
#' ...etc
#' 
#' where <taxon> is some proper taxon name. Use \code{\link{blastDbase16S}} to make such databases.
#' 
#' The identity of each alignment is also computed. This should be close to 1.0 for a classification
#' to be trusted. Identity values below 0.95 could indicate uncertain classifications, but this will
#' vary between taxa.
#' 
#' @return A \code{data.frame} with two columns: Taxon is the predicted taxon for each \code{sequence}
#' and Identity is the corresponding identity-value. If no BLAST hit is seen, the sequence is 
#' \code{"unclassified"}.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{blastDbase16S}}.
#' 
#' @examples 
#' data("small.16S")
#' seq <- small.16S$Sequence
#' tax <- sapply(strsplit(small.16S$Header,split=" "),function(x){x[2]})
#' \dontrun{
#' dbase <- blastDbase16S("test",seq,tax)
#' reads <- amplicon(seq)
#' predicted <- blastClassify16S(reads[nchar(reads)>0],dbase)
#' print(predicted)
#' }
#' 
#' @importFrom utils read.table
#' @export blastClassify16S
#' 
blastClassify16S <- function( sequence, bdb ){
  n <- length( sequence )
  tags <- paste( "Query", 1:n, sep="_" )
  fdta <- data.frame( Header=tags,
                      Sequence=sequence,
                      stringsAsFactors=F )
  writeFasta( fdta, out.file="query.fasta" )
  cmd <- paste( "blastn -query query.fasta -db ", bdb, " -num_alignments 1",
                " -out bres.txt -outfmt \"6 qseqid qlen sseqid length pident bitscore\"",
                sep="")
  system( cmd )
  btab <- read.table( "bres.txt", sep="\t", header=F, stringsAsFactors=F )
  file.remove( c( "query.fasta", "bres.txt" ) )
  btab <- btab[order( btab[,6], decreasing=T ),]
  btab <- btab[which( !duplicated( btab[,1] ) ),]
  tax.hat <- gsub( "_[0-9]+$", "", btab[,3] )
  idty <- (btab[,5]/100) * btab[,4]/btab[,2] + pmax(0,btab[,2]-btab[,4])*0.25
  
  taxon.hat <- rep( "unclassified", n )
  identity <- rep( 0, n )
  idx <- match( btab[,1], tags )
  taxon.hat[idx] <- tax.hat
  identity[idx] <- idty
  
  return( data.frame( Taxon=taxon.hat, Identity=identity, stringsAsFactors=F ) )
}



#' @name blastDbase16S
#' @title Building a BLAST database
#' 
#' @description Building a BLAST database for 16S based classification.
#' 
#' @param name The name of the database (text).
#' @param sequence A character vector with 16S sequence data.
#' @param taxon A character vector with taxon information.
#' 
#' @details This functions builds a database using the \code{makeblastdb} program of the
#' BLAST+ software \url{https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download}.
#' Thus, this software must be available on the system when using this
#' function. If you type \code{system("makeblastdb -help")} in the Console window some meaningful
#' Help-text should be displayed.
#' 
#' This function is most typically used prior to \code{\link{blastClassify16S}} to set up the database
#' before searching and classifying. It can be seen as the 'training step' of a BLAST-based
#' classification procedure.
#' 
#' The \code{sequence} must be a vector of DNA-sequences (16S sequences). The \code{taxon} is a vector of the same 
#' length as \code{sequence}, containing the correpsonding taxon information. 
#' 
#' @return The database files are created, and the name of the database (\code{name}) is returned.
#' 
#' @author Lars Snipen.
#' 
#' @seealso \code{\link{blastClassify16S}}.
#' 
#' @examples # See examples for blastClassify16S.
#' 
#' @export blastDbase16S
#' 
blastDbase16S <- function( name, sequence, taxon ){
  n <- length( taxon )
  fdta <- data.frame( Header=paste( taxon, 1:n, sep="_" ),
                      Sequence=gsub( "[^ACGTNRYSWKMBDHV]", "N", gsub( "U", "T", toupper(sequence) ) ),
                      stringsAsFactors=F )
  writeFasta(fdta,out.file="dbase.fasta")
  cmd <- paste( "makeblastdb -dbtype nucl -in dbase.fasta -out", name )
  system( cmd )
  file.remove( "dbase.fasta" )
  return( name )
}


