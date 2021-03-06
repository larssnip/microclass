% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multinom.R
\name{multinomTrain}
\alias{multinomTrain}
\title{Training multinomial model}
\usage{
multinomTrain(sequence, taxon, K = 8, col.names = FALSE, n.pseudo = 100)
}
\arguments{
\item{sequence}{Character vector of 16S sequences.}

\item{taxon}{Character vector of taxon labels for each sequence.}

\item{K}{Word length (integer).}

\item{col.names}{Logical indicating if column names should be added to the trained model matrix.}

\item{n.pseudo}{Number of pseudo-counts to use (positive numerics, need not be integer). Special case -1
will only return word counts, not log-probabilities.}
}
\value{
A list with two elements. The first element is \code{Method}, which is the text 
\code{"multinom"} in this case. The second element is \code{Fitted}, which is a matrix
of probabilities with one row for each unique \code{taxon} and one column for each possible word of
length\code{K}. The sum of each row is 1.0. No probabilities are 0 if \code{n.pseudo}>0.0.

The matrix \code{Fitted} has an attribute \code{attr("prior",)}, that contains the relative
taxon sizes.
}
\description{
Training the multinomial K-mer method on sequence data.
}
\details{
The training step of the multinomial method (Vinje et al, 2015) means counting K-mers
on all sequences and compute the multinomial probabilities for each K-mer for each unique taxon. 
\code{n.pseudo} pseudo-counts are added, divided equally over all K-mers, before probabilities
are estimated. The optimal choice of \code{n.pseudo} will depend on \code{K} and the 
training data set. The default value \code{n.pseudo=100} has proven good for \code{K=8} and the
\code{\link[microcontax]{contax.trim}} data set (see the \code{microcontax} R-package).

Adding the actual K-mers as column names (\code{col.names=TRUE}) will slow down the
computations.

The relative taxon sizes are also computed, and may be used as an empirical prior in the
classification step (see "prior" below).
}
\examples{
# See examples for multinomClassify

}
\references{
Vinje, H, Liland, KH, Almøy, T, Snipen, L. (2015). Comparing K-mer based methods for
improved classification of 16S sequences. BMC Bioinformatics, 16:205.
}
\seealso{
\code{\link{KmerCount}}, \code{\link{multinomClassify}}.
}
\author{
Kristian Hovde Liland and Lars Snipen.
}
