\name{seqbias.predict}
\alias{seqbias.predict}
\title{Predicting sequencing bias}
\description{Predicts sequencing bias given a fit seqbias model}
\usage{seqbias.predict(sb, I)}
\arguments{
    \item{sb}{a seqbias object}
    \item{I}{a GRanges object}
}
\details{
    Once a seqbias model is fit with 'seqbias.fit', the sequencing bias of any
    region in the reference sequence can be predicted using this function. Given
    the coordinates of a region, this function produces a vector of the same
    length as the sequence. Each position 'i' is given a sequence score 'v_i'.

    A simple procedure is then to normalize read counts given the sequencing
    bias.  The read count of (i.e. the number of reads beginning on) position 'i',
    denoted by 'x_i', can be normalized by computing 'x_i/v_i', giving an
    estimate of abundance that is more accurate in expectation.
}
\value{
    A list of numeric vectors are returned, one for each genomic interval in I.
    The vectors are of equal length to the interval given, and the predicted
    sequencing bias is given for each position.
}
\author{
    Daniel Jones
    \email{dcjones@cs.washington.edu}}
\seealso{
    \code{\link{seqbias.fit}}
}
\examples{
  reads_fn <- system.file( "extra/example.bam", package = "seqbias" )
  ref_fn <- system.file( "extra/example.fa", package = "seqbias" )

  sb <- seqbias.fit( ref_fn, reads_fn )

  I <- GRanges( c('seq1'), IRanges( c(1), c(5000) ), strand = c('-') )

  bias <- seqbias.predict( sb, I )
}

