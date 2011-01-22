\name{count.reads}
\alias{count.reads}
\title{Counting reads across intervals}
\description{Counts the number of reads starting at each position across given
genomic intervals}
\usage{count.reads(reads_fn, I, binary=TRUE)}
\arguments{
  \item{reads_fn}{filename of aligned reads in BAM format}
  \item{I}{a GRanges object giving valid genomic intervals}
  \item{binary}{if \code{TRUE}, return a 0-1 vector, otherwise return a vector
  counting the number of reads mapped to each position}
}
\details{
  Given an indexed BAM file, this function counts the number of reads starting
  at each position of each provided interval. These counts can then be
  normalized by dividing by the predicted bias obtained from 'bias.predict'.

  By default, a 0-1 vector is returned, where positions at which no reads are
  mapped are 0, and those with one or more are 1. If \code{binary} is
  \code{FALSE}, the number of reads mapping to each position is returned.
}
\value{
  A list of numeric vectors is returned, one for each interval provided. Each
  vector gives an integer count of the number of reads beginning on that
  position.
}
\note{
  The BAM file provided should be indexed with 'samtools index'.
}
\author{
    Daniel Jones
    \email{dcjones@cs.washington.edu}
}
\seealso{
    \code{\link{seqbias.predict}}
}
\examples{
  reads_fn <- system.file( "extra/example.bam", package = "seqbias" )

  I <- GRanges( c('seq1'), IRanges( c(1), c(5000) ), strand = c('-') )

  counts <- count.reads( reads_fn, I )
}


