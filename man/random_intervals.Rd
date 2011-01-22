\name{random.intervals}
\alias{random.intervals}
\title{Generating random genomic intervals}
\description{Given a vector of sequence lengths, generate genomic intervals uniformly at
random}
\usage{random.intervals(I, n=1, ms=10000)}
\arguments{
  \item{I}{a GRanges object giving intervals from which to sample from}
  \item{n}{number of intervals to generate}
  \item{ms}{length of the intervals to generate (may be a vector)}
}
\details{
  The function is used to place intervals of fixed sizes at random (possibly
  overlapping) positions across one or more sequences. The input should be a
  \code{GRanges} objects giving the sequence intervals in which the random intervals
  sholud be placed. If they are to be placed anywhere within a reference
  sequence, use the \code{scanFaIndex} function from Rsamtools, to obtain a set
  of intervals.
}
\value{
  Returns a GRanges object giving the generated intervals.
}
\author{
    Daniel Jones
    \email{dcjones@cs.washington.edu}
}
\examples{
  library(Rsamtools)
  ref_fn <- system.file( "extra/example.fa", package = "seqbias" )
  ref_f <- FaFile( ref_fn )
  open.FaFile( ref_f )

  ref_seqs <- scanFaIndex( ref_f )

  I <- random.intervals( ref_seqs, n = 100, ms = 1000 )
}

