\name{kmer.freq}
\alias{kmer.freq}
\title{Measuring positional kmer frequencies}
\description{
  Given a sample of sequences and corresponding read counts, produce a table
  giving the position kmer frequencies relative to read starts
}
\usage{kmer.freq(seqs, counts, L = 50, R = 50, k = 1)}
\arguments{
  \item{seqs}{a list of DNAString objects.}
  \item{counts}{a list of numeric vectors.}
  \item{L}{how many positions to the left of the read start to consider}
  \item{R}{how many positions to the right of the read start to consider}
  \item{k}{the size of each kmer}
}
\details{
  Sequences and read counts are used to produce a table of aggregate kmer
  frequencies for each position relative to the read start. The position on
  which the read starts is numbered 0, positions to the left of the read are
  negative, and those to the right are positive.

  The sequences and counts can be generated with the provided functions
  \code{scanFa} and \code{\link{count.reads}}, respectively. The reverse
  complement of sequences on the negative strand obtained from \code{scanFa}
  should be used. To properly visualize bias a relatively large random sample of
  intervals should be generated.
}
\value{
  A data frame is returned with columns \code{pos}, \code{seq}, and \code{freq}.
  Where \code{pos} gives the position relative to te read start, \code{seq}
  gives the kmer, and \code{freq} gives the frequency of that kmer.
}
\author{
    Daniel Jones
    \email{dcjones@cs.washington.edu}
}
\seealso{
    \code{\link{count.reads}}
}
\examples{
  library(Rsamtools)
  reads_fn <- system.file( "extra/example.bam", package = "seqbias" )
  ref_fn <- system.file( "extra/example.fa", package = "seqbias" )

  I <- GRanges( c('seq1'), IRanges( c(1), c(5000) ), strand = c('-') )

  ref_f <- FaFile( ref_fn )
  open.FaFile( ref_f )

  seqs <- scanFa( ref_f, I )

  neg_idx <- as.logical( I@strand == '-' )
  seqs[neg_idx] <- reverseComplement( seqs[neg_idx] )

  counts <- count.reads( reads_fn, I )

  freqs <- kmer.freq(seqs, counts, L = 30, R = 30, k = 2)
}


