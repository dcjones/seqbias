\name{seqbias.fit}
\alias{seqbias.fit}
\title{Fitting seqbias models}
\description{Fits a seqbias module given a reference sequence and reads in BAM
format}
\usage{seqbias.fit(ref_fn, reads_fn, n = 1e5, L = 15, R = 15)}
\arguments{
    \item{ref_fn}{filename of a reference sequence against which the reads are
    aligned, in FASTA format.}
    \item{reads_fn}{filename of aligned reads in BAM format.}
    \item{n}{train on at most this many reads.}
    \item{L}{consider at most L positions to the left of the read start.}
    \item{R}{consider at most R positions to the right of the read start.}
}
\details{
    A Bayesian network is trained on the first \code{n} unique reads in the provided
    BAM file, predicting the posterior probability of a read beginning at a
    position given the surrounding sequence. This is used to discern the
    sequencing bias: how more or less likely a read is to fall on a particular
    position.

    The abundance of region can be more accurately assessed by normalizing
    (dividing) each position by its predicted bias.
}
\value{A vector of reals giving the predicted sequencing bias for each
position.}
\note{
    Both the BAM file and the FASTA file should be indexed, with,
        'samtools index'
    and,
        'samtools faidx'
    respectively.
}
\author{
    Daniel Jones
    \email{dcjones@cs.washington.edu}}
\seealso{
    \code{\link{seqbias.predict}}
}
\examples{
  reads_fn <- system.file( "extra/example.bam", package = "seqbias" )
  ref_fn <- system.file( "extra/example.fa", package = "seqbias" )

  sb <- seqbias.fit( ref_fn, reads_fn )

  I <- GRanges( c('seq1'), IRanges( c(1), c(5000) ), strand = c('-') )

  bias <- seqbias.predict( sb, I )
}

