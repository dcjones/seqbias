\name{seqbias.load}
\alias{seqbias.load}
\title{Loading seqbias models}
\description{
Loads a seqbias model from a file written with 'seqbias.save'.
}
\usage{seqbias.load(ref_fn, model_fn)}
\arguments{
    \item{ref_fn}{filename of a reference sequence against which the reads are
    aligned in FASTA format.}
    \item{model_fn}{filename of a saved seqbias model}
}
\details{
    A large seqbias model can take some time (several minutes) to fit. It is
    often preferable to do this just once. This function load the model from a
    file in YAML format, having been written with 'seqbias.save'.
}
\value{A seqbias class.}
\author{
    Daniel Jones
    \email{dcjones@cs.washington.edu}}
\seealso{
    \code{\link{seqbias.save}}
}
\examples{
  reads_fn <- system.file( "extra/example.bam", package = "seqbias" )
  ref_fn <- system.file( "extra/example.fa", package = "seqbias" )

  sb <- seqbias.fit( ref_fn, reads_fn )

  seqbias.save( sb, "my_seqbias_model.yml" )

  # load sometime later
  sb <- seqbias.load( ref_fn, "my_seqbias_model.yml" )
}


