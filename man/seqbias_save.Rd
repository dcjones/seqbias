\name{seqbias.save}
\alias{seqbias.save}
\title{Saving seqbias models}
\description{
Writes to a seqbias model to a file suitable to loaded with 'seqbias.load'
}
\usage{seqbias.save(sb, fn)}
\arguments{
    \item{sb}{A seqbias class created with 'seqbias.fit' or 'seqbias.load' }
    \item{fn}{A file name to write the model to}
}
\details{
    A large seqbias model can take some time (several minutes) to fit. It is
    often preferable to do this just once. This function writes the model to a
    file in YAML format, suitable to be read with 'seqbias.load'.
}
\author{
    Daniel Jones
    \email{dcjones@cs.washington.edu}}
\seealso{
    \code{\link{seqbias.load}}
}
\examples{
  reads_fn <- system.file( "extra/example.bam", package = "seqbias" )
  ref_fn <- system.file( "extra/example.fa", package = "seqbias" )

  sb <- seqbias.fit( ref_fn, reads_fn )

  seqbias.save( sb, "my_seqbias_model.yml" )

  # load sometime later
  sb <- seqbias.load( ref_fn, "my_seqbias_model.yml" )
}

