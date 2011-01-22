

setClass( "seqbias",
          representation = representation( ptr = "externalptr" ) )

"seqbias.load" <- function( ref_fn, model_fn )
{
    sb <- .Call( "seqbias_load",
                 path.expand( ref_fn ),
                 path.expand( model_fn ),
                 PACKAGE = "seqbias" )

    new( "seqbias", ptr = sb )
}


"seqbias.save" <- function( sb, fn )
{
    if( class(sb) != "seqbias" ) {
        stop( "The first parameter of seqbias.save must be a seqbias class." )
    }

    .Call( "seqbias_save",
           sb@ptr,
           path.expand(fn),
           PACKAGE = "seqbias" )

    invisible()
}


"seqbias.fit" <- function( ref_fn, reads_fn, n = 1e5, L = 15, R = 15 )
{
    sb <- .Call( "seqbias_fit",
                 path.expand( ref_fn ),
                 path.expand( reads_fn ),
                 as.integer(n),
                 as.integer(L),
                 as.integer(R),
                 PACKAGE = "seqbias" )

    new( "seqbias", ptr = sb )
}


"seqbias.predict" <- function( sb, I )
{
    require(GenomicRanges)

    stopifnot( is( sb, "seqbias" ) )
    stopifnot( is( I, "GRanges" ) )

    tapply( I, INDEX = 1:length(I),
                FUN  = function(x) .Call( "seqbias_predict", sb@ptr,
                                           as.character(seqnames(x)),
                                           start(x), end(x),
                                           as.character(strand(x)),
                                           PACKAGE = "seqbias" ) )
}



"random.intervals" <- function( I, n = 1, ms = 10000 )
{
    require(GenomicRanges)
    stopifnot( is( I, "GRanges" ) )

    seqs <- I@ranges@width
    names(seqs) <- I@seqnames

    # avoid some problems when only one sequence is present
    seqs0 <- c( seqs, 0 )

    # sample a sequence weighted proportionally to the number of length m
    # intervals contained within.
    sample_sequence <- function(m) {
        ps <- pmax( 0, as.numeric( seqs0 - m + 1 ) )
        psum <- sum(ps)
        if( psum == 0 ) {
            stop( paste( 'no sequence is long enough to sample an interval of length ', m, sep='' ) )
        }

        prob <- ps / psum
        sample( seqs0, size = 1, prob = prob )
    }

    xs <- sapply( cbind( 1:n, ms )[,2],
                  FUN    = sample_sequence )


    uniform_ints <- function( min, max )
    {
        apply( as.matrix( cbind( min, max ) ),
               MARGIN = 1,
               FUN    = function(u)
                   as.integer( round(runif( n = 1,  u[[1]], u[[2]] ) ) ) )
    }

    # choose random starts 
    starts <- uniform_ints( min = 0, max = xs - ms )

    strand <- sample( c('+','-'), size = length(starts), replace = TRUE )

    GRanges( seqnames   = names(xs),
             ranges     = IRanges( starts, starts + ms),
             strand     = strand,
             seqlengths = seqs )
}



"count.reads" <- function( reads_fn, I, binary = TRUE )
{
    require(GenomicRanges)
    stopifnot( is( I, "GRanges" ) )

    bam_ptr <- .Call( "seqbias_open_bam", path.expand(reads_fn),
                      PACKAGE = "seqbias" )

    counts <- tapply( I,
                      INDEX = 1:length(I),
                      FUN   = function(x) .Call( "seqbias_count_reads",
                                                 bam_ptr,
                                                 as.character(seqnames(x)),
                                                 start(x), end(x),
                                                 as.character(strand(x)),
                      PACKAGE = "seqbias" ) )

    if( binary ) lapply( counts, FUN = function(c) as.integer( c > 0 ) )
    else counts
}


"kmer.freq" <- function( seqs, counts, L = 50, R = 50, k = 1 )
{
    M <- .Call( "seqbias_alloc_kmer_matrix",
                as.integer(L+1+R), as.integer(k),
                PACKAGE = "seqbias" )

    if( length(seqs) != length(counts) ) {
        stop( "length of seqs and counts must be equal" )
    }

    for( i in 1:length(seqs) ) {
        .Call( "seqbias_tally_kmers",
               M,
               as.character( seqs[[i]] ),
               as.numeric( counts[[i]] ),
               as.integer(L),
               PACKAGE = "seqbias" )
    }


    D <- .Call( "seqbias_dataframe_from_kmer_matrix",
                M, as.integer(L),
                PACKAGE = "seqbias" )


    names(D) <- c( "pos", "seq", "freq" )
    as.data.frame(D)
}




