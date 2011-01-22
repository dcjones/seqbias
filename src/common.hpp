
#ifndef PEAKOLATOR_COMMON
#define PEAKOLATOR_COMMON

#include <algorithm>
#include <cstdlib>
#include "samtools/faidx.h"


typedef long          pos;     /* genomic position */
typedef unsigned long rcount;  /* read count */

const rcount rcount_nan = (rcount)-1;

int nt2num( char c );
void num2nt( int n, char* nt, int k, bool colorspace = false );
void seqlower( char* seq ); /* lowercase */
void seqrc( char* seq, int n );    /* reverse complement */

template <typename T>
void rev( T* xs, int n ) {
    int i = 0;
    int j = n-1;
    while( i < j ) std::swap(xs[i++],xs[j--]);
}

/* Change the behavior of the faidx_fetch_seq function to be more useful. If
 * coordinates are outside the actual sequence, write N's, rather than adjusting
 * the start,end. */
char* faidx_fetch_seq_forced_lower( const faidx_t* fai, const char *c_name, int p_beg_i, int p_end_i );


#endif
