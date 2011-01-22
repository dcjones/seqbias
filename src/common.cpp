
#include "common.hpp"
#include "logger.h"

/* faidx stuff */
#include "samtools/khash.h"


#ifndef _NO_RAZF
#include "samtools/razf.h"
#else
#ifdef _WIN32
#define ftello(fp) ftell(fp)
#define fseeko(fp, offset, whence) fseek(fp, offset, whence)
#else
extern off_t ftello(FILE *stream);
extern int fseeko(FILE *stream, off_t offset, int whence);
#endif
#define RAZF FILE
#define razf_read(fp, buf, size) fread(buf, 1, size, fp)
#define razf_open(fn, mode) fopen(fn, mode)
#define razf_close(fp) fclose(fp)
#define razf_seek(fp, offset, whence) fseeko(fp, offset, whence)
#define razf_tell(fp) ftello(fp)
#endif
#ifdef _USE_KNETFILE
#include "samtools/knetfile.h"
#endif


#include <cctype>
#include <cstring>
#include <cmath>


/*
 * arbitrary numerical encoding of nucleotides
 */
//static const char num2nt[] = { 'a', 'c', 'g', 't', 'n' };

int nt2num( char c ) {
    switch( c ) {
        case 'a': case 'A': return 0;
        case 'c': case 'C': return 1;
        case 'g': case 'G': return 2;
        case 't': case 'T': return 3;
                  /* assume the N case is very rare, and just
                   * return a random nucleotide */
        default: return rand()%4;
    }
}


/* convert a number (n) encoding a k-mer into a string of nucleotides (nt) */
void num2nt( int n, char* nt, int k, bool colorspace )
{
    int i;
    /* read backwards, then reverse */
    for( i = 0; i < k; i++ ) {
        if( colorspace ) {
            switch( n & 0x3 ) {
                case 0: nt[i] = '0'; break;
                case 1: nt[i] = '1'; break;
                case 2: nt[i] = '2'; break;
                case 3: nt[i] = '3'; break;
            }
        }
        else {
            switch( n & 0x3 ) {
                case 0: nt[i] = 'a'; break;
                case 1: nt[i] = 'c'; break;
                case 2: nt[i] = 'g'; break;
                case 3: nt[i] = 't'; break;
            }
        }
        n >>= 2;
    }

    nt[i] = '\0';
    char tmp;
    for( i = 0; i < k; i++ ) {
        tmp = nt[i];
        nt[i] = nt[k-1-i];
        nt[k-1-i] = tmp;
    }
}


void seqlower( char* seq )
{
    char* c;
    for( c = seq; *c; c++ ) *c = tolower(*c);
}



char complement( char c )
{
    switch( c ) {
        case 'a': return 't';
        case 'c': return 'g';
        case 'g': return 'c';
        case 't': return 'a';
        case 'n': return 'n';
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        case 'N': return 'N';
        default:  return 'n';
    }
}

void seqrc( char* seq, int n )
{
    char c;
    int i,j;
    i = 0;
    j = n-1;
    while( i < j ) {
        c = complement(seq[i]);
        seq[i] = complement(seq[j]);
        seq[j] = c;
        i++; j--;
    }

    if( i == j ) seq[i] = complement(seq[i]);
}





/* an array to double encode a sequence as colorspace */
                                     /*   A    C    G    T    N   */
static const char colorspace[5][5] = { { 'a', 'c', 'g', 't', 'n' },   /* A */
                                       { 'c', 'a', 't', 'g', 'n' },   /* C */
                                       { 'g', 't', 'a', 'c', 'n' },   /* G */
                                       { 't', 'g', 'c', 'a', 'n' },   /* T */
                                       { 'n', 'n', 'n', 'n', 'n' } }; /* N */



/* Double encode AB SOLiD data. That is, take the nucleotide sequence, convert
 * in to colorspace, then represent the colorspace sequence with nucleotides. */
#if 0
static void colorspace_encode( char prev, char* seq )
{
    char col;
    uint32_t i;
    for( i = 0; seq[i] != '\0'; i++ ) {
        col    = colorspace[nt2num(prev)][nt2num(seq[i])];
        prev   = seq[i];
        seq[i] = col;
    }
}
#endif




typedef struct {
	uint64_t len:32, line_len:16, line_blen:16;
	uint64_t offset;
} faidx1_t;
KHASH_MAP_INIT_STR(s, faidx1_t)


struct __faidx_t {
	RAZF *rz;
	int n, m;
	char **name;
	khash_t(s) *hash;
};



char* faidx_fetch_seq_forced_lower( const faidx_t* fai, const char *c_name, int p_beg_i, int p_end_i )
{

	int l;
	char c;
    khiter_t iter;
    faidx1_t val;
    char* seq0;
    char* seq = NULL;

    iter = kh_get(s, fai->hash, c_name);
    if(iter == kh_end(fai->hash)) return 0;

    seq0 = seq = (char*)malloc( (p_end_i - p_beg_i + 2) * sizeof(char) );
    if( seq0 == NULL ) fail( "Out of memory.\n" );
    seq0[p_end_i-p_beg_i+1] = '\0';

    val = kh_value(fai->hash, iter);

    /* entirely off the map: all Ns */
    if( p_beg_i >= (int)val.len || p_end_i < 0 ) {
        while( p_beg_i <= p_end_i ) {
            *seq++ ='n';
        }
        return seq0;
    }

    /* beginning is off the map */
    while( p_beg_i < 0 && p_beg_i <= p_end_i ) {
        *seq++ = 'n';
        p_beg_i++;
    }

    /* end is off the map */
    while( p_end_i >= (int)val.len ) {
        seq[p_end_i-p_beg_i] = 'n';
        p_end_i--;
    }

    /* retrieve the sequence */
	l = 0;
	razf_seek(fai->rz, val.offset + p_beg_i / val.line_blen * val.line_len + p_beg_i % val.line_blen, SEEK_SET);
	while (razf_read(fai->rz, &c, 1) == 1 && l < p_end_i - p_beg_i + 1)
		if (isgraph(c)) seq[l++] = tolower(c);
    
    while( p_beg_i+l <= p_end_i ) seq[l++] = 'n';

    return seq0;
}


