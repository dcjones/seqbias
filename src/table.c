
#include "table.h"
#include "logger.h"



#define NUM_PRIMES 28
static const uint32_t primes[NUM_PRIMES] = {
           53U,         97U,        193U,        389U,    
          769U,       1543U,       3079U,       6151U,  
        12289U,      24593U,      49157U,      98317U,  
       196613U,     393241U,     786433U,    1572869U, 
      3145739U,    6291469U,   12582917U,   25165843U,
     50331653U,  100663319U,  201326611U,  402653189U,
    805306457U, 1610612741U, 3221225473U, 4294967291U };


static const double max_load = 0.75;

/* marks a vacant cell */
static const int32_t nilpos = -1;


#ifndef MIN
#define MIN(a,b) ((a)<=(b)?(a):(b))
#endif

#ifndef MAX
#define MAX(a,b) ((a)>=(b)?(a):(b))
#endif


/* From Thomas Wang (http://www.cris.com/~Ttwang/tech/inthash.htm) */
uint32_t hash( uint32_t a)
{
    a = (a ^ 61) ^ (a >> 16);
    a = a + (a << 3);
    a = a ^ (a >> 4);
    a = a * 0x27d4eb2d;
    a = a ^ (a >> 15);
    return a;
}



/* simple quadratic probing */
uint32_t probe( uint32_t h, uint32_t i )
{
    const double c1 = 0.5;
    const double c2 = 0.5;

    return h
           + (uint32_t)(c1 * (double)i)
           + (uint32_t)(c2 * (double)(i*i));
}


void subtable_create( struct subtable* T )
{
    T->n = 0;
    T->m = 0;
    T->A = malloc(sizeof(struct hashed_value)*primes[T->n]);
    size_t i;
    for( i = 0; i < primes[T->n]; i++ ) {
        T->A[i].pos = nilpos;
        T->A[i].count = 0;
    }
    T->max_m = (size_t)(((double)primes[T->n]) * max_load);
}


void subtable_copy( struct subtable* T, const struct subtable* U )
{
    T->n     = U->n;
    T->m     = U->m;
    T->max_m = U->max_m;
    T->A     = malloc(sizeof(struct hashed_value)*primes[T->n]);

    size_t i;
    for( i = 0; i < primes[T->n]; i++ ) {
        T->A[i].pos   = U->A[i].pos;
        T->A[i].count = U->A[i].pos;
    }
}


void subtable_destroy( struct subtable* T )
{
    free( T->A );
    T->A = NULL;
}


void subtable_rehash( struct subtable* T, size_t new_n );


bool subtable_inc( struct subtable* T, int32_t pos )
{
    if( T->m == T->max_m ) subtable_rehash( T, T->n + 1 );

    uint32_t h = hash( (uint32_t)pos );
    uint32_t i = 1;
    uint32_t j = h % primes[T->n];

    while( T->A[j].pos != nilpos && T->A[j].pos != pos ) {
        j = probe( h, ++i ) % primes[T->n];
    }

    if( T->A[j].pos == nilpos ) {
        T->A[j].pos = pos;
        T->A[j].count = 1;
        T->m++;
        return true;
    }
    else {
        T->A[j].count++;
        return false;
    }
}


void subtable_set( struct subtable* T, int32_t pos, uint32_t count )
{
    if( T->m == T->max_m ) subtable_rehash( T, T->n + 1 );

    uint32_t h = hash( (uint32_t)pos );
    uint32_t i = 1;
    uint32_t j = h % primes[T->n];

    while( T->A[j].pos != nilpos && T->A[j].pos != pos ) {
        j = probe( h, ++i ) % primes[T->n];
    }

    if( T->A[j].pos == nilpos ) {
        T->A[j].pos = pos;
        T->A[j].count = count;
    }
    else {
        T->A[j].count = count;
    }
}



void subtable_rehash( struct subtable* T, size_t new_n )
{
    if( new_n >= NUM_PRIMES ) {
        fail( "a table has grown too large" );
    }

    struct subtable U;
    U.n = new_n;
    U.A = malloc( sizeof(struct hashed_value) * primes[U.n] );
    size_t i;
    for( i = 0; i < primes[U.n]; i++ ) {
        U.A[i].pos = nilpos;
        U.A[i].count = 0;
    }

    U.m = 0;
    U.max_m = (size_t)(((double)primes[U.n]) * max_load);


    for( i = 0; i < primes[T->n]; i++ ) {
        if( T->A[i].pos == nilpos ) continue;
        subtable_set( &U, T->A[i].pos, T->A[i].count );
    }

    free(T->A);
    T->A = U.A;
    T->n = U.n;
    T->max_m = U.max_m;
}



uint32_t subtable_count( struct subtable* T, int32_t pos )
{
    uint32_t h = hash( (uint32_t)pos );
    uint32_t i = 1;
    uint32_t j = h % primes[T->n];

    while( T->A[j].pos != nilpos && T->A[j].pos != pos ) {
        j = probe( h, ++i ) % primes[T->n];
    }

    if( T->A[j].pos == pos ) return T->A[j].count;
    else                     return 0;
}


void table_create( struct table* T, size_t n )
{
    T->seq_names = NULL;
    T->n = n;
    T->m = 0;

    T->ts[0] = malloc( n * sizeof(struct subtable) );
    T->ts[1] = malloc( n * sizeof(struct subtable) );

    size_t i, j;
    for( i = 0; i <= 1; i++ ) {
        for( j = 0; j < n; j++ ) {
            subtable_create( &T->ts[i][j] );
        }
    }
}

void table_copy( struct table* T, const struct table* U )
{
    T->seq_names = U->seq_names;
    T->n         = U->n;
    T->m         = U->m;

    T->ts[0] = malloc( T->n * sizeof(struct subtable) );
    T->ts[1] = malloc( T->n * sizeof(struct subtable) );

    size_t i, j;
    for( i = 0; i <= 1; i++ ) {
        for( j = 0; j < T->n; j++ ) {
            subtable_copy( &T->ts[i][j], &U->ts[i][j] );
        }
    }
}


void table_destroy( struct table* T )
{
    size_t i, j;
    for( i = 0; i <= 1; i++ ) {
        for( j = 0; j < T->n; j++ ) {
            subtable_destroy( &T->ts[i][j] );
        }
    }

    free( T->ts[0] );
    free( T->ts[1] );

    T->n = 0;
}



void table_inc( struct table* T, bam1_t* read )
{
    int32_t pos;
    if( bam1_strand(read) ) pos = bam_calend( &read->core, bam1_cigar(read) ) - 1;
    else                    pos = read->core.pos;

    table_inc_pos( T, read->core.tid, pos, bam1_strand(read) );
}




void table_inc_pos( struct table* T, int32_t tid, int32_t pos, uint32_t strand )
{
    if( tid < 0 || tid >= T->n ) return;
    if( subtable_inc( &T->ts[strand][tid], pos ) ) T->m++;
}


uint32_t table_count( struct table* T, bam1_t* read )
{
    int32_t pos;
    if( bam1_strand(read) ) pos = bam_calend( &read->core, bam1_cigar(read) ) - 1;
    else                    pos = read->core.pos;

    return table_count_pos( T, read->core.tid, pos, bam1_strand(read) );
}


uint32_t table_count_pos( struct table* T, int32_t tid, int32_t pos, uint32_t strand )
{
    if( tid < 0 || tid >= T->n ) return 0;
    return subtable_count( &T->ts[strand][tid], pos );
}



int hashed_value_compare( const void* p1, const void* p2 )
{
    return ((struct hashed_value*)p1)->pos - ((struct hashed_value*)p2)->pos;
}


void read_counts_create( struct read_counts* C, const struct table* T )
{
    C->n = T->n;
    C->m = T->m;
    C->seq_names   = T->seq_names;

    C->mss[0] = malloc( C->n * sizeof(size_t) );
    C->mss[1] = malloc( C->n * sizeof(size_t) );
    
    C->xss[0] = malloc( C->n * sizeof(struct hashed_value*) );
    C->xss[1] = malloc( C->n * sizeof(struct hashed_value*) );

    int32_t  tid;
    uint32_t strand;
    size_t i, j;

    size_t m, n;
    struct hashed_value* ys;
    struct hashed_value* xs;


    for( strand = 0; strand <= 1; strand++ ) {
        for( tid = 0; tid < T->n; tid++ ) {
            m  = T->ts[strand][tid].m;
            n  = T->ts[strand][tid].n;
            ys = T->ts[strand][tid].A;
            xs = malloc( m * sizeof(struct hashed_value) );

            for( i = 0, j = 0; j < primes[n]; j++ ) {
                if( ys[j].pos != nilpos ) {
                    xs[i].pos   = ys[j].pos;
                    xs[i].count = ys[j].count;
                    i++;
                }
            }

            qsort( xs, m, sizeof(struct hashed_value), hashed_value_compare );

            C->mss[strand][tid] = m;
            C->xss[strand][tid] = xs;
        }
    }
}

void read_counts_copy( struct read_counts* C, const struct read_counts* B )
{
    C->n = B->n;
    C->m = B->m;
    C->seq_names = B->seq_names;

    int32_t  tid;
    uint32_t strand;
    size_t siz;

    for( strand = 0; strand <= 1; strand++ ) {
        C->mss[strand] = malloc( C->n * sizeof(size_t) );
        C->xss[strand] = malloc( C->n * sizeof(struct hashed_value*) );
        for( tid = 0; tid < C->n; tid++ ) {
            C->mss[strand][tid] = B->mss[strand][tid];
            siz = C->mss[strand][tid] * sizeof(struct hashed_value);
            C->xss[strand][tid] = malloc( siz );
            memcpy( C->xss[strand][tid], B->xss[strand][tid], siz );
        }
    }
}


void read_counts_destroy( struct read_counts* C )
{
    int32_t  tid;
    uint32_t strand;

    for( strand = 0; strand <= 1; strand++ ) {
        for( tid = 0; tid < C->n; tid++ ) {
            free( C->xss[strand][tid] );
            C->xss[strand][tid] = NULL;
        }
    }

    free( C->mss[0] ); C->mss[0] = NULL;
    free( C->mss[1] ); C->mss[1] = NULL;

    free( C->xss[0] ); C->xss[0] = NULL;
    free( C->xss[1] ); C->xss[1] = NULL;
}



/* find an index i, such that
 *      xs[i-1].pos < start <= xs[i].pos
 * using binary search.
 */
size_t bisect( struct hashed_value* xs, size_t m, int32_t start )
{
    size_t a = 0;
    size_t b = m;
    size_t i = 0;

    while( a <= b ) {
        i = a + (b - a) / 2;

        if( xs[i].pos < start )                  a = i + 1;
        else if( i > 0 && start <= xs[i-1].pos ) b = i - 1;

        /* xs[i-1].pos <= start <= xs[i] */
        else break;
    }

    return i;
}


void read_counts_count( const struct read_counts* C,
                        int32_t tid, int32_t start, int32_t end, uint32_t strand,
                        unsigned int* ys )
{
    struct hashed_value* xs = C->xss[strand][tid];
    size_t               m  = C->mss[strand][tid];

    if( m == 0 ) return;

    size_t i = bisect( xs, m, start );

    memset( ys, 0, m*sizeof(unsigned int) );
    while( i < m && xs[i].pos <= end ) {
        ys[ xs[i].pos - start ] = xs[i].count;
        i++;
    }
}

unsigned int read_counts_total( const struct read_counts* C,
                                int32_t tid, int32_t start, int32_t end, uint32_t strand )
{
    struct hashed_value* xs = C->xss[strand][tid];
    size_t               m  = C->mss[strand][tid];

    if( m == 0 ) return 0;

    size_t i = bisect( xs, m, start );

    unsigned int total = 0;
    while( i < m && xs[i].pos <= end ) {
        total += xs[i].count;
        i++;
    }

    return total;
}


void read_count_occurances( const struct read_counts* C,
                            int32_t tid,  int32_t start, int32_t end, uint32_t strand,
                            uint64_t* ks, size_t max_k )
{
    struct hashed_value* xs = C->xss[strand][tid];
    size_t               m  = C->mss[strand][tid];

    if( m == 0 ) return;

    size_t i;

    i = bisect( xs, m, start );
    uint64_t nonzeros = 0;

    while( i < m && xs[i].pos <= end ) {
        if( xs[i].count <= max_k ) ks[xs[i].count]++;
        nonzeros++;
        i++;
    }

    uint64_t zeros = (end - start + 1) - nonzeros;

    /* Ignore and leading or trailing zeros if we are at the start or end of the
     * sequence. Many genome assemblies have several kilobases of N's at the
     * beginning and end. Considering these will lead to deflated statistics. */
    if( start <= xs[0].pos ) {
        zeros -= MIN( end, xs[0].pos ) - start + 1;
    }

    if( end >= xs[m-1].pos ) {
        zeros -= end - MAX( start, xs[m-1].pos ) + 1;
    }

    ks[0] += zeros;
}



void table_dump( struct table* T, struct read_pos** A_, size_t* N_, size_t limit )
{
    struct read_pos* A;
    size_t N = 0;
    size_t i, j;


    for( i = 0; i <= 1; i++ ) {
        for( j = 0; j < T->n; j++ ) {
            N += T->ts[i][j].m;
        }
    }

    if( limit > 0 && N > limit ) N = limit;

    A = malloc( N * sizeof(struct read_pos) );


    size_t u = 0;
    size_t v;
    for( i = 0; i <= 1; i++ ) {
        for( j = 0; j < T->n; j++ ) {
            for( v = 0; v < primes[T->ts[i][j].n]; v++ ) {
                if( T->ts[i][j].A[v].pos != -1 ) {
                    A[u].tid = j;
                    A[u].strand = i;
                    A[u].pos    = T->ts[i][j].A[v].pos;
                    A[u].count  = T->ts[i][j].A[v].count;
                    u++;
                    if( u >= N ) goto table_dump_finish;
                }
            }
        }
    }
table_dump_finish:

    *A_ = A;
    *N_ = N;
}



