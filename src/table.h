
/*
 * A quick little hash table designed to be very good at one thing: hashing the
 * positions of every read in a BAM file.
 *
 *           Daniel Jones <dcjones@cs.washington.edu>
 *           Feb 2011
 *
 */




#ifndef PEAKOLATOR_TABLE
#define PEAKOLATOR_TABLE

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include "samtools/sam.h"


/* The table maps positions to counts. */
struct hashed_value
{
    int32_t  pos;
    uint32_t count;
};


/* Each strand of each sequence is stored in a seperate hash table, since the
 * BAM file already assigns an integer index to sequences */
struct subtable
{
    struct hashed_value* A; /* table proper */
    size_t n;               /* table size (as an index into a prime table) */
    size_t m;               /* hashed items */
    size_t max_m;           /* max hashed items before rehash */
};


struct table
{
    /* an array indexd by strand -> sequence id */
    struct subtable* ts[2];
    size_t m; /* number of unique positions */
    size_t n; /* number of sequences */
    char** seq_names;
};



/* initialize, where n is the number of sequences */
void table_create( struct table* T, size_t n );
void table_copy( struct table* T, const struct table* U );
void table_destroy( struct table* T );

void table_inc( struct table*, bam1_t* read );
void table_inc_pos( struct table*, int32_t tid, int32_t pos, uint32_t strand );

uint32_t table_count( struct table*, bam1_t* read );
uint32_t table_count_pos( struct table*, int32_t tid, int32_t pos, uint32_t strand );





/* A simple structure storing the read count at each nonzero position across the
 * genome, sorted by position to allow for binary search. Using this we can
 * get the read count across any interval very quickly. */
struct read_counts
{
    struct hashed_value** xss[2]; /* read count arrays indexed by strand -> tid */
    size_t* mss[2]; /* lengths of xss arrays */
    size_t m;       /* total number of unique positions */
    size_t n;       /* number of sequences */
    char**  seq_names;
};


void read_counts_create( struct read_counts* C, const struct table* T );
void read_counts_copy( struct read_counts* C, const struct read_counts* B );
void read_counts_destroy( struct read_counts* C );


void read_counts_count( const struct read_counts* C,
                        int32_t tid, int32_t start, int32_t end, uint32_t strand,
                        unsigned int* ys );

unsigned int read_counts_total( const struct read_counts* C,
                                int32_t tid, int32_t start, int32_t end, uint32_t strand );


void read_count_occurances( const struct read_counts* C,
                            int32_t tid,  int32_t start, int32_t end, uint32_t strand,
                            uint64_t* ks, size_t max_k );



/* dump to a flat array (used by sequencing bias) */
struct read_pos
{
    int32_t  tid;
    uint32_t strand;
    int32_t  pos;
    uint32_t count;
};

void table_dump( struct table* T, struct read_pos** A_, size_t* N_, size_t limit );

#ifdef __cplusplus
}
#endif

#endif




