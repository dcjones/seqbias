/*
 *           hash
 *           A cute little hash table implementation, to hash read positions
 *           very quickly.
 *
 *           Daniel Jones <dcjones@cs.washington.edu>
 *           July 2010
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


struct read_pos
{
    int32_t   tid;
    int32_t   pos;
    uint32_t  strand;
};

struct hashed_value
{
    struct   read_pos pos;
    uint32_t count;
    struct   hashed_value* next;
};


/* Hash table structure. */
struct table
{
    struct hashed_value** A; /* table proper */
    size_t n;                /* table size */
    size_t m;                /* hashed items */
    size_t max_m;            /* max hashed items before rehash */
    size_t min_m;            /* min hashed items before rehash */
};


void table_create( struct table* T );
void table_destroy( struct table* T );

void table_inc( struct table*, bam1_t* read );

bool table_member( struct table*, struct read_pos* pos );



void table_sort_by_seq_rand( struct table* T,
                             struct hashed_value*** _S );

void table_sort_by_seq_count( struct table* T,
                                  struct hashed_value*** S );

void table_sort_by_count( struct table* T,
                    struct hashed_value*** S );

void table_sort_by_position( struct table* T,
                       struct hashed_value*** S );

void rehash_tail( struct table* T, int32_t q1, int32_t q2 );
 


#ifdef __cplusplus
}
#endif

#endif
