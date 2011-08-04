
#ifndef ISOLATOR_COMMON_H
#define ISOLATOR_COMMON_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

/* kmers are encoded in (at least) 32 bits, allowing for k <= 16 */
typedef unsigned int kmer;

/* genomic position */
typedef long pos_t;

/* sequence identifier */
typedef int seqid_t;

typedef enum {
    strand_pos = 0,
    strand_neg = 1,
    strand_na  = 2
} strand_t;


/* compare sequence names to sort them in an appealing manner */
int seqname_compare(const char* u, const char* v);


void* malloc_or_die(size_t);
void* realloc_or_die(void*, size_t);

/** reverse complement */
void seqrc(char* seq, int n);


#ifdef __cplusplus
}
#endif

#endif

