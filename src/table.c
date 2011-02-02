/*
 *           hash
 *           A quick and dirty hash table implementation.
 *
 *           Daniel Jones <dcjones@cs.washington.edu>
 *           July 2010
 *
 */

#include "table.h"
#include "superfasthash.h"


#define INITIAL_TABLE_SIZE 10000
#define MAX_LOAD 0.75
#define MIN_LOAD 0.05 /* make sure this is less than MAX_LOAD/2 */


void rehash( struct table* T, size_t new_n );



/* Create an empty hash table. */
void table_create( struct table* T )
{
    T->A = malloc(sizeof(struct hashed_value*)*INITIAL_TABLE_SIZE);
    memset( T->A, 0, sizeof(struct hashed_value*)*INITIAL_TABLE_SIZE );
    T->n = INITIAL_TABLE_SIZE;
    T->m = 0;
    T->max_m = T->n * MAX_LOAD;
    T->min_m = T->n * MIN_LOAD;
}



/* Remove all elements in the table. */
void table_clear( struct table* T )
{
    struct hashed_value* u;
    size_t i;
    for( i = 0; i < T->n; i++ ) {
        while( T->A[i] ){
            u = T->A[i]->next;
            free(T->A[i]);
            T->A[i] = u;
        }
    }
    T->m = 0;
}



/* Free all memory associated with a table. */
void table_destroy( struct table* T )
{
    table_clear(T);
    free(T->A);
}



void table_inc( struct table* T, bam1_t* read )
{
    if( T->m == T->max_m ) rehash( T, T->n*2 );

    struct read_pos pos;
    pos.tid = read->core.tid;
    /*                                      XXX: I don't know why I must
     *                                      subtract one here. It bothers me
     *                                      that I don't, but the results do not
     *                                      come out correctly otherwise.
     *                                      Possibly, the function gives the
     *                                      nucleotide immediately after the
     *                                      read. */
    if( bam1_strand(read) ) pos.pos = bam_calend( &read->core, bam1_cigar(read) ) - 1;
    else                    pos.pos = read->core.pos;
    pos.strand = bam1_strand(read);

    uint32_t h = hash((void*)&pos, sizeof(struct read_pos)) % T->n;

    struct hashed_value* u = T->A[h];

    while(u) {
        if( memcmp( &u->pos, &pos, sizeof(struct read_pos) ) == 0 ) {
            u->count++;
            return;
        }

        u = u->next;
    }

    u = malloc(sizeof(struct hashed_value));
    memcpy( &u->pos, &pos, sizeof(struct read_pos) );

    u->count = 1;

    u->next = T->A[h];
    T->A[h] = u;

    T->m++;
}



/* Insert existing entries without copying sequences. Used for rehashing. */
bool table_insert_without_copy( struct table* T, struct hashed_value* V )
{
    if( T->m == T->max_m ) rehash( T, T->n*2 );

    uint32_t h = hash((void*)&V->pos,sizeof(struct read_pos)) % T->n;

    V->next = T->A[h];
    T->A[h] = V;

    T->m++;

    return true;
}


/* Rezise the table T to new_n. */
void rehash( struct table* T, size_t new_n )
{
    struct table U;
    U.n = new_n;
    U.m = 0;
    U.A = malloc( sizeof(struct hashed_value*) * U.n );
    memset( U.A, 0, sizeof(struct hashed_value*) * U.n );


    struct hashed_value *j,*k;
    size_t i;
    for( i = 0; i < T->n; i++ ) {
        j = T->A[i];
        while( j ) {
            k = j->next;
            table_insert_without_copy( &U, j );
            j = k;
        }
        T->A[i] = NULL;
    }

    free(T->A);
    T->A = U.A;
    T->n = U.n;
    T->max_m = T->n*MAX_LOAD;
    T->min_m = T->n*MIN_LOAD;
}


int compare_seq_hash( const void* x, const void* y )
{
    struct hashed_value* const * a = x;
    struct hashed_value* const * b = y;

    int c = (*a)->pos.tid - (*b)->pos.tid;
    if( c == 0 ) {
        uint32_t ha = hash( (void*)&(*a)->pos, sizeof(struct read_pos) );
        uint32_t hb = hash( (void*)&(*b)->pos, sizeof(struct read_pos) );
        return (int32_t)ha - (int32_t)hb;
    }
    else return c;

}

int compare_seq_count( const void* x, const void* y )
{
    struct hashed_value* const * a = x;
    struct hashed_value* const * b = y;

    int c = (*a)->pos.tid - (*b)->pos.tid;
    if( c == 0 ) {
        if( (*a)->count == (*b)->count ) return 0;
        else return (*a)->count > (*b)->count ? 1 : -1;
    }
    else return c;
}


int compare_count( const void* x, const void* y )
{
    struct hashed_value* const * a = x;
    struct hashed_value* const * b = y;

    return (*a)->count - (*b)->count;
}


int compare_hashed_pos( const void* x, const void* y )
{
    struct hashed_value* const * a = x;
    struct hashed_value* const * b = y;

    int c = (*a)->pos.tid - (*b)->pos.tid;
    if( c == 0 ) {
        if( (*a)->pos.pos == (*b)->pos.pos ) return 0;
        else return (*a)->pos.pos > (*b)->pos.pos ? 1 : -1;
    }
    else return c;
}


void sort_table( struct table* T,
                 struct hashed_value*** S_,
                 int(*compar)(const void *, const void *) )
{
    struct hashed_value** S = malloc( sizeof(struct hashed_value*) * T->m );
    memset( S, 0, sizeof(struct hashed_value*) * T->m );

    struct hashed_value* j;
    size_t i,k;
    for( i=0, k=0; i < T->n; i++ ) {
        j = T->A[i];
        while( j ) {
            S[k] = j;
            k++;
            j = j->next;
        }
    }

    qsort( S, T->m, sizeof(struct hashed_value*), compar );

    *S_ = S;
}

void table_sort_by_seq_rand( struct table* T,
                             struct hashed_value*** S_ )
{
    sort_table( T, S_, compare_seq_hash );
}

void table_sort_by_seq_count( struct table* T,
                              struct hashed_value*** S_ )
{
    sort_table( T, S_, compare_seq_count );
}


void table_sort_by_count( struct table* T,
                    struct hashed_value*** S_ )
{
    sort_table( T, S_, compare_count );
}

void table_sort_by_position( struct table* T,
                    struct hashed_value*** S_ )
{
    sort_table( T, S_, compare_hashed_pos );
}


bool table_member( struct table* T, struct read_pos* pos )
{
    uint32_t h = hash((void*)pos, sizeof(struct read_pos)) % T->n;

    struct hashed_value* u = T->A[h];

    while(u) {
        if( memcmp( &u->pos, pos, sizeof(struct read_pos) ) == 0 ) {
            return true;
        }

        u = u->next;
    }

    return false;
}



void rehash_tail( struct table* T, int32_t q1, int32_t q2 )
{
    struct table U;
    U.n = T->n;
    U.m = 0;
    U.A = malloc( sizeof(struct hashed_value*) * U.n );
    memset( U.A, 0, sizeof(struct hashed_value*) * U.n );

    struct hashed_value** S;
    table_sort_by_count( T, &S );

    int32_t i;
    for( i = q1; i < q2; i++ ) {
        table_insert_without_copy( &U, S[i] );
    }

    /* free the rest */
    for( ; i < T->m; i++ ) free(S[i]);

    free(S);
    free(T->A);
    T->A = U.A;
    T->m = U.m;
}

