
#include "common.h"
#include "logger.h"
#include <stdio.h>
#include <ctype.h>
#include <string.h>


void* malloc_or_die(size_t n)
{
    void* ptr = malloc(n);
    if (ptr == NULL) {
        logger_abort("Falied to allocate %zu bytes. Out of memory.\n", n);
        return NULL;
    }

    return ptr;
}




void* realloc_or_die(void* ptr, size_t n)
{
    ptr = realloc(ptr, n);
    if (ptr == NULL) {
        logger_abort("Falied to (re)allocate %zu bytes. Out of memory.\n", n);
        return NULL;
    }

    return ptr;
}


static int seqname_num(const char* u)
{
    int num = 0;
    while(*u != '\0') {
        if (isdigit(*u)) {
            sscanf(u, "%d", &num);
            break;
        }

        ++u;
    }

    return num;
}


int seqname_compare(const char* u, const char* v)
{
    int x = seqname_num(u);
    int y = seqname_num(v);

    if (x != y) return x - y;
    else        return strcmp(u, v);
}



static char complement(char c)
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

void seqrc(char* seq, int n)
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


