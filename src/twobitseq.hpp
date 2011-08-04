/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */


#ifndef TWOBITSEQ_MATRIX_HPP
#define TWOBITSEQ_MATRIX_HPP

#include "common.h"
#include <cstdlib>

/** a nucleotide sequence encoded packed two-bits per nucleotide.
 *
 * This is useful for a couple reasons:
 *   1. this uses a fourth of the memory as char arrays would use.
 *   2. we need not worry about capitalization, A vs U, N's, or the like.
 */
class twobitseq
{
    public:
        /** create an empty sequence (uninitialized, filled with noise) */
        twobitseq(size_t n);

        /** create a sequence from a DNA/RNA sequence string */
        twobitseq(const char* seq);

        /** copy constructor */
        twobitseq(const twobitseq&);

        ~twobitseq();

        void operator = (const twobitseq&);

        /** extract a kmer.
         *
         * The kmer is made from all the positions with a '1' bit in the 'mask'
         * vector, at the given offset in the sequence.
         */
        int make_kmer(kmer& K, size_t offset, bool* mask, size_t mask_len) const;
        


    private:
        kmer* xs;
        size_t n;

        /** how many nucleotides can we back in a kmer (or, 'k') */
        static const size_t max_kmer;
};

/** Convert a nucleotide charactor to a number, using the same scheme as
 * twobitseq */
kmer nuc_to_num(char c);

/** Convert a number n encoding a kmer into a string of nucleotides */
void num_to_nuc(char* dest, int n, int k);

#endif

