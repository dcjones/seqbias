/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#ifndef ISOLATOR_SEQUENCING_BIAS_HPP
#define ISOLATOR_SEQUENCING_BIAS_HPP


#include "common.h"
#include "motif.hpp"
#include "pos_table.h"
#include "samtools/faidx.h"
#include "yaml-cpp/yaml.h"
#include <string>

/** A representation of sequencing bias for a particular dataset.
 */
class sequencing_bias
{
    public:
        /** Load a model that has bee previously trained. */
        sequencing_bias(const char* ref_fn,
                        const char* model_fn);

        /** Train a new model.
         *
         * \param ref_fn    File name of an indexed FASTA file.
         * \param reads_fn  File name of an indexed BAM file.
         * \param max_reads How many reads (at most) to use.
         * \param L         How many left positions to consider.
         * \param R         How many right positions to consider.
         * \param complexity_penalty Force sparser models by making this number
         *                           larger.
         */
        sequencing_bias(const char* ref_fn,
                        const char* reads_fn,
                        size_t max_reads,
                        pos_t L, pos_t R,
                        double complexity_penalty = 1.0);


        /** Train a new model.
         *
         * \param ref_fn    File name of an indexed FASTA file.
         * \param T         A table of hashed read positions.
         * \param max_reads How many reads (at most) to use.
         * \param L         How many left positions to consider.
         * \param R         How many right positions to consider.
         * \param complexity_penalty Force sparser models by making this number
         *                           larger.
         */
        sequencing_bias(const char* ref_fn,
                        pos_table* T,
                        size_t max_reads,
                        pos_t L, pos_t R,
                        double complexity_penalty = 1.0);

        /** destructor */
        ~sequencing_bias();

        /** Compute the bias across the given region. 
         *
         * The vector returned must be freed with 'delete []'.
         */
        double* get_bias(const char* seqname, pos_t start, pos_t end, strand_t strand);

        /** Serialize the model to a file. */
        void save_to_file(const char* fn) const;

        /** Serialize the model and emit in YAML format. */
        void to_yaml(YAML::Emitter&) const;

        /** Return a string of the model graph in dot format. */
        std::string print_model_graph();


    private:
        sequencing_bias();

        void clear();

        void build(const char* ref_fn,
                   const char* reads_fn,
                   size_t max_reads,
                   pos_t L, pos_t R,
                   double complexity_penalty);

        void build(const char* ref_fn,
                   pos_table* T,
                   size_t max_reads,
                   pos_t L, pos_t R,
                   double complexity_penalty);


        /* left and right sequence context */
        pos_t L, R;

        /* reference sequence */
        faidx_t*    ref_f;
        std::string ref_fn;

        /* trained background (M0) and foreground (M1) models */
        motif* M;
};


#endif
