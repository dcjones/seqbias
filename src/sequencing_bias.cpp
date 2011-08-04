/*
 * This file is part of Isolator.
 *
 * Copyright (c) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 */

#include "sequencing_bias.hpp"
#include "logger.hpp"
#include "twobitseq.hpp"
#include "miscmath.h"
#include "samtools_extra.h"
#include "samtools/faidx.h"
#include <cmath>
#include <cctype>
#include <ctime>
#include <fstream>
#include <algorithm>

using namespace std;


/* round away from zero */
static double round_away(double a)
{
    if (a < 0.0) return floor(a);
    else         return ceil(a);
}


/* simply uniform random numbers */
static double rand_uniform(double a, double b)
{
    return a + b * (double)rand() / (double)RAND_MAX;
}

/* random gaussians (copied from GSL, to avoid dependency) */
static double rand_gauss(double sigma)
{
    double x, y, r2;

    do
    {
        /* choose x,y in uniform square (-1,-1) to (+1,+1) */
        x = -1 + 2 * rand_uniform(0.0,1.0);
        y = -1 + 2 * rand_uniform(0.0,1.0);

        /* see if it is in the unit circle */
        r2 = x * x + y * y;
    }
    while (r2 > 1.0 || r2 == 0);

    /* Box-Muller transform */
    return sigma * y * sqrt (-2.0 * log (r2) / r2);
}


static double rand_trunc_gauss(double sigma, double a, double b)
{
    double x;
    do {
        x = rand_gauss(sigma);
    } while (x < a || x > b);

    return x;
}


double gauss_pdf (const double x, const double sigma)
{
  double u = x / fabs (sigma);
  double p = (1 / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-u * u / 2);
  return p;
}



static int read_pos_tid_compare(const void* p1, const void* p2)
{
    return ((read_pos*)p1)->tid - ((read_pos*)p2)->tid;
}

static int read_pos_count_compare(const void* p1, const void* p2)
{
    return (int)((read_pos*)p2)->count - (int)((read_pos*)p1)->count;
}


/*
static int read_pos_tid_count_compare(const void* p1, const void* p2)
{
    int c = (int)(((read_pos*)p1)->tid - ((read_pos*)p2)->tid);

    if (c == 0) {
        return (int)((read_pos*)p2)->count - (int)((read_pos*)p1)->count;
    }
    else return c;
}
*/



sequencing_bias::sequencing_bias()
    : ref_f(NULL)
    , M(NULL)
{}

sequencing_bias::sequencing_bias(const char* ref_fn,
                                 const char* model_fn)
    : ref_f(NULL)
    , M(NULL)
{
    std::ifstream fin;
    fin.open(model_fn);

    YAML::Parser parser(fin);
    YAML::Node doc;
    parser.GetNextDocument(doc);

    doc["L"] >> L;
    doc["R"] >> R;

    M = new motif(doc["motif"]);

    this->ref_fn = ref_fn;

    ref_f = fai_load(ref_fn);
    if (ref_f == NULL) {
        logger::abort("Can't open indexed FASTA file %s\n", ref_fn);
    }
}


sequencing_bias::sequencing_bias(const char* ref_fn,
                                 const char* reads_fn,
                                 size_t max_reads,
                                 pos_t L, pos_t R,
                                 double complexity_penalty)
    : ref_f(NULL)
    , M(NULL)
{
    build(ref_fn, reads_fn, max_reads, L, R,
          complexity_penalty);
}



sequencing_bias::sequencing_bias(const char* ref_fn,
                                 pos_table* T,
                                 size_t max_reads,
                                 pos_t L, pos_t R,
                                 double complexity_penalty)
    : ref_f(NULL)
    , M(NULL)
{
    build(ref_fn, T, max_reads, L, R, complexity_penalty);
}



void sequencing_bias::to_yaml(YAML::Emitter& out) const
{
    out << YAML::BeginMap;

    out << YAML::Key   << "L";
    out << YAML::Value << L;

    out << YAML::Key   << "R";
    out << YAML::Value << R;

    out << YAML::Key   << "motif";
    out << YAML::Value;
    M->to_yaml(out);

    out << YAML::EndMap;
}


void sequencing_bias::save_to_file(const char* fn) const
{
    FILE* f = fopen(fn, "w");
    if (f == NULL) {
        logger::abort("Can't open file %s for writing.", fn);
    }

    YAML::Emitter out;
    this->to_yaml(out);

    fputs(out.c_str(), f);
    fclose(f);
}


void sequencing_bias::clear()
{
    if (ref_f) {
        fai_destroy(ref_f);
        ref_f = NULL;
    }
    ref_fn.clear();

    delete M;
    M = NULL;
}



void sequencing_bias::build(const char* ref_fn,
                            const char* reads_fn,
                            size_t max_reads,
                            pos_t L, pos_t R,
                            double complexity_penalty)
{
    samfile_t* reads_f = samopen(reads_fn, "rb", NULL);
    if (reads_f == NULL) {
        logger::abort("Can't open bam file '%s'.", reads_fn);
    }

    bam_index_t* reads_index = bam_index_load(reads_fn);
    if (reads_index == NULL) {
        logger::abort("Can't open bam index '%s.bai'.", reads_fn);
    }

    bam_init_header_hash(reads_f->header);

    bam1_t* read = bam_init1();

    pos_table T;
    pos_table_create(&T, reads_f->header->n_targets);
    T.seq_names = reads_f->header->target_name;

    size_t k = 0;

    while (samread(reads_f, read) > 0) {
        if (read->core.n_cigar != 1) continue;
        if (read->core.flag & BAM_FREAD2) continue;
        k++;
        if (k % 1000000 == 0) {
            logger::info("hashed %zu reads.", k);
        }
        pos_table_inc(&T, read);
    }
    logger::info("hashed %zu reads.", k);

    bam_destroy1(read);

    build(ref_fn, &T, max_reads, L, R, complexity_penalty);

    pos_table_destroy(&T);

    bam_index_destroy(reads_index);
    samclose(reads_f);
}


void sequencing_bias::build(const char* ref_fn,
                            pos_table* T, size_t max_reads,
                            pos_t L, pos_t R,
                            double complexity_penalty)
{
    clear();

    this->ref_fn = ref_fn;
    
    this->L = L;
    this->R = R;

    unsigned int i;

    read_pos* S_tmp;
    read_pos* S;
    size_t N;
    const size_t max_dump = 10000000;
    pos_table_dump(T, &S_tmp, &N, max_dump);

    /* sort by count */
    qsort(S_tmp, N, sizeof(read_pos), read_pos_count_compare);

    /* consider only reads with at least one duplicate */
    for (i = 0; i < N; ++i) {
        if (S_tmp[i].count <= 1) break;
    }

    /* (unless there are very few of these reads */
    if (i > 10000) {
        max_reads = std::min<size_t>(max_reads, i);
        logger::info("%zu reads with duplicates.", i);
    }
    else {
        i = N;
    }


    /* ignore the top 1%, as they tend to be vastly higher than anything else,
     * and thus can throw things off when training with small numbers of reads
     * */
    S = S_tmp + i/100;
    max_reads = min<size_t>(max_reads, 99*i/100); 

    /* sort by tid (so we can load one chromosome at a time) */
    qsort(S, std::min<size_t>(max_reads, N), sizeof(read_pos), read_pos_tid_compare);


    /* sample foreground and background kmer frequencies */
    ref_f = fai_load(ref_fn);
    if (ref_f == NULL) {
        logger::abort("Can't open fasta file '%s'.", ref_fn);
    }

    std::deque<twobitseq*> foreground_seqs;
    std::deque<twobitseq*> background_seqs;


    /* background sampling */
    int bg_samples = 2; // make this many samples for each read
    int bg_sample_num;  // keep track of the number of samples made
    pos_t bg_pos;

    
    char*          seqname   = NULL;
    int            seqlen    = 0;
    int            curr_tid  = -1;
    char*          seq       = NULL;

    char* local_seq;
    local_seq = new char[ L + R + 2 ];
    local_seq[L+R+1] = '\0';


    for (i = 0; i < N && i < max_reads; i++) {

        /* Load/switch sequences (chromosomes) as they are encountered in the
         * read stream. The idea here is to avoid thrashing by loading a large
         * sequence, but also avoid overloading memory by only loading one
         * chromosome at a time. */
        if (S[i].tid != curr_tid) {
            seqname = T->seq_names[ S[i].tid ];
            if (seq) free(seq); 

            seq = faidx_fetch_seq(ref_f, seqname, 0, 1000000000, &seqlen);

            logger::info("read sequence %s.", seqname);

            if (seq == NULL) {
                logger::warn("warning: reference sequence not found, skipping.");
            }
            else {
                for (char* c = seq; *c; c++) *c = tolower(*c);
            }

            curr_tid = S[i].tid;
        }

        if (seq == NULL) continue;


        /* add a foreground sequence */
        if (S[i].strand == strand_neg) {
            if (S[i].pos < R || S[i].pos >= seqlen - L) continue;
            memcpy(local_seq, seq + S[i].pos - R, (L+1+R)*sizeof(char));
            seqrc(local_seq, L+1+R);
        }
        else {
            if (S[i].pos < L || S[i].pos >= seqlen - R) continue;
            memcpy(local_seq, seq + (S[i].pos-L), (L+1+R)*sizeof(char));
        }

        if (strchr(local_seq, 'n') != NULL) continue;

        foreground_seqs.push_back(new twobitseq(local_seq));


        /* add a background sequence */
        /* adjust the current read position randomly, and sample */
        for (bg_sample_num = 0; bg_sample_num < bg_samples;) {

            bg_pos = S[i].pos + (pos_t)round_away(rand_trunc_gauss(12, -200, 200));

            if (S[i].strand == strand_neg) {
                if (bg_pos < R || bg_pos >= seqlen - L) continue;
                memcpy(local_seq, seq + bg_pos - R, (L+1+R)*sizeof(char));
                seqrc(local_seq, L+1+R);
            }
            else {
                if (bg_pos < L || bg_pos >= seqlen - R) continue;
                memcpy(local_seq, seq + (bg_pos-L), (L+1+R)*sizeof(char));
            }

            if (strchr(local_seq, 'n') != NULL) continue;

            background_seqs.push_back(new twobitseq(local_seq));
            bg_sample_num++;
        }
    }

    size_t max_parents  = 4;
    size_t max_distance = 10;

    /* A bit of a hack: if we are training on very few reads (a couple thousand,
     * as a opposed to tens of thousands), we tend to end up with too sparse of
     * a model. */
    if (foreground_seqs.size() < 10000) complexity_penalty = 0.25;

    M = new motif(background_seqs,
                  foreground_seqs,
                  L + 1 + R,
                  max_parents,
                  max_distance,
                  complexity_penalty);


    std::deque<twobitseq*>::iterator seqit;
    for (seqit = background_seqs.begin(); seqit != background_seqs.end(); seqit++) {
        delete *seqit;
    }

    for (seqit = foreground_seqs.begin(); seqit != foreground_seqs.end(); seqit++) {
        delete *seqit;
    }

    free(S_tmp);
    free(seq);
    free(local_seq);
}


sequencing_bias::~sequencing_bias()
{
    clear();
}



double* sequencing_bias::get_bias(const char* seqname,
                                  pos_t start, pos_t end,
                                  strand_t strand)
{
    if (strand == strand_na || ref_f == NULL || M == NULL) return NULL;

    pos_t i;
    pos_t seqlen = end - start + 1;

    double* bs = new double [seqlen];
    for (i = 0; i < seqlen; i++) bs[i] = 1.0;

    char* seqstr;

    if (strand == strand_neg) {
        seqstr = faidx_fetch_seq_forced_lower(ref_f, seqname,
                                              start - R, end + L);
        if (seqstr) seqrc(seqstr, seqlen + R + L);
    }
    else {
        seqstr = faidx_fetch_seq_forced_lower(ref_f, seqname,
                                              start - L, end + R);
    }


    if (!seqstr) return bs;

    twobitseq seq(seqstr);

    for (i = 0; i < seqlen; i++) {
        bs[i] = M->eval(seq, i);
    }


    free(seqstr);
    return bs;
}


string sequencing_bias::print_model_graph()
{
    // TODO
    //return M->print_model_graph(L);
    return string();
}



