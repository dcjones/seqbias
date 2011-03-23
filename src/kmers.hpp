
#ifndef PEAKOLATOR_KMERS
#define PEAKOLATOR_KMERS

#include "yaml-cpp/yaml.h"
#include <deque>
#include <set>

/* kmers are encoded in (at least) 16 bits, allowing for k <= 8 */
typedef unsigned short kmer;


/*
 * A class that stores a matrix over kmers, of the form
 *
 *        AAA AAT AAC AAG ATA ...
 * pos 1  a11 a12 a13 ...
 * pos 2  a21 a22 a23 ...
 * pos 3  a31 a32 a33 ...
 *  ...   ...
 */
class kmer_matrix
{
    public:
        kmer_matrix( const YAML::Node& node );
        kmer_matrix( size_t n, size_t k );
        kmer_matrix( const kmer_matrix& );
        ~kmer_matrix();

        void operator=( const kmer_matrix& );
        double& operator()( size_t i, size_t j );

        void to_yaml( YAML::Emitter& out ) const;

        void setall( double x );
        void setrowall( size_t i, double x );

        void getrow( size_t i, double* x );
        void setrow( size_t i, double* x );

        size_t getn() const;
        size_t getm() const;
        size_t getk() const;


        /* normalize to turn each position into a proper distribution over kmers
         * */
        void dist_normalize();
        void dist_normalize_row( size_t i );
        void dist_marginalize( size_t i, size_t j );

        void dist_conditionalize( int effective_k = -1 );
        void dist_conditionalize_row( size_t i, size_t j, int effective_k = -1 );
        void log_transform_row( size_t i, int effective_k = -1 );

    private:

        size_t k; // size of k-mer
        size_t n; // number of positions
        size_t m; // 4^k
        double* A;
};




/*
 * Represent a sequence in 2bit encoding, extract kmers, etc. 
 */
class sequence
{
    public:
        sequence( const char* s, int c = 0, double w = 1.0 );
        sequence( const sequence& );
        void operator=( const sequence& );
        ~sequence();

        kmer get( size_t i ) const;
        bool get( const bool* indexes, size_t maxn, kmer& K, size_t offset = 0 ) const;

        int c; /* class (i.e. foreground (0) or foreground (1)) */
        double w;

    private:
        kmer* xs;
        size_t n;

        static const size_t kmer_max_k;
};





/* A 'bayesian' network representing sequence probability. */


class motif
{
    public:
        motif( const YAML::Node& );
        motif( size_t n, size_t k, int c );
        motif( const motif& );
        ~motif();

        void to_yaml( YAML::Emitter& ) const;

        void add_edge( size_t i, size_t j, const std::deque<sequence*>* data );
        void remove_edge( size_t i, size_t j, const std::deque<sequence*>* data );

        double eval( const sequence&, size_t offset = 0 ) const;
        double eval_node( size_t i, const std::deque<sequence*>* data,
                          size_t offset = 0 ) const;

        size_t num_params() const;

        char* print_model_graph( int offset = 0 );

        int c; /* which subset of the training data to consider */

    private:

        size_t num_parents( size_t i ) const;
        bool has_edge( size_t i, size_t j );
        void set_edge( size_t i, size_t j, bool );

        bool reachable( size_t i, size_t j );
        void compute_reachability();


        void update_likelihood_column( double* L, size_t n, size_t m, size_t j,
                                       const std::deque<sequence*>* training_seqs );

        size_t n; /* number of positions */
        size_t k; /* maximum number of edges */
        kmer_matrix* P;

        bool* parents;
        bool* R; // reachability

        static const double pseudocount;


        friend void train_motifs( motif& M0, motif& M1,
                                  const std::deque<sequence*>* training_seqs,
                                  size_t max_dep_dist, double complexity_penalty );

};

void train_motifs( motif& M0, motif& M1,
                   const std::deque<sequence*>* training_seqs,
                   size_t max_dep_dist = 0, double complexity_penalty = 1.0 );




#endif


