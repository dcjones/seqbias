
#include "kmers.hpp"
#include "common.hpp"
#include "miscmath.hpp"
#include "logger.h"
#include "asprintf.h"

#include <string>
#include <cstring>
#include <cstdio>
#include <cmath>

kmer_matrix::kmer_matrix( const YAML::Node& node )
{
    node["k"] >> k;

    node["n"] >> n;
    node["m"] >> m;

    A = new double[ n * m ];

    const YAML::Node& node_A = node["A"];
    size_t i;
    for( i = 0; i < n*m; i++ ) {
        node_A[i] >> A[i];
    }
}

kmer_matrix::kmer_matrix( size_t n, size_t k )
    : k(k), n(n)
{
    m = 1<<(2*k); /* == 4^k */

    A = new double[ n * m ];
}

kmer_matrix::kmer_matrix( const kmer_matrix& M )
{
    k = M.k;

    n = M.n;
    m = M.m;
    A = new double[ n * m ];
    memcpy( (void*)A, (void*)M.A, n * m * sizeof(double) );
}

size_t kmer_matrix::getn() const { return n; }
size_t kmer_matrix::getm() const { return m; }
size_t kmer_matrix::getk() const { return k; }

void kmer_matrix::operator=( const kmer_matrix& M )
{
    k = M.k;

    if( n != M.n || m != M.m ) {

        n = M.n;
        m = M.m;

        delete[] A;
        A = new double[ n * m ];
    }

    memcpy( (void*)A, (void*)M.A, n * m * sizeof(double) );
}


void kmer_matrix::to_yaml( YAML::Emitter& out ) const
{
    out << YAML::BeginMap;

    out << YAML::Key   << "k";
    out << YAML::Value << (unsigned int)k;
    out << YAML::Key   << "n";
    out << YAML::Value << (unsigned int)n;
    out << YAML::Key   << "m";
    out << YAML::Value << (unsigned int)m;
    out << YAML::Key   << "A";
    out << YAML::Flow;
    out << YAML::Value;
    out << YAML::BeginSeq;
    size_t i;
    for( i = 0; i < n * m; i++ ) {
        out << A[i];
    }
    out << YAML::EndSeq;

    out << YAML::EndMap;
}


kmer_matrix::~kmer_matrix()
{
    delete[] A;
}


double& kmer_matrix::operator()( size_t i, size_t j )
{
    return A[ i * m + j ];
}


void kmer_matrix::setall( double x )
{
    size_t i;
    for( i = 0; i < n * m; i++ ) A[i] = x;
}

void kmer_matrix::setrowall( size_t i, double x )
{
    size_t j;
    for( j = 0; j < m; j++ ) {
        A[ i * m + j ] = x;
    }
}

void kmer_matrix::getrow( size_t i, double* xs )
{
    memcpy( (void*)xs, (void*)(A + i * m), m * sizeof(double) );
}

void kmer_matrix::setrow( size_t i, double* xs )
{
    memcpy( (void*)(A + i * m ), (void*)xs, m * sizeof(double) );
}


void kmer_matrix::dist_normalize()
{
    size_t i;
    for( i = 0; i < n; i++ ) {
        dist_normalize_row( i );
    }
}

void kmer_matrix::dist_normalize_row( size_t i )
{
    double z = 0.0;
    size_t j;
    for( j = 0; j < m; j++ ) z += A[ i * m + j ];
    for( j = 0; j < m; j++ ) A[ i * m + j ] /= z;
}


void kmer_matrix::dist_conditionalize( int effective_k )
{
    size_t i;

    for( i = 0; i < n; i++ ) {
        dist_conditionalize_row( i, effective_k );
    }
}


void kmer_matrix::dist_conditionalize_row( size_t i, size_t j, int effective_k  )
{
    if( effective_k <= 0 ) effective_k = m;
    kmer L, R;
    kmer L_max = 1 << (2*j);
    kmer R_max = 1 << (2*(effective_k-j-1));
    kmer K;
    double z;
    kmer nt;

    for( L = 0; L < L_max; L++ ) {
        for( R = 0; R < R_max; R++ ) { 
            z = 0.0;

            for( nt = 0; nt < 4; nt++ ) {
                K = (L<<(2*(effective_k-j))) | (nt<<(2*(effective_k-j-1))) | R;
                z += A[ i*m + K ];
            }

            for( nt = 0; nt < 4; nt++ ) {
                K = (L<<(2*(effective_k-j))) | (nt<<(2*(effective_k-j-1))) | R;
                A[ i*m + K ] /= z;
            }
        }
    }
}


void kmer_matrix::log_transform_row( size_t i, int effective_k )
{
    if( effective_k <= 0 ) effective_k = m;
    size_t j_max = 1 << (2*effective_k);

    size_t j;
    for( j = 0; j < j_max; j++ ) {
        A[ i * m + j ] = log( A[ i * m + j ] );
    }
}


const size_t sequence::kmer_max_k = 4*sizeof(kmer);


sequence::sequence( const char* s, int c, double w )
    : c(c), w(w), xs(NULL), n(0)
{
    n = strlen(s);
    if( n > 0 ) {
        xs = new kmer[ n/kmer_max_k + 1 ];
        memset( xs, 0, (n/kmer_max_k + 1)*sizeof(kmer) );

        size_t i, block, offset;
        for( i = 0; i < n; i++ ) {
            block  = i / kmer_max_k;
            offset = i % kmer_max_k;
            xs[block] = xs[block] | (nt2num(s[i]) << (2*offset));
        }
    }
}

sequence::sequence( const sequence& s )
    : xs(NULL), n(0)
{
    c = s.c;
    w = s.w;
    n = s.n;
    if( n > 0 ) {
        xs = new kmer[ n/kmer_max_k + 1 ];
        memcpy( xs, s.xs, (n/kmer_max_k + 1)*sizeof(kmer) );
    }
}


void sequence::operator=( const sequence& s )
{
    delete[] xs;
    c = s.c;
    n = s.n;
    xs = new kmer[ n/kmer_max_k + 1 ];
    memcpy( xs, s.xs, (n/kmer_max_k + 1)*sizeof(kmer) );
}


sequence::~sequence()
{
    delete[] xs;
}

kmer sequence::get( size_t i ) const
{
    size_t block  = i / kmer_max_k;
    size_t offset = i % kmer_max_k;
    return (xs[block] >> (2*offset)) & 0x3;
}


bool sequence::get( const bool* indexes, size_t maxn, kmer& K, size_t seq_offset ) const
{
    bool nonempty = false;
    size_t block, offset;
    K = 0;
    size_t i;
    for( i = 0; i < maxn; i++ ) {
        if( indexes[i] ) {
            nonempty = true;
            block  = (i+seq_offset) / kmer_max_k;
            offset = (i+seq_offset) % kmer_max_k;
            K = (K<<2) | ((xs[block] >> (2*offset)) & 0x3);
        }
    }

    return nonempty;
}


const double motif::pseudocount = 1;

motif::motif( const YAML::Node& node )
{
    node["n"] >> n;
    node["k"] >> k;
    node["c"] >> c;

    parents = new bool[n*n];
    const YAML::Node& node_parents = node["parents"];
    size_t i;
    int parents_i;
    for( i = 0; i < n*n; i++ ) {
        node_parents[i] >> parents_i;
        parents[i] = (bool)parents_i;
    }

    P = new kmer_matrix( node["P"] );

    R = new bool[n*n];
    memset( R, 0, n*n*sizeof(bool) );
    compute_reachability();
}

motif::motif( size_t n, size_t k, int c )
    : c(c), n(n), k(k)
{
    P = new kmer_matrix( n, k );
    P->setall( 0.0 );

    parents = new bool[n*n];
    memset( parents, 0, n*n*sizeof(bool) );

    R = new bool[n*n];
    memset( R, 0, n*n*sizeof(bool) );
}


motif::motif( const motif& M )
{
    P = new kmer_matrix( *M.P );
    c = M.c;
    n = M.n;
    k = M.k;

    parents = new bool[n*n];
    memcpy( parents, M.parents, n*n*sizeof(bool) );

    R = new bool[n*n];
    memset( R, 0, n*n*sizeof(bool) );
}



motif::~motif()
{
    delete[] parents;
    delete[] R;
    delete P;
}


void motif::to_yaml( YAML::Emitter& out ) const
{
    out << YAML::BeginMap;

    out << YAML::Key   << "n";
    out << YAML::Value << (unsigned int)n;

    out << YAML::Key   << "k";
    out << YAML::Value << (unsigned int)k;

    out << YAML::Key   << "c";
    out << YAML::Value << (unsigned int)c;

    out << YAML::Key << "parents";
    out << YAML::Value;
    out << YAML::Flow << YAML::BeginSeq;
    size_t i;
    for( i = 0; i < n*n; i++ ) {
        out << (parents[i] ? 1 : 0);
    }
    out << YAML::EndSeq;

    out << YAML::Key   << "P";
    out << YAML::Value;
    P->to_yaml( out );

    out << YAML::EndMap;
}


double motif::eval( const sequence& seq, size_t offset ) const
{
    double p = 0.0;

    const size_t n = P->getn();
    size_t i;
    kmer K;
    for( i = 0; i < n; i++ ) {
        if( !seq.get( parents + i*n, n, K, offset ) ) continue;
        p += (*P)( i, K );
    }

    return p;
}


size_t motif::num_params() const
{
    size_t N = 0;
    size_t i;
    for( i = 0; i < n; i++ ) {
        N += (1 << (2*num_parents(i))) - 1;
    }

    return N;
}

size_t motif::num_parents( size_t i ) const
{
    size_t j;
    size_t M = 0;
    for( j = 0; j < n; j++ ) {
        if( parents[i*n+j] ) M++;
    }

    return M;
}

bool motif::has_edge( size_t i, size_t j )
{
    return parents[j*n+i];
}

void motif::set_edge( size_t i, size_t j, bool x )
{
    parents[j*n+i] = x;
}

bool motif::reachable( size_t i, size_t j )
{
    return R[j*n+i];
}

void motif::compute_reachability()
{
    /* find the transitive closure of the parents matrix via flyod-warshall */

    memcpy( R, parents, n*n*sizeof(bool) );

    size_t k, i, j;
    for( k = 0; k < n; k++ ) {
        for( i = 0; i < n; i++ ) {
            for( j = 0; j < n; j++ ) {
                R[j*n+i] = R[j*n+i] || (R[k*n+i] && R[j*n+k]);
            }
        }
    }
}



char* motif::print_model_graph( int offset )
{
    std::string graph_str;
    char* tmp;
    int r;

    graph_str += "digraph {\n";
    graph_str += "splines=\"true\";\n";
    graph_str += "node [shape=\"box\"];\n";

    int i, j;
    for( j = 0; j < (int)n; j++ ) {
        r = asprintf( &tmp, "n%d [label=\"%d\",pos=\"%d,0\",style=\"%s\"];\n",
                            j, j - offset, j*100,
                            parents[j*n+j] ? "solid" : "dotted" 
                            );
        graph_str += tmp;
        free(tmp);
    }


    for( j = 0; j < (int)n; j++ ) {
        if( !parents[j*n+j] ) continue;

        for( i = 0; i < (int)n; i++ ) {
            if( i == j ) continue;
            if( parents[j*n+i] ) {
                r = asprintf( &tmp, "n%d -> n%d;\n", i, j );
                graph_str += tmp;
                free(tmp);
            }

        }

    }

    graph_str += "}\n";
    return strdup(graph_str.c_str());
}


/* make an edge i --> j
 * That is, condition j on i. */
void motif::add_edge( size_t i, size_t j, const std::deque<sequence*>* data )
{
    set_edge( i, j, true );

    P->setrowall( j, 0.0 );
    size_t n_parents = num_parents(j);
    size_t m = 1 << (2*n_parents);
    kmer K;
    for( K = 0; K < m; K++ ) {
        (*P)( j, K ) = pseudocount;
    }


    std::deque<sequence*>::const_iterator seq;
    for( seq = data->begin(); seq != data->end(); seq++ ) {
        if( (*seq)->c == c && (*seq)->get( parents + j*n, n, K ) ) {
            (*P)( j, K ) += (*seq)->w;
        }
    }

    size_t n_pred_parents = 0;
    size_t u;
    for( u = 0; u < j; u++ ) {
        if( has_edge( u, j ) ) n_pred_parents++;
    }

    P->dist_normalize_row( j );
    P->dist_conditionalize_row( j, n_pred_parents, n_parents );
    P->log_transform_row( j, n_parents );

    compute_reachability();
}



void motif::remove_edge( size_t i, size_t j, const std::deque<sequence*>* data )
{
    set_edge( i, j, false );

    P->setrowall( j, 0.0 );
    size_t n_parents = num_parents(j);
    size_t m = 1 << (2*n_parents);
    kmer K;
    for( K = 0; K < m; K++ ) {
        (*P)( j, K ) = pseudocount;
    }


    std::deque<sequence*>::const_iterator seq;
    for( seq = data->begin(); seq != data->end(); seq++ ) {
        if( (*seq)->c == c && (*seq)->get( parents + j*n, n, K ) ) {
            (*P)( j, K )++;
        }
    }

    size_t n_pred_parents = 0;
    size_t u;
    for( u = 0; u < j; u++ ) {
        if( has_edge( u, j ) ) n_pred_parents++;
    }

    P->dist_normalize_row( j );
    P->dist_conditionalize_row( j, n_pred_parents, n_parents );
    P->log_transform_row( j, n_parents );

    compute_reachability();
}





/* Let x_i be the sequence of training example i, c be the class, and M be the
 * model. Then given l0 and l1, where,
 *      l0[i] = log P( x_i | c = 0, M )
 *      l1[i] = log P( x_i | c = 1, M )
 *
 * This function computes the log marginal likelihood of M under a flat (50/50)
 * prior. That is,
 *
 *     log P( c_1, c_2, ..., c_n | x_1, x_2, ..., x_n, M )
 *                      = sum_i log( P( x_i | c_i, M ) / ( P( x_i | c = 0, M ) + P( x_i | c = 1, M ) )
 *
 */
double conditional_likelihood( size_t n, const double* l0, const double* l1, const int* c,
                               double prior )
{
    double l;
    size_t i;

    double p = log(prior);
    double q = log(1-prior);

    l = 0.0;
    for( i = 0; i < n; i++ ) {
        l += (c[i] == 0 ? (q+l0[i]) : (p+l1[i]))
                    - logaddexp( q+l0[i], p+l1[i] );
    }

    return l;
}


void motif::update_likelihood_column( double* L, size_t n, size_t m, size_t j,
                                      const std::deque<sequence*>* training_seqs )
{
    size_t i;
    kmer K;

    for( i = 0; i < n; i++ ) L[ i * m + j ] = 0.0; /* zero column j */

    for( i = 0; i < n; i++ )  {
        if( (*training_seqs)[i]->get( parents + j * m, m, K ) ) {
            L[ i * m + j ] = (*training_seqs)[i]->w * (*P)( j, K );
        }
    }
}




/* various information criterion to try */

/* Akaike Information Criterion */
double aic( double L, double n_obs, double n_params, double c = 1.0 )
{
    return L - c*n_params;
}

/* Akaike Information Criterion with second order small-sample bias correction.
 * (From Hurvich and Tsai (1989)) */
double aicc( double L, double n_obs, double n_params, double c = 1.0 )
{
    return L - c*n_params
             - (n_params + 1) * (n_params + 2) / (n_obs - n_params - 2);
}

/* Bayesian (Schwarz) Information Criterion */
double bic( double L, double n_obs, double n_params, double c = 1.0 )
{
    return 2.0*L - c*n_params*log(n_obs);
}


void train_motifs( motif& M0, motif& M1,
                   const std::deque<sequence*>* training_seqs,
                   size_t max_dep_dist, double complexity_penalty )
{

    log_puts( LOG_MSG, "training motifs (forwards) ...\n" );
    log_indent();

    if( M0.n != M1.n ) {
        failf( "Motif models of mismatching size. (%zu != %zu)\n", M0.n, M1.n );
    }

    double (*compute_ic)( double, double, double, double ) = aicc;


    int i, j;
    int i_start, i_end;

    const size_t n = training_seqs->size();
    const size_t m = M0.n;

    /* likelihood matrices: 
     * Lk_ij gives the likelihood of training example i on node j in model k.
     */
    double* L0 = new double[ n * m ]; memset( (void*)L0, 0, n * m * sizeof(double) );
    double* L1 = new double[ n * m ]; memset( (void*)L1, 0, n * m * sizeof(double) );

    /* likelihood vectors:
     * summed columns of L0, L1, maintained to minimize redundant computation.
     * */
    double* l0 = new double[ n ]; memset( (void*)l0, 0, n * sizeof(double) );
    double* l1 = new double[ n ]; memset( (void*)l1, 0, n * sizeof(double) );


    /* vectors for saving and restoring parameters, to avoid recomputing
     * frequencies */
    double* M0_row_j = new double[ M0.P->getm() ];
    double* M1_row_j = new double[ M0.P->getm() ];
    double* M0_row_i = new double[ M0.P->getm() ];
    double* M1_row_i = new double[ M0.P->getm() ];


    /* 0-1 vector giving labeling each sequence as foreground or background */
    int* cs = new int[ n ];

    for( i = 0; i < (int)n; i++ ) {
        cs[i] = (*training_seqs)[i]->c;
    }

    /* set prior probability of an example being foreground */
    double ws[2];
    for( i = 0; i < (int)n; i++ ) {
        if( cs[i] <= 1 ) ws[cs[i]] += (*training_seqs)[i]->w;
    }
    double prior = ws[1] / (ws[0] + ws[1]);


    /* backup likelihood matrix columns, to restore state after trying a new edge */
    double* b0_j = new double[ n ];
    double* b1_j = new double[ n ];
    double* b0_i = new double[ n ];
    double* b1_i = new double[ n ];
    

    /* keeping track of the optimal edge */
    double ic, ic_curr;

    double ic_forw_best;
    int j_forw_best, i_forw_best;

    double ic_back_best;
    int j_back_best, i_back_best;

    double ic_rev_best;
    int j_rev_best, i_rev_best;


    /* for cycle detection */
    bool has_ij_path;

    /* parameters to compute information criterion */
    double n_obs    = training_seqs->size();
    double n_params = M0.num_params() + M1.num_params();


    /* log conditional likelihood */
    double l; 

    /* baseline ic */
    l = conditional_likelihood( n, l0, l1, cs, prior );
    ic_curr = compute_ic( l, n_obs, n_params, complexity_penalty );


    /* for pretty output */
    size_t col;
    const char* col_base = "\n%34s";
    const size_t col_max = 30;


    size_t round_num = 0;

    while( true ) {
        round_num++;

        log_printf( LOG_MSG, "round %4zu (ic = %0.4e) ", round_num, ic_curr );
        col = 0;

        ic_forw_best = ic_back_best = ic_rev_best = -HUGE_VAL;
        i_forw_best = i_back_best = i_rev_best = 0;
        j_forw_best = j_back_best = j_rev_best = 0;


        /* phase 1: try all possible edge additions */
        for( j = M0.n-1; j >= 0; j-- ) {

            if( M0.has_edge( j, j ) ) {
                if( max_dep_dist == 0 || j <= (int)max_dep_dist ) i_start = 0;
                else i_start = j - max_dep_dist;

                if( max_dep_dist == 0 ) i_end = M0.n - 1;
                else i_end = std::min( M0.n-1, j + max_dep_dist );
            }
            else i_start = i_end = j;



            for( i = i_start; i <= i_end; i++ ) {

                /* skip existing edges  */
                if( M0.has_edge( i, j ) ) {
                    continue;
                }

                /* skip edges that would introduce cycles */
                if( M0.reachable( j, i ) ) {
                    continue;
                }

                /* skip edges that would exceed the parent limit */
                if( M0.num_parents(j) >= M0.k ) {
                    continue;
                }

                /* skip edges that are equivalent to one already tried */
                if( i > j && M0.num_parents(j) == 1 && M0.num_parents(i) == 1 ) {
                    continue;
                }


                log_puts( LOG_MSG, "+" );
                if( ++col > col_max ) {
                    col = 0;
                    log_printf( LOG_MSG, col_base, "" );
                }


                /* keep track of the old parameters to avoid retraining */
                M0.P->getrow( j, M0_row_j );
                M1.P->getrow( j, M1_row_j );

                /* keep track of the old likelihoods to avoid reevaluating */
                colcpy( b0_j, L0, j, n, m  );
                colcpy( b1_j, L1, j, n, m );

                /* add edge */
                M0.add_edge( i, j, training_seqs );
                M1.add_edge( i, j, training_seqs );

                /* evaluate likelihoods for that column */
                M0.update_likelihood_column( L0, n, m, j, training_seqs );
                M1.update_likelihood_column( L1, n, m, j, training_seqs );


                /* update training example likelihoods */
                vecsub( l0, b0_j, n );
                vecaddcol( l0, L0, n, m, j );

                vecsub( l1, b1_j, n );
                vecaddcol( l1, L1, n, m, j );


                l        = conditional_likelihood( n, l0, l1, cs, prior );
                n_params = M0.num_params() + M1.num_params();
                ic       = compute_ic( l, n_obs, n_params, complexity_penalty );


                if( ic >= ic_forw_best ) {
                    ic_forw_best = ic;
                    i_forw_best = i;
                    j_forw_best = j;
                }
        
                /* TODO: delete me */
                //log_printf( LOG_MSG, "\n(%zu,%zu) : %0.4e\n", i, j, ic );

                /* remove edge */
                M0.set_edge( i, j, false );
                M1.set_edge( i, j, false );

                /* restore previous parameters */
                M0.P->setrow( j, M0_row_j );
                M1.P->setrow( j, M1_row_j );

                /* restore previous likelihoods */
                vecsubcol( l0, L0, n, m, j );
                vecadd( l0, b0_j, n );

                vecsubcol( l1, L1, n, m, j );
                vecadd( l1, b1_j, n );

                matsetcol( L0, b0_j, n, m, j );
                matsetcol( L1, b1_j, n, m, j );
            }
        }

        /* phase 2: try all possible edge removals */
        for( j = 0; j < (int)M0.n; j++ ) {
            for( i = 0; i < (int)M0.n; i++ ) {

                if( !M0.has_edge( i, j ) ) continue;
                if( i == j && M0.num_parents(j) > 1 ) continue;

                log_puts( LOG_MSG, "-" );
                if( ++col > col_max ) {
                    col = 0;
                    log_printf( LOG_MSG, col_base, "" );
                }

                /* keep track of the old parameters to avoid retraining */
                M0.P->getrow( j, M0_row_j );
                M1.P->getrow( j, M1_row_j );

                /* keep track of the old likelihoods to avoid reevaluating */
                colcpy( b0_j, L0, j, n, m  );
                colcpy( b1_j, L1, j, n, m );

                /* remove edge */
                M0.remove_edge( i, j, training_seqs );
                M1.remove_edge( i, j, training_seqs );

                /* evaluate likelihoods for that column */
                M0.update_likelihood_column( L0, n, m, j, training_seqs );
                M1.update_likelihood_column( L1, n, m, j, training_seqs );

                /* update training example likelihoods */
                vecsub( l0, b0_j, n );
                vecaddcol( l0, L0, n, m, j );

                vecsub( l1, b1_j, n );
                vecaddcol( l1, L1, n, m, j );

                /* evaluate likelihood / ic */
                l        = conditional_likelihood( n, l0, l1, cs, prior );
                n_params = M0.num_params() + M1.num_params();
                ic       = compute_ic( l, n_obs, n_params, complexity_penalty );

                if( ic > ic_back_best ) {
                    ic_back_best = ic;
                    i_back_best = i;
                    j_back_best = j;
                }

                /* TODO: delete me */
                //log_printf( LOG_MSG, "\n(%zu,%zu) : %0.4e\n", i, j, ic );

                /* replace edge */
                M0.set_edge( i, j, true );
                M1.set_edge( i, j, true );

                /* restore previous parameters */
                M0.P->setrow( j, M0_row_j );
                M1.P->setrow( j, M1_row_j );

                /* restore previous likelihoods */
                vecsubcol( l0, L0, n, m, j );
                vecadd( l0, b0_j, n );

                vecsubcol( l1, L1, n, m, j );
                vecadd( l1, b1_j, n );

                matsetcol( L0, b0_j, n, m, j );
                matsetcol( L1, b1_j, n, m, j );
            }
        }

        /* phase 3: try all possible edge reversals */
        for( j = 0; j < (int)M0.n; j++ ) {
            for( i = 0; i < (int)M0.n; i++ ) {

                /* skip nonsense reversals */
                if( i == j ) continue;

                /* skip reversals that add parameters */
                if( !(M0.has_edge( i, j ) && M0.has_edge( i, i ) && M0.has_edge( j, j )) ) {
                    continue;
                }

                /* skip reversals that would introduce cycles:
                 * this is determined by cuting the (i,j) edge, then recomputing
                 * reachability to see if an i, j path remains */
                M0.set_edge( i, j, false );
                M0.compute_reachability();

                has_ij_path = M0.reachable( i, j );

                M0.set_edge( i, j, true );
                M0.compute_reachability();

                if( has_ij_path ) continue;



                log_puts( LOG_MSG, "." );
                if( ++col > col_max ) {
                    col = 0;
                    log_printf( LOG_MSG, col_base, "" );
                }

                /* keep track of the old parameters to avoid retraining */
                M0.P->getrow( j, M0_row_j );
                M1.P->getrow( j, M1_row_j );
                M0.P->getrow( i, M0_row_i );
                M1.P->getrow( i, M1_row_i );

                /* keep track of the old likelihoods to avoid reevaluating */
                colcpy( b0_j, L0, j, n, m  );
                colcpy( b1_j, L1, j, n, m );
                colcpy( b0_i, L0, i, n, m  );
                colcpy( b1_i, L1, i, n, m );

                /* reverse edge */
                M0.remove_edge( i, j, training_seqs );
                M1.remove_edge( i, j, training_seqs );
                M0.add_edge( j, i, training_seqs );
                M1.add_edge( j, i, training_seqs );

                /* evaluate likelihoods for those columnn */
                M0.update_likelihood_column( L0, n, m, j, training_seqs );
                M1.update_likelihood_column( L1, n, m, j, training_seqs );
                M0.update_likelihood_column( L0, n, m, i, training_seqs );
                M1.update_likelihood_column( L1, n, m, i, training_seqs );

                /* update training example likelihoods */
                vecsub( l0, b0_j, n );
                vecaddcol( l0, L0, n, m, j );
                vecsub( l1, b1_j, n );
                vecaddcol( l1, L1, n, m, j );

                vecsub( l0, b0_i, n );
                vecaddcol( l0, L0, n, m, i );
                vecsub( l1, b1_i, n );
                vecaddcol( l1, L1, n, m, i );


                /* evaluate likelihood / ic */
                l        = conditional_likelihood( n, l0, l1, cs, prior );
                n_params = M0.num_params() + M1.num_params();
                ic       = compute_ic( l, n_obs, n_params, complexity_penalty );

                if( ic > ic_rev_best ) {
                    ic_rev_best = ic;
                    i_rev_best = i;
                    j_rev_best = j;
                }

                /* TODO: delete me */
                //log_printf( LOG_MSG, "\n(%zu,%zu) : %0.4e\n", i, j, ic );

                /* replace edge */
                M0.set_edge( j, i, false );
                M1.set_edge( j, i, false );
                M0.set_edge( i, j, true );
                M1.set_edge( i, j, true );

                /* restore previous parameters */
                M0.P->setrow( j, M0_row_j );
                M1.P->setrow( j, M1_row_j );
                M0.P->setrow( i, M0_row_i );
                M1.P->setrow( i, M1_row_i );

                /* restore previous likelihoods */
                vecsubcol( l0, L0, n, m, j );
                vecadd( l0, b0_j, n );
                vecsubcol( l1, L1, n, m, j );
                vecadd( l1, b1_j, n );

                vecsubcol( l0, L0, n, m, i );
                vecadd( l0, b0_i, n );
                vecsubcol( l1, L1, n, m, i );
                vecadd( l1, b1_i, n );

                matsetcol( L0, b0_j, n, m, j );
                matsetcol( L1, b1_j, n, m, j );
                matsetcol( L0, b0_i, n, m, i );
                matsetcol( L1, b1_i, n, m, i );
            }
        }



        log_puts( LOG_MSG, "\n" );

        if( std::max( std::max( ic_forw_best, ic_back_best ), ic_rev_best ) <= ic_curr ) break;

        if( ic_forw_best > ic_back_best && ic_forw_best > ic_rev_best ) {
            log_printf( LOG_MSG, " [+] %zu->%zu\n", i_forw_best, j_forw_best );

            ic_curr = ic_forw_best;

            M0.add_edge( i_forw_best, j_forw_best, training_seqs );
            M1.add_edge( i_forw_best, j_forw_best, training_seqs );

            vecsubcol( l0, L0, n, m, j_forw_best );
            vecsubcol( l1, L1, n, m, j_forw_best );

            M0.update_likelihood_column( L0, n, m, j_forw_best, training_seqs );
            M1.update_likelihood_column( L1, n, m, j_forw_best, training_seqs );

            vecaddcol( l0, L0, n, m, j_forw_best );
            vecaddcol( l1, L1, n, m, j_forw_best );
        }
        else if( ic_back_best > ic_rev_best ) {
            log_printf( LOG_MSG, " [-] %zu->%zu\n", i_back_best, j_back_best );
            ic_curr = ic_back_best;

            M0.remove_edge( i_back_best, j_back_best, training_seqs );
            M1.remove_edge( i_back_best, j_back_best, training_seqs );

            vecsubcol( l0, L0, n, m, j_back_best );
            vecsubcol( l1, L1, n, m, j_back_best );

            M0.update_likelihood_column( L0, n, m, j_back_best, training_seqs );
            M1.update_likelihood_column( L1, n, m, j_back_best, training_seqs );

            vecaddcol( l0, L0, n, m, j_back_best );
            vecaddcol( l1, L1, n, m, j_back_best );
        }
        else {
            log_printf( LOG_MSG, " [.] %zu->%zu\n", i_rev_best, j_rev_best );
            ic_curr = ic_rev_best;

            M0.remove_edge( i_rev_best, j_rev_best, training_seqs );
            M1.remove_edge( i_rev_best, j_rev_best, training_seqs );
            M0.add_edge( j_rev_best, i_rev_best, training_seqs );
            M1.add_edge( j_rev_best, i_rev_best, training_seqs );

            vecsubcol( l0, L0, n, m, j_rev_best );
            vecsubcol( l1, L1, n, m, j_rev_best );
            vecsubcol( l0, L0, n, m, i_rev_best );
            vecsubcol( l1, L1, n, m, i_rev_best );

            M0.update_likelihood_column( L0, n, m, j_rev_best, training_seqs );
            M1.update_likelihood_column( L1, n, m, j_rev_best, training_seqs );
            M0.update_likelihood_column( L0, n, m, i_rev_best, training_seqs );
            M1.update_likelihood_column( L1, n, m, i_rev_best, training_seqs );

            vecaddcol( l0, L0, n, m, j_rev_best );
            vecaddcol( l1, L1, n, m, j_rev_best );
            vecaddcol( l0, L0, n, m, i_rev_best );
            vecaddcol( l1, L1, n, m, i_rev_best );
        }
    }



    delete[] cs;
    delete[] L0;
    delete[] L1;
    delete[] l0;
    delete[] l1;
    delete[] M0_row_j;
    delete[] M1_row_j;
    delete[] M0_row_i;
    delete[] M1_row_i;
    delete[] b0_j;
    delete[] b1_j;
    delete[] b0_i;
    delete[] b1_i;
    
    log_unindent();
}



