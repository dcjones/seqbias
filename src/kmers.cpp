
#include "kmers.hpp"
#include "common.hpp"
#include "miscmath.hpp"
#include "logger.h"

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

    stored_row = new double[ m ];
}

kmer_matrix::kmer_matrix( size_t n, size_t k )
    : k(k), n(n)
{
    m = 1<<(2*k); /* == 4^k */

    A = new double[ n * m ];
    stored_row = new double[ m ];
}

kmer_matrix::kmer_matrix( const kmer_matrix& M )
{
    k = M.k;

    n = M.n;
    m = M.m;
    A = new double[ n * m ];
    memcpy( (void*)A, (void*)M.A, n * m * sizeof(double) );


    stored_row = new double[ m ];
    memcpy( (void*)stored_row, (void*)M.stored_row, m * sizeof(double) );
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

        delete[] stored_row;
        stored_row = new double[ m ];
    }

    memcpy( (void*)A, (void*)M.A, n * m * sizeof(double) );
    memcpy( (void*)stored_row, (void*)M.stored_row, m * sizeof(double) );
}


void kmer_matrix::to_yaml( YAML::Emitter& out ) const
{
    out << YAML::BeginMap;

    out << YAML::Key   << "k";
    out << YAML::Value << k;
    out << YAML::Key   << "n";
    out << YAML::Value << n;
    out << YAML::Key   << "m";
    out << YAML::Value << m;
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
    delete[] stored_row;
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

void kmer_matrix::setrow( size_t i, double x )
{
    size_t j;
    for( j = 0; j < m; j++ ) {
        A[ i * m + j ] = x;
    }
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


void kmer_matrix::dist_conditionalize_row( size_t i, int effective_k  )
{
    if( effective_k <= 0 ) effective_k = m;

    kmer L;
    kmer L_max = 1 << (2*(effective_k-1));
    kmer K;
    kmer nt;
    double z;

    for( L = 0; L < L_max; L++ ) {
        K = L<<2;
        z = 0.0;
        for( nt = 0; nt < 4; nt++ ) z += A[ i * m + (nt | K) ];
        for( nt = 0; nt < 4; nt++ ) A[ i * m + (nt | K) ] /= z;
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


void kmer_matrix::store_row( size_t i )
{
    memcpy( (void*)stored_row, (void*)(A + i * m), m * sizeof(double) );
    stored_row_index = i;
}


void kmer_matrix::restore_stored_row()
{
    memcpy( (void*)(A + stored_row_index * m), (void*)stored_row, m * sizeof(double) );
}


const size_t sequence::kmer_max_k = 4*sizeof(kmer);


sequence::sequence( const char* s, int c )
    : c(c), xs(NULL), n(0)
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
}

motif::motif( size_t n, size_t k, int c )
    : c(c), n(n), k(k)
{
    P = new kmer_matrix( n, k );
    P->setall( 0.0 );

    parents = new bool[n*n];
    memset( parents, 0, n*n*sizeof(bool) );
}


motif::motif( const motif& M )
{
    P = new kmer_matrix( *M.P );
    c = M.c;
    n = M.n;
    k = M.k;

    parents = new bool[n*n];
    memcpy( parents, M.parents, n*n*sizeof(bool) );
}



motif::~motif()
{
    delete[] parents;
    delete P;
}


void motif::to_yaml( YAML::Emitter& out ) const
{
    out << YAML::BeginMap;

    out << YAML::Key   << "n";
    out << YAML::Value << n;

    out << YAML::Key   << "k";
    out << YAML::Value << k;

    out << YAML::Key   << "c";
    out << YAML::Value << c;

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
    for( j = 0; j <= i; j++ ) {
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


void motif::store_row( size_t i )
{
    P->store_row( i );
}


void motif::restore_stored_row()
{
    P->restore_stored_row();
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

        for( i = 0; i < j; i++ ) {
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


void motif::add_all_edges( const std::deque<sequence*>* data )
{
    size_t i, j;
    size_t n_parents;
    size_t m;
    kmer K;
    std::deque<sequence*>::const_iterator seq;

    for( j = 0; j < n; j++ ) {
        for( i = k-1 <= j ? j-(k-1) : 0; i <= j; i++ ) {
            set_edge( i, j, true );
        }
        P->setrow( j, 0.0 );
        n_parents = num_parents(j);
        m = 1 << (2*n_parents);

        for( K = 0; K < m; K++ ) {
            (*P)( j, K ) = pseudocount;
        }

        for( seq = data->begin(); seq != data->end(); seq++ ) {
            if( (*seq)->c == c && (*seq)->get( parents + j*n, n, K ) ) {
                (*P)( j, K )++;
            }
        }

        P->dist_normalize_row( j );
        P->dist_conditionalize_row( j, n_parents );
        P->log_transform_row( j, n_parents );
    }
}


/* make an edge i --> j
 * That is, condition j on i. */
void motif::add_edge( size_t i, size_t j, const std::deque<sequence*>* data )
{
    if( i > j ) {
        failf( "Invalid motif edge (%zu, %zu)\n", i, j );
    }

    set_edge( i, j, true );

    P->setrow( j, 0.0 );
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

    P->dist_normalize_row( j );
    P->dist_conditionalize_row( j, n_parents );
    P->log_transform_row( j, n_parents );
}



void motif::remove_edge( size_t i, size_t j, const std::deque<sequence*>* data )
{
    if( i > j ) {
        failf( "Invalid motif edge (%zu, %zu)\n", i, j );
    }

    set_edge( i, j, false );

    P->setrow( j, 0.0 );
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

    P->dist_normalize_row( j );
    P->dist_conditionalize_row( j, n_parents );
    P->log_transform_row( j, n_parents );
}





double motif_log_likelihood( const motif& M0, const motif& M1,
                             const std::deque<sequence*>* training_seqs )
{
    double L, L0, L1;

    L = 0.0;
    
    std::deque<sequence*>::const_iterator i;
    for( i = training_seqs->begin(); i != training_seqs->end(); i++ ) {
        L0 = M0.eval( **i );
        L1 = M1.eval( **i );
        if( (*i)->c == 0 ) {
            L += L0 - logaddexp( L0, L1 );
        }
        else if( (*i)->c == 1 ) {
            L += L1 - logaddexp( L0, L1 );
        }
    }

    return L;
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
double conditional_likelihood( size_t n, const double* l0, const double* l1, const int* c )
{
    double l;
    size_t i;

    l = 0.0;
    for( i = 0; i < n; i++ ) {
        l += (c[i] == 0 ? l0[i] : l1[i])
                    - logaddexp( l0[i], l1[i] );
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
            L[ i * m + j ] = (*P)( j, K );
        }
    }
}




/* various information criterion to try */

/* Akaike Information Criterion */
double aic( double L, double n_obs, double n_params, double c = 1.0 )
{
    return L - c*n_params*n_obs / (n_obs - n_params - 1.0);
}

/* Bayesian (Schwarz) Information Criterion */
double bic( double L, double n_obs, double n_params, double c = 1.0 )
{
    return 2.0*L - c*n_params*log(n_obs);
}


/* Hannan-Quinn Information Criterion */
double hqic( double L, double n_obs, double n_params, double c = 1.0 )
{
    return L - c*n_params*log(log(n_obs));
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

    double (*compute_ic)( double, double, double, double ) = aic;


    size_t i, j;
    size_t i_start;

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


    /* 0-1 vector giving labeling each sequence as foreground or background */
    int* cs = new int[ n ];

    for( i = 0; i < n; i++ ) {
        cs[i] = (*training_seqs)[i]->c;
    }


    /* backup likelihood matrix columns, to restore state after trying a new edge */
    double* b0 = new double[ n ];
    double* b1 = new double[ n ];
    

    /* keeping track of the optimal edge */
    double ic, ic_curr, ic_best;
    size_t j_best, i_best;


    /* parameters to compute information criterion */
    double n_obs    = training_seqs->size();
    double n_params = M0.num_params() + M1.num_params();


    /* log conditional likelihood */
    double l; 


    /* baseline ic */
    l = conditional_likelihood( n, l0, l1, cs );
    ic_curr = compute_ic( l, n_obs, n_params, complexity_penalty );


    size_t round_num = 0;

    while( true ) {
        round_num++;

        log_printf( LOG_MSG, "round %4d (ic = %0.4e) ", round_num, ic_curr );

        ic_best = -HUGE_VAL;
        j_best = i_best = 0;

        for( j = 0; j < M0.n; j++ ) {

            if( M0.has_edge( j, j ) ) {
                if( max_dep_dist == 0 || j <= max_dep_dist ) i_start = 0;
                else i_start = i - max_dep_dist;
            }
            else i_start = j;


            for( i = i_start; i <= j; i++ ) {

                if( M0.has_edge( i, j ) ) {
                    continue;
                }

                if( M0.num_parents(j) >= M0.k ) {
                    continue;
                }

                log_puts( LOG_MSG, "." );

                /* keep track of the old parameters to avoid retraining */
                M0.store_row(j);
                M1.store_row(j);


                /* keep track of the old likelihoods to avoid reevaluating */
                colcpy( b0, L0, j, n, m  );
                colcpy( b1, L1, j, n, m );

                /* add edge */
                M0.add_edge( i, j, training_seqs );
                M1.add_edge( i, j, training_seqs );

                /* evaluate likelihoods for that column */
                M0.update_likelihood_column( L0, n, m, j, training_seqs );
                M1.update_likelihood_column( L1, n, m, j, training_seqs );


                /* update training example likelihoods */
                vecsub( l0, b0, n );
                vecaddcol( l0, L0, n, m, j );

                vecsub( l1, b1, n );
                vecaddcol( l1, L1, n, m, j );


                l        = conditional_likelihood( n, l0, l1, cs );
                n_params = M0.num_params() + M1.num_params();
                ic       = compute_ic( l, n_obs, n_params, complexity_penalty );

                if( ic > ic_best ) {
                    ic_best = ic;
                    i_best = i;
                    j_best = j;
                }
        

                /* remove edge */
                M0.set_edge( i, j, false );
                M1.set_edge( i, j, false );

                /* restore previous parameters */
                M0.restore_stored_row();
                M1.restore_stored_row();

                /* restore previous likelihoods */
                vecsubcol( l0, L0, n, m, j );
                vecadd( l0, b0, n );

                vecsubcol( l1, L1, n, m, j );
                vecadd( l1, b1, n );

                matsetcol( L0, b0, n, m, j );
                matsetcol( L1, b1, n, m, j );
            }
        }

        log_puts( LOG_MSG, "\n" );

        if( ic_best <= ic_curr || ic_best == -HUGE_VAL ) break;

        ic_curr = ic_best;

        M0.add_edge( i_best, j_best, training_seqs );
        M1.add_edge( i_best, j_best, training_seqs );

        vecsubcol( l0, L0, n, m, j_best );
        vecsubcol( l1, L1, n, m, j_best );

        M0.update_likelihood_column( L0, n, m, j_best, training_seqs );
        M1.update_likelihood_column( L1, n, m, j_best, training_seqs );

        vecaddcol( l0, L0, n, m, j_best );
        vecaddcol( l1, L1, n, m, j_best );
    }



    delete[] cs;
    delete[] L0;
    delete[] L1;
    delete[] l0;
    delete[] l1;
    delete[] b0;
    delete[] b1;
    
    log_unindent();
}




void train_motifs_backwards( motif& M0, motif& M1,
                             const std::deque<sequence*>* training_seqs,
                             size_t max_dep_dist, double complexity_penalty )
{
    log_puts( LOG_MSG, "training motifs (backwards) ...\n" );
    log_indent();


    if( M0.n != M1.n ) {
        failf( "Motif models of mismatching size. (%zu != %zu)\n", M0.n, M1.n );
    }

    double (*compute_ic)( double, double, double, double ) = aic;

    size_t i, j;
    size_t i_last;

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


    /* 0-1 vectors giving labeling each sequence as foreground or background */
    int* cs = new int[ training_seqs->size() ];

    for( i = 0; i < training_seqs->size(); i++ ) {
        cs[i] = (*training_seqs)[i]->c;
    }


    /* backup likelihood matrix columns, to restore state after trying a new edge */
    double* b0 = new double[ n ];
    double* b1 = new double[ n ];
    


    /* keeping track of the optimal edge */
    double ic, ic_curr, ic_best;
    size_t j_best, i_best;


    /* start with all edges */
    M0.add_all_edges( training_seqs );
    M1.add_all_edges( training_seqs );


    /* parameters to compute information criterion */
    double n_obs    = training_seqs->size();
    double n_params = M0.num_params() + M1.num_params();


    /* log conditional likelihood */
    double l; 

    /* baseline ic */
    l = conditional_likelihood( n, l0, l1, cs );
    ic_curr = compute_ic( l, n_obs, n_params, complexity_penalty );

    size_t round = 0;

    while( true ) {
        round++;

        log_printf( LOG_MSG, "round %4d (ic = %0.4e) ", round, ic_curr );

        ic_best = -HUGE_VAL;
        j_best = i_best = 0;

        for( j = 0; j < M0.n; j++ ) {

            /* do not remove the position unless it has no edges */
            i_last = M0.num_parents(j) > 1 ? j-1 : j;

            for( i = 0; i <= i_last; i++ ) {

                if( !M0.has_edge( i, j ) ) continue;

                log_puts( LOG_MSG, "." );

                /* keep track of the old parameters to avoid retraining */
                M0.store_row(j);
                M1.store_row(j);

                /* keep track of the old likelihoods to avoid reevaluating */
                colcpy( b0, L0, j, n, m  );
                colcpy( b1, L1, j, n, m );

                /* remove edge */
                M0.remove_edge( i, j, training_seqs );
                M1.remove_edge( i, j, training_seqs );

                /* evaluate likelihoods for that column */
                M0.update_likelihood_column( L0, j, n, m, training_seqs );
                M1.update_likelihood_column( L1, j, n, m, training_seqs );

                /* update training example likelihoods */
                vecsub( l0, b0, n );
                vecaddcol( l0, L0, n, m, j );

                vecsub( l1, b1, n );
                vecaddcol( l1, L1, n, m, j );

                /* evaluate likelihood / ic */
                l        = conditional_likelihood( n, l0, l1, cs );
                n_params = M0.num_params() + M1.num_params();
                ic       = compute_ic( l, n_obs, n_params, complexity_penalty );

                if( ic > ic_best ) {
                    ic_best = ic;
                    i_best = i;
                    j_best = j;
                }

                /* replace edge */
                M0.set_edge( i, j, true );
                M1.set_edge( i, j, true );

                /* restore previous parameters */
                M0.restore_stored_row();
                M1.restore_stored_row();

                /* restore previous likelihoods */
                vecsubcol( l0, L0, n, m, j );
                vecadd( l0, b0, n );

                vecsubcol( l1, L1, n, m, j );
                vecadd( l1, b1, n );

                matsetcol( L0, b0, n, m, j );
                matsetcol( L1, b1, n, m, j );
            }
        }

        log_puts( LOG_MSG, "\n" );


        if( ic_best <= ic_curr || ic_best == -HUGE_VAL ) break;

        ic_curr = ic_best;

        M0.remove_edge( i_best, j_best, training_seqs );
        M1.remove_edge( i_best, j_best, training_seqs );

        vecsubcol( l0, L0, n, m, j_best );
        vecsubcol( l1, L1, n, m, j_best );

        M0.update_likelihood_column( L0, n, m, j_best, training_seqs );
        M1.update_likelihood_column( L1, n, m, j_best, training_seqs );

        vecaddcol( l0, L0, n, m, j_best );
        vecaddcol( l1, L1, n, m, j_best );
    }

    

    delete[] cs;
    delete[] L0;
    delete[] L1;
    delete[] l0;
    delete[] l1;
    delete[] b0;
    delete[] b1;

    log_unindent();
}


