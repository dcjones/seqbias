
#include "sequencing_bias.hpp"
#include "logger.h"
#include "common.hpp"
#include "miscmath.hpp"

#include <cmath>
#include <cctype>
#include <ctime>
#include "samtools/faidx.h"

#include <fstream>
#include <algorithm>
using namespace std;


/* simply uniform random numbers */
static double rand_uniform( double a, double b )
{
    return a + b * (double)rand() / (double)RAND_MAX;
}

/* random gaussians (copied from GSL, to avoid dependency) */
static double rand_gauss( double sigma )
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


static double rand_trunc_gauss( double sigma, double a, double b )
{
    double x;
    do {
        x = rand_gauss( sigma );
    } while( x < a || x > b );

    return x;
}


double gauss_pdf (const double x, const double sigma)
{
  double u = x / fabs (sigma);
  double p = (1 / (sqrt (2 * M_PI) * fabs (sigma))) * exp (-u * u / 2);
  return p;
}


/* pseudocount used when sampling foreground and background nucleotide * frequencies */
const double sequencing_bias::pseudocount = 1;


int read_pos_tid_compare( const void* p1, const void* p2 )
{
    return ((read_pos*)p1)->tid - ((read_pos*)p2)->tid;
}

int read_pos_count_compare( const void* p1, const void* p2 )
{
    return (int)((read_pos*)p2)->count - (int)((read_pos*)p1)->count;
}


int read_pos_tid_count_compare( const void* p1, const void* p2 )
{
    int c = (int)(((read_pos*)p1)->tid - ((read_pos*)p2)->tid);

    if( c == 0 ) {
        return (int)((read_pos*)p2)->count - (int)((read_pos*)p1)->count;
    }
    else return c;
}



sequencing_bias::sequencing_bias()
    : ref_f(NULL)
    , ref_fn(NULL)
    , M0(NULL), M1(NULL)
{}

sequencing_bias::sequencing_bias( const char* ref_fn,
                                  const char* model_fn )
{
    std::ifstream fin;
    fin.open( model_fn );

    YAML::Parser parser( fin );
    YAML::Node doc;
    parser.GetNextDocument( doc );

    doc["L"] >> L;
    doc["R"] >> R;

    M0 = new motif( doc["motif0"] );
    M1 = new motif( doc["motif1"] );

    this->ref_fn = ref_fn ? strdup(ref_fn) : NULL;
    if( ref_fn ) {
        ref_f = fai_load(ref_fn);
        if( ref_f == NULL ) {
            failf( "Can't open fasta file '%s'\n", ref_fn );
        }
    }
    else ref_f = NULL;
}


sequencing_bias::sequencing_bias( const char* ref_fn,
                                  const char* reads_fn,
                                  size_t max_reads, pos L, pos R,
                                  double complexity_penalty,
                                  double offset_std )
    : ref_f(NULL)
    , ref_fn(NULL)
    , M0(NULL), M1(NULL)
{
    build( ref_fn, reads_fn, max_reads, L, R,
           complexity_penalty, offset_std );
}



sequencing_bias::sequencing_bias( const char* ref_fn,
                                  table* T, size_t max_reads,
                                  pos L, pos R,
                                  double complexity_penalty,
                                  double offset_std )
    : ref_f(NULL)
    , ref_fn(NULL)
    , M0(NULL), M1(NULL)
{
    build( ref_fn, T, max_reads, L, R, complexity_penalty, offset_std );
}



sequencing_bias* sequencing_bias::copy() const
{
    sequencing_bias* sb = new sequencing_bias();
    sb->L = L;
    sb->R = R;

    if( M0 && M1 ) {
        sb->M0 = new motif( *M0 );
        sb->M1 = new motif( *M1 );
    }

    sb->ref_fn   = ref_fn   ? strdup(ref_fn)   : NULL;

    if( ref_fn ) {
        sb->ref_f = fai_load(ref_fn);
        if( sb->ref_f == NULL ) {
            failf( "Can't open fasta file '%s'\n", ref_fn );
        }
    }
    else sb->ref_f = NULL;

    return sb;
}


void sequencing_bias::to_yaml( YAML::Emitter& out ) const
{
    out << YAML::BeginMap;

    out << YAML::Key   << "L";
    out << YAML::Value << L;

    out << YAML::Key   << "R";
    out << YAML::Value << R;

    out << YAML::Key   << "motif0";
    out << YAML::Value;
    M0->to_yaml( out );

    out << YAML::Key   << "motif1";
    out << YAML::Value;
    M1->to_yaml( out );

    out << YAML::EndMap;
}


void sequencing_bias::save_to_file( const char* fn ) const
{
    FILE* f = fopen( fn, "w" );
    if( f == NULL ) {
        failf( "Can\'t open file %s\n", fn );
    }

    YAML::Emitter out;
    this->to_yaml( out );

    fputs( out.c_str(), f );
    fclose( f );
}


void sequencing_bias::clear()
{
    if( ref_f ) {
        fai_destroy(ref_f);
        ref_f = NULL;
    }
    free(ref_fn);  ref_fn   = NULL;

    delete M0;
    delete M1;
    M0 = M1 = NULL;
}



void sequencing_bias::build( const char* ref_fn,
                             const char* reads_fn,
                             size_t max_reads, pos L, pos R,
                             double complexity_penalty,
                             double offset_std )
{
    log_printf( LOG_MSG, "hashing reads ... \n" );
    log_indent();

    samfile_t* reads_f = samopen( reads_fn, "rb", NULL );
    if( reads_f == NULL ) {
        failf( "Can't open bam file '%s'.", reads_fn );
    }

    bam_index_t* reads_index = bam_index_load( reads_fn );
    if( reads_index == NULL ) {
        failf( "Can't open bam index '%s.bai'.", reads_fn );
    }

    bam_init_header_hash( reads_f->header );

    bam1_t* read = bam_init1();

    table T;
    table_create( &T, reads_f->header->n_targets );
    T.seq_names = reads_f->header->target_name;

    size_t k = 0;

    while( samread( reads_f, read ) > 0 ) {
        k++;
        if( k % 1000000 == 0 ) {
            log_printf( LOG_MSG, "%zu reads\n", k );
        }
        table_inc( &T, read );
    }
    log_printf( LOG_MSG, "%zu reads\n", k );

    bam_destroy1(read);

    log_puts( LOG_MSG, "done.\n" );
    log_unindent();

    build( ref_fn, &T, max_reads, L, R, complexity_penalty, offset_std );

    table_destroy( &T );

    bam_index_destroy(reads_index);
    samclose( reads_f );
}


void sequencing_bias::build( const char* ref_fn,
                             table* T, size_t max_reads,
                             pos L, pos R,
                             double complexity_penalty,
                             double offset_std )
{
    clear();

    this->ref_fn = strdup(ref_fn);
    
    this->L = L;
    this->R = R;

    unsigned int i;


    log_puts( LOG_MSG, "shuffling ... " );

    read_pos* S;
    size_t N;
    const size_t max_dump = 1000000;
    table_dump( T, &S, &N, max_dump );

    /* shuffle, then sort by position */
    qsort( S, N, sizeof(read_pos), read_pos_count_compare );
    qsort( S, std::min<size_t>( max_reads, N ), sizeof(read_pos), read_pos_tid_compare );

    log_puts( LOG_MSG, "done.\n" );


    /* sample foreground and background kmer frequencies */
    log_puts( LOG_MSG, "sampling sequence bias ...\n" );
    log_indent();


    ref_f = fai_load(ref_fn);
    if( ref_f == NULL ) {
        failf( "Can't open fasta file '%s'\n", ref_fn );
    }

    std::deque<sequence*> training_seqs;


    /* background sampling */
    int bg_samples = 2; // make this many samples for each read
    int bg_sample_num;           // keep track of the number of samples made
    pos bg_pos;

    
    char*          seqname   = NULL;
    int            seqlen    = 0;
    int            curr_tid  = -1;
    char*          seq       = NULL;

    char* local_seq;
    local_seq = new char[ L + R + 2 ];
    local_seq[L+R+1] = '\0';


    for( i = 0; i < N && i < max_reads; i++ ) {

        /* Load/switch sequences (chromosomes) as they are encountered in the
         * read stream. The idea here is to avoid thrashing by loading a large
         * sequence, but also avoid overloading memory by only loading one
         * chromosome at a time. */
        if( S[i].tid != curr_tid ) {
            seqname = T->seq_names[ S[i].tid ];
            if( seq ) free(seq); 

            log_printf( LOG_MSG, "reading sequence '%s' ... ", seqname );

            seq = faidx_fetch_seq( ref_f, seqname, 0, 1000000000, &seqlen );

            if( seq == NULL ) {
                log_puts( LOG_WARN, "warning: reference sequence not found, skipping.\n" );
            }
            else {
                for( char* c = seq; *c; c++ ) *c = tolower(*c);
            }

            curr_tid = S[i].tid;

            log_puts( LOG_MSG, "done.\n" );
        }

        if( seq == NULL ) continue;


        /* add a foreground sequence */
        if( S[i].strand ) {
            if( S[i].pos < R ) continue;
            memcpy( local_seq, seq + S[i].pos - R, (L+1+R)*sizeof(char) );
            seqrc( local_seq, L+1+R );
        }
        else {
            if( S[i].pos < L ) continue;
            memcpy( local_seq, seq + (S[i].pos-L), (L+1+R)*sizeof(char) );
        }


        training_seqs.push_back( new sequence( local_seq, 1 ) );


        /* add a background sequence */
        /* adjust the current read position randomly, and sample */
        for( bg_sample_num = 0; bg_sample_num < bg_samples; bg_sample_num++ ) {

            bg_pos = S[i].pos + (pos)ceil( rand_trunc_gauss( offset_std, -100, 100 ) );

            if( S[i].strand ) {
                if( bg_pos < R ) continue;
                memcpy( local_seq, seq + bg_pos - R, (L+1+R)*sizeof(char) );
                seqrc( local_seq, L+1+R );
            }
            else {
                if( bg_pos < L ) continue;
                memcpy( local_seq, seq + (bg_pos-L), (L+1+R)*sizeof(char) );
            }

            training_seqs.push_back( new sequence( local_seq, 0 ) );
        }
    }

    size_t max_k = 4;
    size_t max_dep_dist = 5;
    M0 = new motif( L+1+R, max_k, 0 );
    M1 = new motif( L+1+R, max_k, 1 );

    train_motifs( *M0, *M1, &training_seqs, max_dep_dist, complexity_penalty );


    std::deque<sequence*>::iterator i_seq;
    for( i_seq = training_seqs.begin(); i_seq != training_seqs.end(); i_seq++ ) {
        delete *i_seq;
    }



    log_unindent();

    free(S);
    free(seq);
    free(local_seq);
}


sequencing_bias::~sequencing_bias()
{
    clear();
}



double* sequencing_bias::get_bias( const char* seqname, pos start, pos end, int strand )
{
    if( strand < 0 || ref_f == NULL || M0 == NULL || M1 == NULL ) return NULL;

    pos i;
    pos seqlen = end - start + 1;

    double* bias = new double[ seqlen ];
    for( i = 0; i < seqlen; i++ ) bias[i] = 1.0;

    char* seqstr;

    if( strand == 1 ) {
        seqstr = faidx_fetch_seq_forced_lower( ref_f, seqname,
                                               start-R, end+L );
        if( seqstr ) seqrc( seqstr, seqlen + R + L );
    }
    else {
        seqstr = faidx_fetch_seq_forced_lower( ref_f, seqname,
                                               start-L, end+R );
    }


    if( !seqstr ) return bias;

    sequence* seq = new sequence( seqstr );
    double L0, L1;

    for( i = 0; i < seqlen; i++ ) {
        L0 = M0->eval( *seq, i );
        L1 = M1->eval( *seq, i );

        bias[i] = exp( L1 - L0 );
        if( !isfinite(bias[i]) ) bias[i] = 1.0;
    }

    free(seqstr);
    delete seq;
    return bias;
}


char* sequencing_bias::print_model_graph()
{
    return M0->print_model_graph( L );
}

