
#include "sequencing_bias.hpp"
#include "logger.h"
#include "common.hpp"

#include <cmath>
#include <cctype>
#include <ctime>
#include "samtools/faidx.h"

#include <fstream>
#include <algorithm>
using namespace std;


/* simply uniform random numbers */
double rand_uniform( double a, double b )
{
    return a + b * (double)rand() / (double)RAND_MAX;
}


/* pseudocount used when sampling foreground and background nucleotide * frequencies */
const double sequencing_bias::pseudocount = 1;


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
                                  size_t n, pos L, pos R,
                                  bool   train_backwards,
                                  double complexity_penalty )
    : ref_f(NULL)
    , ref_fn(NULL)
    , M0(NULL), M1(NULL)
{
    build( ref_fn, reads_fn, n, L, R, train_backwards, complexity_penalty );
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
                             size_t n, pos L, pos R,
                             bool   train_backwards,
                             double complexity_penalty )
{
    log_puts( LOG_MSG, "Determining sequencing bias...\n" );
    log_indent();

    clock_t t0, t1;
    t0 = clock();

    clear();

    this->ref_fn   = strdup(ref_fn);
    
    this->L = L;
    this->R = R;

    unsigned int i;

    samfile_t* reads_f = samopen( reads_fn, "rb", NULL );
    if( reads_f == NULL ) {
        failf( "Can't open bam file '%s'.\n", reads_fn );
    }

    table T;
    hash_reads( &T, reads_f, n );
    n = T.m;

    /* resort the remaining (1-q)*n by position */
    log_puts( LOG_MSG, "sorting by position ... " );
    struct hashed_value** S;
    table_sort_by_position( &T, &S );
    log_puts( LOG_MSG, "done.\n" );


    /* sample foreground and background kmer frequencies */
    log_puts( LOG_MSG, "sampling sequence bias ...\n" );
    log_indent();


    ref_f = fai_load(ref_fn);
    if( ref_f == NULL ) {
        failf( "Can't open fasta file '%s'\n", ref_fn );
    }

    std::deque<sequence*> training_seqs;


    /* sample this many background sequence for each foreground sequence */
    const int bg_count = 1;
    
    int j;
    pos bg_pos;
    pos bg_offset;

    char*          seqname   = NULL;
    int            seqlen    = 0;
    int            curr_tid  = -1;
    char*          seq       = NULL;

    char* local_seq;
    local_seq = new char[ L + R + 2 ];
    local_seq[L+R+1] = '\0';


    for( i = 0; i < n; i++ ) {

        /* Load/switch sequences (chromosomes) as they are encountered in the
         * read stream. The idea here is to avoid thrashing by loading a large
         * sequence, but also avoid overloading memory by only loading one
         * chromosome at a time. */
        if( S[i]->pos.tid != curr_tid ) {
            seqname = reads_f->header->target_name[S[i]->pos.tid];
            seqlen  = reads_f->header->target_len[S[i]->pos.tid];
            if( seq ) free(seq); 

            log_printf( LOG_MSG, "reading sequence '%s' ... ", seqname );

            seq = faidx_fetch_seq( ref_f, seqname, 0, seqlen-1, &seqlen );

            if( seq == NULL ) {
                log_puts( LOG_WARN, "warning: reference sequence not found, skipping.\n" );
            }
            else {
                for( char* c = seq; *c; c++ ) *c = tolower(*c);
            }

            curr_tid = S[i]->pos.tid;

            log_puts( LOG_MSG, "done.\n" );
        }

        if( seq == NULL ) continue;
        
        /* add a foreground sequence */
        if( S[i]->pos.strand ) {
            if( S[i]->pos.pos < R ) continue;
            memcpy( local_seq, seq + S[i]->pos.pos - R, (L+1+R)*sizeof(char) );
            seqrc( local_seq, L+1+R );
        }
        else {
            if( S[i]->pos.pos < L ) continue;
            memcpy( local_seq, seq + (S[i]->pos.pos-L), (L+1+R)*sizeof(char) );
        }

        training_seqs.push_back( new sequence( local_seq, 1 ) );


        /* add a background sequence */
        /* adjust the current read position randomly, and sample */
        for( j = 0; j < bg_count; j++ ) {

            /* make things a bit more robust by sampling the background from a
             * random offset */
            bg_offset = (pos)rand_uniform( 100.0, 200.0 );
            if( rand_uniform( -1.0, 1.0 ) < 0.0 ) bg_offset = -bg_offset;

            bg_pos = S[i]->pos.pos + bg_offset;

            if( S[i]->pos.strand ) {
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


    size_t max_k = 5;
    size_t max_dep_dist = 5;
    M0 = new motif( L+1+R, max_k, 0 );
    M1 = new motif( L+1+R, max_k, 1 );

    if( train_backwards ) {
        train_motifs_backwards( *M0, *M1, &training_seqs, max_dep_dist, complexity_penalty );
    } else {
        train_motifs( *M0, *M1, &training_seqs, max_dep_dist, complexity_penalty );
    }


    std::deque<sequence*>::iterator i_seq;
    for( i_seq = training_seqs.begin(); i_seq != training_seqs.end(); i_seq++ ) {
        delete *i_seq;
    }



    log_unindent();

    free(S);
    free(seq);
    free(local_seq);
    samclose(reads_f);
    table_destroy(&T);

    t1 = clock();
    log_printf( LOG_MSG, "finished in %0.2f seconds\n",
                         (double)(t1-t0)/(double)CLOCKS_PER_SEC );

    log_unindent();
}


sequencing_bias::~sequencing_bias()
{
    clear();
}



void sequencing_bias::hash_reads( table* T, samfile_t* reads_f, size_t limit ) const
{
    log_puts( LOG_MSG, "hashing read positions..." );

    table_create(T);

    bam1_t* read = bam_init1();

    while( (limit == 0 || T->m < limit) && samread( reads_f, read ) > 0 ) {
        table_inc( T, read );
    }

    bam_destroy1(read);

    log_puts( LOG_MSG, "done.\n" );
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

