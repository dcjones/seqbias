
#ifndef PEAKOLATOR_LOGGER
#define PEAKOLATOR_LOGGER

#ifdef __cplusplus
extern "C" {
#endif


typedef enum {
    LOG_ERROR=0,
    LOG_WARN=1,
    LOG_MSG=2,
    LOG_BLAB=3
} log_vl;


void log_puts( int vl, const char* msg );
void log_printf( int vl, const char* fmt, ... );
void log_indent(void);
void log_unindent(void);
void log_verbosity(int);

void fail( const char* msg );
void failf( const char* fmt, ... );


#ifdef __cplusplus
}
#endif


#endif


