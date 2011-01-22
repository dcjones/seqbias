
#ifndef PEAKOLATOR_MISCMATH
#define PEAKOLATOR_MISCMATH

#include <cstdlib>

/* Miscelaneous math functions. */

/* Compute log( exp(a) + exp(b) ) avoiding underflow. */
double logaddexp( double a, double b );



#ifdef HAVE_LIBGSL
/* Compute log( exp(a) - exp(b) ) avoiding underflow. */
double logsubexp( double a, double b );


/* Compute: log( I_x( a, b ) ), where I is the regularized incomplete beta
 * function.  Equivalent to log( gsl_sf_beta_inc( a, b, x ) ), but computed in
 * such a way as to avoid underflow.
 */
double log_beta_inc( double a, double b, double x );


/* The negative binomial log comulative distribution function. */
double lpnbinom( unsigned int q, double r, double p, bool lower_tail = false );
#endif


/* Generalived extreme value cumulative distribution function. */
double pgev( double q, double loc, double scale, double shape, bool lower_tail = false );


/* Generalized extreme value log-density function. */
double ldgev( double x, double loc, double scale, double shape );



/* (very) simple vector/matrix arithmatic */
void colcpy( double* dest, const double* src, size_t j, size_t n, size_t m );
void vecadd( double* u, double* v, size_t n );
void vecsub( double* u, double* v, size_t n );
void vecaddcol( double* u, double* V, size_t n, size_t m, size_t j );
void vecsubcol( double* u, double* V, size_t n, size_t m, size_t j );
void matsetcol( double* U, double* v, size_t n, size_t m, size_t j );

#endif


