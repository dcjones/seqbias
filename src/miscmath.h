/*
 * Miscelaneous and obscure math functions, primarily a supplement to GSL.
 *
 * Copyright (C) 2011 by Daniel C. Jones <dcjones@cs.washington.edu>
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */


#ifndef _MISCMATH_H_
#define _MISCMATH_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "config.h"
#include <stdlib.h>
#include <stdbool.h>

/* is a a valid probability */
bool is_prob(double a);

/* is a a valid log-probability */
bool is_log_prob(double a);

/* Compute log( exp(a) + exp(b) ) avoiding underflow. */
double logaddexp(double a, double b);


/* Compute log( exp(a) - exp(b) ) avoiding underflow. */
double logsubexp(double a, double b);


#ifdef HAVE_LIBGSL

/* geometric log density function (with x >= 0)  */
double ldgeom(unsigned int x, double p);



/* Compute: log( I_x( a, b ) ), where I is the regularized incomplete beta
 * function.  Equivalent to log( gsl_sf_beta_inc( a, b, x ) ), but computed in
 * such a way as to avoid underflow.
 */
double log_beta_inc(double a, double b, double x);


/* The poisson log density function */
double ldpois(unsigned int x, double mu);


/* The binomial log density function */
double ldbinom(unsigned int x, double p, unsigned int n);

/* The binomial log distribution function */
double lpbinom(unsigned int x, double p, unsigned int n, bool lower_tail);


/* The negative binomial log density function */
double ldnbinom(unsigned int x, double r, double p);

/* The negative binomial log comulative distribution function. */
double lpnbinom(unsigned int q, double r, double p, bool lower_tail);



/* The truncated negative binomial log density function */
double ldtnbinom(unsigned int x, double r, double p, unsigned int k);

/* The truncated negative binomial log distribution function. */
double lptnbinom(unsigned int q, double r, double p,
                 unsigned int k, bool lower_tail);



/* binomial coefficient */
double binco(double n, double k);

/* log binomial coefficient */
double lbinco(double n, double k);


/* The log density function of the decapitated negative binomial distribution.
 */
double lddnbinom(unsigned x, double r, double p);


/* The density function of the the distribution over the sum of d i.i.d.
 * decapitated negative binomial variables. */
double ddnbsum(unsigned int x, double r, double p, unsigned int d);

/* The distribution function of the the distribution over the sum of d i.i.d.
 * decapitated negative binomial variables. */
double pdnbsum(unsigned int x, double r, double p,
               unsigned int d, bool lower_tail);


/* The log density function of the the distribution over the sum of d i.i.d.
 * decapitated negative binomial variables. */
double lddnbsum(unsigned int x, double r, double p, unsigned int d);

/* The log distribution function of the the distribution over the sum of d
 * i.i.d.  decapitated negative binomial variables. */
double lpdnbsum(unsigned int x, double r, double p,
                unsigned int d, bool lower_tail);


/* log density of the zero inflated negative binomial distribution Where f( x;
 * r, p ) is the decapitated negative binomial density function, then density
 * function for the zinb distribution is,
 *
 *  g( x; r, p, a ) = a                           if x = 0
 *                    (1.0 - a) * f( x; r, p )    otherwise
 */
double ldzinb(unsigned int x, double r, double p, double a);


/* Log densidy of the distribution over the sum over d i.i.d zero inflated
 * negative binomial distributed variables.
 */
double ldzinbsum(unsigned int x, double r, double p, double a, unsigned int d);


/* Mean/expected value of the sum of d i.i.d. zero inflated negative binomial
 * variables.
 */
double zinb_mean(double r, double p, double a, unsigned int d);


/* The density function for the classical vacancy problem. If k balls are thrown
 * into n bins, landing in each with equal probability, what is the probability
 * that x of the bins are empty?
 */
double dcvp(unsigned int x, unsigned int k, unsigned int n);
double ldcvp(unsigned int x, unsigned int k, unsigned int n);


#endif


/* Generalived extreme value cumulative distribution function. */
double pgev(double q, double loc, double scale, double shape,
            bool lower_tail);


/* Generalized extreme value log-density function. */
double ldgev(double x, double loc, double scale, double shape);


#ifdef __cplusplus
}
#endif

#endif



