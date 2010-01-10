/*******************************************************************************
 * Functions below are implementation of multivariate density and random number
 * functions using the GNU Scientific Library (GSL) functions.
 * Most of the functions are available in R. This program is the reimplementation
 * under GSL standard, mainly designed for the convience of use in pure C code.
 * The program is also a supplement for the existing functions in current GSL
 * library.
 *
 * There is a simlar program written by Ralph dos Santos Silva. I greatly 
 * acknowledge the help from reading his program.  This program is
 * more close to the taste and format of GSL library and have more
 * functions. 
 * 
 * Depends: GSL >= 1.12
 *
 * Copyright (C) 2010 Chunhua Wu
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307, USA.
 *
 *******************************************************************************/

#ifndef __RANDIST_MV_H__
#define __RANDIST_MV_H__
#include <gsl/gsl_rng.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

double ran_invgamma(const gsl_rng *r, const double a, const double b);

int ran_mv_normal(const gsl_rng *r, const gsl_vector *mu,
		  const gsl_matrix *Sigma, gsl_vector *x);

int ran_mv_lognormal(const gsl_rng *r, const gsl_vector *mu,
		     const gsl_matrix *Sigma, gsl_vector *x);

int ran_mv_t(const gsl_rng *r, const gsl_vector *mu,
	     const gsl_matrix *Sigma, const double nu,  gsl_vector *x);

int ran_wishart(const gsl_rng *r, const double nu,
		const gsl_matrix *V, gsl_matrix *X);

int ran_invwishart(const gsl_rng *r, const double nu,
		   const gsl_matrix *V, gsl_matrix *X);

int ran_dirichlet(const gsl_rng *r, const gsl_vector *alpha,
		  gsl_vector *sample);

    
double ran_mv_normal_pdf(const gsl_vector *x, const gsl_vector *mu,
			 const gsl_matrix *Sigma);
  

double ran_mv_t_pdf(const gsl_vector *x, const gsl_vector *mu,
		    const gsl_matrix *Sigma, const double nu);

double ran_wishart_pdf(const gsl_matrix *X, const double nu, const gsl_matrix *V);

double ran_invwishart_pdf(const gsl_matrix *X, const double nu, const gsl_matrix *V);

double sf_mv_gamma(const double x, const int p);
double sf_mv_lngamma(const double x, const int p);

__END_DECLS

#endif /* __RANDIST_MV_H__ */





