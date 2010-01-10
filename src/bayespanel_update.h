/*******************************************************************************
 * Functions below are implementation of updating function for bayespanel
 * regression. These are building blocks for complicated panel models. 
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

#ifndef __BAYESPANEL_UPDATE_H__
#define __BAYESPANEL_UPDATE_H__
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

int bayespanel_update_beta(const gsl_rng *r, const gsl_vector *y, const gsl_matrix *X,
			   const gsl_matrix *W, const gsl_vector_int *ind,
			   const gsl_vector *beta0, const gsl_matrix *B0,
			   gsl_vector *beta, const gsl_vector *eta,
			   const gsl_matrix *D, const gsl_vector *lambda,
			   const gsl_vector *sigma2,
			   const int ar, const int ma, const gsl_vector *phi);

int bayespanel_update_b(const gsl_rng *r, const gsl_vector *y, const gsl_matrix *X,
			const gsl_matrix *W, const gsl_vector_int *ind,
			const gsl_vector *beta, gsl_matrix *b,
			const gsl_vector *eta,
			const gsl_matrix *D, const gsl_vector *lambda,
			const gsl_vector *sigma2,
			const int ar, const int ma, const gsl_vector *phi);

int bayespanel_update_D(const gsl_rng *r, const gsl_matrix *b, const gsl_vector *eta,
			const double rho0, const gsl_matrix *R0,
			gsl_matrix *D);

int bayespanel_update_sigma2(const gsl_rng *r,
			     const gsl_vector *y, const gsl_matrix *X,
			     const gsl_matrix *W, const gsl_vector_int *ind,
			     const gsl_vector *beta, const gsl_matrix *b,
			     const gsl_vector *lambda, const int hetero,
			     const double nu0,
			     const double delta0, gsl_vector *sigma2,
			     const int ar, const int ma, const gsl_vector *phi);

int bayespanel_update_lambda(const gsl_rng *r,
			     const gsl_vector *y, const gsl_matrix *X, 
			     const gsl_matrix *W, const gsl_vector_int *ind,
			     const gsl_vector *beta, const gsl_matrix *b,
			     const double nuG0, gsl_vector *lambda,
			     const gsl_vector *sigma2, 
			     const int ar, const int ma, const gsl_vector *phi);

int bayespanel_update_eta(const gsl_rng *r, const gsl_matrix *b,
			  const gsl_matrix *D, const double nuF0,
			  gsl_vector *eta);

int bayespanel_update_delta0(const gsl_rng *r, const gsl_vector_int *ind,
			     const double nu0,
			     const double nu00,
			     const double delta00, 
			     double *delta0, 
			     const gsl_vector *sigma2);

__END_DECLS

#endif /* __BAYESPANEL_UPDATE_H__ */





