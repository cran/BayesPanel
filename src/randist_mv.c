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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "randist_mv.h"

//Random number generator for multivariate normal distribution

int ran_mv_normal(const gsl_rng *r, const gsl_vector *mu,
		     const gsl_matrix *Sigma, gsl_vector *x)
{
  const int k = mu->size;
  int i;
  gsl_matrix *A = gsl_matrix_alloc(k, k);

  gsl_matrix_memcpy(A, Sigma);
  gsl_linalg_cholesky_decomp(A);

  for (i=0; i<k; i++)
    {
      gsl_vector_set(x, i, gsl_ran_gaussian(r, 1));
    }
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, A, x);
  gsl_vector_add(x, mu);

  gsl_matrix_free(A);

  return 0;
}

//Random number generator for multivariate log-normal distribution
int ran_mv_lognormal(const gsl_rng *r, const gsl_vector *mu,
			const gsl_matrix *Sigma, gsl_vector *x)
{
    const int k = mu->size;
    int i;
    double temp;
    
    ran_mv_normal(r, mu, Sigma, x);
    for(i=0; i<k; i++)
    {
	temp = gsl_vector_get(x, i);
	gsl_vector_set(x, i, exp(temp));
    }

    return 0;
}
    
// random number generator for multivariate t distribution
// nu is degree of freedom, mu is mean value vector.

int ran_mv_t(const gsl_rng *r, const gsl_vector *mu,
		const gsl_matrix *Sigma, const double nu,  gsl_vector *x)
{
  const int k = mu->size;
  gsl_matrix *A = gsl_matrix_alloc(k, k);
  double v;

  v = gsl_ran_chisq(r, nu);
  gsl_matrix_memcpy(A, Sigma);
  gsl_linalg_cholesky_decomp(A);

  ran_mv_normal(r, mu, Sigma, x);
  
  gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, A, x);
  gsl_vector_scale(x, 1/sqrt(v));
  gsl_vector_add(x, mu);

  gsl_matrix_free(A);
  return 0;
  
}


//Wishart distribution random number generator
int ran_wishart(const gsl_rng *r, const double nu,
		   const gsl_matrix *V,	 gsl_matrix *X)
{
  const int k = V->size1;
  int i, j;
  gsl_matrix *A = gsl_matrix_calloc(k, k);
  gsl_matrix *L = gsl_matrix_alloc(k, k);

  for(i=0; i<k; i++)
    {
      gsl_matrix_set(A, i, i, sqrt(gsl_ran_chisq(r, (nu-i))));
      for (j=0; j<i; j++){
       gsl_matrix_set(A, i, j, gsl_ran_gaussian(r, 1));
      }
    }
  gsl_matrix_memcpy(L, V);
  gsl_linalg_cholesky_decomp(L);
  gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, 1.0, L, A);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, A, A, 0.0, X);
  
  gsl_matrix_free(A);
  gsl_matrix_free(L);

  return 0;
}

// Inverse Wishart distribution random number generator
int ran_invwishart(const gsl_rng *r, const double nu,
		      const gsl_matrix *V, gsl_matrix *X)
{
  const int k = V->size1;
  gsl_matrix *Vinv = gsl_matrix_alloc(k, k);
  
  gsl_matrix_memcpy(Vinv, V);
  gsl_linalg_cholesky_decomp(Vinv);
  gsl_linalg_cholesky_invert(Vinv);

  ran_wishart(r, nu, Vinv, X);
  gsl_linalg_cholesky_decomp(X);
  gsl_linalg_cholesky_invert(X);

  gsl_matrix_free(Vinv);
  return 0;
}


int ran_dirichlet(const gsl_rng *r, const gsl_vector *alpha,
		  gsl_vector *sample)
{
    const int k = alpha->size;
    double sum =0.0;
    double alphai =0.0;
    double xi;
    int i;
    
    for(i=0; i<k; i++)
    {
	alphai = gsl_vector_get(alpha, i);
	xi = gsl_ran_gamma(r, alphai, 1.0);
	gsl_vector_set(sample, i, xi);
	sum += xi;
    }
    gsl_vector_scale(sample, 1/sum);

    return 0;
}


    
/* Probability density functions. 
 */


//Multivariate normal distribution

double ran_mv_normal_pdf(const gsl_vector *x, const gsl_vector *mu,
			 const gsl_matrix *Sigma)
{
  const int k = x->size;
  int s;
  double det, den;

  gsl_vector *y = gsl_vector_alloc(k);
  gsl_vector *work_k = gsl_vector_alloc(k);
  
  gsl_matrix *work_k_k = gsl_matrix_alloc(k, k);
  gsl_matrix *Sigmainv = gsl_matrix_alloc(k, k);
  gsl_permutation *p = gsl_permutation_alloc(k);
  
  gsl_vector_memcpy(y, x);
  gsl_vector_sub(y, mu);
  
  gsl_matrix_memcpy(work_k_k, Sigma);
  gsl_linalg_LU_decomp(work_k_k, p, &s);
  gsl_linalg_LU_invert(work_k_k, p, Sigmainv);
  det = gsl_linalg_LU_det(work_k_k, s);

  gsl_blas_dgemv(CblasNoTrans, 1.0, Sigmainv, y, 0.0, work_k);
  gsl_blas_ddot(y, work_k, &den);
  den = exp(-0.5*den) / sqrt(pow((2*M_PI), k)*det);

  gsl_vector_free(y);
  gsl_vector_free(work_k);
  gsl_matrix_free(work_k_k);
  gsl_matrix_free(Sigmainv);
  gsl_permutation_free(p);
  
  return den;
}
  

double ran_mv_t_pdf(const gsl_vector *x, const gsl_vector *mu,
		    const gsl_matrix *Sigma, const double nu)
{
  const int k = x->size;
  int s;
  double det,temp, den;

  gsl_vector *y = gsl_vector_alloc(k);
  gsl_vector *work_k = gsl_vector_alloc(k);

  gsl_matrix *work_k_k = gsl_matrix_alloc(k, k);
  gsl_matrix *Sigmainv = gsl_matrix_alloc(k, k);
  gsl_permutation *p = gsl_permutation_alloc(k);

  gsl_vector_memcpy(y, x);
  gsl_vector_sub(y, mu);

  gsl_matrix_memcpy(work_k_k, Sigma);
  gsl_linalg_LU_decomp(work_k_k, p, &s);
  gsl_linalg_LU_invert(work_k_k, p, Sigmainv);
  det = gsl_linalg_LU_det(work_k_k, s);

  gsl_blas_dgemv(CblasNoTrans, 1.0/k, Sigmainv, y, 0.0, work_k);
  gsl_blas_ddot(y, work_k, &temp);
  temp = pow((1+temp), (nu+ (double) k)/2 );
  temp *= gsl_sf_gamma(nu/2) * pow(nu, k/2) * pow(M_PI, k/2) * sqrt(det);

  den = gsl_sf_gamma((nu+ (double) k)/2);
  den /= temp;

  gsl_vector_free(y);
  gsl_vector_free(work_k);
  gsl_matrix_free(work_k_k);
  gsl_matrix_free(Sigmainv);
  gsl_permutation_free(p);

  return den;
}

  

double ran_wishart_pdf(const gsl_matrix *X, const double nu, const gsl_matrix *V)
{
  const int k = X->size1;
  double detX, detV, den, temp;
  int s, i;

  gsl_matrix *work_k_k = gsl_matrix_alloc(k, k);
  gsl_matrix *Vinv = gsl_matrix_alloc(k, k);
  gsl_permutation *p = gsl_permutation_alloc(k);

  gsl_matrix_memcpy(work_k_k, X);
  gsl_linalg_LU_decomp(work_k_k, p, &s);
  detX = gsl_linalg_LU_det(work_k_k, s);

  gsl_matrix_memcpy(work_k_k, V);
  gsl_linalg_LU_decomp(work_k_k, p, &s);
  detV = gsl_linalg_LU_det(work_k_k, s);
  gsl_linalg_LU_invert(work_k_k, p, Vinv);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, Vinv, X, 0.0, work_k_k);
  den = 0;
  for (i=0; i<k; i++)
    {
      den -= 0.5 * gsl_matrix_get(work_k_k, i, i);
    }
  //den = exp(den);

  temp = 0.5*(nu-k-1)*log(detX) - 0.5*nu*k*log(2) -0.5*nu*log(detV);
  temp -= sf_mv_lngamma(nu/2, k);
  den += temp;

  den = exp(den);
  
  gsl_matrix_free(work_k_k);
  gsl_matrix_free(Vinv);
  gsl_permutation_free(p);

  return den;
}

  
double ran_invwishart_pdf(const gsl_matrix *X, const double nu, const gsl_matrix *V)
{
  const int k = X->size1;
  double den;
  
  gsl_matrix *Xinv = gsl_matrix_alloc(k, k);
  gsl_matrix *Vinv = gsl_matrix_alloc(k, k);

  gsl_matrix_memcpy(Xinv, X);
  gsl_matrix_memcpy(Vinv, V);
  gsl_linalg_cholesky_decomp(Xinv);
  gsl_linalg_cholesky_invert(Xinv);
  gsl_linalg_cholesky_decomp(Vinv);
  gsl_linalg_cholesky_invert(Vinv);

  den = ran_wishart_pdf(Xinv, nu, Vinv);

  gsl_matrix_free(Xinv);
  gsl_matrix_free(Vinv);

  return den;
}

  

double sf_mv_gamma(const double x, const int p)
{
  double v;
  int i;

  v = pow(M_PI, p*(p-1)/4);

  for (i=0; i<p; i++)
    {
      v *= gsl_sf_gamma(x+(1-i)/2);
    }

  return v;
}


double sf_mv_lngamma(const double x, const int p)
{
  double v;
  int i;

  v = p*(p-1)*log(M_PI)/4;

  for (i=0; i<p; i++)
    {
      v += gsl_sf_lngamma(x+(1-i)/2);
    }

  return v;
}






