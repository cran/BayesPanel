/********************************************************************************
 * Modules for BayesPanel updating
 ********************************************************************************/

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "randist_mv.h"
#include "bayespanel_update.h"

int bayespanel_update_beta(const gsl_rng *r, const gsl_vector *y, const gsl_matrix *X,
			   const gsl_matrix *W, const gsl_vector_int *ind,
			   const gsl_vector *beta0, const gsl_matrix *B0,
			   gsl_vector *beta, const gsl_vector *eta,
			   const gsl_matrix *D, const gsl_vector *lambda,
			   const gsl_vector *sigma2,
			   const int ar, const int ma, const gsl_vector *phi)
{
  const int n = ind->size;
  const int nn = X->size1;
  const int k = X->size2;
  const int q = W->size2;

  int i, ni, row = 0;
  double etai, lambdai, sigma2i;

  gsl_matrix *XVX = gsl_matrix_calloc(k, k);
  gsl_vector *XVy = gsl_vector_calloc(k);
  gsl_vector *betahat = gsl_vector_calloc(k);
  gsl_vector *work_k = gsl_vector_calloc(k);
  gsl_matrix *B = gsl_matrix_alloc(k, k);
  gsl_matrix *B0inv = gsl_matrix_alloc(k, k);
  gsl_matrix_memcpy(B0inv, B0);
  gsl_linalg_cholesky_decomp(B0inv);
  gsl_linalg_cholesky_invert(B0inv);
 
  for (i=0; i<n; i++)
    {
      ni = gsl_vector_int_get(ind, i);
      etai = gsl_vector_get(eta, i);
      lambdai = gsl_vector_get(lambda, i);
      sigma2i = gsl_vector_get(sigma2, i);
      gsl_vector_const_view v_yi = gsl_vector_const_subvector(y, row, ni);
      gsl_matrix_const_view m_Xi = gsl_matrix_const_submatrix(X, row, 0, ni, k);
      gsl_matrix_const_view m_Wi = gsl_matrix_const_submatrix(W, row, 0, ni, q);
      gsl_matrix *Vi = gsl_matrix_alloc(ni, ni);
      gsl_matrix *Omegai = gsl_matrix_alloc(ni, ni);
      //construct Omegai based on ar, ma, and corresponding parameters
      if (ar == 0 && ma == 0)
	{
	  gsl_matrix_set_identity(Omegai);
	}
      gsl_matrix *work_ni_q = gsl_matrix_alloc(ni, q);
      gsl_matrix *work_k_ni = gsl_matrix_alloc(k, ni);

      gsl_matrix_memcpy(Vi, Omegai);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, &m_Wi.matrix,
		     D, 0.0, work_ni_q);
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0/etai, work_ni_q,
		     &m_Wi.matrix, sigma2i/lambdai, Vi); //completes Vi
      gsl_linalg_cholesky_decomp(Vi);
      gsl_linalg_cholesky_invert(Vi);
 
      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &m_Xi.matrix, Vi,
		     0.0, work_k_ni); //internal working matrix
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, work_k_ni,
		     &m_Xi.matrix, 1.0, XVX);//completes XVX
      gsl_blas_dgemv(CblasNoTrans, 1.0, work_k_ni, &v_yi.vector,
		     1.0, XVy); //completes XVy
      
      gsl_matrix_free(Vi);
      gsl_matrix_free(Omegai);
      gsl_matrix_free(work_ni_q);
      gsl_matrix_free(work_k_ni);
      row += ni;
    }
  
  gsl_matrix_memcpy(B, B0inv);
  gsl_matrix_add(B, XVX);
  gsl_linalg_cholesky_decomp(B);
  gsl_linalg_cholesky_invert(B);//this compeltes B

  gsl_blas_dgemv(CblasNoTrans, 1.0, B0inv, beta0, 0.0, work_k);
  gsl_vector_add(work_k, XVy);
  gsl_blas_dgemv(CblasNoTrans, 1.0, B, work_k, 0.0, betahat);
  //this completes betahat

  ran_mv_normal(r, betahat, B, beta);//update beta
  gsl_matrix_free(XVX);
  gsl_vector_free(XVy);
  gsl_vector_free(betahat);
  gsl_vector_free(work_k);
  gsl_matrix_free(B);
  gsl_matrix_free(B0inv);

  return 0;
}  


int bayespanel_update_b(const gsl_rng *r, const gsl_vector *y, const gsl_matrix *X,
			const gsl_matrix *W, const gsl_vector_int *ind,
			const gsl_vector *beta, gsl_matrix *b,
			const gsl_vector *eta,
			const gsl_matrix *D, const gsl_vector *lambda,
			const gsl_vector *sigma2,
			const int ar, const int ma, const gsl_vector *phi)
{
  const int n = ind->size;
  const int nn = X->size1;
  const int k = X->size2;
  const int q = W->size2;

  int i, ni, row = 0;
  double etai, lambdai, sigma2i;
  
  gsl_vector *eb = gsl_vector_alloc(nn);
  gsl_vector *bihat = gsl_vector_alloc(q);
  gsl_matrix *Dinv = gsl_matrix_alloc(q, q);
  gsl_matrix *Di = gsl_matrix_alloc(q, q);

  gsl_matrix_memcpy(Dinv, D);
  gsl_linalg_cholesky_decomp(Dinv);
  gsl_linalg_cholesky_invert(Dinv);
  gsl_vector_memcpy(eb, y);
  gsl_blas_dgemv(CblasNoTrans, -1.0, X, beta, 1.0, eb);
  
  for (i=0; i<n; i++)
    {
      ni = gsl_vector_int_get(ind, i);
      etai = gsl_vector_get(eta, i);
      lambdai = gsl_vector_get(lambda, i);
      sigma2i = gsl_vector_get(sigma2, i);
      gsl_vector_view v_ebi = gsl_vector_subvector(eb, row, ni);
      gsl_matrix_const_view m_Wi = gsl_matrix_const_submatrix(W, row, 0, ni, q);
      gsl_vector_view v_bi = gsl_matrix_row(b, i);
      gsl_matrix *Omegai = gsl_matrix_alloc(ni, ni);

      gsl_matrix *work_q_ni = gsl_matrix_alloc(q, ni);
      gsl_vector *work_q = gsl_vector_alloc(q);

      if (ar == 0 && ma == 0)
	{
	  gsl_matrix_set_identity(Omegai);
	}
      gsl_linalg_cholesky_decomp(Omegai);
      gsl_linalg_cholesky_invert(Omegai);
      gsl_matrix_memcpy(Di, Dinv);

      gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0,
		     &m_Wi.matrix, Omegai, 0.0, work_q_ni);
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, lambdai/sigma2i,
		     work_q_ni, &m_Wi.matrix, etai, Di);
      gsl_linalg_cholesky_decomp(Di);
      gsl_linalg_cholesky_invert(Di);

      gsl_blas_dgemv(CblasNoTrans, 1.0, work_q_ni, &v_ebi.vector, 0.0, work_q);
      gsl_blas_dgemv(CblasNoTrans, lambdai/sigma2i, Di, work_q, 0.0, bihat);

      ran_mv_normal(r, bihat, Di, &v_bi.vector);//draw bi
      
      gsl_matrix_free(Omegai);
      gsl_matrix_free(work_q_ni);
      gsl_vector_free(work_q);
      
      row += ni;
    }

  gsl_vector_free(eb);
  gsl_vector_free(bihat);
  gsl_matrix_free(Dinv);
  gsl_matrix_free(Di);

  return 0;
}
  

      
      
int bayespanel_update_D(const gsl_rng *r, const gsl_matrix *b, const gsl_vector *eta,
			const double rho0, const gsl_matrix *R0,
			gsl_matrix *D)
{
  const int n = b->size1;
  const int q = b->size2;

  double rho = rho0 + n;
  double etai;
  gsl_matrix *R = gsl_matrix_alloc(q, q);
  gsl_matrix_memcpy(R, R0);
  gsl_linalg_cholesky_decomp(R);
  gsl_linalg_cholesky_invert(R);

  int i;
  for (i=0; i<n; i++)
    {
      etai = gsl_vector_get(eta, i);
      gsl_vector_const_view v_bi = gsl_matrix_const_row(b, i);
      gsl_blas_dger(etai, &v_bi.vector, &v_bi.vector, R);//rank 1 update
    }
      
  gsl_linalg_cholesky_decomp(R);
  gsl_linalg_cholesky_invert(R);
  
  ran_wishart(r, rho, R, D);
  gsl_linalg_cholesky_decomp(D);
  gsl_linalg_cholesky_invert(D);

  gsl_matrix_free(R);

  return 0;
}


int bayespanel_update_sigma2(const gsl_rng *r,
			     const gsl_vector *y, const gsl_matrix *X,
			     const gsl_matrix *W, const gsl_vector_int *ind,
			     const gsl_vector *beta, const gsl_matrix *b,
			     const gsl_vector *lambda, const int hetero,
			     const double nu0,
			     const double delta0, gsl_vector *sigma2,
			     const int ar, const int ma, const gsl_vector *phi)
{
  const int n = ind->size;
  const int nn = X->size1;
  const int k = X->size2;
  const int q = W->size2;

  int i, ni, row = 0;
  double lambdai, sigma2i;

  double *delta = malloc(sizeof(double));
  double deltasum = 0;
  double s2 = 0;
  
  gsl_vector *e = gsl_vector_alloc(nn);
  gsl_vector_memcpy(e, y);
  gsl_blas_dgemv(CblasNoTrans, -1.0, X, beta, 1.0, e);

  for (i=0; i<n; i++)
    {
      ni = gsl_vector_int_get(ind, i);
      lambdai = gsl_vector_get(lambda, i);
      gsl_vector_view v_ei = gsl_vector_subvector(e, row, ni);
      gsl_vector_const_view v_bi = gsl_matrix_const_row(b, i);
      gsl_matrix_const_view m_Wi = gsl_matrix_const_submatrix(W, row, 0, ni, q);
      gsl_matrix *Omegai = gsl_matrix_alloc(ni, ni);

      gsl_blas_dgemv(CblasNoTrans, -1.0, &m_Wi.matrix, &v_bi.vector,
		     1.0, &v_ei.vector);

      if(ar == 0 && ma == 0)
	{
	  gsl_matrix_set_identity(Omegai);
	}

      gsl_linalg_cholesky_decomp(Omegai);
      gsl_linalg_cholesky_invert(Omegai);
      
      gsl_vector *work_ni = gsl_vector_alloc(ni);
      gsl_blas_dgemv(CblasNoTrans, lambdai, Omegai, &v_ei.vector, 0.0, work_ni);
      gsl_blas_ddot(&v_ei.vector, work_ni, delta);
      
      deltasum += *delta;
      
      sigma2i = 1 / gsl_ran_gamma(r, (nu0+ni)/2.0, 2.0/(delta0+ *delta));
      gsl_vector_set(sigma2, i, sigma2i);

      gsl_matrix_free(Omegai);
      gsl_vector_free(work_ni);
      row += ni;
    }

  if (hetero == 0)
    {
      s2 = 1 / gsl_ran_gamma(r, (nu0+nn)/2.0, 2.0/(delta0+deltasum));
      gsl_vector_set_all(sigma2, s2);
    }

  gsl_vector_free(e);
  free(delta);

  return 0;
}


int bayespanel_update_lambda(const gsl_rng *r,
			     const gsl_vector *y, const gsl_matrix *X, 
			     const gsl_matrix *W, const gsl_vector_int *ind,
			     const gsl_vector *beta, const gsl_matrix *b,
			     const double nuG0, gsl_vector *lambda,
			     const gsl_vector *sigma2, 
			     const int ar, const int ma, const gsl_vector *phi)
{
  const int n = ind->size;
  const int nn = X->size1;
  const int k = X->size2;
  const int q = W->size2;

  int i, ni, row = 0;
  double lambdai, sigma2i;
  double *temp = malloc(sizeof(double));

  gsl_vector *e = gsl_vector_alloc(nn);
  gsl_vector_memcpy(e, y);
  gsl_blas_dgemv(CblasNoTrans, -1.0, X, beta, 1.0, e);

  for (i=0; i<n; i++)
    {
      ni = gsl_vector_int_get(ind, i);
      sigma2i = gsl_vector_get(sigma2, i);

      gsl_vector_view v_ei = gsl_vector_subvector(e, row, ni);
      gsl_vector_const_view v_bi = gsl_matrix_const_row(b, i);
      gsl_matrix_const_view m_Wi = gsl_matrix_const_submatrix(W, row, 0, ni, q);
      gsl_matrix *Omegai = gsl_matrix_alloc(ni, ni);
      gsl_vector *work_ni = gsl_vector_alloc(ni);

      gsl_blas_dgemv(CblasNoTrans, -1.0, &m_Wi.matrix, &v_bi.vector,
		     1.0, &v_ei.vector);//completes ei

      if (ar == 0 && ma == 0)
	{
	  gsl_matrix_set_identity(Omegai);
	}
      gsl_linalg_cholesky_decomp(Omegai);
      gsl_linalg_cholesky_invert(Omegai);

      gsl_blas_dgemv(CblasNoTrans, 1/sigma2i, Omegai, &v_ei.vector, 0.0, work_ni);
      gsl_blas_ddot(&v_ei.vector, work_ni, temp);

      lambdai = gsl_ran_gamma(r, (nuG0+ni)/2.0, (nuG0+ *temp)/2.0);
      gsl_vector_set(lambda, i, lambdai);

      gsl_matrix_free(Omegai);
      gsl_vector_free(work_ni);
      row += ni;
    }

  free(temp);
  gsl_vector_free(e);

  return 0;
}



int bayespanel_update_eta(const gsl_rng *r, const gsl_matrix *b,
			  const gsl_matrix *D, const double nuF0,
			  gsl_vector *eta)
{
  const int n = b->size1;
  const int q = b->size2;

  int i;
  double etai;
  double *temp = malloc(sizeof(double));

  gsl_matrix *Dinv = gsl_matrix_alloc(q, q);
  gsl_matrix_memcpy(Dinv, D);
  gsl_linalg_cholesky_decomp(Dinv);
  gsl_linalg_cholesky_invert(Dinv);
  gsl_vector *work_q = gsl_vector_alloc(q);
  
  for (i=0; i<n; i++)
    {
      gsl_vector_const_view v_bi = gsl_matrix_const_row(b, i);
      gsl_blas_dgemv(CblasNoTrans, 1.0, Dinv, &v_bi.vector, 0.0, work_q);
      gsl_blas_ddot(&v_bi.vector, work_q, temp);

      etai = gsl_ran_gamma(r, (nuF0+q)/2.0, (nuF0+ *temp)/2.0);
    }

  free(temp);
  gsl_matrix_free(Dinv);
  gsl_vector_free(work_q);

  return 0;
}


int bayespanel_update_delta0(const gsl_rng *r, const gsl_vector_int *ind,
			     const double nu0,
			     const double nu00,
			     const double delta00, 
			     double *delta0, 
			     const gsl_vector *sigma2)
{
  const int n = ind->size;

  int i;
  double temp = 0.0;
  for (i=0; i<n; i++)
    {
      temp += 1/gsl_vector_get(sigma2, i);
    }
  *delta0 = gsl_ran_gamma(r, (n*nu0+nu00)/2.0, (temp+delta00)/2.0);

  return 0;
}
	
  

  
  


// Other functions to implement: bayespanel_update_phi for time series



      
      
      
			  
			  


      
      
      

      

      
      

  

  
	

      

      
      

  
  
  
      
      
      

      
      
      
      

  
