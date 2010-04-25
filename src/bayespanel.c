#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "randist_mv.h"
#include "bayespanel_update.h"
#define max(a, b) ((a) > (b) ? (a) : (b))

/* codes for MCMC estimation for panel random effects model
   Allows for unbalanced panel,
   Yi = Xi*beta + Wi*bi + ei
   this is the implementation of alogrithm 2 of Chib and Carlin (1999)
   Variates of the functions are separated for interface purpose
   ,however they share the same building blocks.
 */



void bayespanel (int *control, int *link, int*hetero, int *arma, 
		 double *y, double *X, double *W, int *ind, int *dim,
		 double *beta0, double *B0, double *nuF0, double *rho0, double *R0,
		 double *nuG0, double *nu00, double *delta00,
		 double *nu0, double *delta0, double *phi0, double *P0,
		 double *betam, double *bm, double *Dm, double *sigma2m,
		 double *phim,  double *lnmarglik)
{
    const int n = dim[0];
    const int nn = dim[1];
    const int k = dim[2];
    const int q = dim[3];
    const int M0 = control[0]; // burnin
    const int M = control[1]; // number of iterations
    const int Mo = control[2]; // number of output iterations
    const int G = control[3]; // number of additional iterations used for marginal likelihood
    const unsigned long int seed = control[4];
    const int thin = control[5];
    const int verbose = control[6];
    const int ar = arma[0];
    const int ma = arma[1];
    const int narma = max(1, ar+ma);
    int iter, i, j;
    
    //memory access from input
    gsl_vector_view v_y = gsl_vector_view_array(y, nn);
    gsl_matrix_view m_X = gsl_matrix_view_array(X, nn, k);
    gsl_matrix_view m_W = gsl_matrix_view_array(W, nn, q);
    gsl_vector_int_view v_ind = gsl_vector_int_view_array(ind, n);
    gsl_vector_view v_beta0 = gsl_vector_view_array(beta0, k);
    gsl_matrix_view m_B0 = gsl_matrix_view_array(B0, k, k);
    gsl_matrix_view m_R0 = gsl_matrix_view_array(R0, q, q);
    gsl_vector_view v_phi0 = gsl_vector_view_array(phi0, narma);
    gsl_matrix_view m_P0 = gsl_matrix_view_array(P0, narma, narma);
    
    //memory copy to output
    gsl_matrix_view m_betam = gsl_matrix_view_array(betam, Mo, k);
    gsl_matrix_view m_bm = gsl_matrix_view_array(bm, Mo, n*q);
    gsl_matrix_view m_Dm = gsl_matrix_view_array(Dm, Mo, q*q);
 
    //intermediate storage
    gsl_vector *beta = gsl_vector_alloc(k);
    gsl_matrix *b = gsl_matrix_alloc(n, q);
    gsl_vector_view v_b = gsl_vector_view_array(b->data, n*q);
    gsl_matrix *D = gsl_matrix_alloc(q, q);
    gsl_vector_view v_D = gsl_vector_view_array(D->data, q*q);
    gsl_vector *sigma2 = gsl_vector_alloc(n);
    gsl_vector *eta = gsl_vector_alloc(n);
    gsl_vector *lambda = gsl_vector_alloc(n);
    gsl_vector *phi = gsl_vector_calloc(narma);

    double s2;
    
    //initialize the random variable generation environment
    gsl_rng *r = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(r, seed);

    //initialize values
    gsl_vector_memcpy(beta, &v_beta0.vector);
    ran_wishart(r, *rho0, &m_R0.matrix, D);
    gsl_linalg_cholesky_decomp(D);
    gsl_linalg_cholesky_invert(D);
    s2 = 1/gsl_ran_gamma(r, *nu0/2.0, 2.0/ *delta0);
    gsl_vector_set_all(sigma2, s2);
    gsl_vector_set_all(eta, 1);
    gsl_vector_set_all(lambda, 1);
    gsl_matrix_set_all(b, 0); 

    gsl_vector_set_all(sigma2, 0.00857);

    //begin the iterations
    for (iter=0; iter<M+M0; iter++)
    {
        //Check user interruption
        R_CheckUserInterrupt();
	//beta update
	bayespanel_update_beta(r, &v_y.vector, &m_X.matrix, &m_W.matrix,
			       &v_ind.vector, &v_beta0.vector, &m_B0.matrix,
			       beta, eta, D, lambda, sigma2, ar, ma, phi);
	//b update
	bayespanel_update_b(r, &v_y.vector, &m_X.matrix, &m_W.matrix,
			    &v_ind.vector, beta, b, eta, D, lambda, sigma2,
			    ar, ma, phi);

	if (*link == 1)
	{
	    //lambda update
	  bayespanel_update_lambda(r, &v_y.vector, &m_X.matrix, &m_W.matrix,
				     &v_ind.vector, beta, b, *nuG0, lambda, sigma2,
				     ar, ma, phi);
	    //eta update
	  bayespanel_update_eta(r, b, D, *nuF0, eta);
	}
	//D update
	bayespanel_update_D(r, b, eta, *rho0, &m_R0.matrix, D);
	//sigma2 update
	bayespanel_update_sigma2(r, &v_y.vector, &m_X.matrix, &m_W.matrix,
				 &v_ind.vector, beta, b, lambda, *hetero, *nu0,
				 *delta0, sigma2, ar, ma, phi);
        if (*hetero == 0)
	  {
	    s2 = gsl_vector_get(sigma2, 0);
	  }
	else
	  {
	  //delta0 update
	  bayespanel_update_delta0(r, &v_ind.vector, *nu0, *nu00, *delta00,
				   delta0, sigma2);
	  }
	    

	
        //print out given verbose option
	if(verbose > 0 && (iter - M0) % verbose == 0)
	  {
	    Rprintf("\nBayesPanel iteration %d of %d \n", iter, M+M0);
	    Rprintf("beta = \n");
	    for (i=0; i<k; i++)
	      Rprintf("%10.5f\n",gsl_vector_get(beta, i));
	    Rprintf("D = \n");
	    for (i=0; i<q; i++)
	      {
		for (j=0; j<q; j++)
		  Rprintf("%10.5f ", gsl_matrix_get(D, i, j));
		Rprintf("\n");
	      }
	    if (*hetero == 0)
	      {
		Rprintf("Sigma2 = \n");
		Rprintf("%10.5f\n", s2);
	      }
	  }

	//output matrix
	if(iter >= M0 && ((iter-M0) % thin) == 0 )
	{
	  gsl_matrix_set_row(&m_betam.matrix, (iter-M0) / thin, beta);
	  gsl_matrix_set_row(&m_bm.matrix, (iter-M0) / thin, &v_b.vector);
	  gsl_matrix_set_row(&m_Dm.matrix, (iter-M0) / thin, &v_D.vector);
	  if (*hetero == 0)
	    {
	      gsl_vector_view v_sigma2m = gsl_vector_view_array(sigma2m, Mo);
	      gsl_vector_set(&v_sigma2m.vector, (iter-M0) / thin, s2);
	    }
	  else
	    {
	      gsl_matrix_view m_sigma2m = gsl_matrix_view_array(sigma2m, Mo, n);
	      gsl_matrix_set_row (&m_sigma2m.matrix, (iter-M0) / thin, sigma2);
	    }
	  if (ar+ma>0)
	    {
	      gsl_matrix_view m_phim = gsl_matrix_view_array(phim, Mo, narma);
	      gsl_matrix_set_row(&m_phim.matrix, (iter-M0) / thin, phi);
	    }
	}

    }
    //calculate ln marginal likelihood
    
    
    //free memory
    gsl_vector_free(beta);
    gsl_matrix_free(b);
    gsl_matrix_free(D);
    gsl_vector_free(sigma2);
    gsl_vector_free(eta);
    gsl_vector_free(lambda);
    gsl_vector_free(phi);
    gsl_rng_free(r);
}













 
