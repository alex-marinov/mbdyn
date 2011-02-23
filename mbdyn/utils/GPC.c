/*
 * Copyright (C) 2008
 *
 * Mattia Mattaboni	<mattaboni@aero.polimi.it>
 */

#include <stdlib.h>
#include "GPC.h"
#include "udu.h"
#include "time.h"
#include "sys/time.h"

//#define CTRB
#define ARX_KALMAN

GPC_ResT ARX_Initialize( ARX_Model *ARX, unsigned na, unsigned nb, 
		    unsigned m, unsigned p, double mu,
		    double delta, unsigned FlagSimplyProper){

	unsigned dimension;

	ARX->na = na;
	ARX->nb = nb;
	ARX->m = m;
	ARX->p = p;
	ARX->mu = mu;
	ARX->SimplyProper = FlagSimplyProper;

	dimension = p*na + m*nb + m*FlagSimplyProper;

	matrix_init( &ARX->P, dimension, dimension );
	matrix_eye( &ARX->P, delta );
	
	matrix_init( &ARX->theta, p, dimension );
	
	vector_init( &ARX->phi, dimension );


	vector_init( &ARX->yp, p );
	vector_init( &ARX->eps, p );
	vector_init( &ARX->eps2, p );
	vector_init( &ARX->Pphi, dimension );
	vector_init( &ARX->K, dimension );
	matrix_init( &ARX->DeltaTheta, p, dimension);
	matrix_init( &ARX->DeltaP, dimension, dimension );

	return GPC_OK;

}

GPC_ResT ARX_Destroy( ARX_Model *ARX) {
	
	matrix_destroy(&ARX->P);
	matrix_destroy(&ARX->theta);

	vector_destroy(&ARX->phi);

	vector_destroy(&ARX->yp);
	vector_destroy(&ARX->eps);
	vector_destroy(&ARX->eps2);
	vector_destroy(&ARX->Pphi);
	vector_destroy(&ARX->K);
	matrix_destroy(&ARX->DeltaTheta);
	matrix_destroy(&ARX->DeltaP);

	return GPC_OK;

}

GPC_ResT ARMAX_Initialize( ARMAX_Model *ARMAX, unsigned na, unsigned nb, 
		    unsigned nc, unsigned m, unsigned p, double mu,
		    double delta, unsigned FlagSimplyProper){

	unsigned dimension;
	unsigned i, j;

	ARMAX->na = na;
	ARMAX->nb = nb;
	ARMAX->nc = nc;
	ARMAX->m = m;
	ARMAX->p = p;
	ARMAX->mu = mu;
	ARMAX->SimplyProper = FlagSimplyProper;

	dimension = p*na + m*nb + m*FlagSimplyProper + p*nc;

	matrix_init( &ARMAX->P, dimension, dimension );
	matrix_eye( &ARMAX->P, delta );
	matrix_init( &ARMAX->P_RML, p*dimension, p*dimension );
	matrix_eye( &ARMAX->P_RML, delta );
	
	matrix_init( &ARMAX->theta, p, dimension );
	vector_init( &ARMAX->theta_RML, p*dimension );
	
	vector_init( &ARMAX->phi, dimension );
	matrix_init( &ARMAX->psi, p, p*dimension );

	if( !( ARMAX->ALPHA = (matrix **)calloc( p, sizeof(matrix *) ) ) ){
		return GPC_GEN_ERROR;
	}
	if( !( ARMAX->BETA = (matrix **)calloc( p, sizeof(matrix *) ) ) ){
		return GPC_GEN_ERROR;
	}
	if( !( ARMAX->DELTA = (matrix **)calloc( p, sizeof(matrix *) ) ) ){
		return GPC_GEN_ERROR;
	}
	for( i=0; i<p; i++ ){
		if( !( ARMAX->ALPHA[i] = (matrix *)calloc( nc, sizeof(matrix) ) ) ){
			return GPC_GEN_ERROR;
		}
		if( !( ARMAX->BETA[i] = (matrix *)calloc( nc, sizeof(matrix) ) ) ){
			return GPC_GEN_ERROR;
		}
		if( !( ARMAX->DELTA[i] = (matrix *)calloc( nc, sizeof(matrix) ) ) ){
			return GPC_GEN_ERROR;
		}
	}
		
	for( i=0; i<p; i++ ){ 
		for( j=0; j<nc; j++ ){
			matrix_init( &ARMAX->ALPHA[i][j], p, p*na ); 
			matrix_init( &ARMAX->BETA[i][j], p, m*(nb+FlagSimplyProper) ); 
			matrix_init( &ARMAX->DELTA[i][j], p, p*nc ); 
		}
	}

	vector_init( &ARMAX->yp, p );
	vector_init( &ARMAX->eps, p );
	vector_init( &ARMAX->eps2, p );
	vector_init( &ARMAX->Pphi, dimension );
	vector_init( &ARMAX->K, dimension );
	matrix_init( &ARMAX->DeltaTheta, p, dimension);
	matrix_init( &ARMAX->DeltaP, dimension, dimension );

	matrix_init( &ARMAX->PpsiT, p*dimension, p); 
	matrix_init( &ARMAX->psiPT, p, p*dimension);
	matrix_init( &ARMAX->psiPpsiT, p, p);
	matrix_init( &ARMAX->mu_eye, p, p);
	matrix_eye( &ARMAX->mu_eye, mu);
	matrix_init( &ARMAX->tmp_pp, p, p);
	matrix_init( &ARMAX->tmp_p_pdim, p, p*dimension);
	vector_init( &ARMAX->aa, p*(p+1)/2);
	vector_init( &ARMAX->bb, p);
	matrix_init( &ARMAX->deltaP_RML, p*dimension, p*dimension );
	vector_init( &ARMAX->delta_theta_RML, p*dimension );
 
	return GPC_OK;

}

GPC_ResT ARMAX_Destroy( ARMAX_Model *ARMAX) {

	unsigned i, j;
	
	matrix_destroy(&ARMAX->P);
	matrix_destroy( &ARMAX->P_RML);
	matrix_destroy(&ARMAX->theta);
	vector_destroy( &ARMAX->theta_RML );

	vector_destroy(&ARMAX->phi);
	matrix_destroy( &ARMAX->psi );

	vector_destroy(&ARMAX->yp);
	vector_destroy(&ARMAX->eps);
	vector_destroy(&ARMAX->eps2);
	vector_destroy(&ARMAX->Pphi);
	vector_destroy(&ARMAX->K);
	matrix_destroy(&ARMAX->DeltaTheta);
	matrix_destroy(&ARMAX->DeltaP);

	for( i=0; i<ARMAX->p; i++ ){ 
		for( j=0; j<ARMAX->nc; j++ ){
			matrix_destroy( &ARMAX->ALPHA[i][j] ); 
			matrix_destroy( &ARMAX->BETA[i][j] ); 
			matrix_destroy( &ARMAX->DELTA[i][j] ); 
		}
	}
	for( i=0; i<ARMAX->p; i++ ){
		free(ARMAX->ALPHA[i]);
		free(ARMAX->BETA[i]);
		free(ARMAX->DELTA[i]);
	}
	free(ARMAX->ALPHA);
	free(ARMAX->BETA);
	free(ARMAX->DELTA);

	matrix_destroy( &ARMAX->psiPT);
	matrix_destroy( &ARMAX->psiPpsiT);
	matrix_destroy( &ARMAX->PpsiT); 
	matrix_destroy( &ARMAX->mu_eye);
	matrix_destroy( &ARMAX->tmp_pp);
	matrix_destroy( &ARMAX->tmp_p_pdim);
	vector_destroy( &ARMAX->aa);
	vector_destroy( &ARMAX->bb);
	matrix_destroy( &ARMAX->deltaP_RML );
	vector_destroy( &ARMAX->delta_theta_RML );
 
	return GPC_OK;

}



GPC_ResT ARX_UpdatePhi( ARX_Model *ARX, vector *y, vector *u){

	unsigned i;
	unsigned na, nb, p, m;

	na = ARX->na;
	nb = ARX->nb;
	p = ARX->p;
	m = ARX->m;

	for( i=0; i<p*(na-1); i++ ){
		ARX->phi.vec[p*na-1-i] = ARX->phi.vec[p*(na-1)-1-i];
	}
	for( i=0; i<m*(nb-1+ARX->SimplyProper); i++ ){
		ARX->phi.vec[p*na+m*(nb+ARX->SimplyProper)-1-i] = ARX->phi.vec[p*na+m*(nb-1+ARX->SimplyProper)-1-i];
	}
	
	for( i=0; i<p; i++ ){
		ARX->phi.vec[i] = y->vec[i];
	}
	
	for( i=0; i<m; i++ ){
		ARX->phi.vec[i+p*na] = u->vec[i];
	}
	
	return GPC_OK;
}

GPC_ResT ARMAX_UpdatePhi( ARMAX_Model *ARMAX, vector *y, vector *u, vector *eps){

	unsigned i;
	unsigned na, nb, nc, p, m;
	unsigned FlagSimplyProper;

	na = ARMAX->na;
	nb = ARMAX->nb;
	nc = ARMAX->nc;
	p = ARMAX->p;
	m = ARMAX->m;
	FlagSimplyProper = ARMAX->SimplyProper;

	for( i=0; i<p*(na-1); i++ ){
		ARMAX->phi.vec[ p*na -1 -i ] = ARMAX->phi.vec[ p*(na-1) -1 -i];
	}
	for( i=0; i<m*(nb-1+FlagSimplyProper); i++ ){
		ARMAX->phi.vec[ p*na + m*(nb+FlagSimplyProper) -1 -i ] = ARMAX->phi.vec[ p*na + m*(nb-1+FlagSimplyProper) -1 -i ];
	}
	for( i=0; i<p*(nc-1); i++ ){
		ARMAX->phi.vec[ p*na + m*(nb+FlagSimplyProper) + p*nc -1 -i ] = ARMAX->phi.vec[ p*na + m*(nb+FlagSimplyProper) + p*(nc-1) -1 -i ];
	}
	
	for( i=0; i<p; i++ ){
		ARMAX->phi.vec[i] = y->vec[i];
	}
	for( i=0; i<m; i++ ){
		ARMAX->phi.vec[i+p*na] = u->vec[i];
	}
	for( i=0; i<p; i++ ){
		ARMAX->phi.vec[i+p*na + m*(nb+FlagSimplyProper)] = eps->vec[i];
	}
	
	return GPC_OK;
}
		

GPC_ResT ARX_RLS( ARX_Model *ARX, vector *y){

	double r, a, tr;
	unsigned i;
	
	a = 0.;

	matrix_vector_prod( &ARX->theta, &ARX->phi, &ARX->yp );
	vector_sum( y, &ARX->yp, &ARX->eps, -1. );

	matrix_vector_prod( &ARX->P, &ARX->phi, &ARX->Pphi );

	scalar_prod( &ARX->phi, &ARX->Pphi, &r );

	vector_copy( &ARX->K, &ARX->Pphi, 1./(1.+r));	

	/* IT SHOULD INCREASE ROBUSTNESS */
	for( i=0; i<ARX->p; i++){
		if( ARX->eps.vec[i] >= 0){
			ARX->eps2.vec[i] = ARX->eps.vec[i]/(1.+a*ARX->eps.vec[i]);
		} else {
			ARX->eps2.vec[i] = ARX->eps.vec[i]/(1.-a*ARX->eps.vec[i]);
		}
	}

	//vector_vector_prod( &ARX->eps, &ARX->K, &ARX->DeltaTheta, 1. );
	vector_vector_prod( &ARX->eps2, &ARX->K, &ARX->DeltaTheta, 1. );
	matrix_sum( &ARX->theta, &ARX->DeltaTheta, &ARX->theta, 1. );

	vector_vector_prod( &ARX->Pphi, &ARX->Pphi, &ARX->DeltaP, 1./(r+ARX->mu) );
	matrix_sum( &ARX->P, &ARX->DeltaP, &ARX->P, -1. );
	matrix_copy( &ARX->P, &ARX->P, 1./ARX->mu );

	//tr = matrix_trace( &ARX->P );
	//printf("traccia P %e\n", tr);
	//if ( (matrix_trace(&ARX->P)) < 1.e5) {
		//matrix_copy( &ARX->P, &ARX->P, 1.e5 );
	//	matrix_eye( &ARX->P, 1.e7 );
	//	printf("tr(P2) = %e\n",matrix_trace(&ARX->P));
	//}

	return GPC_OK;

}

GPC_ResT ARMAX_ELS( ARMAX_Model *ARMAX, vector *y){

	double r, a;
	unsigned i;
	
	a = 0.;

	matrix_vector_prod( &ARMAX->theta, &ARMAX->phi, &ARMAX->yp );
	vector_sum( y, &ARMAX->yp, &ARMAX->eps, -1. );

	matrix_vector_prod( &ARMAX->P, &ARMAX->phi, &ARMAX->Pphi );

	scalar_prod( &ARMAX->phi, &ARMAX->Pphi, &r );

	vector_copy( &ARMAX->K, &ARMAX->Pphi, 1./(1.+r));	

	/* IT SHOULD INCREASE ROBUSTNESS */
	for( i=0; i<ARMAX->p; i++){
		if( ARMAX->eps.vec[i] >= 0){
			ARMAX->eps2.vec[i] = ARMAX->eps.vec[i]/(1.+a*ARMAX->eps.vec[i]);
		} else {
			ARMAX->eps2.vec[i] = ARMAX->eps.vec[i]/(1.-a*ARMAX->eps.vec[i]);
		}
	}

	//vector_vector_prod( &ARMAX->eps, &ARMAX->K, &ARMAX->DeltaTheta, 1. );
	vector_vector_prod( &ARMAX->eps2, &ARMAX->K, &ARMAX->DeltaTheta, 1. );
	matrix_sum( &ARMAX->theta, &ARMAX->DeltaTheta, &ARMAX->theta, 1. );

	vector_vector_prod( &ARMAX->Pphi, &ARMAX->Pphi, &ARMAX->DeltaP, 1./(r+ARMAX->mu) );
	//vector_vector_prod( &ARMAX->Pphi, &ARMAX->Pphi, &ARMAX->DeltaP, 1./(r+1.) );
	matrix_sum( &ARMAX->P, &ARMAX->DeltaP, &ARMAX->P, -1. );
	matrix_copy( &ARMAX->P, &ARMAX->P, 1./ARMAX->mu );
	//printf("tr(P) = %e\n",matrix_trace(&ARMAX->P));
/*
	if ( (matrix_trace(&ARMAX->P)) < 1.e3) {
		matrix_copy( &ARMAX->P, &ARMAX->P, 1.e2 );
		printf("tr(P2) = %e\n",matrix_trace(&ARMAX->P));
	}
*/
	return GPC_OK;

}

GPC_ResT ARMAX_RML( ARMAX_Model *ARMAX, vector *y){

	double a, c_rr_kk;
	unsigned r, k, i, j, q;
	double alpha, beta, delta;
	unsigned na, nb, nc, p, m, FlagSimplyProper;
	
	na = ARMAX->na;
	nb = ARMAX->nb;
	nc = ARMAX->nc;
	m = ARMAX->m;
	p = ARMAX->p;
	FlagSimplyProper = ARMAX->SimplyProper;

	a = 0.;


	matrix_vector_prod( &ARMAX->theta, &ARMAX->phi, &ARMAX->yp );
	vector_sum( y, &ARMAX->yp, &ARMAX->eps, -1. );
	for( r=0; r<p; r++ ){
		for( i=0; i<p; i++ ){
			for( j=0; j<p; j++ ){
				for( q=0; q<na; q++ ){
					if( i==r ){
						alpha = -ARMAX->phi.vec[q*p+j];
					} else {
						alpha = 0.;
					} 
					for( k=0; k<nc; k++ ){
						c_rr_kk = ARMAX->theta.mat[r][na*p+(nb+FlagSimplyProper)*m + k*p + r];
						alpha -= c_rr_kk*ARMAX->ALPHA[r][k].mat[i][q*p+j];
					}
					for( k=nc-1; k>0; k-- ){
						ARMAX->ALPHA[r][k].mat[i][q*p+j] = ARMAX->ALPHA[r][k-1].mat[i][q*p+j];
					}
					ARMAX->ALPHA[r][0].mat[i][q*p+j] = alpha;
					ARMAX->psi.mat[r][q*p*p + i*p + j] = alpha; 
				}
			}
			for( j=0; j<m; j++ ){
				for( q=0; q<nb+FlagSimplyProper; q++ ){
					if( i==r ){
						beta = -ARMAX->phi.vec[p*na+q*m+j];
					} else {
						beta = 0.;
					} 
					for( k=0; k<nc; k++ ){
						c_rr_kk = ARMAX->theta.mat[r][na*p+(nb+FlagSimplyProper)*m + k*p + r];
						beta -= c_rr_kk*ARMAX->BETA[r][k].mat[i][q*m+j];
					}
					for( k=nc-1; k>0; k-- ){
						ARMAX->BETA[r][k].mat[i][q*m+j] = ARMAX->BETA[r][k-1].mat[i][q*m+j];
					}
					ARMAX->BETA[r][0].mat[i][q*m+j] = beta;
					ARMAX->psi.mat[r][p*p*na + q*p*m+i*m+j] = beta; 
				}
			}
			for( j=0; j<p; j++ ){
				for( q=0; q<nc; q++ ){
					if( i==r ){
						delta = -ARMAX->phi.vec[na*p+(nb+FlagSimplyProper)*m+q*p+j];
					} else {
						delta = 0.;
					} 
					for( k=0; k<nc; k++ ){
						c_rr_kk = ARMAX->theta.mat[r][na*p+(nb+FlagSimplyProper)*m + k*p + r];
						delta -= c_rr_kk*ARMAX->DELTA[r][k].mat[i][q*p+j];
					}
					for( k=nc-1; k>0; k-- ){
						ARMAX->DELTA[r][k].mat[i][q*p+j] = ARMAX->DELTA[r][k-1].mat[i][q*p+j];
					}
					ARMAX->DELTA[r][0].mat[i][q*p+j] = delta;
					ARMAX->psi.mat[r][p*p*na + p*m*(FlagSimplyProper+nb) + q*p*p+i*p+j] = delta; 
				}
			}

		}
	}	


	matrix_prod_transpose( &ARMAX->P_RML , &ARMAX->psi, &ARMAX->PpsiT, 1. );
	matrix_prod_transpose( &ARMAX->psi , &ARMAX->P_RML, &ARMAX->psiPT, 1. );
	matrix_prod( &ARMAX->psi, &ARMAX->PpsiT, &ARMAX->psiPpsiT, 1.);
	matrix_sum( &ARMAX->psiPpsiT, &ARMAX->mu_eye, &ARMAX->tmp_pp, 1.);
	//matrix_write(&ARMAX->tmp_pp, stdout, W_M_BIN);
	//matrix_write(&ARMAX->psi, stdout, W_M_BIN);
	LinearSystemM( &ARMAX->tmp_pp, &ARMAX->psiPT, &ARMAX->tmp_p_pdim, &ARMAX->aa, &ARMAX->bb);
	matrix_prod( &ARMAX->PpsiT, &ARMAX->tmp_p_pdim, &ARMAX->deltaP_RML, 1.);
	matrix_sum( &ARMAX->P_RML, &ARMAX->deltaP_RML, &ARMAX->P_RML, -1.);
	matrix_copy( &ARMAX->P_RML, &ARMAX->P_RML, 1./ARMAX->mu );
	//matrix_write(&ARMAX->deltaP_RML, stdout, W_M_BIN);
	//vector_write(&ARMAX->eps, stdout, W_M_BIN);

	matrix_prod_transpose( &ARMAX->P_RML , &ARMAX->psi, &ARMAX->PpsiT, 1. );
	matrix_vector_prod( &ARMAX->PpsiT, &ARMAX->eps, &ARMAX->delta_theta_RML);
	vector_sum( &ARMAX->theta_RML, &ARMAX->delta_theta_RML, &ARMAX->theta_RML, -1.); 
	//vector_write(&ARMAX->theta_RML, stdout, W_M_BIN);

	for( i=0; i<p; i++ ){
		for( j=0; j<p; j++ ){
			for( q=0; q<na; q++ ){
				ARMAX->theta.mat[i][q*p+j] = ARMAX->theta_RML.vec[q*p*p+i*p+j];
			}
		}
		for( j=0; j<m; j++ ){
			for( q=0; q<nb+FlagSimplyProper; q++ ){
				ARMAX->theta.mat[i][na*p + q*m + j] = ARMAX->theta_RML.vec[p*p*na + q*p*m+i*m+j];
			}
		}
		for( j=0; j<p; j++ ){
			for( q=0; q<nc; q++ ){
				ARMAX->theta.mat[i][na*p + (nb+FlagSimplyProper)*m + q*p + j] = ARMAX->theta_RML.vec[p*p*na + p*m*(nb+FlagSimplyProper) + q*p*p+i*p+j];
			}
		}
	}

	//matrix_write(&ARMAX->theta, stdout, W_M_BIN);
	//getchar();
	//matrix_vector_prod( &ARMAX->theta, &ARMAX->phi, &ARMAX->yp );
	//vector_sum( y, &ARMAX->yp, &ARMAX->eps, -1. );

}
					



GPC_ResT GPC_Initialize( GPC_Model *GPC, unsigned na, unsigned nb, unsigned nc,
				unsigned m, unsigned p, unsigned s,
				double lambda ){

	unsigned i;

	GPC->na = na;
	GPC->nb = nb;
	GPC->nc = nc;
	GPC->m = m;
	GPC->p = p;
	GPC->s = s;
	GPC->lambda = lambda;

	matrix_init( &GPC->A, s*p, p*na );
	matrix_init( &GPC->B, s*p, m*nb );
	matrix_init( &GPC->C, s*p, s*m );
	matrix_init( &GPC->B2, s*p, m*(nb-1) );
	matrix_init( &GPC->C2, s*p, m*(s+1) );
	if( nc > 0 ){
		matrix_init( &GPC->D, s*p, p*nc );
	}

	vector_init( &GPC->Yp, p*na );
	if( nc > 0 ){
		vector_init( &GPC->Ep, p*nc );
	}
	vector_init( &GPC->Up, m*nb );
	vector_init( &GPC->Up2, m*(nb-1) );
	vector_init( &GPC->Us, s*m );
	vector_init( &GPC->Us2, m*(s+1) );
	vector_init( &GPC->Ydes, s*p );

	matrix_init( &GPC->Wy, s*p, s*p);
	matrix_init( &GPC->Wu, s*m, s*m);
	matrix_eye( &GPC->Wy, 1.);
	matrix_eye( &GPC->Wu, 1.);
#if 0
	for (i=0; i<s; i++){
		GPC->Wy.mat[3+i*5][3+i*5] = 1.e+2;
		GPC->Wy.mat[4+i*5][4+i*5] = 1.e+2;
		GPC->Wu.mat[i*4][i*4] = 1.e-3;
		//GPC->Wu.mat[2+i*3][2+i*3] = 1.e-1;
	}
	for (i=0; i<3*p; i++){
		//GPC->Wy.mat[s*p-1-i][s*p-1-i] = 0.;
		GPC->Wy.mat[i][i] = 10.;
	}
	for (i=0; i<1*m; i++){
		GPC->Wu.mat[i][i] = 100.;
	}
	for (i=0; i<s; i++){
		GPC->Wu.mat[i*m+1][i*m+1] = 10000.;
	}
#endif
	matrix_init( &GPC->WyC, s*p, s*m);

	matrix_init( &GPC->tmp_pp1, p, p);
	matrix_init( &GPC->tmp_pp2, p, p);
	matrix_init( &GPC->tmp_pp3, p, p);
	matrix_init( &GPC->tmp_pp4, p, p);
	matrix_init( &GPC->tmp_pm1, p, m);
	matrix_init( &GPC->tmp_pm2, p, m);
	matrix_init( &GPC->tmp_pm3, p, m);
	matrix_init( &GPC->tmp_pm4, p, m);
	matrix_init( &GPC->tmp_mm, m, m);

	matrix_init( &GPC->tmp_smsm1, s*m, s*m);
	matrix_init( &GPC->tmp_smsm2, s*m, s*m);
	matrix_init( &GPC->tmp_smsm3, s*m, s*m);
	matrix_init( &GPC->tmp_sm1sm11, m*(s+1), m*(s+1));
	matrix_init( &GPC->tmp_sm1sm12, m*(s+1), m*(s+1));
	matrix_init( &GPC->tmp_sm1sm13, m*(s+1), m*(s+1));
	matrix_init( &GPC->tmp_sp_m, s*p, m);
	vector_init( &GPC->tmp_sp1, s*p);
	vector_init( &GPC->tmp_sp2, s*p);
	vector_init( &GPC->tmp_sp3, s*p);
	vector_init( &GPC->tmp_sp4, s*p);
	vector_init( &GPC->tmp_sp5, s*p);
	vector_init( &GPC->tmp_sm, s*m);
	vector_init( &GPC->tmp_sm1, m*(s+1));
	vector_init( &GPC->tmp_aa, ((s*m)*(s*m+1))/2 );
	vector_init( &GPC->tmp_aa1, ((m*(s+1))*(m*(s+1)+1))/2 );

	matrix_init( &GPC->K, m*s, p*s);

	return GPC_OK;
}

GPC_ResT GPC_Destroy( GPC_Model *GPC ){

	matrix_destroy( &GPC->A );
	matrix_destroy( &GPC->B );
	matrix_destroy( &GPC->C );
	matrix_destroy( &GPC->D );
	matrix_destroy( &GPC->B2 );
	matrix_destroy( &GPC->C2 );

	vector_destroy( &GPC->Yp );
	vector_destroy( &GPC->Up );
	vector_destroy( &GPC->Ep );
	vector_destroy( &GPC->Us );
	vector_destroy( &GPC->Up2 );
	vector_destroy( &GPC->Us2 );
	vector_destroy( &GPC->Ydes );

	matrix_destroy( &GPC->Wy);
	matrix_destroy( &GPC->Wu);
	matrix_destroy( &GPC->WyC);

	matrix_destroy(&GPC->tmp_pp1);
	matrix_destroy(&GPC->tmp_pp2);
	matrix_destroy(&GPC->tmp_pp3);
	matrix_destroy(&GPC->tmp_pp4);
	matrix_destroy(&GPC->tmp_pm1);
	matrix_destroy(&GPC->tmp_pm2);
	matrix_destroy(&GPC->tmp_pm3);
	matrix_destroy(&GPC->tmp_pm4);
	matrix_destroy(&GPC->tmp_mm);

	matrix_destroy(&GPC->tmp_smsm1);
	matrix_destroy(&GPC->tmp_smsm2);
	matrix_destroy(&GPC->tmp_smsm3);
	matrix_destroy(&GPC->tmp_sm1sm11);
	matrix_destroy(&GPC->tmp_sm1sm12);
	matrix_destroy(&GPC->tmp_sm1sm13);
	matrix_destroy(&GPC->tmp_sp_m);
	vector_destroy(&GPC->tmp_sp1);
	vector_destroy(&GPC->tmp_sp2);
	vector_destroy(&GPC->tmp_sp3);
	vector_destroy(&GPC->tmp_sp4);
	vector_destroy(&GPC->tmp_sp5);
	vector_destroy(&GPC->tmp_sm);
	vector_destroy(&GPC->tmp_sm1);
	vector_destroy(&GPC->tmp_aa);
	vector_destroy(&GPC->tmp_aa1);

	matrix_destroy(&GPC->K);

	return GPC_OK;
}

GPC_ResT GPC_PredictionFunction( GPC_Model *GPC, matrix *theta, unsigned FlagSimplyProper, unsigned Flag ){

	unsigned i, j;
	unsigned na, nb, p, m, s, nc;

	na = GPC->na;
	nb = GPC->nb;
	nc = GPC->nc;
	p = GPC->p;
	m = GPC->m;
	s = GPC->s;

	for( i=0; i<s; i++ ){
		for( j=0; j<na; j++ ){
			if( i==0 ){
				sub_matrix_extract( theta, &GPC->tmp_pp1, 0, p*j);
				sub_matrix_insert( &GPC->A, &GPC->tmp_pp1, p*(s-1) , p*j);	
	
				//sub_matrix_extract( theta, &GPC->tmp_pm1, 0, p*na + FlagSimplyProper*m + m*j);
				//sub_matrix_insert( &GPC->B, &GPC->tmp_pm1, p*(s-1) , m*j);
			} else {
				if( j<(na-1) ){
					sub_matrix_extract( &GPC->A, &GPC->tmp_pp1, p*(s-i), 0 );
					sub_matrix_extract( &GPC->A, &GPC->tmp_pp2, p*(s-1), p*j );
					sub_matrix_extract( &GPC->A, &GPC->tmp_pp3, p*(s-i), p*(j+1) );
					matrix_prod( &GPC->tmp_pp1, &GPC->tmp_pp2, &GPC->tmp_pp4, 1. );
					matrix_sum( &GPC->tmp_pp4, &GPC->tmp_pp3, &GPC->tmp_pp2, 1. );
					sub_matrix_insert( &GPC->A, &GPC->tmp_pp2, p*(s-i-1), p*j );

					//sub_matrix_extract( &GPC->B, &GPC->tmp_pm2, p*(s-1), m*j );
					//sub_matrix_extract( &GPC->B, &GPC->tmp_pm3, p*(s-i), m*(j+1) );
					//matrix_prod( &GPC->tmp_pp1, &GPC->tmp_pm2, &GPC->tmp_pm4, 1. );
					//matrix_sum( &GPC->tmp_pm4, &GPC->tmp_pm3, &GPC->tmp_pm1, 1. );
					//sub_matrix_insert( &GPC->B, &GPC->tmp_pm1, p*(s-i-1), m*j );
				} else {
					sub_matrix_extract( &GPC->A, &GPC->tmp_pp1, p*(s-i), 0 );
					sub_matrix_extract( &GPC->A, &GPC->tmp_pp2, p*(s-1), p*j );
					matrix_prod( &GPC->tmp_pp1, &GPC->tmp_pp2, &GPC->tmp_pp4, 1. );
					sub_matrix_insert( &GPC->A, &GPC->tmp_pp4, p*(s-i-1), p*j );
				}
			}
		}
		for( j=0; j<nb; j++ ){
			if( i==0 ){
				sub_matrix_extract( theta, &GPC->tmp_pm1, 0, p*na + FlagSimplyProper*m + m*j);
				sub_matrix_insert( &GPC->B, &GPC->tmp_pm1, p*(s-1) , m*j);	
			} else {
				if( j<(nb-1) ){
					sub_matrix_extract( &GPC->A, &GPC->tmp_pp1, p*(s-i), 0 );
					sub_matrix_extract( &GPC->B, &GPC->tmp_pm2, p*(s-1), m*j );
					sub_matrix_extract( &GPC->B, &GPC->tmp_pm3, p*(s-i), m*(j+1) );
					matrix_prod( &GPC->tmp_pp1, &GPC->tmp_pm2, &GPC->tmp_pm4, 1. );
					matrix_sum( &GPC->tmp_pm4, &GPC->tmp_pm3, &GPC->tmp_pm1, 1. );
					sub_matrix_insert( &GPC->B, &GPC->tmp_pm1, p*(s-i-1), m*j );
				} else {
					sub_matrix_extract( &GPC->A, &GPC->tmp_pp1, p*(s-i), 0 );
					sub_matrix_extract( &GPC->B, &GPC->tmp_pm2, p*(s-1), m*j );
					matrix_prod( &GPC->tmp_pp1, &GPC->tmp_pm2, &GPC->tmp_pm4, 1. );
					sub_matrix_insert( &GPC->B, &GPC->tmp_pm4, p*(s-i-1), m*j );
				}
			}
		}
		for( j=0; j<nc; j++ ){
			if( i==0 ){
				sub_matrix_extract( theta, &GPC->tmp_pp1, 0, p*na + FlagSimplyProper*m + m*nb + p*j);
				sub_matrix_insert( &GPC->D, &GPC->tmp_pp1, p*(s-1) , p*j);	
			} else {
				if( j<(nc-1) ){
					sub_matrix_extract( &GPC->A, &GPC->tmp_pp1, p*(s-i), 0 );
					sub_matrix_extract( &GPC->D, &GPC->tmp_pp2, p*(s-1), p*j );
					sub_matrix_extract( &GPC->D, &GPC->tmp_pp3, p*(s-i), p*(j+1) );
					matrix_prod( &GPC->tmp_pp1, &GPC->tmp_pp2, &GPC->tmp_pp4, 1. );
					matrix_sum( &GPC->tmp_pp4, &GPC->tmp_pp3, &GPC->tmp_pp2, 1. );
					sub_matrix_insert( &GPC->D, &GPC->tmp_pp2, p*(s-i-1), p*j );
				} else {
					sub_matrix_extract( &GPC->A, &GPC->tmp_pp1, p*(s-i), 0 );
					sub_matrix_extract( &GPC->D, &GPC->tmp_pp2, p*(s-1), p*j );
					matrix_prod( &GPC->tmp_pp1, &GPC->tmp_pp2, &GPC->tmp_pp4, 1. );
					sub_matrix_insert( &GPC->D, &GPC->tmp_pp4, p*(s-i-1), p*j );
				}
			}
		}

		if( i==0 ){
			if( FlagSimplyProper == 1){
				sub_matrix_extract( theta, &GPC->tmp_pm1, 0, p*na );
			} else {
				matrix_null( &GPC->tmp_pm1 );
			}
		} else {
			sub_matrix_extract( &GPC->A, &GPC->tmp_pp1, p*(s-i), 0 );
			sub_matrix_extract( &GPC->C, &GPC->tmp_pm2, 0, 0 );
			sub_matrix_extract( &GPC->B, &GPC->tmp_pm3, p*(s-i), 0 );
			matrix_prod( &GPC->tmp_pp1, &GPC->tmp_pm2, &GPC->tmp_pm4, 1. );
			matrix_sum( &GPC->tmp_pm4, &GPC->tmp_pm3, &GPC->tmp_pm1, 1. );
		}
		for( j=0; j<s-i; j++ ){
			sub_matrix_insert( &GPC->C, &GPC->tmp_pm1, j*p, j*m+i*m );
		}
			
	}

	/* se voglio calcolare il controllo al passo K sfruttando la misura al passo K
	   per un modello strettamente proprio devo modificare le matrici B e C
	*/
	if( Flag == 1 ){
		sub_matrix_extract( &GPC->B, &GPC->B2, 0, m);
		sub_matrix_insert( &GPC->C2, &GPC->C, 0, 0);
		sub_matrix_extract( &GPC->B, &GPC->tmp_sp_m, 0, 0);
		sub_matrix_insert( &GPC->C2, &GPC->tmp_sp_m, 0, s*m);
	}
	
	return GPC_OK;
}
	

GPC_ResT GPC_VectorUpdate( GPC_Model *GPC, vector *y_k, vector *u_k, vector *e_k ){

	unsigned i;
	unsigned na, nb, p, m, s, nc;

	na = GPC->na;
	nb = GPC->nb;
	nc = GPC->nc;
	p = GPC->p;
	m = GPC->m;
	s = GPC->s;
	
	for( i=0; i<p*(na-1); i++ ){
		GPC->Yp.vec[p*na-1-i] = GPC->Yp.vec[p*(na-1)-1-i];
	}
	for( i=0; i<m*(nb-1); i++ ){
		GPC->Up.vec[m*nb-1-i] = GPC->Up.vec[m*(nb-1)-1-i];
	}
	if( nc > 0 ){
		for( i=0; i<p*(nc-1); i++ ){
			GPC->Ep.vec[p*nc-1-i] = GPC->Ep.vec[p*(nc-1)-1-i];
		}
	}
	
	for( i=0; i<p; i++ ){
		GPC->Yp.vec[i] = y_k->vec[i];
	}
	for( i=0; i<m; i++ ){
		GPC->Up.vec[i] = u_k->vec[i];
	}
	if( nc > 0 ){
		for( i=0; i<p; i++ ){
			GPC->Ep.vec[i] = e_k->vec[i];
		}
	}

	return GPC_OK;	

}

GPC_ResT GPC_VectorUpdate2( GPC_Model *GPC, vector *y_k, vector *u_k, vector *e_k ){

	unsigned i;
	unsigned na, nb, nc, p, m, s;

	na = GPC->na;
	nb = GPC->nb;
	nc = GPC->nc;
	p = GPC->p;
	m = GPC->m;
	s = GPC->s;
	
	for( i=0; i<p*(na-1); i++ ){
		GPC->Yp.vec[p*na-1-i] = GPC->Yp.vec[p*(na-1)-1-i];
	}
	for( i=0; i<m*(nb-2); i++ ){
		GPC->Up2.vec[m*(nb-1)-1-i] = GPC->Up2.vec[m*(nb-2)-1-i];
	}
	if( nc > 0 ){
		for( i=0; i<p*(nc-1); i++ ){
			GPC->Ep.vec[p*nc-1-i] = GPC->Ep.vec[p*(nc-1)-1-i];
		}
	}
	
	for( i=0; i<p; i++ ){
		GPC->Yp.vec[i] = y_k->vec[i];
	}
	for( i=0; i<m; i++ ){
		GPC->Up2.vec[i] = u_k->vec[i];
	}
	if( nc > 0 ){
		for( i=0; i<p; i++ ){
			GPC->Ep.vec[i] = e_k->vec[i];
		}
	}

	return GPC_OK;	

}

GPC_ResT GPC_ControlW( GPC_Model *GPC, vector *Us_ID ){

	//unsigned i;
	unsigned na, nb, nc, p, m, s;
	unsigned i, r, c, j ,k;
	//FILE *fh;

	na = GPC->na;
	nb = GPC->nb;
	nc = GPC->nc;
	p = GPC->p;
	m = GPC->m;
	s = GPC->s;

	matrix_null(&GPC->tmp_pp1);
	matrix_null(&GPC->tmp_pp2);
	matrix_null(&GPC->tmp_pp3);
	matrix_null(&GPC->tmp_pp4);
	matrix_null(&GPC->tmp_pm1);
	matrix_null(&GPC->tmp_pm2);
	matrix_null(&GPC->tmp_pm3);
	matrix_null(&GPC->tmp_pm4);

	matrix_null(&GPC->tmp_smsm1);
	matrix_null(&GPC->tmp_smsm2);
	matrix_null(&GPC->tmp_smsm3);
	vector_null(&GPC->tmp_sp1);
	vector_null(&GPC->tmp_sp2);
	vector_null(&GPC->tmp_sp3);
	vector_null(&GPC->tmp_sm);
	vector_null(&GPC->tmp_aa);

	matrix_prod( &GPC->Wy, &GPC->C, &GPC->WyC, 1.);
	matrix_transpose_prod( &GPC->C, &GPC->WyC, &GPC->tmp_smsm1, 1.);
	//matrix_transpose_prod( &GPC->C, &GPC->C, &GPC->tmp_smsm1, 1.);

	matrix_copy( &GPC->tmp_smsm2, &GPC->Wu, GPC->lambda );
	matrix_null(&GPC->tmp_smsm3);
	matrix_sum( &GPC->tmp_smsm1, &GPC->tmp_smsm2, &GPC->tmp_smsm3, 1. );

	matrix_vector_prod( &GPC->A, &GPC->Yp, &GPC->tmp_sp1);
	matrix_vector_prod( &GPC->B, &GPC->Up, &GPC->tmp_sp2);
	matrix_vector_prod( &GPC->C, Us_ID, &GPC->tmp_sp3);
	if( nc > 0 ){
		matrix_vector_prod( &GPC->D, &GPC->Ep, &GPC->tmp_sp5);
	}
	vector_null(&GPC->tmp_sp4);
	vector_sum( &GPC->Ydes, &GPC->tmp_sp1, &GPC->tmp_sp4, -1. );
	vector_null(&GPC->tmp_sp1);
	vector_sum( &GPC->tmp_sp4, &GPC->tmp_sp2, &GPC->tmp_sp1, -1. );
	vector_null(&GPC->tmp_sp2);
	vector_sum( &GPC->tmp_sp1, &GPC->tmp_sp3, &GPC->tmp_sp2, -1. );
	if( nc > 0 ){
		vector_null(&GPC->tmp_sp4);
		vector_sum( &GPC->tmp_sp2, &GPC->tmp_sp5, &GPC->tmp_sp4, -1. );
	}

	//printf("WARNING!!!!\n");
	//matrix_vector_prod( &GPC->K, &GPC->tmp_sp1, &GPC->Us );
	
	if( nc > 0 ){
		vector_null(&GPC->tmp_sp2);
		matrix_vector_prod( &GPC->Wy, &GPC->tmp_sp4, &GPC->tmp_sp2);
		matrixT_vector_prod( &GPC->C, &GPC->tmp_sp2, &GPC->tmp_sm);
	} else {
		vector_null(&GPC->tmp_sp4);
		matrix_vector_prod( &GPC->Wy, &GPC->tmp_sp2, &GPC->tmp_sp4);
		matrixT_vector_prod( &GPC->C, &GPC->tmp_sp4, &GPC->tmp_sm);
	}
	vector_null(&GPC->Us);
	vector_null(&GPC->tmp_aa);
	LinearSystemV( &GPC->tmp_smsm3, &GPC->tmp_sm, &GPC->Us, &GPC->tmp_aa );

	return GPC_OK;

}


GPC_ResT GPC_Control( GPC_Model *GPC, vector *Us_ID ){

	//unsigned i;
	unsigned na, nb, nc, p, m, s;
	unsigned i, r, c, j ,k;
	//FILE *fh;

	na = GPC->na;
	nb = GPC->nb;
	nc = GPC->nc;
	p = GPC->p;
	m = GPC->m;
	s = GPC->s;

	matrix_null(&GPC->tmp_pp1);
	matrix_null(&GPC->tmp_pp2);
	matrix_null(&GPC->tmp_pp3);
	matrix_null(&GPC->tmp_pp4);
	matrix_null(&GPC->tmp_pm1);
	matrix_null(&GPC->tmp_pm2);
	matrix_null(&GPC->tmp_pm3);
	matrix_null(&GPC->tmp_pm4);

	matrix_null(&GPC->tmp_smsm1);
	matrix_null(&GPC->tmp_smsm2);
	matrix_null(&GPC->tmp_smsm3);
	vector_null(&GPC->tmp_sp1);
	vector_null(&GPC->tmp_sp2);
	vector_null(&GPC->tmp_sp3);
	vector_null(&GPC->tmp_sm);
	vector_null(&GPC->tmp_aa);

	//matrix_prod( &GPC->Wy, &GPC->C, &GPC->WyC, 1.);
	//matrix_transpose_prod( &GPC->C, &GPC->WyC, &GPC->tmp_smsm1, 1.);
	//matrix_transpose_prod( &GPC->C, &GPC->C, &GPC->tmp_smsm1, 1.);
#if 1
	/* eseguo l'operazione C^T*C sfruttando il fatto che la matrice C2 ha una struttura
	  particolare, in questo modo rendo il codice più efficiente */
	for( i=0; i<s; i++){
		/* blocco 1 */
		for(r=0; r<p; r++){
			for(c=0; c<m; c++){
				GPC->tmp_pm1.mat[r][c] = GPC->C.mat[r][c+(i)*m];
			}
		}
		for(j=i; j<s;j++){
			/* blocco 2 */
			for(r=0; r<p; r++){
				for(c=0; c<m; c++){
					GPC->tmp_pm2.mat[r][c] = GPC->C.mat[r][c+(j)*m];
				}
			}
			matrix_transpose_prod( &GPC->tmp_pm1, &GPC->tmp_pm2, &GPC->tmp_mm, 1.);
			if (j!=i){
				for(k=j;k<s;k++){
					for(r=0; r<m; r++){
						for(c=0; c<m; c++){
							GPC->tmp_smsm1.mat[m*(k-j)+i*m+r][m*(k)+c] += GPC->tmp_mm.mat[r][c];
						}
					}
				}
				matrix_transpose_prod( &GPC->tmp_pm2, &GPC->tmp_pm1, &GPC->tmp_mm, 1.);
				for(k=j;k<s;k++){
					for(r=0; r<m; r++){
						for(c=0; c<m; c++){
							GPC->tmp_smsm1.mat[m*(k)+r][m*(k-j)+i*m+c] += GPC->tmp_mm.mat[r][c];
						}
					}
				}
			} else {
				for(k=j;k<s;k++){
					for(r=0; r<m; r++){
						for(c=0; c<m; c++){
							GPC->tmp_smsm1.mat[m*(k)+r][m*(k)+c] += GPC->tmp_mm.mat[r][c];
						}
					}
				}
			}
		}
	}
#endif

	matrix_eye( &GPC->tmp_smsm2, GPC->lambda );
	//matrix_copy( &GPC->tmp_smsm2, &GPC->Wu, GPC->lambda );
	matrix_null(&GPC->tmp_smsm3);
	matrix_sum( &GPC->tmp_smsm1, &GPC->tmp_smsm2, &GPC->tmp_smsm3, 1. );

	matrix_vector_prod( &GPC->A, &GPC->Yp, &GPC->tmp_sp1);
	matrix_vector_prod( &GPC->B, &GPC->Up, &GPC->tmp_sp2);
	matrix_vector_prod( &GPC->C, Us_ID, &GPC->tmp_sp3);
	if( nc > 0 ){
		matrix_vector_prod( &GPC->D, &GPC->Ep, &GPC->tmp_sp5);
	}
	vector_null(&GPC->tmp_sp4);
	vector_sum( &GPC->Ydes, &GPC->tmp_sp1, &GPC->tmp_sp4, -1. );
	vector_null(&GPC->tmp_sp1);
	vector_sum( &GPC->tmp_sp4, &GPC->tmp_sp2, &GPC->tmp_sp1, -1. );
	vector_null(&GPC->tmp_sp2);
	vector_sum( &GPC->tmp_sp1, &GPC->tmp_sp3, &GPC->tmp_sp2, -1. );
	if( nc > 0 ){
		vector_null(&GPC->tmp_sp4);
		vector_sum( &GPC->tmp_sp2, &GPC->tmp_sp5, &GPC->tmp_sp4, -1. );
	}

	//printf("WARNING!!!!\n");
	//matrix_vector_prod( &GPC->K, &GPC->tmp_sp1, &GPC->Us );
	
	if( nc > 0 ){
		//vector_null(&GPC->tmp_sp2);
		//matrix_vector_prod( &GPC->Wy, &GPC->tmp_sp4, &GPC->tmp_sp2);
		//matrixT_vector_prod( &GPC->C, &GPC->tmp_sp2, &GPC->tmp_sm);
		matrixT_vector_prod( &GPC->C, &GPC->tmp_sp4, &GPC->tmp_sm);
	} else {
		//vector_null(&GPC->tmp_sp4);
		//matrix_vector_prod( &GPC->Wy, &GPC->tmp_sp2, &GPC->tmp_sp4);
		//matrixT_vector_prod( &GPC->C, &GPC->tmp_sp4, &GPC->tmp_sm);
		matrixT_vector_prod( &GPC->C, &GPC->tmp_sp2, &GPC->tmp_sm);
	}
	vector_null(&GPC->Us);
	vector_null(&GPC->tmp_aa);
	LinearSystemV( &GPC->tmp_smsm3, &GPC->tmp_sm, &GPC->Us, &GPC->tmp_aa );

	return GPC_OK;

}

//typedef struct timeval time_eval;
GPC_ResT GPC_Control2( GPC_Model *GPC, vector *Us_ID ){

	//unsigned i;
	unsigned na, nb, nc, p, m, s;
	unsigned i, r, c, j ,k;
	/*
	FILE *fh;
	double time1, time2, time3, delta_time1, delta_time2;
	struct timeval t1, t2 ,t3;
	
	
	fh = fopen("Time1.txt","a");
	*/
	
	na = GPC->na;
	nb = GPC->nb;
	nc = GPC->nc;
	p = GPC->p;
	m = GPC->m;
	s = GPC->s;

	matrix_null(&GPC->tmp_pp1);
	matrix_null(&GPC->tmp_pp2);
	matrix_null(&GPC->tmp_pp3);
	matrix_null(&GPC->tmp_pp4);
	matrix_null(&GPC->tmp_pm1);
	matrix_null(&GPC->tmp_pm2);
	matrix_null(&GPC->tmp_pm3);
	matrix_null(&GPC->tmp_pm4);

	matrix_null(&GPC->tmp_sm1sm11);
	matrix_null(&GPC->tmp_sm1sm12);
	matrix_null(&GPC->tmp_sm1sm13);
	vector_null(&GPC->tmp_sp1);
	vector_null(&GPC->tmp_sp2);
	vector_null(&GPC->tmp_sp3);
	vector_null(&GPC->tmp_sm1);
	vector_null(&GPC->tmp_aa1);
	/*
	gettimeofday( &t1, 0 );
    	time1 = ( t1.tv_sec + t1.tv_usec*1e-6 );
	*/
	/* matrix_transpose_prod( &GPC->C2, &GPC->C2, &GPC->tmp_sm1sm11, 1.);*/
	/* eseguo l'operazione C2^T*C2 sfruttando il fatto che la matrice C2 ha una struttura
	  particolare, in questo modo rendo il codice più efficiente */
	for( i=0; i<s; i++){
		/* blocco 1 */
		for(r=0; r<p; r++){
			for(c=0; c<m; c++){
				GPC->tmp_pm1.mat[r][c] = GPC->C2.mat[r][c+(i+1)*m];
			}
		}
		for(j=i; j<s;j++){
			/* blocco 2 */
			for(r=0; r<p; r++){
				for(c=0; c<m; c++){
					GPC->tmp_pm2.mat[r][c] = GPC->C2.mat[r][c+(j+1)*m];
				}
			}
			matrix_transpose_prod( &GPC->tmp_pm1, &GPC->tmp_pm2, &GPC->tmp_mm, 1.);
			if (j!=i){
				for(k=j;k<s;k++){
					for(r=0; r<m; r++){
						for(c=0; c<m; c++){
							GPC->tmp_sm1sm11.mat[m*(k-j+1)+i*m+r][m*(k+1)+c] += GPC->tmp_mm.mat[r][c];
						}
					}
				}
				matrix_transpose_prod( &GPC->tmp_pm2, &GPC->tmp_pm1, &GPC->tmp_mm, 1.);
				for(k=j;k<s;k++){
					for(r=0; r<m; r++){
						for(c=0; c<m; c++){
							GPC->tmp_sm1sm11.mat[m*(k+1)+r][m*(k-j+1)+i*m+c] += GPC->tmp_mm.mat[r][c];
						}
					}
				}
			} else {
				for(k=j;k<s;k++){
					for(r=0; r<m; r++){
						for(c=0; c<m; c++){
							GPC->tmp_sm1sm11.mat[m*(k+1)+r][m*(k+1)+c] += GPC->tmp_mm.mat[r][c];
						}
					}
				}
			}
		}
	}
	/*
	gettimeofday( &t2, 0 );
    	time2 = ( t2.tv_sec + t2.tv_usec*1e-6 );
	*/
	matrix_eye( &GPC->tmp_sm1sm12, GPC->lambda );
	matrix_null(&GPC->tmp_sm1sm13);
	matrix_sum( &GPC->tmp_sm1sm11, &GPC->tmp_sm1sm12, &GPC->tmp_sm1sm13, 1. );
	
	matrix_vector_prod( &GPC->A, &GPC->Yp, &GPC->tmp_sp1);
	matrix_vector_prod( &GPC->B2, &GPC->Up2, &GPC->tmp_sp2);
	matrix_vector_prod( &GPC->C2, Us_ID, &GPC->tmp_sp3);
	if( nc > 0 ){
		matrix_vector_prod( &GPC->D, &GPC->Ep, &GPC->tmp_sp5);
	}
	vector_null(&GPC->tmp_sp4);
	vector_sum( &GPC->Ydes, &GPC->tmp_sp1, &GPC->tmp_sp4, -1. );
	vector_null(&GPC->tmp_sp1);
	vector_sum( &GPC->tmp_sp4, &GPC->tmp_sp2, &GPC->tmp_sp1, -1. );
	vector_null(&GPC->tmp_sp2);
	vector_sum( &GPC->tmp_sp1, &GPC->tmp_sp3, &GPC->tmp_sp2, -1. );
	if( nc > 0 ){
		vector_null(&GPC->tmp_sp4);
		vector_sum( &GPC->tmp_sp2, &GPC->tmp_sp5, &GPC->tmp_sp4, -1. );
	}

	//printf("WARNING!!!!\n");
	//matrix_vector_prod( &GPC->K, &GPC->tmp_sp1, &GPC->Us );

	if( nc > 0 ){
		matrixT_vector_prod( &GPC->C2, &GPC->tmp_sp4, &GPC->tmp_sm1);
	} else {
		matrixT_vector_prod( &GPC->C2, &GPC->tmp_sp2, &GPC->tmp_sm1);
	}
	vector_null(&GPC->Us2);
	vector_null(&GPC->tmp_aa1);

	LinearSystemV( &GPC->tmp_sm1sm13, &GPC->tmp_sm1, &GPC->Us2, &GPC->tmp_aa1 );
	/*
	gettimeofday( &t3, 0 );
    	time3 = ( t3.tv_sec + t3.tv_usec*1e-6 );
	delta_time1 = time2 - time1;
	delta_time2 = time3 - time2;
	fprintf(fh, "%le %le\n", delta_time1, delta_time2);
	fclose(fh);
	*/
	return GPC_OK;

}

int kalman_initialize( kalman_t *kalman, unsigned int n, unsigned int m, unsigned int p ){

	kalman->n_states = n;
	kalman->n_input = m;
	kalman->n_output = p;

	vector_init(&kalman->xp, n);
	vector_null(&kalman->xp);
	vector_init(&kalman->xm, n);
	vector_null(&kalman->xm);

	vector_init(&kalman->zp, p);
	vector_null(&kalman->zp);
	vector_init(&kalman->z, p);
	vector_null(&kalman->z);
	vector_init(&kalman->z_zp, p);
	vector_null(&kalman->z_zp);

	vector_init(&kalman->u, m);
	vector_null(&kalman->u);

	matrix_init(&kalman->Pm, n, n);
	matrix_null(&kalman->Pm);
	matrix_init(&kalman->Pp, n, n);
	matrix_null(&kalman->Pp);

	matrix_init(&kalman->K, n, p);
	matrix_null(&kalman->K);

	matrix_init(&kalman->Q, n, n);
	matrix_null(&kalman->Q);
	matrix_init(&kalman->R, p, p);
	matrix_null(&kalman->R);

	matrix_init(&kalman->PHI, n, n);
	matrix_null(&kalman->PHI);
	matrix_init(&kalman->H, p, n);
	matrix_null(&kalman->H);

	matrix_init(&kalman->HPH_R, p, p);
	matrix_null(&kalman->HPH_R);


	matrix_init(&kalman->mat_nn_1, n, n);
	matrix_null(&kalman->mat_nn_1);
	matrix_init(&kalman->mat_nn_2, n, n);
	matrix_null(&kalman->mat_nn_2);
	matrix_init(&kalman->mat_nn_3, n, n);
	matrix_null(&kalman->mat_nn_3);
	
	matrix_init(&kalman->mat_np, n, p);
	matrix_null(&kalman->mat_np);
	matrix_init(&kalman->mat_pn_1, p, n);
	matrix_null(&kalman->mat_pn_1);
	matrix_init(&kalman->mat_pn_2, p, n);
	matrix_null(&kalman->mat_pn_2);

	matrix_init(&kalman->mat_pp, p, p);
	matrix_null(&kalman->mat_pp);
	matrix_init(&kalman->mat_pp2, p, p);
	matrix_null(&kalman->mat_pp2);

	vector_init(&kalman->vec_p, p);
	vector_null(&kalman->vec_p);
	vector_init(&kalman->vec_n, n);
	vector_null(&kalman->vec_n);

	vector_init(&kalman->aa, p*(p+1)/2 );
	vector_null(&kalman->aa);
	vector_init(&kalman->bb, p );
	vector_null(&kalman->bb);

	return 0;
	
}

int kalman_initialize2( kalman_t *kalman, unsigned int n, unsigned int m, unsigned int p ){

	unsigned i;

	kalman->n_states = n;
	kalman->n_input = m;
	kalman->n_output = p;

	printf("AAAAA");
        if( !( kalman->XP = (vector *)calloc( p, sizeof(vector) ) ) ){
		fprintf( stdout, "KALMAN INITIALIZATION: MEMORY ERROR\n" );
                return 1;
        }
        if( !( kalman->XM = (vector *)calloc( p, sizeof(vector) ) ) ){
		fprintf( stdout, "KALMAN INITIALIZATION: MEMORY ERROR\n" );
                return 1;
        }
        if( !( kalman->PP = (matrix *)calloc( p, sizeof(matrix) ) ) ){
		fprintf( stdout, "KALMAN INITIALIZATION: MEMORY ERROR\n" );
                return 1;
        }
        if( !( kalman->PM = (matrix *)calloc( p, sizeof(matrix) ) ) ){
		fprintf( stdout, "KALMAN INITIALIZATION: MEMORY ERROR\n" );
                return 1;
        }
        if( !( kalman->QQ = (matrix *)calloc( p, sizeof(matrix) ) ) ){
		fprintf( stdout, "KALMAN INITIALIZATION: MEMORY ERROR\n" );
                return 1;
	}        
        if( !( kalman->RR = (matrix *)calloc( p, sizeof(matrix) ) ) ){
		fprintf( stdout, "KALMAN INITIALIZATION: MEMORY ERROR\n" );
                return 1;
        }
	for( i=0; i<p; i++){
		vector_init(&(kalman->XP[i]), n);
		vector_null(&(kalman->XP[i]));
		vector_init(&(kalman->XM[i]), n);
		vector_null(&(kalman->XM[i]));
		matrix_init(&(kalman->PP[i]), n, n);
		matrix_null(&(kalman->PP[i]));
		matrix_init(&(kalman->PM[i]), n, n);
		matrix_null(&(kalman->PM[i]));
		matrix_init(&(kalman->RR[i]), 1, 1);
		matrix_null(&(kalman->RR[i]));
		matrix_init(&(kalman->QQ[i]), n, n);
		matrix_null(&(kalman->QQ[i]));
	}		
	
	vector_init(&kalman->zp, 1);
	vector_null(&kalman->zp);
	vector_init(&kalman->Z, 1);
	vector_null(&kalman->Z);
	vector_init(&kalman->z, p);
	vector_null(&kalman->z);
	vector_init(&kalman->z_zp, 1);
	vector_null(&kalman->z_zp);

	vector_init(&kalman->u, m);
	vector_null(&kalman->u);

	matrix_init(&kalman->K, n, 1);
	matrix_null(&kalman->K);

	matrix_init(&kalman->PHI, 1, 1);
	matrix_null(&kalman->PHI);
	matrix_init(&kalman->H, 1, n);
	matrix_null(&kalman->H);

	matrix_init(&kalman->HPH_R, 1, 1);
	matrix_null(&kalman->HPH_R);

	matrix_init(&kalman->mat_nn_1, n, n);
	matrix_null(&kalman->mat_nn_1);
	matrix_init(&kalman->mat_nn_2, n, n);
	matrix_null(&kalman->mat_nn_2);
	matrix_init(&kalman->mat_nn_3, n, n);
	matrix_null(&kalman->mat_nn_3);
	
	matrix_init(&kalman->mat_np, n, 1);
	matrix_null(&kalman->mat_np);
	matrix_init(&kalman->mat_pn_1, 1, n);
	matrix_null(&kalman->mat_pn_1);
	matrix_init(&kalman->mat_pn_2, 1, n);
	matrix_null(&kalman->mat_pn_2);

	matrix_init(&kalman->mat_pp, 1, 1);
	matrix_null(&kalman->mat_pp);
	matrix_init(&kalman->mat_pp2, 1, 1);
	matrix_null(&kalman->mat_pp2);

	vector_init(&kalman->vec_p, 1);
	vector_null(&kalman->vec_p);
	vector_init(&kalman->vec_n, n);
	vector_null(&kalman->vec_n);

	vector_init(&kalman->aa, p*(p+1)/2 );
	vector_null(&kalman->aa);
	vector_init(&kalman->bb, p );
	vector_null(&kalman->bb);

	return 0;
	
}
int kalman_destroy( kalman_t *kalman ){

	vector_destroy(&kalman->xp);
	vector_destroy(&kalman->xm);

	vector_destroy(&kalman->zp);
	vector_destroy(&kalman->z);
	vector_destroy(&kalman->z_zp);

	vector_destroy(&kalman->u);

	matrix_destroy(&kalman->Pm);
	matrix_destroy(&kalman->Pp);

	matrix_destroy(&kalman->K);

	matrix_destroy(&kalman->Q);
	matrix_destroy(&kalman->R);

	matrix_destroy(&kalman->PHI);
	matrix_destroy(&kalman->H);

	matrix_destroy(&kalman->HPH_R);


	matrix_destroy(&kalman->mat_nn_1);
	matrix_destroy(&kalman->mat_nn_2);
	matrix_destroy(&kalman->mat_nn_3);
	
	matrix_destroy(&kalman->mat_np);
	matrix_destroy(&kalman->mat_pn_1);
	matrix_destroy(&kalman->mat_pn_2);

	matrix_destroy(&kalman->mat_pp);
	matrix_destroy(&kalman->mat_pp2);

	vector_destroy(&kalman->vec_p);
	vector_destroy(&kalman->vec_n);

	vector_destroy(&kalman->aa);
	vector_destroy(&kalman->bb);

	return 0;
	
}

int kalman_destroy2( kalman_t *kalman ){

	unsigned i;

	for( i=0; i<kalman->n_output; i++){
		vector_destroy(&(kalman->XP[i]));
		vector_destroy(&(kalman->XM[i]));
		matrix_destroy(&(kalman->PP[i]));
		matrix_destroy(&(kalman->PM[i]));
		matrix_destroy(&(kalman->QQ[i]));
		matrix_destroy(&(kalman->RR[i]));
	}
	free(kalman->XP);
	free(kalman->XM);
	free(kalman->PP);
	free(kalman->PM);
	free(kalman->QQ);
	free(kalman->RR);
		
	vector_destroy(&kalman->zp);
	vector_destroy(&kalman->z);
	vector_destroy(&kalman->Z);
	vector_destroy(&kalman->z_zp);

	vector_destroy(&kalman->u);

	matrix_destroy(&kalman->K);

	matrix_destroy(&kalman->R);

	matrix_destroy(&kalman->PHI);
	matrix_destroy(&kalman->H);

	matrix_destroy(&kalman->HPH_R);


	matrix_destroy(&kalman->mat_nn_1);
	matrix_destroy(&kalman->mat_nn_2);
	matrix_destroy(&kalman->mat_nn_3);
	
	matrix_destroy(&kalman->mat_np);
	matrix_destroy(&kalman->mat_pn_1);
	matrix_destroy(&kalman->mat_pn_2);

	matrix_destroy(&kalman->mat_pp);
	matrix_destroy(&kalman->mat_pp2);

	vector_destroy(&kalman->vec_p);
	vector_destroy(&kalman->vec_n);

	vector_destroy(&kalman->aa);
	vector_destroy(&kalman->bb);

	return 0;
	
}
int KalmanFilter( kalman_t *kalman, SysEq_t *data ){

	double det;
	FILE *fh;
	//double time1, time2, time3, delta_time1, delta_time2;
	//double time4, time5, time6, delta_time3, delta_time4;
	//double time7, time8, time9, delta_time5, delta_time6, delta_time7, delta_time8;
	//struct timeval t1, t2 ,t3, t4, t5, t6, t7, t8, t9;
	
	
	//fh = fopen("Time1.txt","a");

	/* computing the predicted state estimate
	 *      x(-)(k+1) = f( x(+)(k) )          */
	system_f( &kalman->xm, &kalman->xp, &kalman->u, data );

	/* cmputing the predicted measurment
	 *      z(k+1) = h( x(-)(k+1) )          */
	system_h( &kalman->xm, &kalman->zp, &kalman->u, data );
	//system_h( &kalman->xp, &kalman->zp, &kalman->u, data );

	/* computing the system linearization    */
	system_PHI( &kalman->PHI, &kalman->xm, &kalman->u, data );
	system_H( &kalman->H, &kalman->xm, &kalman->u, data );

	/* computing the a priori covariance matrix */
#ifdef CTRB
	//gettimeofday( &t1, 0 );
    	//time1 = ( t1.tv_sec + t1.tv_usec*1e-6 );
	matrix_prod_transpose( &kalman->Pp, &kalman->PHI, &kalman->mat_nn_1, 1. );
	//gettimeofday( &t2, 0 );
    	//time2 = ( t2.tv_sec + t2.tv_usec*1e-6 );
	matrix_prod( &kalman->PHI, &kalman->mat_nn_1, &kalman->mat_nn_2, 1. );
	//gettimeofday( &t3, 0 );
    	//time3 = ( t3.tv_sec + t3.tv_usec*1e-6 );
	matrix_sum( &kalman->mat_nn_2, &kalman->Q, &kalman->Pm, 1. );
	//gettimeofday( &t4, 0 );
    	//time4 = ( t4.tv_sec + t4.tv_usec*1e-6 );
#endif
#ifdef ARX_KALMAN
	matrix_sum( &kalman->Pp, &kalman->Q, &kalman->Pm, 1. );
#endif

	matrix_prod_transpose( &kalman->Pm, &kalman->H, &kalman->mat_np, 1. );
	matrix_prod( &kalman->H, &kalman->mat_np, &kalman->mat_pp, 1. );
	matrix_sum( &kalman->mat_pp, &kalman->R, &kalman->HPH_R, 1. );
	vector_sum( &kalman->z, &kalman->zp, &kalman->z_zp, -1. );
	if( kalman->adaptive == 1 ){
		matrix_vector_prod( &kalman->HPH_R, &kalman->z_zp, &kalman->vec_p );
		scalar_prod( &kalman->z_zp, &kalman->vec_p, &kalman->eps );
		if( kalman->eps > kalman->eps_max ){
			matrix_copy( &kalman->Pm, &kalman->Pm, 1+kalman->alpha_add );
		}
		else{
			matrix_copy( &kalman->Pm, &kalman->Pm, 1-kalman->alpha_sub );
		}
		matrix_prod_transpose( &kalman->Pm, &kalman->H, &kalman->mat_np, 1. );
		matrix_prod( &kalman->H, &kalman->mat_np, &kalman->mat_pp, 1. );
		matrix_sum( &kalman->mat_pp, &kalman->R, &kalman->HPH_R, 1. );
	}

	/* computing Kalman Gain */
	matrix_transpose( &kalman->mat_pp, &kalman->HPH_R );
	matrix_prod_transpose( &kalman->H, &kalman->Pm, &kalman->mat_pn_1, 1. );
	//matrix_write(&kalman->mat_pp, stdout, W_M_BIN);
	//getchar();
	LinearSystemM(&kalman->mat_pp, &kalman->mat_pn_1, &kalman->mat_pn_2, &kalman->aa, &kalman->bb);
	//det = kalman->mat_pp.mat[0][0]*kalman->mat_pp.mat[1][1] - kalman->mat_pp.mat[0][1]*kalman->mat_pp.mat[1][0];
	//kalman->mat_pp2.mat[0][0] = kalman->mat_pp.mat[1][1]/det;
	//kalman->mat_pp2.mat[0][1] = -kalman->mat_pp.mat[0][1]/det;
	//kalman->mat_pp2.mat[1][0] = -kalman->mat_pp.mat[1][0]/det;
	//kalman->mat_pp2.mat[1][1] = kalman->mat_pp.mat[0][0]/det;
	//matrix_prod(&kalman->mat_pp2, &kalman->mat_pn_1, &kalman->mat_pn_2, 1.);
	matrix_transpose( &kalman->K, &kalman->mat_pn_2 );
	fh = fopen("MAT.txt","w");
	matrix_write( &kalman->Pp, fh, W_M_BIN);
	fclose(fh);

	/* conditioning the predicted estimate on the measurement */
	matrix_vector_prod( &kalman->K, &kalman->z_zp, &kalman->vec_n );
	vector_sum( &kalman->xm, &kalman->vec_n, &kalman->xp, 1. );

	/* conditioning the a posteriori covariance matrix */
	matrix_eye( &kalman->mat_nn_1, 1. );
	matrix_prod( &kalman->K, &kalman->H, &kalman->mat_nn_2, 1. );
	matrix_sum( &kalman->mat_nn_1, &kalman->mat_nn_2, &kalman->mat_nn_3, -1. );
	//matrix_prod( &kalman->mat_nn_3, &kalman->Pm, &kalman->Pp, 1. );
	matrix_prod_sym( &kalman->mat_nn_3, &kalman->Pm, &kalman->Pp, 1. );

	//fh = fopen("MAT.txt","w");
	//matrix_write( &kalman->mat_nn_3, fh, W_M_BIN);
	//fclose(fh);
	//gettimeofday( &t9, 0 );
    	//time9 = ( t9.tv_sec + t9.tv_usec*1e-6 );
	//delta_time1 = time2 - time1;
	//delta_time2 = time3 - time2;
	//fprintf(fh, "%le %le\n", delta_time1, delta_time2);
	//fclose(fh);

	return 0;

}
int KalmanFilter2( kalman_t *kalman, SysEq_t *data ){

	unsigned i, j, k;
	unsigned n, m, p;
	
	n = data->n;
	m = data->m;
	p = data->p;

	for( i=0; i<p; i++){	

		/* computing the predicted state estimate
		 *      x(-)(k+1) = f( x(+)(k) )          */
		system_f( &(kalman->XM[i]), &(kalman->XP[i]), &kalman->u, data );

		/* cmputing the predicted measurment
		 *      z(k+1) = h( x(-)(k+1) )          */
		//system_h( &(kalman->XM[i]), &kalman->zp, &kalman->u, data );
		//system_h( &kalman->xp, &kalman->zp, &kalman->u, data );
		matrix_vector_prod(&data->C_ARX2, &(kalman->XM[i]), &kalman->zp);

		//system_H( &kalman->H, &kalman->xm, &kalman->u, data );
		matrix_copy(&kalman->H, &data->C_ARX2, 1.);

		matrix_sum( &(kalman->PP[i]), &(kalman->QQ[i]), &(kalman->PM[i]), 1. );

		matrix_prod_transpose( &(kalman->PM[i]), &kalman->H, &kalman->mat_np, 1. );
		matrix_prod( &kalman->H, &kalman->mat_np, &kalman->mat_pp, 1. );
		matrix_sum( &kalman->mat_pp, &(kalman->RR[i]), &kalman->HPH_R, 1. );
		kalman->Z.vec[0] = kalman->z.vec[i];
		
		vector_sum( &kalman->Z, &kalman->zp, &kalman->z_zp, -1. );
		//vector_write( &kalman->z_zp, stdout, W_M_BIN );

		/* computing Kalman Gain */
		matrix_prod_transpose( &kalman->H, &(kalman->PM[i]), &kalman->mat_pn_1, 1. );
		matrix_copy( &kalman->mat_pn_2, &kalman->mat_pn_1, 1./kalman->HPH_R.mat[0][0]);	
		//matrix_transpose( &kalman->mat_pp, &kalman->HPH_R );
		//LinearSystemM(&kalman->mat_pp, &kalman->mat_pn_1, &kalman->mat_pn_2, &kalman->aa, &kalman->bb);
		matrix_transpose( &kalman->K, &kalman->mat_pn_2 );
		//matrix_write( &kalman->K, stdout, W_M_BIN);

		/* conditioning the predicted estimate on the measurement */
		
		matrix_vector_prod( &kalman->K, &kalman->z_zp, &kalman->vec_n );
		vector_sum( &(kalman->XM[i]), &kalman->vec_n, &(kalman->XP[i]), 1. );

		/* conditioning the a posteriori covariance matrix */
		matrix_eye( &kalman->mat_nn_1, 1. );
		matrix_prod( &kalman->K, &kalman->H, &kalman->mat_nn_2, 1. );
		matrix_sum( &kalman->mat_nn_1, &kalman->mat_nn_2, &kalman->mat_nn_3, -1. );
		matrix_prod( &kalman->mat_nn_3, &(kalman->PM[i]), &(kalman->PP[i]), 1. );
		//matrix_prod_sym( &kalman->mat_nn_3, &(kalman->PM[i]), &(kalman->PP[i]), 1. );
	}
	//getchar();


	return 0;

}

#ifdef 	CTRB
int SysEq_initialize( SysEq_t *sys, unsigned int n, unsigned int m, unsigned int p, vector *lam ){

	unsigned int i;

	sys->n = n;
	sys->m = m;
	sys->p = p;
	sys->N = n + n*(m+p);

	vector_init(&sys->lambda, m);
	vector_null(&sys->lambda);

	vector_init(&sys->states, n);
	vector_null(&sys->states);
	vector_init(&sys->parameters, n*m+n*p);
	vector_null(&sys->parameters);

	vector_init(&sys->vec_n_1, n);
	vector_null(&sys->vec_n_1);
	vector_init(&sys->vec_n_2, n);
	vector_null(&sys->vec_n_2);
	vector_init(&sys->vec_p_1, p);
	vector_null(&sys->vec_p_1);
	vector_init(&sys->vec_p_2, p);
	vector_null(&sys->vec_p_2);

	matrix_init(&sys->A, n, n);
	matrix_null(&sys->A);
	matrix_init(&sys->B, n, m);
	matrix_null(&sys->B);
	matrix_init(&sys->C, p, n);
	matrix_null(&sys->C);
	matrix_init(&sys->D, p, m);
	matrix_null(&sys->D);

	vector_copy(&sys->lambda, lam, 1. );

	for( i=1; i<m; i++){
		sys->lambda.vec[i] = sys->lambda.vec[i-1] + sys->lambda.vec[i];
	}

	return 0;
	
}

int SysEq_destroy( SysEq_t *sys ){

	vector_destroy(&sys->lambda);

	vector_destroy(&sys->states);
	vector_destroy(&sys->parameters);

	vector_destroy(&sys->vec_n_1);
	vector_destroy(&sys->vec_n_2);
	vector_destroy(&sys->vec_p_1);
	vector_destroy(&sys->vec_p_2);

	matrix_destroy(&sys->A);
	matrix_destroy(&sys->B);
	matrix_destroy(&sys->C);
	matrix_destroy(&sys->D);

	return 0;
	
}
#endif

int CanonicalForm_ctrb( SysEq_t *data, vector *x ){

	unsigned int i, j, ind;

	matrix_null( &data->A );
	matrix_null( &data->B );
	matrix_null( &data->C );
	matrix_null( &data->D );

	ind = 0;
	data->B.mat[0][0] = 1.;
	for( i=0; i<data->n; i++){
	       //if( local_find(&data->lambda2, i)==1 ){
	       if( local_find(&data->lambda, i)==0 ){
	       		data->A.mat[i+1][i] = 1.;
	       }
	       else{
			if( ind < data->m-1 ){
				data->B.mat[i+1][ind+1] = 1.;
			}
			for( j=0; j<data->n; j++){
				data->A.mat[j][i] = x->vec[ ind*data->n + j];
			}
			ind++;
	       }
	}

	for( i=0; i<data->p; i++){
		for( j=0; j<data->n; j++){
			data->C.mat[i][j] = x->vec[ data->m*data->n + i*(data->n) + j ];
		}
	}
	/*
	vector_write(x, stdout, W_M_BIN);
	matrix_write(&data->A, stdout, W_M_BIN);
	matrix_write(&data->B, stdout, W_M_BIN);
	matrix_write(&data->C, stdout, W_M_BIN);
	matrix_write(&data->D, stdout, W_M_BIN);
	getchar();
	*/

	return 0;
}

int local_find( vector *VEC, unsigned int a ){

	unsigned int i;

	for( i=0; i<VEC->dimension; i++ ){
		if(VEC->vec[i] == a){
			return 1;
		}
	}

	return 0;
}

#ifdef CTRB
int system_f( vector *xm, vector *xp, vector *u, SysEq_t *data){

	unsigned int i;


	for( i=0; i<data->n; i++ ){
		data->states.vec[i] = xp->vec[i];
	}
	for( i=0; i<data->n*(data->p+data->m); i++ ){
		data->parameters.vec[i] = xp->vec[data->n+i];
	}

	CanonicalForm_ctrb( data, &data->parameters );

	matrix_vector_prod( &data->A, &data->states, &data->vec_n_1 );
	matrix_vector_prod( &data->B, u, &data->vec_n_2 );
	vector_sum( &data->vec_n_1, &data->vec_n_2, &data->states, 1. );

	for( i=0; i<data->n; i++ ){
		xm->vec[i] = data->states.vec[i];
	}
	for( i=0; i<data->n*(data->p+data->m); i++ ){
		xm->vec[i+data->n] = data->parameters.vec[i];
	}

	return 0;
}


int system_h( vector *x, vector *y, vector *u, SysEq_t *data){

	unsigned int i;

	for( i=0; i<data->n; i++ ){
		data->states.vec[i] = x->vec[i];
	}
	for( i=0; i<data->n*(data->p+data->m); i++ ){
		data->parameters.vec[i] = x->vec[data->n+i];
	}

	CanonicalForm_ctrb( data, &data->parameters );

	matrix_vector_prod( &data->C, &data->states, &data->vec_p_1 );
	matrix_vector_prod( &data->D, u, &data->vec_p_2 );
	vector_sum( &data->vec_p_1, &data->vec_p_2, y, 1. );

	return 0;
}

int system_PHI( matrix *PHI, vector *x, vector *u, SysEq_t *data){

	unsigned int i, j;

	for( i=0; i<data->n; i++ ){
		data->states.vec[i] = x->vec[i];
	}
	for( i=0; i<data->n*(data->p+data->m); i++ ){
		data->parameters.vec[i] = x->vec[data->n+i];
	}

	CanonicalForm_ctrb( data, &data->parameters );

	matrix_null(PHI);
	for( i=0; i<data->m; i++ ){
		for( j=0; j<data->n; j++ ){
			PHI->mat[j][data->n+j+i*data->n] = x->vec[(int)data->lambda.vec[i]];
		}
	}

	for( i=0; i<data->n; i++ ){
		for( j=0; j<data->n; j++ ){
			PHI->mat[i][j] = data->A.mat[i][j];
		}
	}

	for( i=0; i<data->n*(data->m+data->p); i++ ){
		PHI->mat[data->n+i][data->n+i] = 1.;
	}

	return 0;
}


int system_H( matrix *H, vector *x, vector *u, SysEq_t *data){

	unsigned int i, j;

	for( i=0; i<data->n; i++ ){
		data->states.vec[i] = x->vec[i];
	}
	for( i=0; i<data->n*(data->p+data->m); i++ ){
		data->parameters.vec[i] = x->vec[data->n+i];
	}

	CanonicalForm_ctrb( data, &data->parameters );

	matrix_null(H);
	for( i=0; i<data->n; i++ ){
		for( j=0; j<data->p; j++ ){
			H->mat[j][data->n + data->n*data->m + j + i*data->p] = x->vec[i];
		}
	}

	for( i=0; i<data->p; i++ ){
		for( j=0; j<data->n; j++ ){
			H->mat[i][j] = data->C.mat[i][j];
		}
	}

	return 0;
}
#endif

#ifdef ARX_KALMAN

int SysEq_initialize( SysEq_t *sys, unsigned int n, unsigned int m, unsigned int p, vector *lam ){

	sys->n = n;
	sys->m = m;
	sys->p = p;
	sys->N = n*p*p + n*p*m;

	matrix_init(&sys->Cy, p, p*p*(n-1));
	matrix_null(&sys->Cy);
	matrix_init(&sys->Cu, p, p*m*(n-1));
	matrix_null(&sys->Cu);
	matrix_init(&sys->Cy2, 1, p*(n-1));
	matrix_null(&sys->Cy2);
	matrix_init(&sys->Cu2, 1, m*(n-1));
	matrix_null(&sys->Cu2);

	matrix_init(&sys->C_ARX, p, sys->N);
	matrix_null(&sys->C_ARX);
	matrix_init(&sys->C_ARX2, 1, n*p+m*n);
	matrix_null(&sys->C_ARX2);


	return 0;
	
}

int SysEq_destroy( SysEq_t *sys ){

	matrix_destroy(&sys->Cy);
	matrix_destroy(&sys->Cu);

	matrix_destroy(&sys->C_ARX);

	return 0;
	
}


int system_f( vector *xm, vector *xp, vector *u, SysEq_t *data){

	vector_copy(xm, xp, 1.);

	return 0;
}

int system_h( vector *x, vector *y, vector *u, SysEq_t *data){

	matrix_vector_prod(&data->C_ARX, x, y);

	return 0;
}
	
int system_PHI( matrix *PHI, vector *x, vector *u, SysEq_t *data){

	matrix_eye(PHI, 1.);

	return 0;
}

int system_H( matrix *H, vector *x, vector *u, SysEq_t *data){

	matrix_copy(H, &data->C_ARX, 1.);
	
	return 0;
}
#endif


int UpdateC_ARX( SysEq_t *data  ,vector *y, vector *u){

	unsigned n, p, m, i, j;

	n = data->n;
	p = data->p;
	m = data->m;

	sub_matrix_extract( &data->C_ARX, &data->Cy, 0, 0);	
	sub_matrix_extract( &data->C_ARX, &data->Cu, 0, p*p*n);	

	sub_matrix_insert( &data->C_ARX, &data->Cy, 0, p*p );
	sub_matrix_insert( &data->C_ARX, &data->Cu, 0, p*p*n+p*m );

	for( i=0; i<p; i++){
		for( j=0; j<p; j++){
			data->C_ARX.mat[i][i*p+j] = y->vec[j];
		}
		for( j=0; j<m; j++){
			data->C_ARX.mat[i][n*p*p+i*m+j] = u->vec[j];
		}
	}

	return 0;
}

int UpdateC_ARX2( SysEq_t *data  ,vector *y, vector *u){

	unsigned n, p, m, i, j;

	n = data->n;
	p = data->p;
	m = data->m;

	sub_matrix_extract( &data->C_ARX2, &data->Cy2, 0, 0);	
	sub_matrix_extract( &data->C_ARX2, &data->Cu2, 0, p*n);	

	sub_matrix_insert( &data->C_ARX2, &data->Cy2, 0, p );
	sub_matrix_insert( &data->C_ARX2, &data->Cu2, 0, p*n+m );

	for( j=0; j<p; j++){
		data->C_ARX2.mat[0][j] = y->vec[j];
	}
	for( j=0; j<m; j++){
		data->C_ARX2.mat[0][n*p+j] = u->vec[j];
	}

	return 0;
}

int x2theta( vector * x, matrix * theta, unsigned n, unsigned p, unsigned m){

	unsigned i, j, k;

	for( k=0; k<n; k++){
		for( i=0; i<p; i++){
			for( j=0; j<p; j++){
				theta->mat[i][j+k*p] = x->vec[j+i*p+k*p*p];
			}
			for( j=0; j<m; j++){
				theta->mat[i][p*n+j+k*m] = x->vec[p*p*n+j+i*m+k*p*m];
			}
		}
	}
	//matrix_write(theta, stdout, W_M_BIN);

	return 0;
}
#if 0		
int X2theta( vector **X, matrix * theta, vector * v, unsigned n, unsigned p, unsigned m){

	unsigned i, j, k;

	for( i=0; i<p; i++){
		vector_copy( v, X[i], 1.);
		for( j=0; j<p*n+m*n; j++){
			//theta->mat[i][j] = X[i]->vec[j];
			theta->mat[i][j] = v->vec[j];
		}
	}
	matrix_write(theta, stdout, W_M_BIN);
	getchar();

	return 0;
}
#endif
#if 0
GPC_ResT GPC( ARX_Model *ARX, GPC_Model *GPC, vector *y_k, vector *y_k_1, vector *u, vector *u_kp1, unsigned FlagID, unsigned FlagCON){

	unsigned i;

	ARX_UpdatePhi( ARX, y_k_1, u);

	if( FlagID == 1 ){
		ARX_RLS( ARX, y_k);
		//matrix_write(&ARX->theta, stdout, W_M_BIN);
		//vector_write(&ARX->phi, stdout, W_M_BIN);
	}
	
	vector_null(u_kp1);
	if( FlagCON == 1 ){
		GPC_PredictionFunction( GPC, &ARX->theta, ARX->SimplyProper, 0 );
		//GPC_Control( GPC, &ARX->phi, ARX->SimplyProper );
		GPC_Control( GPC,  );
		for( i=0; i<GPC->m; i++) {
			u_kp1->vec[i] = GPC->Us.vec[GPC->m*(GPC->s-2)+i];
		}
	}

	return GPC_OK;

}
		
#endif		
	
	
#if 0
int main (){

	ARX_Model ARX;
	GPC_Model GPC;
	vector y, u;
	matrix theta, A, A0,b, x;
	vector aa, bb;
	unsigned i;
	unsigned na, nb, m, p, flag;
	FILE *fh, *fh_A, *fh_b, *fh_x;

	fh = fopen("theta.txt","w");
	fh_A = fopen("A.txt","w");
	fh_b = fopen("b.txt","w");
	fh_x = fopen("x.txt","w");

	na = 3;
	nb = 4;
	m = 3;
	p = 2;
	flag = 1;

	ARX_Initialize( &ARX, na, nb, m, p, 0.1, 10., flag);
	GPC_Initialize( &GPC, na, nb, m, p, 3, 1.);
	vector_init( &y, p);
	vector_init( &u, m);
	matrix_init( &theta, p, p*na + m*nb + m*flag);
	matrix_random( &theta, -1, 1 );
	matrix_write( &theta, stdout, W_M_BIN);
	matrix_write( &theta, fh, W_M_BIN);

	matrix_init( &A0, 4,4);
	matrix_random( &A0, -1, 1 );
	matrix_init( &A, 4,4);
	matrix_prod_transpose( &A0, &A0, &A, 1.);
	matrix_write( &A, fh_A, W_M_BIN);
	matrix_init( &b, 4,3);
	matrix_random( &b, -1, 1 );
	matrix_write( &b, fh_b, W_M_BIN);
	matrix_init( &x, 4,3);
	vector_init( &aa, 4*5/2);
	vector_init( &bb, 4);
	LinearSystemM( &A, &b, &x, &aa, &bb);
	matrix_write( &x, fh_x, W_M_BIN);

	GPC_PredictionFunction( &GPC, &theta, flag );
	matrix_write( &GPC.A, stdout, W_M_BIN);
	matrix_write( &GPC.B, stdout, W_M_BIN);
	matrix_write( &GPC.C, stdout, W_M_BIN);

	for( i=0; i<0; i++ ){
		vector_random( &y, -1., 1.);
		//vector_write( &y, stdout, W_M_BIN);
		vector_random( &u, -1., 1.);
		//vector_write( &u, stdout, W_M_BIN);

		ARX_UpdatePhi( &ARX, &y, &u);
		ARX_RLS( &ARX, &y );
		matrix_write( &ARX.P, stdout, W_M_BIN);
		matrix_write( &ARX.theta, stdout, W_M_BIN);
		//vector_write( &ARX.phi, stdout, W_M_BIN);

	}

	fclose(fh);
	fclose(fh_A);
	fclose(fh_b);
	fclose(fh_x);
	vector_destroy( &y );
	vector_destroy( &u );
	ARX_Destroy( &ARX );
	GPC_Destroy( &GPC );
	matrix_destroy(&theta);
	matrix_destroy(&A);
	matrix_destroy(&A0);
	matrix_destroy(&b);
	matrix_destroy(&x);
	vector_destroy(&aa);
	vector_destroy(&bb);
		
	return 0;

}

#endif	
