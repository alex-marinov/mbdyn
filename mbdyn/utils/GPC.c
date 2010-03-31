/*
 * Copyright (C) 2008
 *
 * Mattia Mattaboni	<mattaboni@aero.polimi.it>
 */

#include "GPC.h"
#include "udu.h"

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
	
	matrix_init( &ARMAX->theta, p, dimension );
	
	vector_init( &ARMAX->phi, dimension );


	vector_init( &ARMAX->yp, p );
	vector_init( &ARMAX->eps, p );
	vector_init( &ARMAX->eps2, p );
	vector_init( &ARMAX->Pphi, dimension );
	vector_init( &ARMAX->K, dimension );
	matrix_init( &ARMAX->DeltaTheta, p, dimension);
	matrix_init( &ARMAX->DeltaP, dimension, dimension );

	return GPC_OK;

}

GPC_ResT ARMAX_Destroy( ARMAX_Model *ARMAX) {
	
	matrix_destroy(&ARMAX->P);
	matrix_destroy(&ARMAX->theta);

	vector_destroy(&ARMAX->phi);

	vector_destroy(&ARMAX->yp);
	vector_destroy(&ARMAX->eps);
	vector_destroy(&ARMAX->eps2);
	vector_destroy(&ARMAX->Pphi);
	vector_destroy(&ARMAX->K);
	matrix_destroy(&ARMAX->DeltaTheta);
	matrix_destroy(&ARMAX->DeltaP);

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

	tr = matrix_trace( &ARX->P );
	//printf("traccia P %e\n", tr/ARX->P.Nrow);

	return GPC_OK;

}

GPC_ResT ARMAX_RLS( ARMAX_Model *ARMAX, vector *y){

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
	matrix_sum( &ARMAX->P, &ARMAX->DeltaP, &ARMAX->P, -1. );
	matrix_copy( &ARMAX->P, &ARMAX->P, 1./ARMAX->mu );

	return GPC_OK;

}

GPC_ResT GPC_Initialize( GPC_Model *GPC, unsigned na, unsigned nb,
				unsigned m, unsigned p, unsigned s,
				double lambda ){

	GPC->na = na;
	GPC->nb = nb;
	GPC->m = m;
	GPC->p = p;
	GPC->s = s;
	GPC->lambda = lambda;

	matrix_init( &GPC->A, s*p, p*na );
	matrix_init( &GPC->B, s*p, m*nb );
	matrix_init( &GPC->C, s*p, s*m );
	matrix_init( &GPC->B2, s*p, m*(nb-1) );
	matrix_init( &GPC->C2, s*p, m*(s+1) );

	vector_init( &GPC->Yp, p*na );
	vector_init( &GPC->Up, m*nb );
	vector_init( &GPC->Up2, m*(nb-1) );
	vector_init( &GPC->Us, s*m );
	vector_init( &GPC->Us2, m*(s+1) );
	vector_init( &GPC->Ydes, s*p );

	matrix_init( &GPC->tmp_pp1, p, p);
	matrix_init( &GPC->tmp_pp2, p, p);
	matrix_init( &GPC->tmp_pp3, p, p);
	matrix_init( &GPC->tmp_pp4, p, p);
	matrix_init( &GPC->tmp_pm1, p, m);
	matrix_init( &GPC->tmp_pm2, p, m);
	matrix_init( &GPC->tmp_pm3, p, m);
	matrix_init( &GPC->tmp_pm4, p, m);

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
	matrix_destroy( &GPC->B2 );
	matrix_destroy( &GPC->C2 );

	vector_destroy( &GPC->Yp );
	vector_destroy( &GPC->Up );
	vector_destroy( &GPC->Us );
	vector_destroy( &GPC->Up2 );
	vector_destroy( &GPC->Us2 );
	vector_destroy( &GPC->Ydes );

	matrix_destroy(&GPC->tmp_pp1);
	matrix_destroy(&GPC->tmp_pp2);
	matrix_destroy(&GPC->tmp_pp3);
	matrix_destroy(&GPC->tmp_pp4);
	matrix_destroy(&GPC->tmp_pm1);
	matrix_destroy(&GPC->tmp_pm2);
	matrix_destroy(&GPC->tmp_pm3);
	matrix_destroy(&GPC->tmp_pm4);

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
	vector_destroy(&GPC->tmp_sm);
	vector_destroy(&GPC->tmp_sm1);
	vector_destroy(&GPC->tmp_aa);
	vector_destroy(&GPC->tmp_aa1);

	matrix_destroy(&GPC->K);

	return GPC_OK;
}

GPC_ResT GPC_PredictionFunction( GPC_Model *GPC, matrix *theta, unsigned FlagSimplyProper, unsigned Flag ){

	unsigned i, j;
	unsigned na, nb, p, m, s;

	na = GPC->na;
	nb = GPC->nb;
	p = GPC->p;
	m = GPC->m;
	s = GPC->s;

	for( i=0; i<s; i++ ){
		for( j=0; j<na; j++ ){
			if( i==0 ){
				sub_matrix_extract( theta, &GPC->tmp_pp1, 0, p*j);
				sub_matrix_insert( &GPC->A, &GPC->tmp_pp1, p*(s-1) , p*j);	
	
				sub_matrix_extract( theta, &GPC->tmp_pm1, 0, p*na + FlagSimplyProper*m + m*j);
				sub_matrix_insert( &GPC->B, &GPC->tmp_pm1, p*(s-1) , m*j);
			} else {
				if( j<(na-1) ){
					sub_matrix_extract( &GPC->A, &GPC->tmp_pp1, p*(s-i), 0 );
					sub_matrix_extract( &GPC->A, &GPC->tmp_pp2, p*(s-1), p*j );
					sub_matrix_extract( &GPC->A, &GPC->tmp_pp3, p*(s-i), p*(j+1) );
					matrix_prod( &GPC->tmp_pp1, &GPC->tmp_pp2, &GPC->tmp_pp4, 1. );
					matrix_sum( &GPC->tmp_pp4, &GPC->tmp_pp3, &GPC->tmp_pp2, 1. );
					sub_matrix_insert( &GPC->A, &GPC->tmp_pp2, p*(s-i-1), p*j );

					sub_matrix_extract( &GPC->B, &GPC->tmp_pm2, p*(s-1), m*j );
					sub_matrix_extract( &GPC->B, &GPC->tmp_pm3, p*(s-i), m*(j+1) );
					matrix_prod( &GPC->tmp_pp1, &GPC->tmp_pm2, &GPC->tmp_pm4, 1. );
					matrix_sum( &GPC->tmp_pm4, &GPC->tmp_pm3, &GPC->tmp_pm1, 1. );
					sub_matrix_insert( &GPC->B, &GPC->tmp_pm1, p*(s-i-1), m*j );
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
	

GPC_ResT GPC_VectorUpdate( GPC_Model *GPC, vector *y_k, vector *u_k ){

	unsigned i;
	unsigned na, nb, p, m, s;

	na = GPC->na;
	nb = GPC->nb;
	p = GPC->p;
	m = GPC->m;
	s = GPC->s;
	
	for( i=0; i<p*(na-1); i++ ){
		GPC->Yp.vec[p*na-1-i] = GPC->Yp.vec[p*(na-1)-1-i];
	}
	for( i=0; i<m*(nb-1); i++ ){
		GPC->Up.vec[m*nb-1-i] = GPC->Up.vec[m*(nb-1)-1-i];
	}
	
	for( i=0; i<p; i++ ){
		GPC->Yp.vec[i] = y_k->vec[i];
	}
	
	for( i=0; i<m; i++ ){
		GPC->Up.vec[i] = u_k->vec[i];
	}

	return GPC_OK;	

}

GPC_ResT GPC_VectorUpdate2( GPC_Model *GPC, vector *y_k, vector *u_k ){

	unsigned i;
	unsigned na, nb, p, m, s;

	na = GPC->na;
	nb = GPC->nb;
	p = GPC->p;
	m = GPC->m;
	s = GPC->s;
	
	for( i=0; i<p*(na-1); i++ ){
		GPC->Yp.vec[p*na-1-i] = GPC->Yp.vec[p*(na-1)-1-i];
	}
	for( i=0; i<m*(nb-2); i++ ){
		GPC->Up2.vec[m*(nb-1)-1-i] = GPC->Up2.vec[m*(nb-2)-1-i];
	}
	
	for( i=0; i<p; i++ ){
		GPC->Yp.vec[i] = y_k->vec[i];
	}
	
	for( i=0; i<m; i++ ){
		GPC->Up2.vec[i] = u_k->vec[i];
	}

	return GPC_OK;	

}

GPC_ResT GPC_Control( GPC_Model *GPC, vector *Us_ID ){

	//unsigned i;
	unsigned na, nb, p, m, s;
	//FILE *fh;

	na = GPC->na;
	nb = GPC->nb;
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

	matrix_transpose_prod( &GPC->C, &GPC->C, &GPC->tmp_smsm1, 1.);
	matrix_eye( &GPC->tmp_smsm2, GPC->lambda );
	matrix_null(&GPC->tmp_smsm3);
	matrix_sum( &GPC->tmp_smsm1, &GPC->tmp_smsm2, &GPC->tmp_smsm3, 1. );

	matrix_vector_prod( &GPC->A, &GPC->Yp, &GPC->tmp_sp1);
	matrix_vector_prod( &GPC->B, &GPC->Up, &GPC->tmp_sp2);
	matrix_vector_prod( &GPC->C, Us_ID, &GPC->tmp_sp3);
	vector_null(&GPC->tmp_sp4);
	vector_sum( &GPC->Ydes, &GPC->tmp_sp1, &GPC->tmp_sp4, -1. );
	vector_null(&GPC->tmp_sp1);
	vector_sum( &GPC->tmp_sp4, &GPC->tmp_sp2, &GPC->tmp_sp1, -1. );
	vector_null(&GPC->tmp_sp2);
	vector_sum( &GPC->tmp_sp1, &GPC->tmp_sp3, &GPC->tmp_sp2, -1. );

	//printf("WARNING!!!!\n");
	//matrix_vector_prod( &GPC->K, &GPC->tmp_sp1, &GPC->Us );

	matrixT_vector_prod( &GPC->C, &GPC->tmp_sp2, &GPC->tmp_sm);
	vector_null(&GPC->Us);
	vector_null(&GPC->tmp_aa);
	LinearSystemV( &GPC->tmp_smsm3, &GPC->tmp_sm, &GPC->Us, &GPC->tmp_aa );

	return GPC_OK;

}

GPC_ResT GPC_Control2( GPC_Model *GPC, vector *Us_ID ){

	//unsigned i;
	unsigned na, nb, p, m, s;
	//FILE *fh;

	na = GPC->na;
	nb = GPC->nb;
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

	matrix_transpose_prod( &GPC->C2, &GPC->C2, &GPC->tmp_sm1sm11, 1.);
	matrix_eye( &GPC->tmp_sm1sm12, GPC->lambda );
	matrix_null(&GPC->tmp_sm1sm13);
	matrix_sum( &GPC->tmp_sm1sm11, &GPC->tmp_sm1sm12, &GPC->tmp_sm1sm13, 1. );
	
	matrix_vector_prod( &GPC->A, &GPC->Yp, &GPC->tmp_sp1);
	matrix_vector_prod( &GPC->B2, &GPC->Up2, &GPC->tmp_sp2);
	matrix_vector_prod( &GPC->C2, Us_ID, &GPC->tmp_sp3);
	vector_null(&GPC->tmp_sp4);
	vector_sum( &GPC->Ydes, &GPC->tmp_sp1, &GPC->tmp_sp4, -1. );
	vector_null(&GPC->tmp_sp1);
	vector_sum( &GPC->tmp_sp4, &GPC->tmp_sp2, &GPC->tmp_sp1, -1. );
	vector_null(&GPC->tmp_sp2);
	vector_sum( &GPC->tmp_sp1, &GPC->tmp_sp3, &GPC->tmp_sp2, -1. );

	//printf("WARNING!!!!\n");
	//matrix_vector_prod( &GPC->K, &GPC->tmp_sp1, &GPC->Us );

	matrixT_vector_prod( &GPC->C2, &GPC->tmp_sp2, &GPC->tmp_sm1);
	vector_null(&GPC->Us2);
	vector_null(&GPC->tmp_aa1);
	LinearSystemV( &GPC->tmp_sm1sm13, &GPC->tmp_sm1, &GPC->Us2, &GPC->tmp_aa1 );

	return GPC_OK;

}

	
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
