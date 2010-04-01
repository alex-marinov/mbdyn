/*
 * Copyright (C) 2008
 *
 * Mattia Mattaboni	<mattaboni@aero.polimi.it>
 */

#ifndef GPC_H
#define GPC_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <matrix.h>

/* codici errore */
typedef enum {
        GPC_OK = 0,
	GPC_GEN_ERROR
} GPC_ResT;

/* ARX MODEL
        na               nb            
y(k) =  sum Ai*y(k-i) +  sum Bi*u(k-i)
        i=1              i=0          
*/
typedef struct ARX_Model {
	unsigned na;	/* model order */
	unsigned nb;	/* model order */
	unsigned m;	/* input number */
	unsigned p;	/* output number */
	double mu;	/* forgetting factor */
	double delta;	/* covariance matrix initialization: P = deltaI */
	unsigned SimplyProper;
	vector phi;
	matrix theta;
	matrix P;
	/* workign vectors and matrices: */
	vector yp, eps, eps2, Pphi, K;
	matrix DeltaTheta, DeltaP;
} ARX_Model;

/* ARMAX MODEL
        na               nb              nc         
y(k) =  sum Ai*y(k-i) +  sum Bi*u(k-i) + sum Ci*eps(k-i)
        i=1              i=0             i=1
*/
typedef struct ARMAX_Model {
	unsigned na;	/* model order */
	unsigned nb;	/* model order */
	unsigned nc;	/* model order */
	unsigned m;	/* input number */
	unsigned p;	/* output number */
	double mu;	/* forgetting factor */
	double delta;	/* covariance matrix initialization: P = deltaI */
	unsigned SimplyProper;
	vector phi;
	matrix theta;
	matrix P;
	/* workign vectors and matrices: */
	vector yp, eps, eps2, Pphi, K;
	matrix DeltaTheta, DeltaP;
} ARMAX_Model;


GPC_ResT ARX_Initialize( ARX_Model *ARX, unsigned na, unsigned nb, 
		    unsigned m, unsigned p, double mu,
		    double delta, unsigned FlagProper);

GPC_ResT ARX_Destroy( ARX_Model *ARX ); 

GPC_ResT ARMAX_Initialize( ARMAX_Model *ARMAX, unsigned na, unsigned nb, 
		    unsigned nc, unsigned m, unsigned p, double mu,
		    double delta, unsigned FlagProper);

GPC_ResT ARMAX_Destroy( ARMAX_Model *ARMAX ); 

GPC_ResT ARX_UpdatePhi( ARX_Model *ARX, vector *y, vector *u);

GPC_ResT ARMAX_UpdatePhi( ARMAX_Model *ARMAX, vector *y, vector *u, vector *eps);

GPC_ResT ARX_RLS( ARX_Model *ARX, vector *y);	

GPC_ResT ARMAX_RLS( ARMAX_Model *ARMAX, vector *y);	

typedef struct GPC_Model {
	unsigned na;	/* model order */
	unsigned nb;	/* model order */
	unsigned m;	/* input number */
	unsigned p;	/* output number */
	unsigned s;
	double lambda;
	matrix A, B, C, B2, C2;
	vector Yp, Up, Us, Ydes, Up2, Us2;
	matrix tmp_pp1, tmp_pp2, tmp_pp3, tmp_pp4;
	matrix tmp_pm1, tmp_pm2, tmp_pm3, tmp_pm4;
	matrix tmp_smsm1, tmp_smsm2, tmp_smsm3;
	matrix tmp_sm1sm11, tmp_sm1sm12, tmp_sm1sm13;
	matrix tmp_sp_m;
	vector tmp_sp1, tmp_sp2, tmp_sp3, tmp_sp4, tmp_sm, tmp_aa, tmp_sm1, tmp_aa1;
	matrix K;
} GPC_Model;

GPC_ResT GPC_Initialize( GPC_Model *GPC, unsigned na, unsigned nb,
				unsigned m, unsigned p, unsigned s, double lambda );

GPC_ResT GPC_Destroy( GPC_Model *GPC );

GPC_ResT GPC_PredictionFunction( GPC_Model *GPC, matrix *theta, unsigned FlagSimplyProper, unsigned FlagK );

GPC_ResT GPC_VectorUpdate( GPC_Model *GPC, vector *y_k, vector *u_k );
GPC_ResT GPC_VectorUpdate2( GPC_Model *GPC, vector *y_k, vector *u_k );

GPC_ResT GPC_Control( GPC_Model *GPC, vector *Us_ID  );
GPC_ResT GPC_Control2( GPC_Model *GPC, vector *Us_ID  );
//GPC_ResT GPC_Control( GPC_Model *GPC, vector *phi, unsigned FlagSimplyProper );

GPC_ResT GPC( ARX_Model *ARX, GPC_Model *GPC, vector *y_k, vector *y_k_1, vector *u, vector *u_kp1, unsigned FlagID, unsigned FlagCON);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* GPC_H */

