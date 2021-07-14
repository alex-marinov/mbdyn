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

#include "matrix.h"

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
	matrix psi;
	matrix theta;
	vector theta_RML;
	matrix P, P_RML;
	matrix **ALPHA, **BETA, **DELTA;
	/* workign vectors and matrices: */
	vector yp, eps, eps2, Pphi, K;
	matrix DeltaTheta, DeltaP;
	matrix PpsiT, psiPT, psiPpsiT, mu_eye, tmp_pp, tmp_p_pdim, deltaP_RML;
	vector aa, bb, delta_theta_RML;
 
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

GPC_ResT ARMAX_ELS( ARMAX_Model *ARMAX, vector *y);	
GPC_ResT ARMAX_RML( ARMAX_Model *ARMAX, vector *y);	

typedef struct GPC_Model {
	unsigned na;	/* model order */
	unsigned nb;	/* model order */
	unsigned nc;	/* model order */
	unsigned m;	/* input number */
	unsigned p;	/* output number */
	unsigned s;
	double lambda;
	matrix A, B, C, B2, C2, D;
	vector Yp, Up, Us, Ydes, Up2, Us2, Ep;
	matrix Wy, Wu, WyC;
	matrix tmp_pp1, tmp_pp2, tmp_pp3, tmp_pp4;
	matrix tmp_pm1, tmp_pm2, tmp_pm3, tmp_pm4;
	matrix tmp_smsm1, tmp_smsm2, tmp_smsm3;
	matrix tmp_sm1sm11, tmp_sm1sm12, tmp_sm1sm13;
	matrix tmp_sp_m, tmp_mm;
	vector tmp_sp1, tmp_sp2, tmp_sp3, tmp_sp4, tmp_sp5, tmp_sm, tmp_aa, tmp_sm1, tmp_aa1;
	matrix K;
} GPC_Model;

GPC_ResT GPC_Initialize( GPC_Model *GPC, unsigned na, unsigned nb, unsigned nc,
				unsigned m, unsigned p, unsigned s, double lambda );

GPC_ResT GPC_Destroy( GPC_Model *GPC );

GPC_ResT GPC_PredictionFunction( GPC_Model *GPC, matrix *theta, unsigned FlagSimplyProper, unsigned FlagK );

GPC_ResT GPC_VectorUpdate( GPC_Model *GPC, vector *y_k, vector *u_k, vector *e_k );
GPC_ResT GPC_VectorUpdate2( GPC_Model *GPC, vector *y_k, vector *u_k, vector *e_k );

GPC_ResT GPC_Control( GPC_Model *GPC, vector *Us_ID  );
GPC_ResT GPC_Control2( GPC_Model *GPC, vector *Us_ID  );
GPC_ResT GPC_ControlW( GPC_Model *GPC, vector *Us_ID  );
//GPC_ResT GPC_Control( GPC_Model *GPC, vector *phi, unsigned FlagSimplyProper );

GPC_ResT GPC( ARX_Model *ARX, GPC_Model *GPC, vector *y_k, vector *y_k_1, vector *u, vector *u_kp1, unsigned FlagID, unsigned FlagCON);

typedef struct kalman_t {

	/* kalman filter states number */
	unsigned int n_states;

	/* kalman filter deterministic input number */
	unsigned int n_input;

	/* kalman filter output number */
	unsigned int n_output;

	vector xp, xm, u, zp, z, z_zp;
	vector *XP, *XM, Z;
	matrix PHI, H, Pp, Pm, K, Q, R, HPH_R;
	matrix *PP, *PM, *QQ, *RR;
	
	unsigned int adaptive;
	double eps, eps_max, alpha_add, alpha_sub;

	matrix mat_nn_1, mat_nn_2, mat_nn_3, mat_np, mat_pn_1, mat_pn_2, mat_pp, mat_pp2;
	vector vec_n, vec_p, aa, bb;

}kalman_t;

typedef struct SysEq_t {

	/* states number */
	unsigned int n;

	/* deterministic input number */
	unsigned int m;

	/* output number */
	unsigned int p;

	/* states total number */
	unsigned int N;

	vector lambda;
	matrix A, B, C, D;
	vector parameters, states;
	
	matrix C_ARX, Cy, Cu, C_ARX2, Cy2, Cu2;

	vector vec_n_1, vec_n_2, vec_p_1, vec_p_2;
	
}SysEq_t;

int SysEq_initialize( SysEq_t *, unsigned int, unsigned int, unsigned int, vector *);
int SysEq_destroy( SysEq_t *);
int CanonicalForm_ctrb( SysEq_t *, vector * );
int local_find( vector *, unsigned int );
int system_f( vector *, vector *, vector *, SysEq_t * );
int system_h( vector *, vector *, vector *, SysEq_t * );
int system_PHI( matrix *, vector *, vector *, SysEq_t * );
int system_H( matrix *, vector *, vector *, SysEq_t * );
int UpdateC_ARX( SysEq_t *, vector *, vector *);
int UpdateC_ARX2( SysEq_t *, vector *, vector *);
int x2theta( vector *, matrix *, unsigned, unsigned, unsigned);
int X2theta( vector **, matrix *, vector *, unsigned, unsigned, unsigned);


int kalman_initialize( kalman_t *, unsigned int, unsigned int, unsigned int);
int kalman_destroy( kalman_t *);
int KalmanFilter( kalman_t *, SysEq_t *);
int kalman_initialize2( kalman_t *, unsigned int, unsigned int, unsigned int);
int kalman_destroy2( kalman_t *);
int KalmanFilter2( kalman_t *, SysEq_t *);


#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* GPC_H */

