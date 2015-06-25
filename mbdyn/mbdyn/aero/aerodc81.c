/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ac/f2c.h"
#include "aerodc81.h"
#include "bisec.h"
#include "c81data.h"

/*
cC*************************          AEROD          *********************
c                        MODIFICATA NEL CALCOLO TNG
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE AEROD2(W, VAM, TNG, OUTA, INST, RSPEED, JPRO)
C   W: vettore lungo 6; velocita' del corpo aerodinamico nel sistema locale,
C      riferita ad un punto arbitrario
C   VAM:    vettore lungo 6, contiene dati del problema
C   TNG: vettore lungo 12, presumubilmente il termine noto
C        (ora l'ho ridotto a 6)
C   OUTA: vettore lungo 20, usato in i/o con subroutines di aerod
C   INST: flag di instazionario (>0)
C   RSPEED -> MODOMEGA: modulo della velocita' di rotazione
C   JPRO: tipo di profilo      
C     DEFINIZIONI VETTORE VAM
C
C VAM(1): densita' dell'aria
C VAM(2): celerita' del suono
C VAM(3): corda
C VAM(4): 1/4 corda
C VAM(5): 3/4 corda
C VAM(6): svergolamento
C
C W = [v1 v2 v3 w1 w2 w3 ]
C con:
C 1 direzione di avanzamento,
C 2 direzione normale,
C 3 direzione lungo l'apertura,
C DM*W da' la velocita' nel punto a 3/4 della corda      
*/

static int
get_coef(int nm, doublereal* m, int na, doublereal* a,
		doublereal alpha, doublereal mach,
		doublereal* c, doublereal* c0);

static doublereal
get_dcla(int nm, doublereal* m, doublereal* s, doublereal mach);

#ifdef USE_GET_STALL
static int
get_stall(int nm, doublereal* m, doublereal* s, doublereal mach,
          doublereal *dcpa, doublereal *dasp, doublereal *dasm);
#endif /* USE_GET_STALL */

const outa_t outa_Zero;

doublereal
c81_data_get_coef(int nm, doublereal* m, int na, doublereal* a, doublereal alpha, doublereal mach)
{
	doublereal c;
	
	get_coef(nm, m, na, a, alpha, mach, &c, NULL);

	return c;
}

int 
c81_aerod2(doublereal* W, const vam_t *VAM, doublereal* TNG, outa_t* OUTA, c81_data* data)
{
   	/* 
	 * velocita' del punto in cui sono calcolate le condizioni al contorno
	 */
   	doublereal v[3];
   	doublereal vp, vp2, vtot;
	doublereal rho = VAM->density;
	doublereal cs = VAM->sound_celerity;
	doublereal chord = VAM->chord;
	
	doublereal cl = 0., cl0 = 0., cd = 0., cd0 = 0., cm = 0.;
	doublereal alpha, gamma, cosgam, mach, q;
	doublereal dcla;
	
	doublereal ca = VAM->force_position;
	doublereal c34 = VAM->bc_position;
	
	const doublereal RAD2DEG = 180.*M_1_PI;
	const doublereal M_PI_3 = M_PI/3.;

	enum { V_X = 0, V_Y = 1, V_Z = 2, W_X = 3, W_Y = 4, W_Z = 5 };
	
	/* 
	 * porta la velocita' al punto di calcolo delle boundary conditions 
	 */
	v[V_X] = W[V_X];
	v[V_Y] = W[V_Y] + c34*W[W_Z];
	v[V_Z] = W[V_Z] - c34*W[W_Y];
	
	vp2 = v[V_X]*v[V_X] + v[V_Y]*v[V_Y];
	vp = sqrt(vp2);
	
	vtot = sqrt(vp2 + v[V_Z]*v[V_Z]);
	
	/*
	 * non considera velocita' al di sotto di 1.e-6
	 * FIXME: rendere parametrico?
	 */
	if (vp/cs < 1.e-6) {
		TNG[V_X] = 0.;
		TNG[V_Y] = 0.;
		TNG[V_Z] = 0.;
		TNG[W_X] = 0.;
		TNG[W_Y] = 0.;
		TNG[W_Z] = 0.;
		
		OUTA->alpha = 0.;
		OUTA->gamma = 0.;
		OUTA->mach = 0.;
		OUTA->cl = 0.;
		OUTA->cd = 0.;
		OUTA->cm = 0.;
		OUTA->clalpha = 0.;
		
		return 0;
	}
	
	/*
	 * FIXME: gestire quali angoli indicano cosa
	 *
	 * Idea di base: capire in base ad un  criterio di linearita'
	 * a quale angolo di incidenza il profilo stalla, quindi
	 * oltre quell'angolo prendere la correzione per flusso trasverso
	 * in un certo modo, altrimenti in un altro ?!?
	 */
	alpha = atan2(-v[V_Y], v[V_X]);
	OUTA->alpha = alpha*RAD2DEG;  
	gamma = atan2(-v[V_Z], fabs(v[V_X]));	/* come in COE0 (aerod2.f) */
	/* gamma = atan2(-v[V_Z], vp); */		/* secondo me (?!?) */
	OUTA->gamma = gamma*RAD2DEG;
	
	if (fabs(gamma) > M_PI_3) {
		/* tanto ne viene preso il coseno ... */
		gamma = M_PI_3;
	}
	
	cosgam = cos(gamma);
	mach = (vtot*sqrt(cosgam))/cs;
	OUTA->mach = mach;
	
	/*
	 * Note: all angles in c81 files MUST be in degrees
	 */
	get_coef(data->NML, data->ml, data->NAL, data->al,
			OUTA->alpha, mach, &cl, &cl0);
	get_coef(data->NMD, data->md, data->NAD, data->ad,
			OUTA->alpha, mach, &cd, &cd0);
	get_coef(data->NMM, data->mm, data->NAM, data->am,
			OUTA->alpha, mach, &cm, NULL);

	dcla = get_dcla(data->NML, data->ml, data->stall, mach);
	
/*
 * da COE0 (aerod2.f):
 * 
	ASLRF = ASLOP0
	IF(DABS(ALFA).LT.1.D-6) GOTO 10
	ASLRF = CLIFT/(ALFA*COSGAM)
	IF(ASLRF.GT.ASLOP0) ASLRF = ASLOP0
     10 CLIFT = ASLRF*ALFA
 *
 */

	/*
	 * in soldoni: se si e' oltre il tratto lineare, prende la
	 * secante con l'angolo corretto per la freccia (fino a 60 
	 * gradi) e poi ricalcola il cl con l'angolo vero; in questo
	 * modo il cl viene piu' grande di circa 1/cos(gamma) ma solo
	 * fuori del tratto lineare.
	 */
	dcla *= RAD2DEG;
	if (fabs(alpha) > 1.e-6) {
		doublereal dclatmp = (cl - cl0)/(alpha*cosgam);
		if (dclatmp < dcla) {
			dcla = dclatmp;
		}
	}
	cl = cl0 + dcla*alpha;
	
	OUTA->cl = cl;
	OUTA->cd = cd;
	OUTA->cm = cm;
	OUTA->clalpha = dcla;

	q = .5*rho*chord*vp2;

	TNG[V_X] = -q*(cl*v[V_Y] + cd*v[V_X])/vp;
	TNG[V_Y] = q*(cl*v[V_X] - cd*v[V_Y])/vp;
	TNG[V_Z] = -q*cd0*v[V_Z]/vp;
	TNG[W_X] = 0.;
	TNG[W_Y] = -ca*TNG[V_Z];
	TNG[W_Z] = q*chord*cm + ca*TNG[V_Y];
	
	/* 
	 * Radial drag (TNG[V_Z]) consistent with Harris, JAHS 1970
	 * and with CAMRAD strip theory section forces 
	 */

	return 0;
}

int 
c81_aerod2_u(doublereal* W, const vam_t *VAM, doublereal* TNG, outa_t* OUTA, 
		c81_data* data, long unsteadyflag)
{
   	/* 
	 * velocita' del punto in cui sono calcolate le condizioni al contorno
	 */
   	doublereal v[3];
   	doublereal vp, vp2, vtot;
	doublereal rho = VAM->density;
	doublereal cs = VAM->sound_celerity;
	doublereal chord = VAM->chord;
	
	doublereal cl = 0., cl0 = 0., cd = 0., cd0 = 0., cm = 0.;
	doublereal alpha, gamma, cosgam, mach, q;
	doublereal dcla;
	
	doublereal ca = VAM->force_position;
	doublereal c34 = VAM->bc_position;

	const doublereal RAD2DEG = 180.*M_1_PI;
	const doublereal M_PI_3 = M_PI/3.;

	enum { V_X = 0, V_Y = 1, V_Z = 2, W_X = 3, W_Y = 4, W_Z = 5 };
	
	/* 
	 * porta la velocita' al punto di calcolo delle boundary conditions 
	 */
	v[V_X] = W[V_X];
	v[V_Y] = W[V_Y] + c34*W[W_Z];
	v[V_Z] = W[V_Z] - c34*W[W_Y];
	
	vp2 = v[V_X]*v[V_X] + v[V_Y]*v[V_Y];
	vp = sqrt(vp2);
	
	vtot = sqrt(vp2 + v[V_Z]*v[V_Z]);
	
	/*
	 * non considera velocita' al di sotto di 1.e-6
	 * FIXME: rendere parametrico?
	 */

	if (vp/cs < 1.e-6) {
		TNG[V_X] = 0.;
		TNG[V_Y] = 0.;
		TNG[V_Z] = 0.;
		TNG[W_X] = 0.;
		TNG[W_Y] = 0.;
		TNG[W_Z] = 0.;
		
		OUTA->alpha = 0.;
		OUTA->gamma = 0.;
		OUTA->mach = 0.;
		OUTA->cl = 0.;
		OUTA->cd = 0.;
		OUTA->cm = 0.;
		OUTA->clalpha = 0.;
		
		return 0;
	}
	
	/*
	 * FIXME: gestire quali angoli indicano cosa
	 *
	 * Idea di base: capire in base ad un  criterio di linearita'
	 * a quale angolo di incidenza il profilo stalla, quindi
	 * oltre quell'angolo prendere la correzione per flusso trasverso
	 * in un certo modo, altrimenti in un altro ?!?
	 */
	alpha = atan2(-v[V_Y], v[V_X]);
	OUTA->alpha = alpha*RAD2DEG;  
	gamma = atan2(-v[V_Z], fabs(v[V_X]));	/* come in COE0 (aerod2.f) */
	/* gamma = atan2(-v[V_Z], vp); */		/* secondo me (?!?) */
	OUTA->gamma = gamma*RAD2DEG;
	
	if (fabs(gamma) > M_PI_3) {
		/* tanto ne viene preso il coseno ... */
		gamma = M_PI_3;
	}
	
	cosgam = cos(gamma);
	mach = (vtot*sqrt(cosgam))/cs;
	OUTA->mach = mach;

	/*
	 * mach cannot be more than .99 (see aerod.f)
	 */
	if (mach > .99) {
		mach = .99;
	}

	/*
	 * Compute cl, cd, cm based on selected theory
	 */
	switch (unsteadyflag) {
	case 0: 

		/*
		 * Note: all angles in c81 files MUST be in degrees
		 */
		get_coef(data->NML, data->ml, data->NAL, data->al, 
				OUTA->alpha, mach, &cl, &cl0);
		get_coef(data->NMD, data->md, data->NAD, data->ad, 
				OUTA->alpha, mach, &cd, &cd0);
		get_coef(data->NMM, data->mm, data->NAM, data->am, 
				OUTA->alpha, mach, &cm, NULL);

		dcla = get_dcla(data->NML, data->ml, data->stall, mach);
	
/*
 * da COE0 (aerod2.f):
 * 
	ASLRF = ASLOP0
	IF(DABS(ALFA).LT.1.D-6) GOTO 10
	ASLRF = CLIFT/(ALFA*COSGAM)
	IF(ASLRF.GT.ASLOP0) ASLRF = ASLOP0
     10 CLIFT = ASLRF*ALFA
 *
 */

		/*
		 * in soldoni: se si e' oltre il tratto lineare, prende la
		 * secante con l'angolo corretto per la freccia (fino a 60 
		 * gradi) e poi ricalcola il cl con l'angolo vero; in questo
		 * modo il cl viene piu' grande di circa 1/cos(gamma) ma solo
		 * fuori del tratto lineare.
		 */
		dcla *= RAD2DEG;
		if (fabs(alpha) > 1.e-6) {
			doublereal dclatmp = (cl - cl0)/(alpha*cosgam);
			if (dclatmp < dcla) {
				dcla = dclatmp;
				cl = cl0 + dcla*alpha;
			}
		}
		break;

	case 1:
		return -1;

	case 2: {

		/*
		 * Constants from unsteady theory
		 * synthetized by Richard L. Bielawa,
		 * 31th A.H.S. Forum Washington D.C. 
		 * May 1975
		 */
		doublereal A, B, A2, B2, ETA, ASN, ASM, 
			SGN, SGM, SGMAX, 
			DAN, DCN, DAM, DCM, 
			S2, alphaN, alphaM, C1,
			dcma, dclatan, ALF1, ALF2,
			cn;
		
		const doublereal PN[] = { 
			-3.464003e-1, 
			-1.549076e+0, 
			4.306330e+1, 
			-5.397529e+1,
			5.781402e+0,
			-3.233003e+1,
			-2.162257e+1,
			1.866347e+1,
			4.198390e+1,
			3.295461e+2,
		};
	
		const doublereal QN[] = {
			1.533717e+0,
			6.977203e+0,
			1.749010e+3,
			1.694829e+3,
			-1.771899e+3,
			-3.291665e+4,
			2.969051e+0,
			-3.632448e+1,
			-2.268578e+3,
			6.601995e+3,
			-9.654208e+3, 
			8.533930e+4,
			-1.492624e+0,
			1.163661e+1
		};

		const doublereal PM[] = {
			1.970065e+1,
			-6.751639e+1,
			7.265269e+2,
			4.865945e+4,
			2.086279e+4,
			6.024672e+3,
			1.446334e+2,
			8.586896e+2,
			-7.550329e+2,
			-1.021613e+1,
			2.247664e+1,
		};
	
		const doublereal QM[] = {
			-2.322808e+0,
			-1.322257e+0,
			-2.633891e+0,
			-2.180321e-1,
			4.580014e+0,
			3.125497e-1,
			-2.828806e+1,
			-4.396734e+0,
			2.565870e+2,
			-1.204976e+1,
			-1.157802e+2,
			8.612138e+0,
		};

		enum {
			U_1 = 0,
			U_2 = 1,
			U_3 = 2,
			U_4 = 3,
			U_5 = 4,
			U_6 = 5,
			U_7 = 6,
			U_8 = 7,
			U_9 = 8,
			U10 = 9,
			U11 = 10,
			U12 = 11,
			U13 = 12,
			U14 = 13
		};

		/*
		 * This is the static stall angle for Mach = 0
		 * (here a symmetric airfoil is assumed; the real
		 * _SIGNED_ static stall should be considered ...)
		 */
		const doublereal ASN0 = .22689, ASM0 = .22689;

		ALF1 = OUTA->alf1;
		ALF2 = OUTA->alf2;

		A = .5*chord*ALF1/vp;
		B = .25*chord*chord*ALF2/vp2;

		ETA = sqrt(pow(A/.048, 2) + pow(B/.016, 2));

		if (alpha < 0.) {
			A = -A;
			B = -B;
		}

		if (ETA > 1.) {
			A /= ETA;
			B /= ETA;
		}

		A2 = A*A;
		B2 = B*B;

		ASN = ASN0*(1. - mach);
		ASM = ASM0*(1. - mach);

		SGN = fabs(alpha/ASN);
		SGM = fabs(alpha/ASM);
		SGMAX = 1.839 - 70.33*fabs(B);
		if (SGMAX > 1.86) {
			SGMAX = 1.86;
		}
		if (SGN > SGMAX) {
			SGN = SGMAX;
		}
		if (SGM > SGMAX) {
			SGM = SGMAX;
		}

		DAN = (A*(PN[U_1] + PN[U_5]*SGN) + B*(PN[U_2] + PN[U_6]*SGN)
			+ exp(-1072.52*A2)*(A*(PN[U_3] + PN[U_7]*SGN)
				+A2*(PN[U_9] + PN[U10]*SGN))
			+ exp(-40316.42*B2)*B*(PN[U_4] + PN[U_8]*SGN))*ASN;
		DCN = A*(QN[U_1] + QN[U_3]*A2 + SGN*(QN[U_7] + QN[U_9]*A2 + QN[U13]*SGN)
				+ B2*(QN[U_5] + QN[U11]*SGN))
			+ B*(QN[U_2] + QN[U_4]*A2
					+ SGN*(QN[U_8] + QN[U10]*A2 + QN[U14]*SGN)
					+ B2*(QN[U_6] + QN[U12]*SGN));
		DAM = (A*(PM[U_1] + PM[U_3]*A2 + PM[U_5]*B2 + PM[U10]*SGM + PM[U_7]*A)
				+ B*(PM[U_2] + PM[U_4]*B2 + PM[U_6]*A2
					+ PM[U11]*SGM + PM[U_8]*B + PM[U_9]*A))*ASM;

		S2 = SGM*SGM;

		DCM = A*(QM[U_2] + QM[U_8]*A + SGM*(QM[U_4] + QM[U10]*A)
				+ S2*(QM[U_6] + QM[U12]*A))
			+ B*(QM[U_1] + QM[U_7]*B + SGM*(QM[U_3] + QM[U_9]*B)
					+ S2*(QM[U_5] + QM[U11]*B));

		OUTA->dan = DAN*RAD2DEG;
		OUTA->dam = DAM*RAD2DEG;
		OUTA->dcn = DCN;
		OUTA->dcm = DCM;

		/* 
		 * I think I need to apply this contribution 
		 * with the sign of alpha, because otherwise 
		 * it gets discontinuous as alpha changes sign;
		 * I definitely need the original reference :(
		 */
		if (alpha < 0.) {
			DAN = -DAN;
			DCN = -DCN;
			DAM = -DAM;
			DCM = -DCM;
		}

		alphaN = (alpha - DAN)*RAD2DEG;
		get_coef(data->NML, data->ml, data->NAL, data->al, 
				alphaN, mach, &cl, &cl0);
		get_coef(data->NMD, data->md, data->NAD, data->ad, 
				alphaN, mach, &cd, &cd0);

		alphaM = (alpha - DAM)*RAD2DEG;
		get_coef(data->NMM, data->mm, data->NAM, data->am, 
				alphaM, mach, &cm, NULL);

		dcla = get_dcla(data->NML, data->ml, data->stall, mach);
		dcma = get_dcla(data->NMM, data->mm, data->mstall, mach);

		/* note: cl/alpha in 1/deg */
		dclatan = dcla;
		if (fabs(alphaN) > 1.e-6) {
			dclatan = (cl - cl0)/(alphaN*cosgam);
		}
		cl = cl0 + dclatan*alphaN;

		/* back to 1/rad */
		dcla *= RAD2DEG;
		dcma *= RAD2DEG;

		C1 = .9457/sqrt(1. - mach*mach);

		/* 
		 * the unsteady correction is "cn", 
		 * so split it in "cl" and "cd"
		 * (note: if "vp" is too small the routine exits
		 * earlier without computing forces)
		 */
		cn = dcla*DAN + DCN*C1;
		cl += cn*v[V_X]/vp;	/* cos(x) */
		cd -= cn*v[V_Y]/vp;	/* sin(x) */
		cm += dcma*DAM + DCM*C1;

		break;
	}
	}

	/*
	 * Save cl, cd, cm for output purposes
	 */
	OUTA->cl = cl;
	OUTA->cd = cd;
	OUTA->cm = cm;
	OUTA->clalpha = dcla;

	/*
	 * Local dynamic pressure
	 */
	q = .5*rho*chord*vp2;

	/*
	 * airfoil forces and moments in the airfoil frame
	 */
	TNG[V_X] = -q*(cl*v[V_Y] + cd*v[V_X])/vp;
	TNG[V_Y] = q*(cl*v[V_X] - cd*v[V_Y])/vp;
	TNG[V_Z] = -q*cd0*v[V_Z]/vp;
	TNG[W_X] = 0.;
	TNG[W_Y] = -ca*TNG[V_Z];
	TNG[W_Z] = q*chord*cm + ca*TNG[V_Y];
	
	/* 
	 * Radial drag (TNG[V_Z]) consistent with Harris, JAHS 1970
	 * and with CAMRAD strip theory section forces 
	 */

	return 0;
}

/*
 * trova un coefficiente dato l'angolo ed il numero di Mach
 *
 * il numero di Mach mach viene cercato nell'array m di lunghezza nm
 * ed interpolato linearmente; quindi il coefficiente alpha viene cercato
 * nella prima colonna della matrice a di dimensioni na x nm + 1 ed interpolato
 * linearmente; infine il coefficiente corrispondente alla combinazione di
 * mach e alpha viene restituito.
 */
static int
get_coef(int nm, doublereal* m, int na, doublereal* a, doublereal alpha, doublereal mach,
		doublereal* c, doublereal* c0)
{
   	int im;
   	int ia, ia0 = -1;
	
	while (alpha < -180.) {
		alpha += 360.;
	}
	
	while (alpha >= 180.) {
		alpha -= 360.;
	}
	
	mach = fabs(mach);
	
	/*
	 * im e' l'indice di m in cui si trova
	 * l'approssimazione per difetto di mach
	 */
	im = bisec_d(m, mach, 0, nm - 1);
	
	/*
	 * ia e' l'indice della prima colonna di a in cui si trova
	 * l'approssimazione per difetto di alpha
	 */
	if (c0 != NULL) {
		ia0 = bisec_d(a, 0., 0, na - 1);
	}

	ia = bisec_d(a, alpha, 0, na - 1);

	if (im == nm - 1) {
		if (c0 != NULL) {
			if (ia0 == na - 1) {
				*c0 = a[na*(nm + 1) - 1];

			} else if (ia0 == -1) {
				*c0 = a[na*nm];

			} else {
				doublereal da;

				ia0++;
				da = -a[ia0 - 1]/(a[ia0] - a[ia0 - 1]);
				*c0 = (1. - da)*a[na*nm + ia0 - 1] + da*a[na*nm + ia0];
			}
		}
		
		if (ia == na - 1) {
			*c = a[na*(nm + 1) - 1];

		} else if (ia == -1) {
			*c = a[na*nm];

		} else {
			doublereal da;

			ia++;
			da = (alpha - a[ia - 1])/(a[ia] - a[ia - 1]);
			*c = (1. - da)*a[na*nm + ia - 1] + da*a[na*nm + ia];
		}

	} else if (im == -1) {
		if (c0 != NULL) {
			if (ia0 == na - 1) {
				*c0 = a[na*2 - 1];

			} else if (ia0 == -1) {
				*c0 = a[na];

			} else {
				doublereal da;

				ia0++;
				da = -a[ia0 - 1]/(a[ia0] - a[ia0 - 1]);
				*c0 = (1. - da)*a[na + ia0 - 1] + da*a[na + ia0];
			}
		}
		
		if (ia == na - 1) {
			*c = a[na*2 - 1];

		} else if (ia == -1) {
			*c = a[na];

		} else {
			doublereal da;

			ia++;
			da = (alpha - a[ia - 1])/(a[ia] - a[ia - 1]);
			*c = (1. - da)*a[na + ia - 1] + da*a[na + ia];
		}

	} else {
		doublereal d;

		im++;
		d = (mach - m[im - 1])/(m[im] - m[im - 1]);

		if (c0 != NULL) {
			if (ia0 == na) {
				*c0 = (1. - d)*a[na*(im + 1) - 1] + d*a[na*(im + 2) - 1];
			} else {
				doublereal a1, a2, da;
				a1 = (1. - d)*a[na*im + ia0 - 1] + d*a[na*(im + 1) + ia0 - 1];
				a2 = (1. - d)*a[na*im + ia0] + d*a[na*(im + 1) + ia0];
				da = -a[ia0 - 1]/(a[ia0] - a[ia0 - 1]);
				*c0 = (1. - da)*a1 + da*a2;
			}
		}

		if (ia == na - 1) {
			*c = (1. - d)*a[na*(im + 1) - 1] + d*a[na*(im + 2) - 1];

		} else if (ia == -1) {
			*c = (1. - d)*a[na*im] + d*a[na*(im + 1)];

		} else {
			doublereal a1, a2, da;

			ia++;
			a1 = (1. - d)*a[na*im + ia - 1] + d*a[na*(im + 1) + ia - 1];
			a2 = (1. - d)*a[na*im + ia] + d*a[na*(im + 1) + ia];
			da = (alpha - a[ia - 1])/(a[ia] - a[ia - 1]);
			*c = (1. - da)*a1 + da*a2;
		}
	}

	return 0;
}

static doublereal
get_dcla(int nm, doublereal* m, doublereal* s, doublereal mach)
{
	int im;
	
	mach = fabs(mach);
	
	/*
	 * im e' l'indice di m in cui si trova l'approssimazione per eccesso
	 * di mach
	 */
	im = bisec_d(m, mach, 0, nm - 1);
	
	if (im == nm - 1) {
		return s[3*nm - 1];

	} else if (im == -1) {
		return s[2*nm];

	} else {
		doublereal d;

		im++;
		d = (mach - m[im - 1])/(m[im] - m[im - 1]);

		return (1. - d)*s[2*nm + im - 1] + d*s[2*nm + im];
	}
}

#ifdef USE_GET_STALL
static int
get_stall(int nm, doublereal* m, doublereal* s, doublereal mach,
	  doublereal *dcpa, doublereal *dasp, doublereal *dasm)
{
	int im;

	mach = fabs(mach);

	/*
	 * im e' l'indice di m in cui si trova l'approssimazione per eccesso
	 * di mach
	 */
	im = bisec_d(m, mach, 0, nm - 1);
	if (im != nm) {
		im++;
	}
	
	if (im == nm) {
		*dcpa = s[3*nm - 1];
		*dasp = s[nm - 1];
		*dasm = s[2*nm - 1];
	} else {
		doublereal d = (mach - m[im - 1])/(m[im] - m[im - 1]);

		*dcpa = (1. - d)*s[2*nm + im - 1] + d*s[2*nm + im];
		*dasp = (1. - d)*s[im - 1] + d*s[im];
		*dasm = (1. - d)*s[nm + im - 1] + d*s[nm + im];
	}
	
	return 0;
}
#endif /* USE_GET_STALL */

