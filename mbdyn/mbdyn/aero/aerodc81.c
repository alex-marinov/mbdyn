/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <stdio.h>
#include <mymath.h>

#include <aerodc81.h>

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

extern c81_data *get_c81_data(int jpro);

static int bisec(double* v, double val, int lb, int ub);
static int
get_coef(int nm, double* m, int na, double* a, double alpha, double mach,
		double* c, double* c0);
static double
get_dcla(int nm, double* m, double* s, double mach);
static int
get_stall(int nm, double* m, double* s, double mach,
          double *dcpa, double *dasp, double *dasm);

double
get_c81_coef(int nm, double* m, int na, double* a, double alpha, double mach)
{
	double c;
	
	get_coef(nm, m, na, a, alpha, mach, &c, NULL);

	return c;
}

int 
c81_aerod2(double* W, double* VAM, double* TNG, double* OUTA, c81_data* data)
{
   	/* 
	 * velocita' del punto in cui sono calcolate le condizioni al contorno
	 */
   	double v[3];
   	double vp, vp2, vtot;
	double rho = VAM[0];
	double cs = VAM[1];
	double chord = VAM[2];
	
	double cl = 0., cl0 = 0., cd = 0., cd0 = 0., cm = 0.;
	double alpha, gamma, cosgam, mach, q;
	double dcla;
	
	double ca = VAM[3];
	double c34 = VAM[4];
	
	const double RAD2DEG = 180.*M_1_PI;
	const double M_PI_3 = M_PI/3.;
	
	/* 
	 * porta la velocita' al punto di calcolo delle boundary conditions 
	 */
	v[0] = W[0];
	v[1] = W[1]+c34*W[5];
	v[2] = W[2]-c34*W[4];
	
	vp2 = v[0]*v[0]+v[1]*v[1];
	vp = sqrt(vp2);
	
	vtot = sqrt(vp2+v[2]*v[2]);
	
	/*
	 * non considera velocita' al di sotto di 1.e-6
	 * FIXME: rendere parametrico?
	 */
	if (vp/cs < 1.e-6) {
		TNG[0] = 0.;
		TNG[1] = 0.;
		TNG[2] = 0.;
		TNG[3] = 0.;
		TNG[4] = 0.;
		TNG[5] = 0.;
		
		OUTA[1] = 0.;
		OUTA[2] = 0.;
		OUTA[3] = 0.;
		OUTA[4] = 0.;
		OUTA[5] = 0.;
		OUTA[6] = 0.;
		
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
	alpha = atan2(-v[1], v[0]);
	OUTA[1] = alpha*RAD2DEG;  
	gamma = atan2(-v[2], fabs(v[0]));	/* come in COE0 (aerod2.f) */
	/* gamma = atan2(-v[2], vp); */		/* secondo me (?!?) */
	OUTA[2] = gamma*RAD2DEG;
	
	if (fabs(gamma) > M_PI_3) {
		/* tanto ne viene preso il coseno ... */
		gamma = M_PI_3;
	}
	
	cosgam = cos(gamma);
	mach = (vtot*sqrt(cosgam))/cs;
	OUTA[3] = mach;
	
	/*
	 * Note: all angles in c81 files MUST be in degrees
	 */
	get_coef(data->NML, data->ml, data->NAL, data->al, OUTA[1], mach, 
			&cl, &cl0);
	get_coef(data->NMD, data->md, data->NAD, data->ad, OUTA[1], mach,
			&cd, &cd0);
	get_coef(data->NMM, data->mm, data->NAM, data->am, OUTA[1], mach,
			&cm, NULL);

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
		double dclatmp = (cl-cl0)/(alpha*cosgam);
		if (dclatmp < dcla) {
			dcla = dclatmp;
		}
	}
	cl = cl0+dcla*alpha;
	
	OUTA[4] = cl;
	OUTA[5] = cd;
	OUTA[6] = cm;

	q = .5*rho*chord*vp2;

	TNG[0] = -q*(cl*v[1]+cd*v[0])/vp;
	TNG[1] = q*(cl*v[0]-cd*v[1])/vp;
	TNG[2] = -q*cd0*v[2]/vp;
	TNG[3] = 0.;
	TNG[4] = -ca*TNG[2];
	TNG[5] = -q*chord*cm-ca*TNG[1];
	
	return 0;
}

/*
 * algoritmo di bisezione per ricerca efficiente degli angoli
 *
 * v:	array dei valori
 * val:	valore da cercare
 * lb:	indice inferiore
 * ub:	indice superiore
 *
 * restituisce l'indice corrispondente al valore di v
 * che approssima per eccesso val;
 * se val > v[ub] restituisce ub+1;
 * se val < v[lb] restituisce lb.
 * quindi v[ret-1] <= val <= v[ret]
 *
 * Nota: ovviamente si presuppone che i valori di v siano ordinati in modo
 * crescente e strettamente monotono, ovvero v[i] < v[i+1].
 */
static int 
bisec(double* v, double val, int lb, int ub)
{
	if (val < v[lb]) {
		return lb;
	}
	
	if (val > v[ub]) {
		return ub+1;
	}
	
	while (ub > lb+1) {
		int b = (lb+ub)/2;
		if (v[b] > val) {
			ub = b;
		} else if (v[b] < val) {
			lb = b;
		} else {
			return b;
		}
	}
   
   	return lb;
}

/*
 * trova un coefficiente dato l'angolo ed il numero di Mach
 *
 * il numero di Mach mach viene cercato nell'array m di lunghezza nm
 * ed interpolato linearmente; quindi il coefficiente alpha viene cercato
 * nella prima colonna della matrice a di dimensioni na x nm+1 ed interpolato
 * linearmente; infine il coefficiente corrispondente alla combinazione di
 * mach e alpha viene restituito.
 */
static int
get_coef(int nm, double* m, int na, double* a, double alpha, double mach,
		double* c, double* c0)
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
	 * im e' l'indice di m in cui si trova l'approssimazione per eccesso 
	 * di mach
	 */
	im = bisec(m, mach, 0, nm-1);
	if (im != nm) {
		im++;
	}
	
	/*
	 * ia e' l'indice della prima colonna di a in cui si trova
	 * l'approssimazione per eccesso di alpha
	 */
	if (c0 != NULL) {
		ia0 = bisec(a, 0., 0, na-1);
		if (ia0 != na) {
			ia0++;
		}
	}

	ia = bisec(a, alpha, 0, na-1);
	if (ia != na) {
		ia++;
	}

	/*
	 * nota: sono stati scartati i casi im == 0 e ia == 0 perche'
	 * impossibili per vari motivi
	 */
	if (im == nm) {
		if (c0 != NULL) {
			if (ia0 == na) {
				*c0 = a[na*(nm+1)-1];
			} else {
				double da = -a[ia0-1]/(a[ia0]-a[ia0-1]);
				*c0 = (1.-da)*a[na*nm+ia0-1]+da*a[na*nm+ia0];
			}
		}
		
		if (ia == na) {
			*c = a[na*(nm+1)-1];
		} else {
			double da = (alpha-a[ia-1])/(a[ia]-a[ia-1]);
			*c = (1.-da)*a[na*nm+ia-1]+da*a[na*nm+ia];
		}
	} else {
		double d;
		d = (mach-m[im-1])/(m[im]-m[im-1]);

		if (c0 != NULL) {
			if (ia0 == na) {
				*c0 = (1.-d)*a[na*(im+1)-1]+d*a[na*(im+2)-1];
			} else {
				double a1, a2, da;
				a1 = (1.-d)*a[na*im+ia0-1]+d*a[na*(im+1)+ia0-1];
				a2 = (1.-d)*a[na*im+ia0]+d*a[na*(im+1)+ia0];
				da = -a[ia0-1]/(a[ia0]-a[ia0-1]);
				*c0 = (1.-da)*a1+da*a2;
			}
		}

		if (ia == na) {
			*c = (1.-d)*a[na*(im+1)-1]+d*a[na*(im+2)-1];
		} else {
			double a1, a2, da;
			a1 = (1.-d)*a[na*im+ia-1]+d*a[na*(im+1)+ia-1];
			a2 = (1.-d)*a[na*im+ia]+d*a[na*(im+1)+ia];
			da = (alpha-a[ia-1])/(a[ia]-a[ia-1]);
			*c = (1.-da)*a1+da*a2;
		}
	}

	return 0;
}

static double
get_dcla(int nm, double* m, double* s, double mach)
{
	int im;
	
	mach = fabs(mach);
	
	/*
	 * im e' l'indice di m in cui si trova l'approssimazione per eccesso
	 * di mach
	 */
	im = bisec(m, mach, 0, nm-1);
	if (im != nm) {
		im++;
	}
	
	if (im == nm) {
		return s[3*nm-1];
	} else {
		double d = (mach-m[im-1])/(m[im]-m[im-1]);

		return (1.-d)*s[2*nm+im-1]+d*s[2*nm+im];
	}
}

static int
get_stall(int nm, double* m, double* s, double mach,
	  double *dcpa, double *dasp, double *dasm)
{
	int im;

	mach = fabs(mach);

	/*
	 * im e' l'indice di m in cui si trova l'approssimazione per eccesso
	 * di mach
	 */
	im = bisec(m, mach, 0, nm-1);
	if (im != nm) {
		im++;
	}
	
	if (im == nm) {
		*dcpa = s[3*nm-1];
		*dasp = s[nm-1];
		*dasm = s[2*nm-1];
	} else {
		double d = (mach-m[im-1])/(m[im]-m[im-1]);

		*dcpa = (1.-d)*s[2*nm+im-1]+d*s[2*nm+im];
		*dasp = (1.-d)*s[im-1]+d*s[im];
		*dasm = (1.-d)*s[nm+im-1]+d*s[nm+im];
	}
	
	return 0;
}

