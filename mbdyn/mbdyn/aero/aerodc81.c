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

extern c81_data* get_c81_data(int jpro);
static int bisec(double* v, double val, int lb, int ub);
double get_coef(int nm, double* m, int na, double* a, double alpha, double mach);

int 
c81_aerod2(double* W, double* VAM, double* TNG, double* OUTA, c81_data* data)
{
   	/* 
	 * velocita' del punto in cui sono calcolate le condizioni al contorno
	 */
   	double v[3];
   	double vp, vp2, vtot, vtot2;
	double rho = VAM[0];
	double cs = VAM[1];
	double chord = VAM[2];
	
	double cl, cd, cd0, cm;
	double alpha, gamma, cosgam, mach, q;
	
	double ca = VAM[3];
	double c34 = VAM[4];
	
	const double RAD2DEG = 180.*M_1_PI;
	const double M_PI_3 = M_PI/3.;
	
	/* porta la velocita' al punto di calcolo delle boundary conditions */
	v[0] = W[0];
	v[1] = W[1]+c34*W[5];
	v[2] = W[2]-c34*W[4];
	
	vp2 = v[0]*v[0]+v[1]*v[1];
	vp = sqrt(vp2);
	
	vtot2 = vp2+v[2]*v[2];
	vtot = sqrt(vtot2);
	
	/*
	 * non considera velocita' al di sotto di 1.e-3
	 * FIXME: rendere parametrico?
	 */
	if (vp/cs < 1.e-3) {
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
	/* gamma = atan2(-v[2], fabs(v[0])); */
	gamma = atan2(-v[2], vp);
	OUTA[2] = gamma*RAD2DEG;
	
	if (fabs(gamma) > M_PI_3) {
		gamma = M_PI_3;
	}
	
	cosgam = cos(gamma);
	mach = vtot/cs*sqrt(cosgam);
	OUTA[3] = mach;
	
	cl = get_coef(data->NML, data->ml, data->NAL, data->al, OUTA[1], mach);
	cd = get_coef(data->NMD, data->md, data->NAD, data->ad, OUTA[1], mach);
	cd0 = get_coef(data->NMD, data->md, data->NAD, data->ad, 0., mach);
	cm = get_coef(data->NMM, data->mm, data->NAM, data->am, OUTA[1], mach);
	
	OUTA[4] = cl;
	OUTA[5] = cd;
	OUTA[6] = cm;
	
	q = .5*rho*chord*vp2;
	
	TNG[0] = -q*(cl/cosgam*v[1]+cd*v[0])/vp;
	TNG[1] = q*(cl/cosgam*v[0]-cd*v[1])/vp;
	TNG[2] = -q*cd0*v[2]/vp;
	TNG[3] = 0.;
	TNG[4] = -ca*TNG[2];
	TNG[5] = -q*chord*cm-ca*TNG[1];
	
	return 0;
}

/* algoritmo di bisezione per ricerca efficiente di angoli */
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

/* trova un coefficiente dato l'angolo ed il numero di Mach */
double
get_coef(int nm, double* m, int na, double* a, double alpha, double mach)
{
   	int im;
   	int ia;
  	double c = 0.;
	
	while (alpha < -180.) {
		alpha += 360.;
	}
	
	while (alpha >= 180.) {
		alpha -= 360.;
	}
	
	mach = fabs(mach);
	
	im = bisec(m, mach, 0, nm-1);
	if (im != nm) {
		im++;
	}
	
	ia = bisec(a, alpha, 0, na-1);
	if (ia != na) {
		ia++;
	}
	
	if (im == 0) {
		if (ia == 0) {
			c = a[na];
		} else if (ia == na) {
			c = a[2*na-1];
		} else {
			double da = (alpha-a[ia-1])/(a[ia]-a[ia-1]);
			c = (1.-da)*a[na+ia-1]+da*a[na+ia];
		}
	} else if (im == nm) {
		if (ia == 0) {
			c = a[na*nm];
		} else if (ia == na) {
			c = a[na*(nm+1)-1];
		} else {
			double da = (alpha-a[ia-1])/(a[ia]-a[ia-1]);
			c = (1.-da)*a[na*nm+ia-1]+da*a[na*nm+ia];
		}
	} else {
		double d;
		d = (mach-m[im-1])/(m[im]-m[im-1]);      
		if (ia == 0) {
			c = (1.-d)*a[na*im]+d*a[na*(im+1)];
		} else if (ia == na) {
			c = (1.-d)*a[na*(im+1)-1]+d*a[na*(im+2)-1];
		} else {
			double a1, a2, da;
			a1 = (1.-d)*a[na*im+ia-1]+d*a[na*(im+1)+ia-1];
			a2 = (1.-d)*a[na*im+ia]+d*a[na*(im+1)+ia];
			da = (alpha-a[ia-1])/(a[ia]-a[ia-1]);
			c = (1.-da)*a1+da*a2;
		}
	}
	
	return c;
}

