/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#ifndef AERODC81_H
#define AERODC81_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

enum {
	// 0 unused?
	OUTA_ALPHA	= 1,
	OUTA_GAMMA	= 2,
	OUTA_MACH	= 3,
	OUTA_CL		= 4,
	OUTA_CD		= 5,
	OUTA_CM		= 6,
	// 7 unused?
	OUTA_CLALPHA	= 7,
	OUTA_ALF1	= 8,
	OUTA_ALF2	= 9,
	OUTA_DAN	= 10,
	OUTA_DAM	= 11,
	OUTA_DCN	= 12,
	OUTA_DCM	= 13,
	// 14-19 unused?

	OUTA_LAST	= 20
};

typedef struct outa_t {
	doublereal outa_00;
	doublereal alpha;
	doublereal gamma;
	doublereal mach;
	doublereal cl;
	doublereal cd;
	doublereal cm;
	// doublereal outa_07;
	doublereal clalpha;
	doublereal alf1;
	doublereal alf2;
	doublereal dan;
	doublereal dam;
	doublereal dcn;
	doublereal dcm;
	doublereal outa_14;
	doublereal outa_15;
	doublereal outa_16;
	doublereal outa_17;
	doublereal outa_18;
	doublereal outa_19;
} outa_t;

extern const outa_t outa_Zero;

/*
C     DEFINIZIONI VETTORE VAM
C
C VAM(1): densita' dell'aria
C VAM(2): celerita' del suono
C VAM(3): corda
C VAM(4): 1/4 corda
C VAM(5): 3/4 corda
C VAM(6): svergolamento
 */
enum {
	VAM_DENSITY = 0,
	VAM_SOUND_CELERITY,
	VAM_CHORD,
	VAM_FORCE_POSITION,
	VAM_BC_POSITION,
	VAM_TWIST,

	VAM_LAST
};

typedef struct vam_t {
	doublereal density;
	doublereal sound_celerity;
	doublereal chord;
	doublereal force_position;
	doublereal bc_position;
	doublereal twist;
} vam_t;

/* 
 * dati in formato c81;
 * le array mX contengono gli NMX numeri di mach;
 * le array aX contengono gli NMX*NAX coefficienti, 
 * preceduti dall'angolo di incidenza e row-oriented.
 * 
 * Quindi, il numero di Mach j-esimo (a base 0) e'
 * 
 *     mX[j]
 * 
 * l'angolo di incidenza i-esimo (a base 0) e'
 * 
 *     aX[i]
 * 
 * e per avere il coefficiente relativo all'angolo
 * di incidenza i-esimo, al numero di mach j-esimo
 * 
 *     aX[NAX*(j+1)+i]
 *
 * MODIFICA 2000/10/02:
 * angolo di stallo o fine linearita' per correzione flusso trasverso
 *
 * Ora si aggiunge:
 *   - un vettore di angoli ai quali si perde linearita' (+)
 *   - un vettore di angoli ai quali si perde linearita' (-)
 *   - un vettore con il Cp/alpha
 */
typedef struct c81_data {
   	char header[31];
   
   	int NML;
   	int NAL;
   	doublereal *ml;
   	doublereal *al;

	/*
	 * matrice dei dati di stallo:
	 *  0 -> NAL-1 : 	angoli di stallo positivi
	 *  NAL -> 2*NAL-1 :	angoli di stallo negativi
	 *  2*NAL -> 3*NAL-1 :	Cp/alpha
	 */
	doublereal *stall;
	doublereal *mstall;
   
   	int NMD;
   	int NAD;
   	doublereal *md;
   	doublereal *ad;
   
   	int NMM;
   	int NAM;
   	doublereal *mm;
   	doublereal *am;
} c81_data;

extern int 
c81_aerod2(doublereal* W, const vam_t *VAM, doublereal* TNG, outa_t* OUTA, c81_data* data);

extern int 
c81_aerod2_u(doublereal* W, const vam_t *VAM, doublereal* TNG, outa_t* OUTA, 
		c81_data* data, long unsteadyflag);

#ifdef __cplusplus
}
#endif /* __cplusplus */
   
#endif /* AERODC81_H */

