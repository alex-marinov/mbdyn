/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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
   	double *ml;
   	double *al;

	/*
	 * matrice dei dati di stallo:
	 *  0 -> NAL-1 : 	angoli di stallo positivi
	 *  NAL -> 2*NAL-1 :	angoli di stallo negativi
	 *  2*NAL -> 3*NAL-1 :	Cp/alpha
	 */
	double *stall;
	double *mstall;
   
   	int NMD;
   	int NAD;
   	double *md;
   	double *ad;
   
   	int NMM;
   	int NAM;
   	double *mm;
   	double *am;
} c81_data;

extern int 
c81_aerod2(double* W, double* VAM, double* TNG, double* OUTA, c81_data* data);

extern int 
c81_aerod2_u(double* W, double* VAM, double* TNG, double* OUTA, 
		c81_data* data, long unsteadyflag);

#ifdef __cplusplus
}
#endif /* __cplusplus */
   
#endif /* AERODC81_H */

