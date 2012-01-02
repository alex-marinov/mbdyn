/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/* Elementi aerodinamici
 * 
 * - Proprieta' dell'aria:
 * elemento unico, contiene:
 *   - direzione e modulo del vento relativo. Il modulo e' associato
 *     ad un driver che ne consente la variazione, per consentire la 
 *     simulazione di un transitorio di galleria.
 *     Verra' aggiunta una variazione della direzione per simulare la raffica.
 *   - densita' dell'aria
 *   - celerita' del suono
 *   - ...
 * 
 * - Elemento di rotore:
 * classe virtuale che possiede alcuni dati topologici e geometrici di rotore
 * ed i metodi per calcolare grandezze utili agli elementi aerodinamici che
 * fanno parte di un rotore
 * 
 * - Elemento aerodinamico:
 * stazionario o quasi-stazionario, con metodo p-k di ordine 0, 1 o 2,
 * associato ad un corpo rigido, basato sulla strip theory.
 * 
 * - Elemento aerodinamico:
 * analogo al precedente, ma associato alla trave a Volumi Finiti a tre nodi
 * 
 * - Elemento aerodinamico instazionario:
 * in fase di sviluppo, modella dinamicamente alcuni stati a dare
 * il comportamento instazionario di una superficie aerodinamica
 * modellata con la strip theory
 *
 * - Elemento Aerodinamico modale:
 * in fase di sviluppo, contiene una rappresentazione agli stati delle forze
 * aerodinamiche generalizzate associate ai modi di vibrare e alle raffiche
 * deve essere associato ad un elemento modale su cui il modello si appoggia.
 *
 * -Elemento aerodinamico external:
 * in fase di sviluppo, elemento che permette di svolgere simulazioni 
 * integrate di interazioni fluido-struttura interfacciando MBDyn con 
 * codici aerodynamici esterni. L'elemento manda all'esterno le informazioni 
 * riguardo alla posizione dei nodi associati e ottiene i carichi che il codice 
 * aerodinamico esterno genera 
 * 
 * 
 */

#ifndef AERODYN__H
#define AERODYN__H

#include "aerodyn.h"

/* BasicAirProperties - begin */

class BasicAirProperties 
: virtual public Elem, public AirProperties {
protected:
	DriveOwner AirDensity;
	doublereal dSoundSpeed;
   
public:
	BasicAirProperties(const TplDriveCaller<Vec3>* pDC,
		DriveCaller *pRho, doublereal dSS, std::vector<Gust *>& g,
		const RigidBodyKinematics *pRBK,
		flag fOut);
   
	virtual ~BasicAirProperties(void);

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;
   
	/*
	 * Deprecated; use GetAirProps instead
	 */
	virtual doublereal dGetAirDensity(const Vec3& /* X */ ) const;
   
	virtual doublereal dGetAirPressure(const Vec3& /* X */ ) const;

	virtual doublereal dGetAirTemperature(const Vec3& /* X */ ) const;

	virtual doublereal dGetSoundSpeed(const Vec3& /* X */ ) const;

	/*
	 * End of deprecated; use GetAirProps instead
	 */

	virtual bool GetAirProps(const Vec3& X, doublereal& rho,
		doublereal& c, doublereal& p, doublereal& T) const;
};

/* BasicAirProperties - end */


/* StdAirProperties - begin */

class StdAirProperties 
: virtual public Elem, public AirProperties {
protected:
	doublereal PRef;
	DriveCaller *RhoRef;
	doublereal TRef;
	doublereal a;
	doublereal R;
	doublereal g0;
	doublereal z0;
	doublereal z1;
	doublereal z2;
   
public:
	StdAirProperties(const TplDriveCaller<Vec3>* pDC,
		doublereal PRef_, DriveCaller *RhoRef_, doublereal TRef_,
		doublereal a_, doublereal R_, doublereal g0_,
		doublereal z0_, doublereal z1_, doublereal z2_,
		std::vector<Gust *>& g,
		const RigidBodyKinematics *pRBK,
		flag fOut);
   
	virtual ~StdAirProperties(void);

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;
   
	/*
	 * Deprecated; use GetAirProps instead
	 */
	virtual doublereal dGetAirDensity(const Vec3& /* X */ ) const;
   
	virtual doublereal dGetAirPressure(const Vec3& /* X */ ) const;

	virtual doublereal dGetAirTemperature(const Vec3& /* X */ ) const;

	virtual doublereal dGetSoundSpeed(const Vec3& /* X */ ) const;

	/*
	 * End of deprecated; use GetAirProps instead
	 */

	virtual bool GetAirProps(const Vec3& X, doublereal& rho,
		doublereal& c, doublereal& p, doublereal& T) const;
};

/* StdAirProperties - end */

#endif /* AERODYN__H */

