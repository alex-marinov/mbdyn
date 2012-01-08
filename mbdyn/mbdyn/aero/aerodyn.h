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

#ifndef AERODYN_H
#define AERODYN_H

#include "elem.h"
#include "tpldrive.h"
#include "rbk.h"
#include "gust.h"

extern const char* psAeroNames[];

/* AirProperties - begin */

class AirProperties 
: virtual public Elem, public InitialAssemblyElem, public TplDriveOwner<Vec3> {
protected:
	mutable Vec3 Velocity;
	std::vector<const Gust *> gust;

	// rigid body kinematics
	const RigidBodyKinematics *pRBK;
   
public:
	AirProperties(const TplDriveCaller<Vec3>* pDC,
		std::vector<const Gust *>& pg,
		const RigidBodyKinematics *pRBK,
		flag fOut);
	virtual ~AirProperties(void);

	virtual void AddGust(const Gust *pG);

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Tipo dell'elemento (usato per debug ecc.) */
	virtual Elem::Type GetElemType(void) const { 
		return Elem::AIRPROPERTIES; 
	};
      
	/* funzioni di servizio */

	/* Il metodo iGetNumDof() serve a ritornare il numero di gradi di liberta'
	 * propri che l'elemento definisce. Non e' virtuale in quanto serve a 
	 * ritornare 0 per gli elementi che non possiedono gradi di liberta'.
	 * Viene usato nella costruzione dei DofOwner e quindi deve essere 
	 * indipendente da essi. In genere non comporta overhead in quanto il 
	 * numero di dof aggiunti da un tipo e' una costante e non richede dati 
	 * propri.
	 * Il metodo pGetDofOwner() ritorna il puntatore al DofOwner dell'oggetto.
	 * E' usato da tutti quelli che agiscono direttamente sui DofOwner.
	 * Non e' virtuale in quanto ritorna NULL per tutti i tipi che non hanno
	 * dof propri.
	 * Il metodo GetDofType() ritorna, per ogni dof dell'elemento, l'ordine.
	 * E' usato per completare i singoli Dof relativi all'elemento.
	 */
   
	/* funzioni proprie */
   
	/* Dimensioni del workspace */
	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat, doublereal /* dCoef */ , 
			const VectorHandler& /* XCurr */ ,
			const VectorHandler& /* XPrimeCurr */ );

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec, doublereal /* dCoef */ ,
			const VectorHandler& /* XCurr */ , 
			const VectorHandler& /* XPrimeCurr */ );
   
	/* Numero di GDL iniziali */
	virtual unsigned int iGetInitialNumDof(void) const;
     
	/* Dimensioni initiali del workspace */
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
			const VectorHandler& /* XCurr */ );

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

	virtual void Output(OutputHandler&) const;

	/*
	 * Deprecated; use GetAirProps instead
	 */
	virtual Vec3 GetVelocity(const Vec3& /* X */ ) const;
	virtual bool GetVelocity(const Vec3& /* X */ , Vec3& V) const;
	virtual doublereal dGetAirDensity(const Vec3& /* X */ ) const = 0;
	virtual doublereal dGetAirPressure(const Vec3& /* X */ ) const = 0;
	virtual doublereal dGetAirTemperature(const Vec3& /* X */ ) const = 0;
	virtual doublereal dGetSoundSpeed(const Vec3& /* X */ ) const = 0;

	/*
	 * End of deprecated; use GetAirProps instead
	 */

	virtual bool
	GetAirProps(const Vec3& X, doublereal& rho, doublereal& c,
			doublereal& p, doublereal& T) const = 0;

	/* *******PER IL SOLUTORE BLOCK JACOBI-BROYDEN******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual int GetNumConnectedNodes(void) const {
		return 0;
	};

	/* Dati privati */
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;
};

class DataManager;

extern Elem *
ReadAirProperties(DataManager* pDM, MBDynParser& HP);

/* AirProperties - end */


/* AirPropOwner - begin */

class AirPropOwner {
protected:
	const AirProperties* pAirProperties;
public:
	AirPropOwner(void);
	virtual ~AirPropOwner(void);
   
	virtual void PutAirProperties(const AirProperties* pAP);

	/*
	 * Deprecated; use GetAirProps instead
	 */
	virtual flag fGetAirVelocity(Vec3& Velocity, const Vec3& X) const;
	virtual doublereal dGetAirDensity(const Vec3& X) const;
	virtual doublereal dGetAirPressure(const Vec3& X) const;
	virtual doublereal dGetAirTemperature(const Vec3& X) const;
	virtual doublereal dGetSoundSpeed(const Vec3& X) const;

	/*
	 * End of deprecated; use GetAirProps instead
	 */

	virtual bool
	GetAirProps(const Vec3& X, doublereal& rho, doublereal& c,
			doublereal& p, doublereal& T) const;
};

/* AirPropOwner - end */


/* AerodynamicElem - begin */

class AerodynamicElem :
	virtual public Elem,
	public ElemWithDofs,
	public AirPropOwner
{
public:
	/* Tipi di elementi aerodinamici */
	enum Type {
		UNKNOWN = -1,

		INDUCEDVELOCITY = 0,
		AEROMODAL,
		AERODYNAMICBODY,
		AERODYNAMICBEAM,
		AERODYNAMICEXTERNAL,
		AERODYNAMICEXTERNALMODAL,
		AERODYNAMICLOADABLE,
		AIRCRAFTINSTRUMENTS,
		GENERICFORCE,
		
		LASTAEROTYPE
	};

protected:
 
public:
	AerodynamicElem(unsigned int uL, const DofOwner *pDO, flag fOut);
	virtual ~AerodynamicElem(void);

	/* Tipo di elemento aerodinamico */
	virtual AerodynamicElem::Type GetAerodynamicElemType(void) const = 0;

	virtual bool NeedsAirProperties(void) const;
	virtual const InducedVelocity *pGetInducedVelocity(void) const;
};

/* AerodynamicElem - end */

#endif // AERODYN_H

