/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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

/* Cyclocopter inflow models */

#ifndef CYCLOCOPTER_H
#define CYCLOCOPTER_H

#include "indvel.h"

/* CyclocopterInflow - begin */

/* Classe base per i modelli di influsso del ciclocottero:
 * tutti gli altri modelli sono ottenuti ereditando questa 
 * classe
*/

class CyclocopterInflow
: virtual public Elem, public InducedVelocityElem {
protected:
	const StructNode* pRotor;

	unsigned int iFlagAverage;			// == 1: usa la media delle forze sul giro
					//	per calcolare la velocità indotta
					//	per il giro successivo
					// == 0: usa il valore istantaneo (eventualmente
					// 	filtrato)


	doublereal dRadius;		// Rotor radius
	doublereal dSpan;		// Blade length
	doublereal dArea;		// Cylinder longitudinal area

	DriveOwner Weight;
	doublereal dWeight;

	Mat3x3 RRot;
	/* dati che servono per il calcolo dell'output */
	Mat3x3 RRotorTranspose;
	doublereal dUindMean;

	// Coefficienti del filtro di butterworth del second'ordine
	doublereal a1, a2, b0, b1, b2;	

public:
	CyclocopterInflow(unsigned int uL, const DofOwner* pDO,
		const StructNode* pC, const Mat3x3& rrot,
		const StructNode* pR, ResForceSet **ppres, 
		flag fOut);
	virtual ~CyclocopterInflow(void);

	virtual InducedVelocity::Type GetInducedVelocityType(void) const;

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	// Relativo ai ...WithDofs
	virtual void SetInitialValue(VectorHandler& X);

	// *******PER IL SOLUTORE PARALLELO********
	// Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	// utile per l'assemblaggio della matrice di connessione fra i dofs
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	// ************************************************
	
	// Coefficienti del filtro di butterworth del second'ordine
	void 
	SetFilterCoefficients( const doublereal dOmegaFilter, 
		const doublereal dDeltaT );
};

/* CyclocopterInflow - end */

/* CyclocopterNoInflow - begin */

class CyclocopterNoInflow
: virtual public Elem, public CyclocopterInflow {
public:
	CyclocopterNoInflow(unsigned int uL, const DofOwner* pDO,
		const StructNode* pC, const Mat3x3& rrot,
		const StructNode* pR, ResForceSet **ppres, 
		flag fOut);
	virtual ~CyclocopterNoInflow(void);

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// Restituisce ad un elemento la velocità indotta
	// in base alla posizione azimuthale (usata per
	// iterare tra calcolo della valocità indotta e 
	// calcolo delle forze aerodinamiche - solo
	// per il modello di influsso KARI per il ciclocottero)
	virtual void GetInducedVelocityIter(const Vec3& X, const Vec3& T, doublereal *UindM, doublereal dTn0, doublereal dTn_dUindM) {};

	// Restituisce la velocità indotta dalla metà superiore del rotore
	// ( solo per KARI ciclocottero) 
	virtual doublereal GetW( const Vec3& X) const {
		return 0;
	}
	virtual doublereal GetPsi( const Vec3& X) const {
		return 0;
	}
	virtual Mat3x3 GetRRotor( const Vec3& X) const {
		return Zero3x3;
	}

};

/* CyclocopterNoInflow - end */

/* CyclocopterUniform1D - begin */

/*

From Moble: the induced velocity is opposite to the
the force generated by the rotor in the direction 3
of the rotor reference (the direction 1 of this 
reference must be align with the rotation axis of 
the rotor). It should make sense just for a multiblade
rotor (not for the one-blade rotor) when the generated 
force is mainly in one direction.

*/

class CyclocopterUniform1D
: virtual public Elem, public CyclocopterInflow {
protected:

	Vec3 RRot3;
	Mat3x3 RRotor;

	mutable doublereal dUindMeanPrev;

	bool bFlagIsFirstBlade;

	doublereal dAzimuth, dAzimuthPrev;

	doublereal dTz, dTzMean;
	Vec3 F, FMean, FMeanOut;

	unsigned int iStepCounter;

	/* dati per il filtraggio delle forze */
	doublereal Uk, Uk_1, Uk_2, Yk, Yk_1, Yk_2;

public:
	CyclocopterUniform1D(unsigned int uL, const DofOwner* pDO,
		const StructNode* pC, const Mat3x3& rrot,
		const StructNode* pR, ResForceSet **ppres, 
		const unsigned int& iFlagAve, const doublereal& dR,
		const doublereal& dL, const doublereal& dOmegaFilter,
		const doublereal& dDeltaT, DriveCaller *pdW,
		flag fOut);
	virtual ~CyclocopterUniform1D(void);

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// Restituisce ad un elemento la velocità indotta
	// in base alla posizione azimuthale (usata per
	// iterare tra calcolo della valocità indotta e 
	// calcolo delle forze aerodinamiche - solo
	// per il modello di influsso KARI per il ciclocottero)
	virtual void GetInducedVelocityIter(const Vec3& X, const Vec3& T, doublereal *UindM, doublereal dTn0, doublereal dTn_dUindM) {};

	// Restituisce la velocità indotta dalla metà superiore del rotore
	// ( solo per KARI ciclocottero) 
	virtual doublereal GetW( const Vec3& X) const {
		return 0;
	}
	virtual doublereal GetPsi( const Vec3& X) const {
		return 0;
	}
	virtual Mat3x3 GetRRotor( const Vec3& X) const {
		return Zero3x3;
	}
};

/* CyclocopterUnifor1D - end */

/* CyclocopterUniform2D - begin */

/*

Uniform inflow: the indiced velocity is opposite to
the force generated by the rotor in the plane 
perdendicular to the rotor rotation axis. The
rotor reference must have the direction 1 aligned 
with the rotor rotation axis!

*/

class CyclocopterUniform2D
: virtual public Elem, public CyclocopterInflow {
protected:

	Mat3x3 RRotor;
	Vec3 dUind;
	mutable Vec3 dUindPrev;

	bool bFlagIsFirstBlade;

	doublereal dAzimuth, dAzimuthPrev;

	Vec3 F, FMean, FMeanOut;

	unsigned int iStepCounter;

	/* dati per il filtraggio delle forze */
	Vec3 Uk, Uk_1, Uk_2, Yk, Yk_1, Yk_2;

public:
	CyclocopterUniform2D(unsigned int uL, const DofOwner* pDO,
		const StructNode* pC, const Mat3x3& rrot,
		const StructNode* pR, ResForceSet **ppres, 
		const unsigned int& iFlagAve, const doublereal& dR,
		const doublereal& dL, const doublereal& dOmegaFilter,
		const doublereal& dDeltaT, DriveCaller *pdW,
		flag fOut);
	virtual ~CyclocopterUniform2D(void);

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// Restituisce ad un elemento la velocità indotta
	// in base alla posizione azimuthale (usata per
	// iterare tra calcolo della valocità indotta e 
	// calcolo delle forze aerodinamiche - solo
	// per il modello di influsso KARI per il ciclocottero)
	virtual void GetInducedVelocityIter(const Vec3& X, const Vec3& T, doublereal *UindM, doublereal dTn0, doublereal dTn_dUindM) {};

	// Restituisce la velocità indotta dalla metà superiore del rotore
	// ( solo per KARI ciclocottero) 
	virtual doublereal GetW( const Vec3& X) const {
		return 0;
	}
	virtual doublereal GetPsi( const Vec3& X) const {
		return 0;
	}
	virtual Mat3x3 GetRRotor( const Vec3& X) const {
		return Zero3x3;
	}

};

/* CyclocopterUnifor2D - end */

/* CyclocopterPolimi - begin */

/*

The indiced velocity is opposite to
the force generated by the rotor in the plane 
perdendicular to the rotor rotation axis. The
rotor reference must have the direction 1 aligned 
with the rotor rotation axis!
The inflow velocity distribution is constan along
the cylinder span and is function of r:
Vi(r) = Kc cos( (pi/2)*(r/R) ) + Ks sin( pi*(r/R) )
where the cofficients Kc and Ks are based on the
momentum theory

*/

class CyclocopterPolimi
: virtual public Elem, public CyclocopterInflow {
protected:

	Mat3x3 RRotor;
	Vec3 dUind;
	mutable Vec3 dUindPrev;

	doublereal dXi;

	bool bFlagIsFirstBlade;

	doublereal dAzimuth, dAzimuthPrev;

	Vec3 F, FMean, FMeanOut;

	unsigned int iStepCounter;

	/* dati per il filtraggio delle forze */
	Vec3 Uk, Uk_1, Uk_2, Yk, Yk_1, Yk_2;

	unsigned int iCounter;
	unsigned int iRotationCounter;

	doublereal dpPrev, dp;
	
	bool flag_print;

public:
	CyclocopterPolimi(unsigned int uL, const DofOwner* pDO,
		const StructNode* pC, const Mat3x3& rrot,
		const StructNode* pR, ResForceSet **ppres, 
		const unsigned int& iFlagAve, const doublereal& dR,
		const doublereal& dL, const doublereal& dOmegaFilter,
		const doublereal& dDeltaT, DriveCaller *pdW,
		flag fOut);
	virtual ~CyclocopterPolimi(void);

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// Restituisce ad un elemento la velocità indotta
	// in base alla posizione azimuthale (usata per
	// iterare tra calcolo della valocità indotta e 
	// calcolo delle forze aerodinamiche - solo
	// per il modello di influsso KARI per il ciclocottero)
	virtual void GetInducedVelocityIter(const Vec3& X, const Vec3& T, doublereal *UindM, doublereal dTn0, doublereal dTn_dUindM) {};

	// Restituisce la velocità indotta dalla metà superiore del rotore
	// ( solo per KARI ciclocottero) 
	virtual doublereal GetW( const Vec3& X) const {
		return 0;
	}
	virtual doublereal GetPsi( const Vec3& X) const {
		return 0;
	}
	virtual Mat3x3 GetRRotor( const Vec3& X) const {
		return Zero3x3;
	}

};

/* CyclocopterPolimi - end */

/* CyclocopterKARI - begin */

/*

From:

"A New VTOL UAV Cyclocopter with Cycloidal Blades System",

Chul Yong Yun, Illkyung Park,
Ho Yong Lee, Jai Sang Jung, In Seong Hwang,
Seung Jo Kim,
Sung Nam Jung

Presented at the American Helicopter Society 60th Annual Forum,
Baltimore, MD, June 7-10, 2004.
*/

class CyclocopterKARI
: virtual public Elem, public InducedVelocityElem {
protected:
	const StructNode* pRotor;
	Mat3x3 RRot;

public:
	CyclocopterKARI(unsigned int uL, const DofOwner* pDO,
		const StructNode* pC, const Mat3x3& rrot,
		const StructNode* pR, ResForceSet **ppres, flag fOut);
	virtual ~CyclocopterKARI(void);

	virtual InducedVelocity::Type GetInducedVelocityType(void) const;

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	// Relativo ai ...WithDofs
	virtual void SetInitialValue(VectorHandler& X);

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

#if 0
	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);
#endif

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// *******PER IL SOLUTORE PARALLELO********
	// Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	// utile per l'assemblaggio della matrice di connessione fra i dofs
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	// ************************************************
};

/* CyclocopterKARI - end */


class DataManager;
class MBDynParser;

extern Elem*
ReadCyclocopter(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner* pDO, 
	unsigned int uLabel,
	const StructNode* pC,
	const Mat3x3& rrot,
	const StructNode* pR);

#endif /* CYCLOCOPTER_H */

