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

/* Elementi di rotore */

#ifndef ROTOR_H
#define ROTOR_H

#include "indvel.h"

extern const char* psRotorNames[];

/* Rotor - begin */

class Rotor
: virtual public Elem, public InducedVelocityElem {
protected:
	const StructNode* pRotor;
	const StructNode* pGround;

	doublereal dOmegaRef;		// Velocita' di rotazione di riferimento

	doublereal dRadius;		// Raggio del rotore
	doublereal dVTipRef;		// Raggio del rotore
	doublereal dArea;		// Area del disco
	doublereal dUMean;		// Velocita' indotta media
	doublereal dUMeanRef;		// Velocita' indotta media (nominale)
	mutable doublereal dUMeanPrev;	// Vel. indotta media al passo prec.

	// iterations for dUMeanRef
	unsigned int iMaxIter;
	unsigned int iCurrIter;
	doublereal dTolerance;
	doublereal dEta;
	bool bUMeanRefConverged;

	DriveOwner Weight;
	// Peso della velocita' indotta media
	// (peso della V al passo precedente, def = 0.)
	doublereal dWeight;
	// Correzione H (scala la velocita' indotta)
	doublereal dHoverCorrection;
	// FF
	doublereal dForwardFlightCorrection;

	// Trasposta della matrice di rotazione rotore
	Mat3x3 RRotTranspose;
	Mat3x3 RRot;
	Vec3 RRot3;		// Direzione dell'asse del rotore
	Vec3 VCraft;		// Velocita' di traslazione del velivolo
	doublereal dPsi0;	// Angolo di azimuth del rotore
	doublereal dSinAlphad;	// Angolo di incidenza del disco?
	doublereal dCosAlphad;	// ???
	doublereal dMu;		// Parametro di'avanzamento
	doublereal dLambda;	// Parametro di influsso
	doublereal dChi;	// ???

	doublereal dVelocity;	// Velocita' di riferimento
	doublereal dOmega;	// Velocita' di rotazione di riferimento



	// temporaneo
	mutable int iNumSteps;

	// Questa funzione non viene usata, da Rotor, ma da altre classi derivate
	// con modellazioni piu' sofisticate della velocita' indotta
	virtual doublereal dGetPsi(const Vec3& X) const;

	// Calcola la distanza di un punto dall'asse di rotazione
	// in coordinate adimensionali
	virtual doublereal dGetPos(const Vec3& X) const;

	// Combina i due ...
	virtual void GetPos(const Vec3& X, doublereal& dr, doublereal& dp) const;

	// Calcola la velocita' di traslazione del rotore
	virtual void InitParam(bool bComputeMeanInducedVelocity = true);

	virtual void Init(const StructNode* pC, const Mat3x3& rrot,
		const StructNode* pR, const StructNode* pG,
		ResForceSet **ppres,
		const doublereal& dR,
		unsigned int iMaxIt,
		const doublereal& dTol,
		const doublereal& dE,
		flag fOut);

public:
	Rotor(unsigned int uL, const DofOwner* pDO);
	Rotor(unsigned int uL, const DofOwner* pDO,
		const StructNode* pC, const Mat3x3& rrot,
		const StructNode* pR, const StructNode* pG,
		ResForceSet **ppres,
		const doublereal& dR,
		unsigned int iMaxIt,
		const doublereal& dTol,
		const doublereal& dE,
		flag fOut);
	virtual ~Rotor(void);

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	// Relativo ai ...WithDofs
	virtual void SetInitialValue(VectorHandler& /* X */ ) {
		NO_OP;
	};

	virtual inline const Vec3& GetXCurr(void) const {
		return pRotor->GetXCurr();
	};

	// accesso a dati
	virtual inline doublereal dGetOmega(void) const {
#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
		Wait();
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES
		return dOmega;
	};

	virtual inline doublereal dGetRadius(void) const {
#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
		Wait();
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES
		return dRadius;
	};

	virtual inline doublereal dGetMu(void) const {
#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
		Wait();
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES
		return dMu;
	};

	virtual inline const Vec3& GetForces(void) const {
#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
		Wait();
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES
		return Res.Force();
	};

	virtual inline const Vec3& GetMoments(void) const {
#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
		Wait();
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES
		return Res.Moment();
	};

	// *******PER IL SOLUTORE PARALLELO********
	// Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	// utile per l'assemblaggio della matrice di connessione fra i dofs
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		if (pGround != 0) {
			connectedNodes.resize(3);
			connectedNodes[2] = pGround;
		} else {
			connectedNodes.resize(2);
		}

		connectedNodes[0] = pCraft;
		connectedNodes[1] = pRotor;
	};
	// ************************************************
};

/* Rotor - end */


/* NoRotor - begin */

class NoRotor : virtual public Elem, public Rotor {
protected:
	virtual void Init(const StructNode* pCraft,
		const Mat3x3& rrot,
		const StructNode* pRotor,
		ResForceSet **ppres,
		const doublereal& dR,
		flag fOut);

public:
	NoRotor(unsigned int uLabel, const DofOwner* pDO);
	NoRotor(unsigned int uLabel,
		const DofOwner* pDO,
		const StructNode* pCraft,
		const Mat3x3& rrot,
		const StructNode* pRotor,
		ResForceSet **ppres,
		const doublereal& dR,
		flag fOut);
	virtual ~NoRotor(void);

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual InducedVelocity::Type GetInducedVelocityType(void) const {
		return InducedVelocity::NO;
	};

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;
};

/* NoRotor - end */


/* UniformRotor - begin */

class UniformRotor : virtual public Elem, public Rotor {
protected:
	virtual void Init(const StructNode* pCraft,
	   	const Mat3x3& rrot,
		const StructNode* pRotor,
		const StructNode* pGround,
		ResForceSet **ppres,
		const doublereal& dOR,
		const doublereal& dR,
		DriveCaller *pdW,
		unsigned int iMaxIt,
		const doublereal& dTol,
		const doublereal& dE,
		const doublereal& dCH,
		const doublereal& dCFF,
		flag fOut);

public:
	UniformRotor(unsigned int uLabel, const DofOwner* pDO);
	UniformRotor(unsigned int uLabel,
		const DofOwner* pDO,
		const StructNode* pCraft,
	   	const Mat3x3& rrot,
		const StructNode* pRotor,
		const StructNode* pGround,
		ResForceSet **ppres,
		const doublereal& dOR,
		const doublereal& dR,
		DriveCaller *pdW,
		unsigned int iMaxIt,
		const doublereal& dTol,
		const doublereal& dE,
		const doublereal& dCH,
		const doublereal& dCFF,
		flag fOut);
	virtual ~UniformRotor(void);

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	InducedVelocity::Type GetInducedVelocityType(void) const {
		return InducedVelocity::UNIFORM;
	};

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;
};

class UniformRotor2 : virtual public Elem, public UniformRotor {
protected:
#if 0 // not needed because identical to UniformRotor::Init()
	virtual void Init(const StructNode* pCraft,
	   	const Mat3x3& rrot,
		const StructNode* pRotor,
		const StructNode* pGround,
		ResForceSet **ppres,
		const doublereal& dOR,
		const doublereal& dR,
		DriveCaller *pdW,
		unsigned int iMaxIt,
		const doublereal& dTol,
		const doublereal& dE,
		const doublereal& dCH,
		const doublereal& dCFF,
		flag fOut);
#endif

public:
	UniformRotor2(unsigned int uLabel, const DofOwner* pDO);
	UniformRotor2(unsigned int uLabel,
		const DofOwner* pDO,
		const StructNode* pCraft,
	   	const Mat3x3& rrot,
		const StructNode* pRotor,
		const StructNode* pGround,
		ResForceSet **ppres,
		const doublereal& dOR,
		const doublereal& dR,
		DriveCaller *pdW,
		unsigned int iMaxIt,
		const doublereal& dTol,
		const doublereal& dE,
		const doublereal& dCH,
		const doublereal& dCFF,
		flag fOut);
	virtual ~UniformRotor2(void);

	virtual bool bSectionalForces(void) const;
	virtual void AddSectionalForce(Elem::Type type,
		const Elem *pEl, unsigned uPnt,
		const Vec3& F, const Vec3& M, doublereal dW,
		const Vec3& X, const Mat3x3& R,
		const Vec3& V, const Vec3& W);
};

/* UniformRotor - end */


/* GlauertRotor - begin */

class GlauertRotor : virtual public Elem, public Rotor {
public:
	enum Type {
		UNKNOWN = -1,

		// from (?)
		GLAUERT,

		// from J. Gordon Leishman, Principles of Helicopter Aerodynamics, 2000
		COLEMAN_ET_AL,
		DREES_1,
		PAYNE,
		WHITE_AND_BLAKE,
		PITT_AND_PETERS,
		HOWLETT,

		// from Massimo Gennaretti, Roma Tre
		DREES_2,

		LAST
	};

protected:
	Type gtype;

	virtual void Init(const StructNode* pCraft,
		const Mat3x3& rrot,
		const StructNode* pRotor,
		const StructNode* pGround,
		ResForceSet **ppres,
		const doublereal& dOR,
		const doublereal& dR,
		DriveCaller *pdW,
		unsigned int iMaxIt,
		const doublereal& dTol,
		const doublereal& dE,
		const doublereal& dCH,
		const doublereal& dCFF,
		GlauertRotor::Type gtype,
		flag fOut);

public:
	GlauertRotor(unsigned int uLabel, const DofOwner* pDO);
	GlauertRotor(unsigned int uLabel,
		const DofOwner* pDO,
		const StructNode* pCraft,
		const Mat3x3& rrot,
		const StructNode* pRotor,
		const StructNode* pGround,
		ResForceSet **ppres,
		const doublereal& dOR,
		const doublereal& dR,
		DriveCaller *pdW,
		unsigned int iMaxIt,
		const doublereal& dTol,
		const doublereal& dE,
		const doublereal& dCH,
		const doublereal& dCFF,
		GlauertRotor::Type gtype,
		flag fOut);
	virtual ~GlauertRotor(void);

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	InducedVelocity::Type GetInducedVelocityType(void) const {
		return InducedVelocity::GLAUERT;
	};

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;
};

/* GlauertRotor - end */


/* ManglerRotor - begin */

class ManglerRotor : virtual public Elem, public Rotor {
protected:
	virtual void Init(const StructNode* pCraft,
		const Mat3x3& rrot,
		const StructNode* pRotor,
		const StructNode* pGround,
		ResForceSet **ppres,
		const doublereal& dOR,
		const doublereal& dR,
		DriveCaller *pdW,
		unsigned int iMaxIt,
		const doublereal& dTol,
		const doublereal& dE,
		const doublereal& dCH,
		const doublereal& dCFF,
		flag fOut);

public:
	ManglerRotor(unsigned int uLabel, const DofOwner* pDO);
	ManglerRotor(unsigned int uLabel,
		const DofOwner* pDO,
		const StructNode* pCraft,
		const Mat3x3& rrot,
		const StructNode* pRotor,
		const StructNode* pGround,
		ResForceSet **ppres,
		const doublereal& dOR,
		const doublereal& dR,
		DriveCaller *pdW,
		unsigned int iMaxIt,
		const doublereal& dTol,
		const doublereal& dE,
		const doublereal& dCH,
		const doublereal& dCFF,
		flag fOut);
	virtual ~ManglerRotor(void);

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	InducedVelocity::Type GetInducedVelocityType(void) const {
		return InducedVelocity::MANGLER;
	};

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;
};

/* ManglerRotor - end */


/* DynamicInflowRotor - begin */

/*
 * Based on the 3 state dynamic inflow by Pitt-Peters:
 *
 * D. M. Pitt,  D. A. Peters,
 * "Theoretical Prediction of Dynamic Inflow Derivatives",
 * Vertica, Vol. 5, pp.21-34, 1981
 *
 * as discussed by Chen in:
 *
 * R. T. N. Chen,
 * "A Survey of Nonuniform Inflow Models for Rotorcraft
 * Flight Dynamics and Control Applications"
 * Vertica, Vol 14, No. 2, pp.147-184, 1990
 */

class DynamicInflowRotor : virtual public Elem, public Rotor {
protected:
	doublereal dVConst;
	doublereal dVSine;
	doublereal dVCosine;

	doublereal dL11;
	doublereal dL13;
	doublereal dL22;
	doublereal dL31;
	doublereal dL33;

	virtual void Init(const StructNode* pCraft,
		const Mat3x3& rrot,
		const StructNode* pRotor,
      		const StructNode* pGround,
		ResForceSet **ppres,
		const doublereal& dOR,
		const doublereal& dR,
		unsigned int iMaxIt,
		const doublereal& dTol,
		const doublereal& dE,
		const doublereal& dCH,
		const doublereal& dCFF,
		const doublereal& dVConstTmp,
		const doublereal& dVSineTmp,
		const doublereal& dVCosineTmp,
		flag fOut);

public:
	DynamicInflowRotor(unsigned int uLabel, const DofOwner* pDO);
	DynamicInflowRotor(unsigned int uLabel,
		const DofOwner* pDO,
		const StructNode* pCraft,
		const Mat3x3& rrot,
		const StructNode* pRotor,
      		const StructNode* pGround,
		ResForceSet **ppres,
		const doublereal& dOR,
		const doublereal& dR,
		unsigned int iMaxIt,
		const doublereal& dTol,
		const doublereal& dE,
		const doublereal& dCH,
		const doublereal& dCFF,
		const doublereal& dVConstTmp,
		const doublereal& dVSineTmp,
		const doublereal& dVCosineTmp,
		flag fOut);
	virtual ~DynamicInflowRotor(void);

	// ritorna il numero di Dofs per gli elementi che sono anche DofOwner
	virtual unsigned int iGetNumDof(void) const {
		return 3;
	};

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// Dimensioni del workspace
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 3;
		*piNumCols = 3;
	};

	// assemblaggio jacobiano
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	// Relativo ai ...WithDofs
	virtual void SetInitialValue(VectorHandler& X);

	// Relativo ai ...WithDofs
	virtual void
	SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	InducedVelocity::Type GetInducedVelocityType(void) const {
		return InducedVelocity::DYNAMICINFLOW;
	};

	// Somma alla trazione il contributo di un elemento */
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;
};

/* DynamicInflowRotor - end */

class DataManager;
class MBDynParser;

extern Elem*
ReadRotor(DataManager* pDM,
	MBDynParser& HP,
	const DofOwner* pDO,
	unsigned int uLabel);

#endif /* ROTOR_H */

