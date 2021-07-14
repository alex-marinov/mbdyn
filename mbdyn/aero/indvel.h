/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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

#ifndef INDVEL_H
#define INDVEL_H

#include <cfloat>

#include "ac/pthread.h"
#ifdef USE_MPI
#include "ac/mpi.h"
#endif // USE_MPI

#include "aerodyn.h"
#include "strnode.h"
#include "resforces.h"

/* InducedVelocity - begin */

class InducedVelocity
: virtual public Elem {
public:
	enum Type {
		UNKNOWN = -1,

		// non-rotating...

		USER_DEFINED	= 0x01000000U,
		ROTOR		= 0x10000000U,

		// rotating...

		NO		= (0U | ROTOR),
		UNIFORM		= (1U | ROTOR),
		GLAUERT		= (2U | ROTOR),
		MANGLER		= (3U | ROTOR),
		DYNAMICINFLOW	= (4U | ROTOR),
		PETERS_HE	= (5U | ROTOR),

		CYCLOCOPTER	= (11U | ROTOR),

		LASTROTORTYPE
	};

public:
	class ErrInfiniteMeanInducedVelocity : public MBDynErrBase {
	public:
		ErrInfiniteMeanInducedVelocity(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

protected:
#ifdef USE_MPI
	// Communicator per il calcolo della trazione totale
	bool is_parallel;
	// Gruppo di macchine su cui si scambiano dati relativi al rotore
	MPI::Intracomm IndVelComm;
	int* pBlockLenght;
	// vettore di indirizzi di dati
	MPI::Aint* pDispl;
	// request per le comunicazioni Send Receive
	MPI::Request ReqV;
	// dimensioni vettore forze scambiato fra i processi
	integer iForcesVecDim;
	// vettori temporanei per scambio forze
	doublereal*  pTmpVecR;
	doublereal*  pTmpVecS;
	// datatype che contiene le posizioni
	// dei dati da scambiare ad ogni passo
	MPI::Datatype* pIndVelDataType;
#endif // USE_MPI

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
	mutable pthread_mutex_t forces_mutex;

	mutable pthread_mutex_t induced_velocity_mutex;
	mutable pthread_cond_t induced_velocity_cond;
	mutable bool bDone;

	void Wait(void) const;
	void Done(void) const;
#endif // USE_MULTITHREAD && MBDYN_X_MT_ASSRES

	const StructNode* pCraft;

	// force, couple and pole for resultants
	ExternResForces Res;
	// extra forces
	ResForceSet **ppRes;

public:
	InducedVelocity(unsigned int uL,
		const StructNode* pCraft,
		ResForceSet **ppres, flag fOut);
	virtual ~InducedVelocity(void);

	// funzioni di servizio

	/* Metodi per l'estrazione di dati "privati".
	 * Si suppone che l'estrattore li sappia interpretare.
	 * Come default non ci sono dati privati estraibili */
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;

	// Return "true" if sectional forces are needed
	virtual bool bSectionalForces(void) const;

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

	// ritorna il numero di Dofs per gli elementi che sono anche DofOwners
	virtual unsigned int iGetNumDof(void) const {
		return 0;
	};

	// Type
	virtual InducedVelocity::Type GetInducedVelocityType(void) const = 0;

	virtual inline const Vec3& GetXCurr(void) const {
		return pCraft->GetXCurr();
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

	// Elements use this function to pass induced velocity models
	// their contribution to forces and moments at a specific location X.
	// Each element can call this function multiple times to provide
	// multiple contributions in one iteration.
	// Elements must not call this function when bSectionalForces() 
	// returns "true"
	virtual void AddForce(const Elem *pEl, const StructNode *pNode,
		const Vec3& F, const Vec3& M, const Vec3& X);

	// Elements use this function to pass induced velocity models
	// their contribution to sectional forces and moments at a specific
	// location X.
	// Each element can call this function multiple times to provide
	// multiple contributions in one iteration.
	// Elements should not call this function when bSectionalForces() 
	// returns "false"
	virtual void AddSectionalForce(Elem::Type type,
		const Elem *pEl, unsigned uPnt,
		const Vec3& F, const Vec3& M, doublereal dW,
		const Vec3& X, const Mat3x3& R,
		const Vec3& V, const Vec3& W);

	virtual void ResetForce(void);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const = 0;

	// Dimensioni del workspace
	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 0;
		*piNumCols = 0;
	};

	// assemblaggio jacobiano
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// Relativo ai ...WithDofs
	virtual void SetInitialValue(VectorHandler& /* X */ ) {
		NO_OP;
	};

	// Relativo ai ...WithDofs
	virtual void
	SetValue(DataManager *pDM,
		VectorHandler& /* X */ , VectorHandler& /* XP */ ,
		SimulationEntity::Hints *ph = 0)
	{
		NO_OP;
	};

	// *******PER IL SOLUTORE PARALLELO********
	// Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	// utile per l'assemblaggio della matrice di connessione fra i dofs
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(1);
		connectedNodes[0] = pCraft;
	};
	// ************************************************

#ifdef USE_MPI
	void ExchangeLoads(flag fWhat);	// ExchangeLoads
	void InitializeIndVelComm(MPI::Intracomm* Rot);	// InitializeIndVelComm
	void ExchangeVelocity(void);
#endif /* USE_MPI */
};

/* InducedVelocity - end */

/* InducedVelocityElem - begin */

class InducedVelocityElem
: virtual public Elem, public AerodynamicElem, public InducedVelocity {
public:
	InducedVelocityElem(unsigned int uL, const DofOwner* pDO,
		const StructNode* pCraft,
		ResForceSet **ppres, flag fOut);
	virtual ~InducedVelocityElem(void);

	// funzioni di servizio

	// Tipo dell'elemento (usato per debug ecc.)
	virtual Elem::Type GetElemType(void) const;
	virtual AerodynamicElem::Type GetAerodynamicElemType(void) const;

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

	// esegue operazioni sui dof di proprieta' dell'elemento
	virtual DofOrder::Order GetDofType(unsigned int i) const {
		ASSERT(i >= 0 && i < this->iGetNumDof());
		return DofOrder::DIFFERENTIAL;
	};
};

/* InducedVelocityElem - end */

#endif // INDVEL_H

