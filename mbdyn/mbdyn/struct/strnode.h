/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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

#ifndef STRNODE_H
#define STRNODE_H

#include "node.h"
#include "matvec3.h"
#include "rbk.h"
#include "invdyn.h"

extern const char* psStructNodeNames[];


/* StructNode - begin */

/* Nodo strutturale: possiede i gradi di liberta' di:
 *  - spostamento assoluto,
 *  - parametri di rotazione incrementali,
 * inoltre, se dinamico, i gdl di:
 *  - quantita' di moto,
 *  - momento della quantita' di moto rispetto al polo mobile.
 * Il nodo di per se' non ha caratteristiche inerziali, che gli vengono date
 * dagli elementi ad esso collegati. In particolare, gli elementi Body, corpo
 * rigido, sono responsabili dell'attribuzione di inerzia ai nodi. Altri
 * contributi possono giungere da travi con matrice di inerzia consistente
 * (non ancora implementate). */

const int iStructNodeInitialDofNumber = 12;

class StructNode : public Node, public RigidBodyKinematics {
public:
	class ErrGeneric : public MBDynErrBase {
  	public:
 		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

	enum Type {
		UNKNOWN = -1,

		DYNAMIC = 0,
		STATIC,
		MODAL,
		DUMMY,

		LASTSTRUCTNODETYPE
	};

	enum Output {
		OUTPUT_ACCELERATIONS = (ToBeOutput::OUTPUT_PRIVATE << 0),
		OUTPUT_INERTIA = (ToBeOutput::OUTPUT_PRIVATE << 1)
	};

protected:
	mutable Mat3x3 RPrev;   /* Matrice di rotazione da zero al passo prec. */
	Mat3x3 RRef;            /* Matrice di rotazione predetta al passo corr. */
	mutable Mat3x3 RCurr;   /* Matrice di rotazione all'iterazione corrente */

	mutable Vec3 gRef;
	mutable Vec3 gCurr;     /* parametri e derivate correnti */
	mutable Vec3 gPRef;
	mutable Vec3 gPCurr;

	/* Valgono le relazioni:
	 *        RCurr = RDelta*RRef                (1)
	 *        RDelta = RCurr*RRef^T              (2)
	 * In base a questo, dal momento che la matrice RDelta e' richiesta
	 * solo in fase di aggiornamento ed e' usata dal nodo stesso, conviene
	 * non conservarla e calcolarla in base alla relazione (2).
	 */

	mutable Vec3 XPrev;   /* Posizione al passo precedente */
	mutable Vec3 XCurr;   /* Posizione corrente */

	mutable Vec3 VPrev;   /* Velocita' al passo precedente */
	mutable Vec3 VCurr;   /* Velocita' corrente */

	mutable Vec3 WPrev;   /* Velocita' angolare al passo precedente */
	Vec3 WRef;            /* Velocita' angolare predetta al passo corrente */
	mutable Vec3 WCurr;   /* Velocita' angolare corrente */

	mutable Vec3 XPPCurr;   /* Accelerazione lineare  corrente */
	mutable Vec3 WPCurr;    /* Accelerazione angolare corrente */
	mutable Vec3 XPPPrev;   /* Accelerazione lineare  al passo prec. */
	mutable Vec3 WPPrev;    /* Accelerazione angolare al passo prec. */
	
	const StructNode *pRefNode;	/* Reference node for relative prediction */

#ifdef USE_NETCDF
	NcVar	*Var_X,
		*Var_Phi,
		*Var_XP,
		*Var_Omega,
		*Var_XPP,
		*Var_OmegaP;
#endif /* USE_NETCDF */

	OrientationDescription od;

	/* Rigidezze fittizie usate nell'assemblaggio dei vincoli */
	doublereal dPositionStiffness;
	doublereal dVelocityStiffness;
	bool bOmegaRot;       /* Flag di velocita' angolare solidale col nodo */

	// reference motion, for relative kinematics
	const RigidBodyKinematics *pRefRBK;

	// makes sense also for dummy nodes, as they may inherit
	// accelerations from the parent node
	bool bOutputAccels;

public:
	/* Costruttore definitivo */
	StructNode(unsigned int uL,
		const DofOwner* pDO,
		const Vec3& X0,
		const Mat3x3& R0,
		const Vec3& V0,
		const Vec3& W0,
		const StructNode *pRN,
		const RigidBodyKinematics *pRBK,
		doublereal dPosStiff,
		doublereal dVelStiff,
		bool bOmRot,
		OrientationDescription ood,
		flag fOut);

	/* Distruttore (per ora e' banale) */
	virtual ~StructNode(void);

	/* Tipo di nodo */
	virtual Node::Type GetNodeType(void) const;

	/* rigid-body kinematics */
	const RigidBodyKinematics *pGetRBK(void) const;

	const Vec3& GetX(void) const;
	const Mat3x3& GetR(void) const;
	const Vec3& GetV(void) const;
	const Vec3& GetW(void) const;
	const Vec3& GetXPP(void) const;
	const Vec3& GetWP(void) const;

	/* Contributo del nodo strutturale al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual std::ostream& DescribeDof(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;

	virtual void DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;

	virtual std::ostream& DescribeEq(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;

	virtual void DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;

	/* Restituisce il valore del dof iDof;
	 * se differenziale, iOrder puo' essere = 1 per la derivata */
	virtual const doublereal& dGetDofValue(int iDof, int iOrder = 0) const;

	/* Restituisce il valore del dof iDof al passo precedente;
	 * se differenziale, iOrder puo' essere = 1 per la derivata */
	virtual const doublereal& dGetDofValuePrev(int iDof, int iOrder = 0) const;

	/* Setta il valore del dof iDof a dValue;
	 * se differenziale, iOrder puo' essere = 1 per la derivata */
	virtual void SetDofValue(const doublereal& dValue,
		unsigned int iDof, unsigned int iOrder = 0);

	virtual DofOrder::Order GetDofType(unsigned int) const;

	/* Tipo di nodo strutturale */
	virtual StructNode::Type GetStructNodeType(void) const = 0;

	/* Ritorna il numero di dofs usato nell'assemblaggio iniziale */
	virtual inline unsigned int iGetInitialNumDof(void) const;

	/* Ritorna il primo indice (-1) di posizione */
	virtual inline integer iGetFirstPositionIndex(void) const;

	/* Ritorna il primo indice (-1) di Quantita' di moto */
	virtual integer iGetFirstMomentumIndex(void) const = 0;

	/* Restituiscono i dati correnti */
	/* Attenzione: restituiscono un reference ai dati veri, per limitare
	 * l'overhead. Tuttavia, una loro modifica e' permanente. Valutare
	 * quindi la possibilita' di far passare un const Mat3x3& ecc., in modo
	 * da obbligare il chiamante a farsi una copia dei dati */
	virtual inline const Mat3x3& GetRPrev(void) const;
	virtual inline const Mat3x3& GetRRef(void) const;
	virtual inline const Mat3x3& GetRCurr(void) const;

	virtual inline const Vec3& GetgRef(void) const;
	virtual inline const Vec3& GetgCurr(void) const;

	virtual inline const Vec3& GetgPRef(void) const;
	virtual inline const Vec3& GetgPCurr(void) const;

	virtual inline const Vec3& GetXPrev(void) const;
	virtual inline const Vec3& GetXCurr(void) const;

	virtual inline const Vec3& GetVPrev(void) const;
	virtual inline const Vec3& GetVCurr(void) const;

	virtual inline const Vec3& GetWPrev(void) const;
	virtual inline const Vec3& GetWRef(void) const;
	virtual inline const Vec3& GetWCurr(void) const;

	virtual inline const Vec3& GetXPPCurr(void) const;
	virtual inline const Vec3& GetXPPPrev(void) const;

	virtual inline const Vec3& GetWPCurr(void) const;
	virtual inline const Vec3& GetWPPrev(void) const;

	virtual inline const doublereal& dGetPositionStiffness(void) const;
	virtual inline const doublereal& dGetVelocityStiffness(void) const;
	virtual inline bool bOmegaRotates(void) const;

	virtual bool ComputeAccelerations(bool b);
	virtual inline bool bComputeAccelerations(void) const;
	virtual inline bool bOutputAccelerations(void) const;
	virtual void OutputAccelerations(bool bOut);

	virtual void OutputPrepare(OutputHandler &OH);

	/* Output del nodo strutturale (da mettere a punto) */
	virtual void Output(OutputHandler& OH) const;

	/* Output della soluzione perturbata (modi ...) */
	virtual void Output(OutputHandler& OH,
		const VectorHandler& X, const VectorHandler& XP) const;

	/* Aggiorna dati in base alla soluzione */
	virtual void Update(const VectorHandler& X,
		const VectorHandler& XP);

	/* Aggiorna dati durante l'iterazione fittizia iniziale */
	virtual void DerivativesUpdate(const VectorHandler& X,
		const VectorHandler& XP);

	/* Aggiorna dati in base alla soluzione durante l'assemblaggio iniziale */
	virtual void InitialUpdate(const VectorHandler& X);

	/* Inverse Dynamics: */
	/* Do Update on node position, velocity or acceleration 
	 * depending on iOrder */
	void Update(const VectorHandler& X, InverseDynamics::Order iOrder);

	/* Funzioni di inizializzazione, ereditate da DofOwnerOwner */
	virtual void SetInitialValue(VectorHandler& X);
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	/* Elaborazione vettori e dati prima e dopo la predizione
	 * per MultiStepIntegrator */
	virtual void BeforePredict(VectorHandler& X, VectorHandler& XP,
		VectorHandler& XPrev,
		VectorHandler& XPPrev) const;
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
	
	/*Inverse Dynamics: reset orientation parameters*/
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP, 
			const VectorHandler& XPP);

	/*
	 * Metodi per l'estrazione di dati "privati".
	 * Si suppone che l'estrattore li sappia interpretare.
	 * Come default non ci sono dati privati estraibili
	 */
	virtual unsigned int iGetNumPrivData(void) const;

	/*
	 * Maps a string (possibly with substrings) to a private data;
	 * returns a valid index ( > 0 && <= iGetNumPrivData()) or 0 
	 * in case of unrecognized data; error must be handled by caller
	 */
	virtual unsigned int iGetPrivDataIdx(const char *s) const;

	/*
	 * Returns the current value of a private data
	 * with 0 < i <= iGetNumPrivData()
	 */
	virtual doublereal dGetPrivData(unsigned int i) const;
}; /* End class StructNode */

/* Ritorna il numero di dofs usato nell'assemblaggio iniziale */
inline unsigned int
StructNode::iGetInitialNumDof(void) const
{
	return iStructNodeInitialDofNumber;
}


/* Ritorna il primo indice (-1) di posizione */
inline integer
StructNode::iGetFirstPositionIndex(void) const
{
	return DofOwnerOwner::iGetFirstIndex();
}


/* Restituiscono i dati correnti */
/* Attenzione: restituiscono un reference ai dati veri, per limitare
 * l'overhead. Tuttavia, una loro modifica e' permanente. Valutare
 * quindi la possibilita' di far passare un const Mat3x3& ecc., in modo
 * da obbligare il chiamante a farsi una copia dei dati */
inline const Mat3x3&
StructNode::GetRPrev(void) const
{
	return RPrev;
}

inline const Mat3x3&
StructNode::GetRRef(void) const
{
	return RRef;
}

inline const Mat3x3&
StructNode::GetRCurr(void) const
{
	return RCurr;
}

inline const Vec3&
StructNode::GetgRef(void) const
{
	return gRef;
}

inline const Vec3&
StructNode::GetgCurr(void) const
{
	return gCurr;
}

inline const Vec3&
StructNode::GetgPRef(void) const
{
	return gPRef;
}

inline const Vec3&
StructNode::GetgPCurr(void) const
{
	return gPCurr;
}

inline const Vec3&
StructNode::GetXPrev(void) const
{
	return XPrev;
}

inline const Vec3&
StructNode::GetXCurr(void) const
{
	return XCurr;
}

inline const Vec3&
StructNode::GetVPrev(void) const
{
	return VPrev;
}

inline const Vec3&
StructNode::GetVCurr(void) const
{
	return VCurr;
}

inline const Vec3&
StructNode::GetWPrev(void) const
{
	return WPrev;
}

inline const Vec3&
StructNode::GetWRef(void) const
{
	return WRef;
}

inline const Vec3&
StructNode::GetWCurr(void) const
{
	return WCurr;
}

inline const Vec3&
StructNode::GetXPPPrev(void) const
{
	return XPPPrev;
}

inline const Vec3&
StructNode::GetXPPCurr(void) const
{
	return XPPCurr;
}

inline const Vec3&
StructNode::GetWPPrev(void) const
{
	return WPPrev;
}

inline const Vec3&
StructNode::GetWPCurr(void) const
{
	return WPCurr;
}

inline const doublereal&
StructNode::dGetPositionStiffness(void) const
{
	return dPositionStiffness;
}

inline const doublereal&
StructNode::dGetVelocityStiffness(void) const
{
	return dVelocityStiffness;
}

inline bool
StructNode::bOmegaRotates(void) const
{
	return bOmegaRot;
}

inline bool
StructNode::bComputeAccelerations(void) const
{
	return false;
}

inline bool
StructNode::bOutputAccelerations(void) const
{
	return bOutputAccels;
}

inline void
StructNode::OutputAccelerations(bool bOut)
{
	bOutputAccels = bOut;
}

/* StructNode - end */


/* DynamicStructNode - begin */

/* Nodo strutturale per problemi dinamici: possiede i gradi di liberta' di:
 *  - spostamento assoluto,
 *  - parametri di rotazione incrementali,
 *  - quantita' di moto,
 *  - momento della quantita' di moto rispetto al polo mobile.
 * Il nodo di per se' non ha caratteristiche inerziali, che gli vengono date
 * dagli elementi ad esso collegati. In particolare, gli elementi Body, corpo
 * rigido, sono responsabili dell'attribuzione di inerzia ai nodi. Altri
 * contributi possono giungere da travi con matrice di inerzia consistente
 * (non ancora implementate). Fa eccezione il caso di un nodo incastrato.
 * In questo caso non e' necessario attribuirgli inerzia perche' il vincolo
 * di incastro, ClampJoint, si occupa di rendere non singolare la matrice
 * jacobiana. */


/* Numero di dof del tipo di nodo */
const int iDynamicStructNodeDofNumber = 12;

/* Forward declaration */
class AutomaticStructElem;

class DynamicStructNode : public StructNode {
	friend class AutomaticStructElem;

protected:
	/* Acceleration and angular acceleration; DynamicStructNode uses them
	 * only for output; ModalNode uses them to store actual unknowns */

	// "mutable" because can be set in iGetPrivDataIdx
	mutable bool bComputeAccels;
	mutable AutomaticStructElem *pAutoStr;

	virtual inline void SetAutoStr(const AutomaticStructElem *p);
	virtual inline const AutomaticStructElem *pGetAutoStr(void) const;

public:
	/* Costruttore definitivo (da mettere a punto) */
	/* I dati sono passati a mezzo di reference, quindi i relativi oggetti
	 * devono essere creati da chi costruisce il nodo, ovvero la funzione
	 * DataManager::ReadStructNode(). Non e' il modo piu' efficiente ma e'
	 * comodo e sicuro */
	DynamicStructNode(unsigned int uL,
		const DofOwner* pDO,
		const Vec3& X0,
		const Mat3x3& R0,
		const Vec3& V0,
		const Vec3& W0,
	        const StructNode *pRN,
		const RigidBodyKinematics *pRBK,
		doublereal dPosStiff,
		doublereal dVelStiff,
		bool bOmRot,
		OrientationDescription ood,
		flag fOut);

	/* Distruttore (per ora e' banale) */
	virtual ~DynamicStructNode(void);

	/* Tipo di nodo strutturale */
	virtual StructNode::Type GetStructNodeType(void) const;

	/* rigid-body kinematics */
	const Vec3& GetXPP(void) const;
	const Vec3& GetWP(void) const;

	/* Ritorna il numero di dofs (comune a tutto cio' che possiede dof) */
	virtual inline unsigned int iGetNumDof(void) const;

	virtual std::ostream& DescribeDof(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;

	virtual void DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;

	virtual std::ostream& DescribeEq(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;

	virtual void DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;

	/* Ritorna il primo indice (-1) di quantita' di moto */
	virtual inline integer iGetFirstMomentumIndex(void) const;

	/* Usato dalle forze astratte, dai bulk ecc., per assemblare le forze
	 * al posto giusto */
	virtual integer iGetFirstRowIndex(void) const;

	virtual void AddInertia(const doublereal& dm, const Vec3& dS,
		const Mat3x3& dJ) const;

   	/* Accesso ai suoi dati */
	virtual const Vec3& GetBCurr(void) const;
	virtual const Vec3& GetGCurr(void) const;     
	virtual const Vec3& GetBPCurr(void) const;
	virtual const Vec3& GetGPCurr(void) const;   

	virtual void AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP);

	/* Elaborazione vettori e dati prima e dopo la predizione
	 * per MultiStepIntegrator */
	virtual void BeforePredict(VectorHandler& X, VectorHandler& XP,
		VectorHandler& XPrev,
		VectorHandler& XPPrev) const;
	
	/* Restituisce il valore del dof iDof;
	 * se differenziale, iOrder puo' essere = 1 per la derivata */
	virtual const doublereal& dGetDofValue(int iDof, int iOrder = 0) const;

	/* Restituisce il valore del dof iDof al passo precedente;
	 * se differenziale, iOrder puo' essere = 1 per la derivata */
	virtual const doublereal& dGetDofValuePrev(int iDof, int iOrder = 0) const;

	/* Setta il valore del dof iDof a dValue;
	 * se differenziale, iOrder puo' essere = 1 per la derivata */
	virtual void SetDofValue(const doublereal& dValue,
		unsigned int iDof, unsigned int iOrder = 0);

	virtual inline bool bComputeAccelerations(void) const;
	virtual bool ComputeAccelerations(bool b);
	virtual void SetOutputFlag(flag f = flag(1));
};

inline void
DynamicStructNode::SetAutoStr(const AutomaticStructElem *p)
{
	pAutoStr = const_cast<AutomaticStructElem *>(p);
}

inline bool
DynamicStructNode::bComputeAccelerations(void) const
{
	return bComputeAccels;
}

inline const AutomaticStructElem *
DynamicStructNode::pGetAutoStr(void) const
{
	return pAutoStr;
}

/* Ritorna il numero di dofs (comune a tutto cio' che possiede dof) */
inline unsigned int
DynamicStructNode::iGetNumDof(void) const
{
	return iDynamicStructNodeDofNumber;
}

/* Ritorna il primo indice (-1) di quantita' di moto */
inline integer
DynamicStructNode::iGetFirstMomentumIndex(void) const
{
	return DofOwnerOwner::iGetFirstIndex() + 6;
}

/* DynamicStructNode - end */


/* StaticStructNode - begin */

/* Nodo strutturale per problemi statici: possiede i gradi di liberta' di:
 *  - spostamento assoluto,
 *  - parametri di rotazione incrementali
 * Il nodo puo' essere usato:
 * - in problemi statici e quasi statici
 * - quando e' vincolato da un incastro
 * - per punti geometrici statici, di cui si intende trascurare la dinamica,
 *   la cui non-singolarita' sia garantita da elementi elastici
 *   o da vincoli */

/* Numero di dof del tipo di nodo - usato anche dal DofManager (?) */
const int iStaticStructNodeDofNumber = 6;

class StaticStructNode : public StructNode {
protected:
public:
	/* Costruttore definitivo */
	StaticStructNode(unsigned int uL,
		const DofOwner* pDO,
		const Vec3& X0,
		const Mat3x3& R0,
		const Vec3& V0,
		const Vec3& W0,
	        const StructNode *pRN,
		const RigidBodyKinematics *pRBK,
		doublereal dPosStiff,
		doublereal dVelStiff,
		bool bOmRot,
		OrientationDescription ood,
		flag fOut);

	/* Distruttore (per ora e' banale) */
	virtual ~StaticStructNode(void);

	/* Tipo di nodo strutturale */
	virtual StructNode::Type GetStructNodeType(void) const;

	/* Ritorna il numero di dofs (comune a tutto cio' che possiede dof) */
	virtual inline unsigned int iGetNumDof(void) const;

	/* Ritorna il primo indice (-1) di quantita' di moto */
	virtual inline integer iGetFirstMomentumIndex(void) const;
};

/* Ritorna il numero di dofs (comune a tutto cio' che possiede dof) */
inline unsigned int
StaticStructNode::iGetNumDof(void) const
{
	return iStaticStructNodeDofNumber;
}


/* Ritorna il primo indice (-1) di quantita' di moto */
inline integer
StaticStructNode::iGetFirstMomentumIndex(void) const
{
	return DofOwnerOwner::iGetFirstIndex();
}

/* StaticStructNode - end */


/* classe ModalNode derivato da Dynamic */

class ModalNode : public DynamicStructNode {
protected:
public:
	/* Costruttore definitivo (da mettere a punto) */
	/* I dati sono passati a mezzo di reference, quindi i relativi oggetti
	 * devono essere creati da chi costruisce il nodo, ovvero la funzione
	 * DataManager::ReadStructNode(). Non e' il modo piu' efficiente ma e'
	 * comodo e sicuro */
	ModalNode(unsigned int uL,
		const DofOwner* pDO,
		const Vec3& X0,
		const Mat3x3& R0,
		const Vec3& V0,
		const Vec3& W0,
		const RigidBodyKinematics *pRBK,
		doublereal dPosStiff,
		doublereal dVelStiff,
		bool bOmRot,
		OrientationDescription ood,
		flag fOut);

	/* Distruttore (per ora e' banale) */
	virtual ~ModalNode(void);

	/* Tipo di nodo strutturale */
	virtual StructNode::Type GetStructNodeType(void) const;

	/* Ritorna il numero di dofs (comune a tutto cio' che possiede dof) */
	virtual inline unsigned int iGetNumDof(void) const;

	virtual std::ostream& DescribeDof(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;

	virtual void DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;

	virtual std::ostream& DescribeEq(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;

	virtual void DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;

	/* Ritorna il primo indice (-1) di quantita' di moto */
	virtual inline integer iGetFirstMomentumIndex(void) const;

	/* Usato dalle forze astratte, dai bulk ecc., per assemblare le forze
	 * al posto giusto */
	virtual integer iGetFirstRowIndex(void) const;

	/* Aggiorna dati in base alla soluzione */
	virtual void Update(const VectorHandler& X,
		const VectorHandler& XP);
};


/* Ritorna il numero di dofs (comune a tutto cio' che possiede dof) */
inline unsigned int
ModalNode::iGetNumDof(void) const
{
	return iDynamicStructNodeDofNumber;
}

/* Ritorna il primo indice (-1) di quantita' di moto */
inline integer
ModalNode::iGetFirstMomentumIndex(void) const
{
	return DofOwnerOwner::iGetFirstIndex() + 6;
}

/* ModalNode - end */


/* DummyStructNode - begin */

class DummyStructNode : public StructNode {
public:
	enum Type {
		UNKNOWN = -1,

		OFFSET = 0,
		RELATIVEFRAME,
		PIVOTRELATIVEFRAME,

		LASTTYPE
	};

protected:
	const StructNode* pNode;

	virtual void Update_int(void) = 0;

public:
	/* Costruttore definitivo */
	DummyStructNode(unsigned int uL,
		const DofOwner* pDO,
		const StructNode* pNode,
		OrientationDescription ood,
		flag fOut);

	/* Distruttore (per ora e' banale) */
	virtual ~DummyStructNode(void);

	/* tipo */
	virtual DummyStructNode::Type GetDummyType(void) const = 0;

	/* Ritorna il numero di dofs (comune a tutto cio' che possiede dof) */
	virtual inline unsigned int iGetNumDof(void) const;

	/* Restituisce il valore del dof iDof;
	 * se differenziale, iOrder puo' essere = 1 per la derivata */
	virtual const doublereal& dGetDofValue(int iDof, int iOrder = 0) const;

	/* Restituisce il valore del dof iDof al passo precedente;
	 * se differenziale, iOrder puo' essere = 1 per la derivata */
	virtual const doublereal& dGetDofValuePrev(int iDof, int iOrder = 0) const;

	/* Setta il valore del dof iDof a dValue;
	 * se differenziale, iOrder puo' essere = 1 per la derivata */
	virtual void SetDofValue(const doublereal& dValue,
		unsigned int iDof, unsigned int iOrder = 0);

	/* Tipo di nodo strutturale */
	virtual StructNode::Type GetStructNodeType(void) const;

	/* Ritorna il numero di dofs usato nell'assemblaggio iniziale */
	virtual inline unsigned int iGetInitialNumDof(void) const;

	/* Ritorna il primo indice (-1) di posizione */
	virtual inline integer iGetFirstPositionIndex(void) const;

	/* Ritorna il primo indice (-1) di Quantita' di moto */
	virtual inline integer iGetFirstMomentumIndex(void) const;

	/* Aggiorna dati durante l'iterazione fittizia iniziale */
	virtual void DerivativesUpdate(const VectorHandler& X,
		const VectorHandler& XP);

	/* Aggiorna dati in base alla soluzione durante l'assemblaggio iniziale */
	virtual void InitialUpdate(const VectorHandler& X);

	/* Funzioni di inizializzazione, ereditate da DofOwnerOwner */
	virtual void SetInitialValue(VectorHandler& X);
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	/* Elaborazione vettori e dati prima e dopo la predizione
	 * per MultiStepIntegrator */
	virtual void BeforePredict(VectorHandler& X, VectorHandler& XP,
		VectorHandler& XPrev,
		VectorHandler& XPPrev) const;
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);

	virtual inline bool bComputeAccelerations(void) const;
	virtual bool ComputeAccelerations(bool b);
};

/* Ritorna il numero di dofs usato nell'assemblaggio iniziale */
inline unsigned int
DummyStructNode::iGetInitialNumDof(void) const
{
	return 0;
}

/* Ritorna il primo indice (-1) di posizione */
inline integer
DummyStructNode::iGetFirstPositionIndex(void) const
{
	silent_cerr("DummyStructNode(" << GetLabel() << ") has no dofs"
		<< std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* Ritorna il primo indice (-1) di Quantita' di moto */
inline integer
DummyStructNode::iGetFirstMomentumIndex(void) const
{
	silent_cerr("DummyStructNode(" << GetLabel() << ") has no dofs"
		<< std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* Ritorna il numero di dofs (comune a tutto cio' che possiede dof) */
inline unsigned int
DummyStructNode::iGetNumDof(void) const
{
	return 0;
}

inline bool
DummyStructNode::bComputeAccelerations(void) const
{
	return pNode->bComputeAccelerations();
}

/* DummyStructNode - end */


/* OffsetDummyStructNode - begin */

class OffsetDummyStructNode : public DummyStructNode {
protected:
	Vec3 f;
	Mat3x3 R;

	void Update_int(void);

public:
	/* Costruttore definitivo */
	OffsetDummyStructNode(unsigned int uL,
		const DofOwner* pDO,
		const StructNode* pNode,
		const Vec3& f,
		const Mat3x3& R,
		OrientationDescription ood,
		flag fOut);

	/* Distruttore (per ora e' banale) */
	virtual ~OffsetDummyStructNode(void);

	/* tipo */
	virtual DummyStructNode::Type GetDummyType(void) const;

	/* Aggiorna dati in base alla soluzione */
	virtual void Update(const VectorHandler& X,
		const VectorHandler& XP);
};

/* OffsetDummyStrNode - end */


/* RelFrameDummyStructNode - begin */

class RelFrameDummyStructNode : public DummyStructNode {
protected:
	const StructNode* pNodeRef;
	const Mat3x3 RhT;
	const Vec3 fhT;

	void Update_int(void);

public:
	/* Costruttore definitivo */
	RelFrameDummyStructNode(unsigned int uL,
		const DofOwner* pDO,
		const StructNode* pNode,
		const StructNode* pNodeRef,
		const Vec3& fh,
		const Mat3x3& Rh,
		OrientationDescription ood,
		flag fOut);

	/* Distruttore (per ora e' banale) */
	virtual ~RelFrameDummyStructNode(void);

	/* tipo */
	virtual DummyStructNode::Type GetDummyType(void) const;

	/* Aggiorna dati in base alla soluzione */
	virtual void Update(const VectorHandler& X,
		const VectorHandler& XP);

	virtual inline bool bComputeAccelerations(void) const;
	virtual bool ComputeAccelerations(bool b);
};

inline bool
RelFrameDummyStructNode::bComputeAccelerations(void) const
{
	return pNode->bComputeAccelerations() && pNodeRef->bComputeAccelerations();
}

/* RelFrameDummyStructNode - end */

/* PivotRelFrameDummyStructNode - begin */

class PivotRelFrameDummyStructNode : public RelFrameDummyStructNode {
protected:
	const StructNode* pNodeRef2;
	const Mat3x3 Rh2;
	const Vec3 fh2;

	void Update_int(void);

public:
	/* Costruttore definitivo */
	PivotRelFrameDummyStructNode(unsigned int uL,
		const DofOwner* pDO,
		const StructNode* pNode,
		const StructNode* pNodeRef,
		const Vec3& fh,
		const Mat3x3& Rh,
		const StructNode* pNodeRef2,
		const Vec3& fh2,
		const Mat3x3& Rh2,
		OrientationDescription ood,
		flag fOut);

	/* Distruttore (per ora e' banale) */
	virtual ~PivotRelFrameDummyStructNode(void);

	/* tipo */
	virtual DummyStructNode::Type GetDummyType(void) const;

	/* Aggiorna dati in base alla soluzione */
	virtual void Update(const VectorHandler& X,
		const VectorHandler& XP);

	virtual inline bool bComputeAccelerations(void) const;
	virtual bool ComputeAccelerations(bool b);
};

inline bool
PivotRelFrameDummyStructNode::bComputeAccelerations(void) const
{
	return pNode->bComputeAccelerations()
		&& pNodeRef->bComputeAccelerations()
		&& pNodeRef2->bComputeAccelerations();
}

/* PivotRelFrameDummyStructNode - end */

class DataManager;
class MBDynParser;

extern Node*
ReadStructNode(DataManager* pDM,
	MBDynParser& HP,
	DofOwner* pDO,
	unsigned int uLabel);

#endif /* STRNODE_H */
