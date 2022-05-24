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

#ifndef STRNODE_H
#define STRNODE_H

#include "node.h"
#include "matvec3.h"
#include "rbk.h"
#include "invdyn.h"
#include "output.h"

#ifdef USE_AUTODIFF
#include "gradient.h"
#include "matvec.h"
#endif

#ifdef USE_SPARSE_AUTODIFF
#include "sp_gradient.h"
#include "sp_gradient_op.h"
#include "sp_matrix_base.h"
#endif

extern const char* psStructNodeNames[];


class StructDispNode;

class StructNode;

/* StructDispNode - begin */

/* Nodo strutturale: possiede i gradi di liberta' di:
 *  - spostamento assoluto,
 * inoltre, se dinamico, i gdl di:
 *  - quantita' di moto,
 * Il nodo di per se' non ha caratteristiche inerziali, che gli vengono date
 * dagli elementi ad esso collegati. In particolare, gli elementi Body, corpo
 * rigido, sono responsabili dell'attribuzione di inerzia ai nodi. Altri
 * contributi possono giungere da travi con matrice di inerzia consistente
 * (non ancora implementate). */

class StructDispNode : public Node, public RigidBodyKinematics {
public:
	class ErrGeneric : public MBDynErrBase {
  	public:
 		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

	enum Type {
		UNKNOWN = -1,

		DYNAMIC = 0,
		STATIC,

		LASTSTRUCTDISPNODETYPE
	};

	enum Output {
		OUTPUT_ACCELERATIONS = (ToBeOutput::OUTPUT_PRIVATE << 0),
		OUTPUT_INERTIA = (ToBeOutput::OUTPUT_PRIVATE << 1)
	};

protected:
	mutable Vec3 XPrev;   /* Posizione al passo precedente */
	mutable Vec3 XCurr;   /* Posizione corrente */

	mutable Vec3 VPrev;   /* Velocita' al passo precedente */
	mutable Vec3 VCurr;   /* Velocita' corrente */

	mutable Vec3 XPPCurr;   /* Accelerazione lineare  corrente */
	mutable Vec3 XPPPrev;   /* Accelerazione lineare  al passo prec. */

#ifdef USE_SPARSE_AUTODIFF
        mutable Vec3 XY;
#endif
	const StructNode *pRefNode;	/* Reference node for relative prediction
					WARNING: used only if the relative macro is
					active (not default, see configuration options!) */

#ifdef USE_NETCDF
	size_t Var_X,
		Var_Phi,
		Var_XP,
		Var_Omega,
		Var_XPP,
		Var_OmegaP;
#endif /* USE_NETCDF */

	OrientationDescription od;

	/* Rigidezze fittizie usate nell'assemblaggio dei vincoli */
	doublereal dPositionStiffness;
	doublereal dVelocityStiffness;

	// FIXME: reference motion, for relative kinematics
	const RigidBodyKinematics *pRefRBK;

	// makes sense also for dummy nodes, as they may inherit
	// accelerations from the parent node
	bool bOutputAccels;

#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
	/*
	 * Returns the dof index -1 of VCurr during initial assembly
	 */
	virtual integer iGetInitialFirstIndexPrime() const=0;
#endif

public:
	/* Costruttore definitivo */
	StructDispNode(unsigned int uL,
		const DofOwner* pDO,
		const Vec3& X0,
		const Vec3& V0,
		const StructNode *pRN,
		const RigidBodyKinematics *pRBK,
		doublereal dPosStiff,
		doublereal dVelStiff,
		OrientationDescription od,
		flag fOut);

	/* Distruttore (per ora e' banale) */
	virtual ~StructDispNode(void);

	/* Tipo di nodo */
	virtual Node::Type GetNodeType(void) const;

	/* FIXME: rigid-body kinematics */
	const RigidBodyKinematics *pGetRBK(void) const;

	// RBK
	const Vec3& GetX(void) const;
	const Mat3x3& GetR(void) const;
	const Vec3& GetV(void) const;
	const Vec3& GetW(void) const;
	const Vec3& GetXPP(void) const;
	const Vec3& GetWP(void) const;

	/* Ritorna il primo indice (-1) di posizione */
	virtual inline integer iGetFirstPositionIndex(void) const;

	/* Ritorna il primo indice (-1) di Quantita' di moto */
	virtual integer iGetFirstMomentumIndex(void) const = 0;

	/* Tipo di nodo strutturale */
	virtual StructDispNode::Type GetStructDispNodeType(void) const = 0;

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

	virtual inline const Vec3& GetXPrev(void) const;
	virtual inline const Vec3& GetXCurr(void) const;

	virtual inline const Vec3& GetVPrev(void) const;
	virtual inline const Vec3& GetVCurr(void) const;

	virtual inline const Vec3& GetXPPPrev(void) const;
	virtual inline const Vec3& GetXPPCurr(void) const;

#ifdef USE_AUTODIFF
	inline void GetXCurr(grad::Vector<doublereal, 3>& X, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;

	template <grad::index_type N_SIZE>
	inline void GetXCurr(grad::Vector<grad::Gradient<N_SIZE>, 3>& X, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;

	inline void GetVCurr(grad::Vector<doublereal, 3>& V, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;

	template <grad::index_type N_SIZE>
	inline void GetVCurr(grad::Vector<grad::Gradient<N_SIZE>, 3>& V, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;
#endif

#ifdef USE_SPARSE_AUTODIFF
     	inline void GetXCurr(sp_grad::SpColVector<doublereal, 3>& X, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void GetXCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& X, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void GetVCurr(sp_grad::SpColVector<doublereal, 3>& V, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void GetVCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& V, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void GetXCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& X, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void GetVCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& V, doublereal dCoef, sp_grad::SpFunctionCall func) const;

        virtual void UpdateJac(const VectorHandler& Y, doublereal dCoef) override;
#endif
     
	virtual inline const doublereal& dGetPositionStiffness(void) const;
	virtual inline const doublereal& dGetVelocityStiffness(void) const;

	virtual bool ComputeAccelerations(bool b);
	virtual inline bool bComputeAccelerations(void) const;
	virtual inline bool bOutputAccelerations(void) const;
	virtual void OutputAccelerations(bool bOut);

	virtual void OutputPrepare(OutputHandler &OH);

	/* Output del nodo strutturale (da mettere a punto) */
	virtual void Output(OutputHandler& OH) const;

	virtual const OrientationDescription& GetOrientationDescription(void) const;

#if 0
	/* Output della soluzione perturbata (modi ...) */
	virtual void Output(OutputHandler& OH,
		const VectorHandler& X, const VectorHandler& XP) const;
#endif

	/* Aggiorna dati in base alla soluzione */
	virtual void Update(const VectorHandler& X,
		const VectorHandler& XP);

	/* Aggiorna dati durante l'iterazione fittizia iniziale */
	virtual void DerivativesUpdate(const VectorHandler& X,
		const VectorHandler& XP);

	/* Ritorna il numero di dofs usato nell'assemblaggio iniziale */
	virtual inline unsigned int iGetInitialNumDof(void) const;

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
	
	/* Inverse Dynamics: reset orientation parameters */
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP, 
			const VectorHandler& XPP);

	/* Metodi per l'estrazione di dati "privati".
	 * Si suppone che l'estrattore li sappia interpretare.
	 * Come default non ci sono dati privati estraibili */
	virtual unsigned int iGetNumPrivData(void) const;

	/* Maps a string (possibly with substrings) to a private data;
	 * returns a valid index ( > 0 && <= iGetNumPrivData()) or 0 
	 * in case of unrecognized data; error must be handled by caller */
	virtual unsigned int iGetPrivDataIdx(const char *s) const;

	/* Returns the current value of a private data
	 * with 0 < i <= iGetNumPrivData() */
	virtual doublereal dGetPrivData(unsigned int i) const;

	/* test code for getting dimension of components */
	const virtual MBUnits::Dimensions GetEquationDimension(integer index) const;
};

/* Ritorna il numero di dofs usato nell'assemblaggio iniziale */
inline unsigned int
StructDispNode::iGetInitialNumDof(void) const
{
	return 6;
}

inline const Vec3&
StructDispNode::GetXPrev(void) const
{
	return XPrev;
}

inline const Vec3&
StructDispNode::GetXCurr(void) const
{
	return XCurr;
}

inline const Vec3&
StructDispNode::GetVPrev(void) const
{
	return VPrev;
}

inline const Vec3&
StructDispNode::GetVCurr(void) const
{
	return VCurr;
}

inline const Vec3&
StructDispNode::GetXPPPrev(void) const
{
	return XPPPrev;
}

inline const Vec3&
StructDispNode::GetXPPCurr(void) const
{
	return XPPCurr;
}

#ifdef USE_AUTODIFF
inline void
StructDispNode::GetXCurr(grad::Vector<doublereal, 3>& X, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const
{
	X = XCurr;
}

template <grad::index_type N_SIZE>
inline void
StructDispNode::GetXCurr(grad::Vector<grad::Gradient<N_SIZE>, 3>& X, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const
{
	using namespace grad;

	index_type iFirstDofIndex = -1;

	switch (func) {
	case INITIAL_ASS_JAC:
		GRADIENT_ASSERT(dCoef == 1.);
	case INITIAL_DER_JAC:
	case REGULAR_JAC:
		iFirstDofIndex = iGetFirstIndex();
		break;

	default:
		GRADIENT_ASSERT(false);
	}

	for (index_type i = 1; i <= 3; ++i) {
		Gradient<N_SIZE>& g = X(i);
		g.SetValuePreserve(XCurr(i));
		g.DerivativeResizeReset(pDofMap,
					iFirstDofIndex + 1,
					iFirstDofIndex + 4,
					MapVectorBase::GLOBAL,
					0.);
		g.SetDerivativeGlobal(iFirstDofIndex + i, -dCoef);
	}
}

inline void
StructDispNode::GetVCurr(grad::Vector<doublereal, 3>& V, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const
{
	V = VCurr;
}

template <grad::index_type N_SIZE>
inline void
StructDispNode::GetVCurr(grad::Vector<grad::Gradient<N_SIZE>, 3>& V, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const
{
	using namespace grad;

	index_type iFirstDofIndex = -1;

	switch (func) {
	case INITIAL_ASS_JAC:
		GRADIENT_ASSERT(dCoef == 1.);
		iFirstDofIndex = iGetInitialFirstIndexPrime();
		break;

	case INITIAL_DER_JAC:
	case REGULAR_JAC:
		iFirstDofIndex = iGetFirstIndex();
		break;

	default:
		GRADIENT_ASSERT(false);
	}

	for (index_type i = 1; i <= 3; ++i) {
		Gradient<N_SIZE>& g = V(i);
		g.SetValuePreserve(VCurr(i));
		g.DerivativeResizeReset(pDofMap,
					iFirstDofIndex + 1,
					iFirstDofIndex + 4,
					MapVectorBase::GLOBAL,
					0.);
		g.SetDerivativeGlobal(iFirstDofIndex + i, -1.);
	}
}
#endif

#ifdef USE_SPARSE_AUTODIFF
inline void StructDispNode::GetXCurr(sp_grad::SpColVector<doublereal, 3>& X, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     X = XCurr;
}

inline void StructDispNode::GetXCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& X, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     sp_grad::index_type iFirstDofIndex = -1;

     switch (func) {
     case sp_grad::SpFunctionCall::INITIAL_ASS_JAC:
	  SP_GRAD_ASSERT(dCoef == 1.);
     case sp_grad::SpFunctionCall::INITIAL_DER_JAC:
     case sp_grad::SpFunctionCall::REGULAR_JAC:
	  iFirstDofIndex = iGetFirstIndex();
	  break;

     default:
	  SP_GRAD_ASSERT(false);
     }
     
     X.ResizeReset(3, 1);

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
	  X(i).Reset(XCurr(i), iFirstDofIndex + i, -dCoef);
     }
}

inline void StructDispNode::GetXCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& X, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     X.ResizeReset(3, 1);

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          X(i).Reset(XCurr(i), -dCoef * XY(i));
     }
}

inline void StructDispNode::GetVCurr(sp_grad::SpColVector<doublereal, 3>& V, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     V = VCurr;
}

inline void StructDispNode::GetVCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& V, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     sp_grad::index_type iFirstDofIndex = -1;

     switch (func) {
     case sp_grad::SpFunctionCall::INITIAL_ASS_JAC:
	  SP_GRAD_ASSERT(dCoef == 1.);
	  iFirstDofIndex = iGetInitialFirstIndexPrime();
	  break;

     case sp_grad::SpFunctionCall::INITIAL_DER_JAC:
     case sp_grad::SpFunctionCall::REGULAR_JAC:
	  iFirstDofIndex = iGetFirstIndex();
	  break;

     default:
	  SP_GRAD_ASSERT(false);
     }
     
     V.ResizeReset(3, 1);

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
	  V(i).Reset(VCurr(i), iFirstDofIndex + i, -1.);
     }
}

inline void StructDispNode::GetVCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& V, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     V.ResizeReset(3, 1);

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          V(i).Reset(VCurr(i), -XY(i));
     }
}
#endif

inline const doublereal&
StructDispNode::dGetPositionStiffness(void) const
{
	return dPositionStiffness;
}

inline const doublereal&
StructDispNode::dGetVelocityStiffness(void) const
{
	return dVelocityStiffness;
}

inline bool
StructDispNode::bComputeAccelerations(void) const
{
	return false;
}

inline bool
StructDispNode::bOutputAccelerations(void) const
{
	return bOutputAccels;
}

inline void
StructDispNode::OutputAccelerations(bool bOut)
{
	bOutputAccels = bOut;
}

/* Ritorna il primo indice (-1) di posizione */
inline integer
StructDispNode::iGetFirstPositionIndex(void) const
{
	return DofOwnerOwner::iGetFirstIndex();
}

/* StructDispNode - end */

/* DynamicStructDispNode - begin */

/* Nodo strutturale per problemi dinamici: possiede i gradi di liberta' di:
 *  - spostamento assoluto,
 *  - quantita' di moto,
 * Il nodo di per se' non ha caratteristiche inerziali, che gli vengono date
 * dagli elementi ad esso collegati. In particolare, gli elementi Body, corpo
 * rigido, sono responsabili dell'attribuzione di inerzia ai nodi. Altri
 * contributi possono giungere da travi con matrice di inerzia consistente
 * (non ancora implementate). Fa eccezione il caso di un nodo incastrato.
 * In questo caso non e' necessario attribuirgli inerzia perche' il vincolo
 * di incastro, ClampJoint, si occupa di rendere non singolare la matrice
 * jacobiana. */


/* Forward declaration */
class AutomaticStructDispElem;

class DynamicStructDispNode : virtual public StructDispNode {
protected:
	/* Acceleration and angular acceleration; DynamicStructNode uses them
	 * only for output; ModalNode uses them to store actual unknowns */

	// "mutable" because it can be set in iGetPrivDataIdx
	mutable bool bComputeAccels;
	mutable AutomaticStructDispElem *pAutoStr;

#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
	virtual inline integer iGetInitialFirstIndexPrime() const;
#endif

public:
	/* Costruttore definitivo (da mettere a punto) */
	/* I dati sono passati a mezzo di reference, quindi i relativi oggetti
	 * devono essere creati da chi costruisce il nodo, ovvero la funzione
	 * DataManager::ReadStructNode(). Non e' il modo piu' efficiente ma e'
	 * comodo e sicuro */
	DynamicStructDispNode(unsigned int uL,
		const DofOwner* pDO,
		const Vec3& X0,
		const Vec3& V0,
	        const StructNode *pRN,
		const RigidBodyKinematics *pRBK,
		doublereal dPosStiff,
		doublereal dVelStiff,
		OrientationDescription od,
		flag fOut);

	/* Distruttore (per ora e' banale) */
	virtual ~DynamicStructDispNode(void);

	/* Tipo di nodo strutturale */
	virtual StructDispNode::Type GetStructDispNodeType(void) const;

	virtual inline void SetAutoStr(const AutomaticStructDispElem *p);

	/* rigid-body kinematics */
	const Vec3& GetXPP(void) const;

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

	virtual void AddInertia(const doublereal& dm) const;

   	/* Accesso ai suoi dati */
	virtual const Vec3& GetBCurr(void) const;
	virtual const Vec3& GetBPCurr(void) const;

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

	/* Aggiorna dati in base alla soluzione */
	virtual void Update(const VectorHandler& X,
		const VectorHandler& XP);

	virtual inline bool bComputeAccelerations(void) const;
	virtual bool ComputeAccelerations(bool b);
	virtual void SetOutputFlag(flag f = flag(1));

	/* for getting dimension of equations */
	const virtual MBUnits::Dimensions GetEquationDimension (integer index) const;
};

inline void
DynamicStructDispNode::SetAutoStr(const AutomaticStructDispElem *p)
{
	pAutoStr = const_cast<AutomaticStructDispElem *>(p);
}

inline bool
DynamicStructDispNode::bComputeAccelerations(void) const
{
	return bComputeAccels;
}

/* Ritorna il numero di dofs (comune a tutto cio' che possiede dof) */
inline unsigned int
DynamicStructDispNode::iGetNumDof(void) const
{
	return 6;
}

/* Ritorna il primo indice (-1) di quantita' di moto */
inline integer
DynamicStructDispNode::iGetFirstMomentumIndex(void) const
{
	return DofOwnerOwner::iGetFirstIndex() + 3;
}

#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
inline integer
DynamicStructDispNode::iGetInitialFirstIndexPrime() const
{
	// FIXME: Is it correct this way?
	return iGetFirstIndex() + 3;
}
#endif

/* DynamicStructDispNode - end */


/* StaticStructDispNode - begin */

/* Nodo strutturale per problemi statici: possiede i gradi di liberta' di:
 *  - spostamento assoluto,
 * Il nodo puo' essere usato:
 * - in problemi statici e quasi statici
 * - quando e' vincolato da un incastro
 * - per punti geometrici statici, di cui si intende trascurare la dinamica,
 *   la cui non-singolarita' sia garantita da elementi elastici
 *   o da vincoli */

/* Numero di dof del tipo di nodo - usato anche dal DofManager (?) */
class StaticStructDispNode : virtual public StructDispNode {
protected:

#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
	virtual inline integer iGetInitialFirstIndexPrime() const;
#endif

public:
	/* Costruttore definitivo */
	StaticStructDispNode(unsigned int uL,
		const DofOwner* pDO,
		const Vec3& X0,
		const Vec3& V0,
	        const StructNode *pRN,
		const RigidBodyKinematics *pRBK,
		doublereal dPosStiff,
		doublereal dVelStiff,
		OrientationDescription od,
		flag fOut);

	/* Distruttore (per ora e' banale) */
	virtual ~StaticStructDispNode(void);

	/* Tipo di nodo strutturale */
	virtual StructDispNode::Type GetStructDispNodeType(void) const;

	/* Ritorna il numero di dofs (comune a tutto cio' che possiede dof) */
	virtual inline unsigned int iGetNumDof(void) const;

	/* Ritorna il primo indice (-1) di quantita' di moto */
	virtual inline integer iGetFirstMomentumIndex(void) const;
};

/* Ritorna il numero di dofs (comune a tutto cio' che possiede dof) */
inline unsigned int
StaticStructDispNode::iGetNumDof(void) const
{
	return 3;
}


/* Ritorna il primo indice (-1) di quantita' di moto */
inline integer
StaticStructDispNode::iGetFirstMomentumIndex(void) const
{
	return DofOwnerOwner::iGetFirstIndex();
}
#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
inline integer
StaticStructDispNode::iGetInitialFirstIndexPrime() const
{
	// FIXME: Is it correct this way?
	return iGetFirstIndex() + 3;
}
#endif

/* StaticStructDispNode - end */


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

class StructNode : virtual public StructDispNode
#if defined(USE_AUTODIFF) && !defined(USE_SPARSE_AUTODIFF)
, public grad::AlignedAlloc
#endif
{
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

	mutable Vec3 WPrev;   /* Velocita' angolare al passo precedente */
	Vec3 WRef;            /* Velocita' angolare predetta al passo corrente */
	mutable Vec3 WCurr;   /* Velocita' angolare corrente */

	mutable Vec3 WPCurr;    /* Accelerazione angolare corrente */
	mutable Vec3 WPPrev;    /* Accelerazione angolare al passo prec. */

	// FIXME: qui o in StructDispNode	
	// const StructNode *pRefNode;	/* Reference node for relative prediction */

	/* Rigidezze fittizie usate nell'assemblaggio dei vincoli */
	bool bOmegaRot;       /* Flag di velocita' angolare solidale col nodo */

	// reference motion, for relative kinematics
	// FIXME: qui o in StructDispNode	
	// const RigidBodyKinematics *pRefRBK;

	// makes sense also for dummy nodes, as they may inherit
	// accelerations from the parent node
	// FIXME: qui o in StructDispNode	
	// bool bOutputAccels;

#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
private:
       mutable bool bNeedRotation, bUpdateRotation, bUpdateRotationGradProd;
#endif

#if defined(USE_AUTODIFF) && !defined(USE_SPARSE_AUTODIFF)
	grad::FunctionCall eCurrFunc;
        void UpdateRotation(doublereal dCoef) const;
		
	static const grad::index_type iNumADVars = 3; // Account for the initial assembly phase
	mutable grad::Matrix<grad::Gradient<iNumADVars>, 3, 3> RCurr_grad;
	mutable grad::Vector<grad::Gradient<iNumADVars>, 3> WCurr_grad;

	template <typename T>
	inline void UpdateRotation(const Mat3x3& RRef, const Vec3& WRef, const grad::Vector<T, 3>& g, const grad::Vector<T, 3>& gP, grad::Matrix<T, 3, 3>& RCurr, grad::Vector<T, 3>& WCurr, doublereal dCoef, enum grad::FunctionCall func) const;

	inline void GetgCurrInt(grad::Vector<grad::Gradient<iNumADVars>, 3>& g, doublereal dCoef, enum grad::FunctionCall func) const;

	inline void GetgPCurrInt(grad::Vector<grad::Gradient<iNumADVars>, 3>& gP, doublereal dCoef, enum grad::FunctionCall func) const;
        inline void GetWCurrInt(grad::Vector<grad::Gradient<iNumADVars>, 3>& gP, doublereal dCoef, enum grad::FunctionCall func) const;
#endif

#ifdef USE_SPARSE_AUTODIFF
        sp_grad::SpFunctionCall eCurrFunc;
        mutable sp_grad::SpMatrixA<sp_grad::SpGradient, 3, 3, 3> RCurr_grad;
        mutable sp_grad::SpColVectorA<sp_grad::SpGradient, 3, 3> WCurr_grad;
        mutable sp_grad::SpMatrixA<sp_grad::GpGradProd, 3, 3> RCurr_gradp;
        mutable sp_grad::SpColVectorA<sp_grad::GpGradProd, 3> WCurr_gradp;
        mutable Vec3 gY;
        void UpdateRotation(doublereal dCoef) const;
        void UpdateRotation(const VectorHandler& Y, doublereal dCoef) const;
        template <typename T>
        inline void UpdateRotation(const Mat3x3& RRef, const Vec3& WRef, const sp_grad::SpColVector<T, 3>& g, const sp_grad::SpColVector<T, 3>& gP, sp_grad::SpMatrix<T, 3, 3>& RCurr, sp_grad::SpColVector<T, 3>& WCurr, doublereal dCoef, sp_grad::SpFunctionCall func) const;
        inline void GetWCurrInitAss(sp_grad::SpColVector<doublereal, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const;
        inline void GetWCurrInitAss(sp_grad::SpColVector<sp_grad::SpGradient, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const;
        inline void GetWCurrInitAss(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const;
#endif

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

	// RBK
	const Mat3x3& GetR(void) const;
	const Vec3& GetW(void) const;
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

	/* Tipo di nodo strutturale */
	virtual StructNode::Type GetStructNodeType(void) const = 0;

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

	/* Ritorna il numero di dofs usato nell'assemblaggio iniziale */
	virtual inline unsigned int iGetInitialNumDof(void) const;

	/* Restituiscono i dati correnti */
	/* Attenzione: restituiscono un reference ai dati veri, per limitare
	 * l'overhead. Tuttavia, una loro modifica e' permanente. Valutare
	 * quindi la possibilita' di far passare un const Mat3x3& ecc., in modo
	 * da obbligare il chiamante a farsi una copia dei dati */
	virtual inline const Vec3& GetgRef(void) const;
	virtual inline const Vec3& GetgCurr(void) const;

	virtual inline const Vec3& GetgPRef(void) const;
	virtual inline const Vec3& GetgPCurr(void) const;

	virtual inline const Mat3x3& GetRPrev(void) const;
	virtual inline const Mat3x3& GetRRef(void) const;
	virtual inline const Mat3x3& GetRCurr(void) const;

	virtual inline const Vec3& GetWPrev(void) const;
	virtual inline const Vec3& GetWRef(void) const;
	virtual inline const Vec3& GetWCurr(void) const;

	virtual inline const Vec3& GetWPCurr(void) const;
	virtual inline const Vec3& GetWPPrev(void) const;

#ifdef USE_AUTODIFF
	inline void GetgCurr(grad::Vector<doublereal, 3>& g, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;

	template <grad::index_type N_SIZE>
	inline void GetgCurr(grad::Vector<grad::Gradient<N_SIZE>, 3>& g, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;

        inline void GetgPCurr(grad::Vector<doublereal, 3>& gP, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;
                
	template <grad::index_type N_SIZE>
	inline void GetgPCurr(grad::Vector<grad::Gradient<N_SIZE>, 3>& gP, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;

	inline void GetRCurr(grad::Matrix<doublereal, 3, 3>& R, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;

	template <grad::index_type N_SIZE>
	inline void GetRCurr(grad::Matrix<grad::Gradient<N_SIZE>, 3, 3>& R, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;

	inline void GetWCurr(grad::Vector<doublereal, 3>& W, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;

	template <grad::index_type N_SIZE>
	inline void GetWCurr(grad::Vector<grad::Gradient<N_SIZE>, 3>& W, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;
#endif

#ifdef USE_SPARSE_AUTODIFF
     	inline void GetgCurr(sp_grad::SpColVector<doublereal, 3>& g, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void GetgCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& g, doublereal dCoef, sp_grad::SpFunctionCall func) const;

        inline void GetgPCurr(sp_grad::SpColVector<doublereal, 3>& gP, doublereal dCoef, sp_grad::SpFunctionCall func) const;
                
	inline void GetgPCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& gP, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void GetRCurr(sp_grad::SpMatrix<doublereal, 3, 3>& R, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void GetRCurr(sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& R, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void GetWCurr(sp_grad::SpColVector<doublereal, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void GetWCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void GetgCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& g, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void GetgPCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& gP, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void GetRCurr(sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& R, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void GetWCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const;     
#endif

#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
        virtual void UpdateJac(doublereal dCoef) override;
#endif
#ifdef USE_SPARSE_AUTODIFF
        virtual void UpdateJac(const VectorHandler& Y, doublereal dCoef) override;
#endif
	virtual inline bool bOmegaRotates(void) const;

	virtual void OutputPrepare(OutputHandler &OH);

	/* Output del nodo strutturale (da mettere a punto) */
	virtual void Output(OutputHandler& OH) const;

#if 0
	/* Output della soluzione perturbata (modi ...) */
	virtual void Output(OutputHandler& OH,
		const VectorHandler& X, const VectorHandler& XP) const;
#endif

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

	/* test code for getting dimension of components */
	const virtual MBUnits::Dimensions GetEquationDimension(integer index) const;
}; /* End class StructNode */

/* Ritorna il numero di dofs usato nell'assemblaggio iniziale */
inline unsigned int
StructNode::iGetInitialNumDof(void) const
{
	return 12;
}


/* Restituiscono i dati correnti */
/* Attenzione: restituiscono un reference ai dati veri, per limitare
 * l'overhead. Tuttavia, una loro modifica e' permanente. Valutare
 * quindi la possibilita' di far passare un const Mat3x3& ecc., in modo
 * da obbligare il chiamante a farsi una copia dei dati */
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
StructNode::GetWPPrev(void) const
{
	return WPPrev;
}

inline const Vec3&
StructNode::GetWPCurr(void) const
{
	return WPCurr;
}

#ifdef USE_AUTODIFF
inline void StructNode::GetgCurr(grad::Vector<doublereal, 3>& g, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const
{
	g = gCurr;
}

template <grad::index_type N_SIZE>
inline void StructNode::GetgCurr(grad::Vector<grad::Gradient<N_SIZE>, 3>& g, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const
{
	using namespace grad;

	index_type iFirstDofIndex = -1;

	switch (func) {
	case INITIAL_ASS_JAC:
		GRADIENT_ASSERT(dCoef == 1.);

	case INITIAL_DER_JAC:
	case REGULAR_JAC:
		iFirstDofIndex = iGetFirstIndex();
		break;

	default:
		GRADIENT_ASSERT(false);
	}

	for (index_type i = 1; i <= 3; ++i) {
		Gradient<N_SIZE>& g_i = g(i);
		g_i.SetValuePreserve(gCurr(i));
		g_i.DerivativeResizeReset(pDofMap,
					  iFirstDofIndex + 4,
					  iFirstDofIndex + 7,
					  MapVectorBase::GLOBAL,
					  0.);
		g_i.SetDerivativeGlobal(iFirstDofIndex + i + 3, -dCoef);
	}
}

inline void StructNode::GetgPCurr(grad::Vector<doublereal, 3>& gP, doublereal, enum grad::FunctionCall func, grad::LocalDofMap*) const
{
        GRADIENT_ASSERT(func == grad::INITIAL_DER_RES || func == grad::REGULAR_RES);
        
        gP = gPCurr;
}
        
template <grad::index_type N_SIZE>
inline void StructNode::GetgPCurr(grad::Vector<grad::Gradient<N_SIZE>, 3>& gP, doublereal, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const
{
	using namespace grad;

	index_type iFirstDofIndex = -1;

	switch (func) {
	case INITIAL_ASS_JAC:
	     gP = gPCurr;
	     return;
	     
	case INITIAL_DER_JAC:
	case REGULAR_JAC:
	     iFirstDofIndex = iGetFirstIndex() + 3;
	     break;

	default:
		GRADIENT_ASSERT(false);
	}

	for (index_type i = 1; i <= 3; ++i) {
		Gradient<N_SIZE>& g_i = gP(i);
		g_i.SetValuePreserve(gPCurr(i));
		g_i.DerivativeResizeReset(pDofMap,
                                          iFirstDofIndex + 1,
                                          iFirstDofIndex + 3,
                                          MapVectorBase::GLOBAL,
                                          0.);
		g_i.SetDerivativeGlobal(iFirstDofIndex + i, -1.);
	}
}

inline void StructNode::GetRCurr(grad::Matrix<doublereal, 3, 3>& R, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const
{
	bNeedRotation = true;
	R = RCurr;
}

template <grad::index_type N_SIZE>
inline void StructNode::GetRCurr(grad::Matrix<grad::Gradient<N_SIZE>, 3, 3>& R, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const
{
     GRADIENT_ASSERT(bNeedRotation);
     GRADIENT_ASSERT(!bUpdateRotation);
     
     using namespace grad;

#if !defined(USE_SPARSE_AUTODIFF)
     index_type iFirstDofIndex = -1;

     switch (func) {
     case INITIAL_ASS_JAC:
	  GRADIENT_ASSERT(dCoef == 1.);
     case INITIAL_DER_JAC:
     case REGULAR_JAC:
	  iFirstDofIndex = iGetFirstIndex() + 4;
	  break;

     default:
	  GRADIENT_ASSERT(false);
     }
#endif
     
     for (index_type i = 1; i <= 3; ++i) {
	  for (index_type j = 1; j <= 3; ++j) {
	       const auto& RCurr_ij = RCurr_grad(i, j);	       
	       Gradient<N_SIZE>& R_ij = R(i, j);

	       R_ij.SetValuePreserve(RCurr_ij.dGetValue());
	       
#if !defined(USE_SPARSE_AUTODIFF)
	       R_ij.DerivativeResizeReset(pDofMap,
					  iFirstDofIndex + RCurr_ij.iGetStartIndexLocal(),
					  iFirstDofIndex + RCurr_ij.iGetEndIndexLocal(),
					  MapVectorBase::GLOBAL,
					  0.);
	       
	       for (index_type k = RCurr_ij.iGetStartIndexLocal(); k < RCurr_ij.iGetEndIndexLocal(); ++k) {
		    R_ij.SetDerivativeGlobal(iFirstDofIndex + k, RCurr_ij.dGetDerivativeLocal(k));
	       }
#else
	       sp_grad::SpGradDofStat oDofStat;

	       RCurr_ij.GetDofStat(oDofStat);
	       
	       R_ij.DerivativeResizeReset(pDofMap,
					  oDofStat.iMinDof,
					  oDofStat.iMaxDof + 1,
					  MapVectorBase::GLOBAL,
					  0.);
	       
	       for (const auto& r:RCurr_ij) {
		    R_ij.SetDerivativeGlobal(r.iDof, r.dDer);
	       }
#endif
	  }
     }
}

inline void StructNode::GetWCurr(grad::Vector<doublereal, 3>& W, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const
{
	bNeedRotation = true;
	W = WCurr;
}

template <grad::index_type N_SIZE>
inline void StructNode::GetWCurr(grad::Vector<grad::Gradient<N_SIZE>, 3>& W, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const
{
	using namespace grad;

	GRADIENT_ASSERT(bNeedRotation);
	GRADIENT_ASSERT(!bUpdateRotation);

#if !defined(USE_SPARSE_AUTODIFF)
	index_type iFirstDofIndex = -1;

	switch (func) {
	case INITIAL_ASS_JAC:
		GRADIENT_ASSERT(dCoef == 1.);
		/*
		 * XCurr = [X1,		#1
		 * 			X2,		#2
		 * 			X3,		#3
		 * 			g1,		#4
		 * 			g2,		#5
		 * 			g3,		#6
		 * 			V1,		#7
		 * 			V2,		#8
		 * 			V3,		#9
		 * 			W1,		#10
		 * 			W2,		#11
		 * 			W3];	#12
		 */
		iFirstDofIndex = iGetFirstIndex() + 10;
		break;

	case REGULAR_JAC:
		/*
		 * XCurr = [X1,		#1
		 * 			X2,		#2
		 * 			X3,		#3
		 * 			g1		#4
		 * 			g2		#5
		 * 			g3];	#6
		 */
		iFirstDofIndex = iGetFirstIndex() + 4;
		break;

	default:
		GRADIENT_ASSERT(false);
	}
#endif
	
    for (int i = 1; i <= 3; ++i) {
    	Gradient<N_SIZE>& W_i = W(i);
    	const auto& WCurr_i = WCurr_grad(i);

        W_i.SetValuePreserve(WCurr_i.dGetValue());

#if !defined(USE_SPARSE_AUTODIFF)
        W_i.DerivativeResizeReset(pDofMap,
				  iFirstDofIndex + WCurr_i.iGetStartIndexLocal(),
				  iFirstDofIndex + WCurr_i.iGetEndIndexLocal(),
				  MapVectorBase::GLOBAL,
				  0.);

        for (index_type j = WCurr_i.iGetStartIndexLocal(); j < WCurr_i.iGetEndIndexLocal(); ++j) {
	     W_i.SetDerivativeGlobal(iFirstDofIndex + j, WCurr_i.dGetDerivativeLocal(j));
        }
#else
	sp_grad::SpGradDofStat oDofStat;

	WCurr_i.GetDofStat(oDofStat);
	       
	W_i.DerivativeResizeReset(pDofMap,
				  oDofStat.iMinDof,
				  oDofStat.iMaxDof + 1,
				  MapVectorBase::GLOBAL,
				  0.);
	       
	for (const auto& r:WCurr_i) {
	     W_i.SetDerivativeGlobal(r.iDof, r.dDer);
	}
#endif
    }
}
#endif

#ifdef USE_SPARSE_AUTODIFF
inline void StructNode::GetgCurr(sp_grad::SpColVector<doublereal, 3>& g, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     g = gCurr;
}
     
inline void StructNode::GetgCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& g, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     sp_grad::index_type iFirstDofIndex = -1;

     switch (func) {
     case sp_grad::SpFunctionCall::INITIAL_ASS_JAC:
	  SP_GRAD_ASSERT(dCoef == 1.);

     case sp_grad::SpFunctionCall::INITIAL_DER_JAC:
     case sp_grad::SpFunctionCall::REGULAR_JAC:
	  iFirstDofIndex = iGetFirstIndex();
	  break;

     default:
	  SP_GRAD_ASSERT(false);
     }

     g.ResizeReset(3, 1);

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
	  g(i).Reset(gCurr(i), iFirstDofIndex + i + 3, -dCoef);
     }
}

inline void StructNode::GetgPCurr(sp_grad::SpColVector<doublereal, 3>& gP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     gP = gPCurr;
}
                
inline void StructNode::GetgPCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& gP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     sp_grad::index_type iFirstDofIndex = -1;

     switch (func) {
     case sp_grad::SpFunctionCall::INITIAL_ASS_JAC:
	  gP = gPCurr;
	  return;
	  
     case sp_grad::SpFunctionCall::INITIAL_DER_JAC:
     case sp_grad::SpFunctionCall::REGULAR_JAC:
	  iFirstDofIndex = iGetFirstIndex() + 3;
	  break;

     default:
	  SP_GRAD_ASSERT(false);
     }

     gP.ResizeReset(3, 1);

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
	  gP(i).Reset(gPCurr(i), iFirstDofIndex + i, -1.);
     }
}


inline void StructNode::GetgCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& g, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     g.ResizeReset(3, 1);
     
     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          g(i).Reset(gCurr(i), -dCoef * gY(i));
     }
}

inline void StructNode::GetgPCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& gP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     gP.ResizeReset(3, 1);

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          gP(i).Reset(gPCurr(i), -gY(i));
     }
}

inline void StructNode::GetRCurr(sp_grad::SpMatrix<sp_grad::GpGradProd, 3, 3>& R, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     SP_GRAD_ASSERT(bNeedRotation);
     SP_GRAD_ASSERT(!bUpdateRotationGradProd);
     
     R = RCurr_gradp;
}

inline void StructNode::GetWCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     SP_GRAD_ASSERT(bNeedRotation);
     SP_GRAD_ASSERT(!bUpdateRotationGradProd);
     
     W = WCurr_gradp;
}

inline void StructNode::GetRCurr(sp_grad::SpMatrix<doublereal, 3, 3>& R, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     bNeedRotation = true;
     R = RCurr;
}

inline void StructNode::GetRCurr(sp_grad::SpMatrix<sp_grad::SpGradient, 3, 3>& R, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     SP_GRAD_ASSERT(bNeedRotation);
     SP_GRAD_ASSERT(!bUpdateRotation);
     
     R = RCurr_grad;
}

inline void StructNode::GetWCurr(sp_grad::SpColVector<doublereal, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     bNeedRotation = true;
     W = WCurr;
}

inline void StructNode::GetWCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     SP_GRAD_ASSERT(bNeedRotation);
     SP_GRAD_ASSERT(!bUpdateRotation);

     W = WCurr_grad;
}

inline void StructNode::GetWCurrInitAss(sp_grad::SpColVector<doublereal, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     W = WCurr;
}

inline void StructNode::GetWCurrInitAss(sp_grad::SpColVector<sp_grad::SpGradient, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     sp_grad::index_type iFirstDofIndex = -1;

     switch (func) {
     case sp_grad::SpFunctionCall::INITIAL_ASS_JAC:
	  iFirstDofIndex = iGetFirstIndex() + 9;
	  break;

     default:
	  SP_GRAD_ASSERT(false);
     }

     W.ResizeReset(3, 1);

     for (sp_grad::index_type i = 1; i <= 3; ++i) {
	  W(i).Reset(WCurr(i), iFirstDofIndex + i, -1.);
     }     
}

inline void StructNode::GetWCurrInitAss(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& W, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
}
#endif

inline bool
StructNode::bOmegaRotates(void) const
{
	return bOmegaRot;
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


/* Forward declaration */
class AutomaticStructElem;

class DynamicStructNode
: virtual public StructDispNode,
public DynamicStructDispNode,
public StructNode
{
protected:

#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
	virtual inline integer iGetInitialFirstIndexPrime() const;
#endif

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
	virtual const Vec3& GetGCurr(void) const;     
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

	/* Aggiorna dati in base alla soluzione */
	virtual void Update(const VectorHandler& X,
		const VectorHandler& XP);

	/* to get dimensions of equations */
	const virtual MBUnits::Dimensions GetEquationDimension(integer index) const;
};

/* Ritorna il numero di dofs (comune a tutto cio' che possiede dof) */
inline unsigned int
DynamicStructNode::iGetNumDof(void) const
{
	return 12;
}

/* Ritorna il primo indice (-1) di quantita' di moto */
inline integer
DynamicStructNode::iGetFirstMomentumIndex(void) const
{
	return DofOwnerOwner::iGetFirstIndex() + 6;
}

#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
inline integer
DynamicStructNode::iGetInitialFirstIndexPrime() const
{
	return iGetFirstIndex() + 6;
}
#endif

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

class StaticStructNode
: virtual public StructDispNode,
public StaticStructDispNode,
public StructNode
{
protected:
#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
	virtual inline integer iGetInitialFirstIndexPrime() const;
#endif

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
	return 6;
}


/* Ritorna il primo indice (-1) di quantita' di moto */
inline integer
StaticStructNode::iGetFirstMomentumIndex(void) const
{
	return DofOwnerOwner::iGetFirstIndex();
}

#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
inline integer
StaticStructNode::iGetInitialFirstIndexPrime() const
{
	return iGetFirstIndex() + 6;
}
#endif

/* StaticStructNode - end */


/* classe ModalNode derivato da Dynamic */

class ModalNode : virtual public StructDispNode, public DynamicStructNode {
protected:
#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
	virtual inline integer iGetInitialFirstIndexPrime() const;
#endif

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

#if 0
	/* Ritorna il primo indice (-1) di quantita' di moto */
	virtual inline integer iGetFirstMomentumIndex(void) const;
#endif

	/* Usato dalle forze astratte, dai bulk ecc., per assemblare le forze
	 * al posto giusto */
	virtual integer iGetFirstRowIndex(void) const;

	/* Aggiorna dati in base alla soluzione */
	virtual void Update(const VectorHandler& X,
		const VectorHandler& XP);

	virtual void DerivativesUpdate(const VectorHandler& X,
                                       const VectorHandler& XP);
     
	virtual void AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP);

#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
       using StructNode::GetXPPCurr;
       using StructNode::GetWPCurr;
#endif

#ifdef USE_AUTODIFF                
       inline void
       GetXPPCurr(grad::Vector<doublereal, 3>& XPP, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;
                
       template <grad::index_type N_SIZE>
       inline void
       GetXPPCurr(grad::Vector<grad::Gradient<N_SIZE>, 3>& XPP, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;

       inline void
       GetWPCurr(grad::Vector<doublereal, 3>& WP, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;
                
       template <grad::index_type N_SIZE>
       inline void
       GetWPCurr(grad::Vector<grad::Gradient<N_SIZE>, 3>& WP, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const;
#endif

#ifdef USE_SPARSE_AUTODIFF                
	inline void
	GetXPPCurr(sp_grad::SpColVector<doublereal, 3>& XPP, doublereal dCoef, sp_grad::SpFunctionCall func) const;
                
	inline void
	GetXPPCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& XPP, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void
	GetWPCurr(sp_grad::SpColVector<doublereal, 3>& XPP, doublereal dCoef, sp_grad::SpFunctionCall func) const;
                
	inline void
	GetWPCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& WP, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void
	GetXPPCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& XPP, doublereal dCoef, sp_grad::SpFunctionCall func) const;

	inline void
	GetWPCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& WP, doublereal dCoef, sp_grad::SpFunctionCall func) const;

        virtual void UpdateJac(const VectorHandler& Y, doublereal dCoef) override;
#endif

	/* to get dimension of equations */
	const virtual MBUnits::Dimensions GetEquationDimension (integer index) const;
private:
#ifdef USE_SPARSE_AUTODIFF
        Vec3 XPPY, WPY;
#endif
};


#if 0
/* Ritorna il primo indice (-1) di quantita' di moto */
inline integer
ModalNode::iGetFirstMomentumIndex(void) const
{
	return DofOwnerOwner::iGetFirstIndex() + 6;
}
#endif

#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
inline integer
ModalNode::iGetInitialFirstIndexPrime() const
{
	// FIXME: Don't know how it should be implemented!
	silent_cerr("ModalNode::iGetInitialFirstIndexPrime() not supported yet!" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}
#endif

#if defined(USE_AUTODIFF)
inline void
ModalNode::GetXPPCurr(grad::Vector<doublereal, 3>& XPP, doublereal, enum grad::FunctionCall func, grad::LocalDofMap*) const {
        GRADIENT_ASSERT(func == grad::INITIAL_DER_RES || func == grad::REGULAR_RES);
        XPP = XPPCurr;
}

template <grad::index_type N_SIZE>
inline void
ModalNode::GetXPPCurr(grad::Vector<grad::Gradient<N_SIZE>, 3>& XPP, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const {
	using namespace grad;

	index_type iFirstDofIndex = -1;

	switch (func) {
	case INITIAL_DER_JAC:
	case REGULAR_JAC:
		iFirstDofIndex = iGetFirstIndex();
		break;

	default:
		GRADIENT_ASSERT(false);
	}

	for (index_type i = 1; i <= 3; ++i) {
		Gradient<N_SIZE>& g = XPP(i);
		g.SetValuePreserve(XPPCurr(i));
		g.DerivativeResizeReset(pDofMap,
                                        iFirstDofIndex + 7,
                                        iFirstDofIndex + 10,
                                        MapVectorBase::GLOBAL,
                                        0.);
		g.SetDerivativeGlobal(iFirstDofIndex + i + 6, -1.);
	}
}

inline void
ModalNode::GetWPCurr(grad::Vector<doublereal, 3>& WP, doublereal, enum grad::FunctionCall func, grad::LocalDofMap*) const {
        GRADIENT_ASSERT(func == grad::INITIAL_DER_RES || func == grad::REGULAR_RES);
        WP = WPCurr;
}

template <grad::index_type N_SIZE>
inline void
ModalNode::GetWPCurr(grad::Vector<grad::Gradient<N_SIZE>, 3>& WP, doublereal dCoef, enum grad::FunctionCall func, grad::LocalDofMap* pDofMap) const {
	using namespace grad;

	index_type iFirstDofIndex = -1;

	switch (func) {
	case INITIAL_DER_JAC:
	case REGULAR_JAC:
		iFirstDofIndex = iGetFirstIndex();
		break;

	default:
		GRADIENT_ASSERT(false);
	}

	for (index_type i = 1; i <= 3; ++i) {
		Gradient<N_SIZE>& g = WP(i);
		g.SetValuePreserve(WPCurr(i));
		g.DerivativeResizeReset(pDofMap,
                                        iFirstDofIndex + 10,
                                        iFirstDofIndex + 13,
                                        MapVectorBase::GLOBAL,
                                        0.);
		g.SetDerivativeGlobal(iFirstDofIndex + i + 9, -1.);
	}        
}
        
#endif

#ifdef USE_SPARSE_AUTODIFF                
inline void
ModalNode::GetXPPCurr(sp_grad::SpColVector<doublereal, 3>& XPP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
	XPP = XPPCurr;
}
                
inline void
ModalNode::GetXPPCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& XPP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
	sp_grad::index_type iFirstDofIndex = -1;

	switch (func) {
	case sp_grad::SpFunctionCall::INITIAL_DER_JAC:
	case sp_grad::SpFunctionCall::REGULAR_JAC:
		iFirstDofIndex = iGetFirstIndex();
		break;

	default:
		SP_GRAD_ASSERT(false);
	}

	XPP.ResizeReset(3, 1);
	
	for (sp_grad::index_type i = 1; i <= 3; ++i) {
	     XPP(i).Reset(XPPCurr(i), iFirstDofIndex + i + 6, -1.);
	}		
}

inline void
ModalNode::GetXPPCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& XPP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          XPP(i).Reset(XPPCurr(i), -XPPY(i));
     }
}
                

inline void
ModalNode::GetWPCurr(sp_grad::SpColVector<doublereal, 3>& WP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
	WP = WPCurr;
}
                
inline void
ModalNode::GetWPCurr(sp_grad::SpColVector<sp_grad::SpGradient, 3>& WP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
	sp_grad::index_type iFirstDofIndex = -1;

	switch (func) {
	case sp_grad::SpFunctionCall::INITIAL_DER_JAC:
	case sp_grad::SpFunctionCall::REGULAR_JAC:
		iFirstDofIndex = iGetFirstIndex();
		break;

	default:
		SP_GRAD_ASSERT(false);
	}

	WP.ResizeReset(3, 1);
	
	for (sp_grad::index_type i = 1; i <= 3; ++i) {
	     WP(i).Reset(WPCurr(i), iFirstDofIndex + i + 9, -1.);
	}
}

inline void
ModalNode::GetWPCurr(sp_grad::SpColVector<sp_grad::GpGradProd, 3>& WP, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     for (sp_grad::index_type i = 1; i <= 3; ++i) {
          WP(i).Reset(WPCurr(i), -WPY(i));
     }
}
#endif
	
/* ModalNode - end */


/* DummyStructNode - begin */

class DummyStructNode : virtual public StructDispNode, public StructNode {
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

#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
	virtual inline integer iGetInitialFirstIndexPrime() const;
#endif

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
	virtual StructDispNode::Type GetStructDispNodeType(void) const;

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

     	virtual inline integer iGetFirstIndex(void) const;

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

#if defined(USE_AUTODIFF) || defined(USE_SPARSE_AUTODIFF)
        virtual void UpdateJac(doublereal dCoef) override;
#endif
#ifdef USE_SPARSE_AUTODIFF
        virtual void UpdateJac(const VectorHandler& Y, doublereal dCoef) override;
#endif     
};

/* DummyStructNode - end */


/* OffsetDummyStructNode - begin */

class OffsetDummyStructNode : virtual public StructDispNode, public DummyStructNode {
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

class RelFrameDummyStructNode : virtual public StructDispNode, public DummyStructNode {
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

class PivotRelFrameDummyStructNode : virtual public StructDispNode, public RelFrameDummyStructNode {
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
