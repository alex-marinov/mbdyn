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

/* Elemento strutturale automatico */

#ifndef AUTOSTR_H
#define AUTOSTR_H

#include "elem.h"
#include "strnode.h"
#include "matvec3.h"

#ifdef USE_SPARSE_AUTODIFF
#include "sp_gradient.h"
#include "sp_matrix_base.h"
#include "sp_matvecass.h"
#endif

/* AutomaticStructDispElem - begin */

class AutomaticStructDispElem : virtual public Elem {
	friend class DynamicStructDispNode;
	friend class DynamicStructNode;

protected:
	mutable DynamicStructDispNode* pNode;
	Vec3 B;
	Vec3 BP;

	/* Accesso ai suoi dati */
	// FIXME: return 0?
	virtual inline const Vec3& GetBCurr(void) const { return B; };
	virtual inline const Vec3& GetGCurr(void) const { return ::Zero3; };
	virtual inline const Vec3& GetBPCurr(void) const { return BP; };
	virtual inline const Vec3& GetGPCurr(void) const { return ::Zero3; };

	mutable doublereal m;

	virtual void ComputeAccelerations(Vec3& XPP) const;

#ifdef USE_NETCDF
	MBDynNcVar Var_B,
		Var_G,
		Var_BP,
		Var_GP;
#endif /* USE_NETCDF */

public:
	AutomaticStructDispElem(const DynamicStructDispNode* pN);

	virtual ~AutomaticStructDispElem(void) {
		NO_OP;
	};

	/* inizializza i dati */
	void Init(const Vec3& b, const Vec3& bp);

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Tipo dell'elemento (usato per debug ecc.) */
	virtual Elem::Type GetElemType(void) const {
		return Elem::AUTOMATICSTRUCTURAL;
	};

#if 0
	// FIXME: ?
	virtual void AddInertia(const doublereal& dm, const Vec3& dS,
		const Mat3x3& dJ);
#endif
	virtual void AddInertia(const doublereal& dm);

	virtual doublereal dGetM(void) const {
		return m;
	};

	// FIXME:
	virtual const Vec3& GetS(void) const {
		return ::Zero3;
	};

	// FIXME: ?
	virtual const Mat3x3& GetJ(void) const {
		return ::Zero3x3;
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
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 6;
		*piNumCols = 3;
	};

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* assemblaggio eig */
	void AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* output; si assume che ogni tipo di elemento sappia, attraverso
	 * l'OutputHandler, dove scrivere il proprio output */
	void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	/* Setta i valori iniziali delle variabili (e fa altre cose)
	 * prima di iniziare l'integrazione */
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(1);
		connectedNodes[0] = pNode;
	};
	/* ************************************************ */

	/* Dati privati */
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;
};

/* AutomaticStructDispElem - end */

/* AutomaticStructElem - begin */

class AutomaticStructElem : virtual public Elem, public AutomaticStructDispElem {
	friend class DynamicStructNode;

protected:
	Vec3 G;
	Vec3 GP;

	/* Accesso ai suoi dati */
	virtual inline const Vec3& GetGCurr(void) const { return G; };
	virtual inline const Vec3& GetGPCurr(void) const { return GP; };

	mutable Vec3 S;
	mutable Mat3x3 J;

	virtual void ComputeAccelerations(Vec3& XPP, Vec3& WP) const;

public:
	AutomaticStructElem(const DynamicStructNode* pN);

	virtual ~AutomaticStructElem(void) {
		NO_OP;
	};

	/* inizializza i dati */
	void Init(const Vec3& b, const Vec3& g, const Vec3& bp, const Vec3& gp);

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void AddInertia(const doublereal& dm, const Vec3& dS,
		const Mat3x3& dJ);

	virtual const Vec3& GetS(void) const {
		return S;
	};

	virtual const Mat3x3& GetJ(void) const {
		return J;
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
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 12;
		*piNumCols = 6;
	};

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* assemblaggio eig */
	void AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

#ifdef USE_SPARSE_AUTODIFF
     template <typename T>
     void
     AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
            doublereal dCoef,
            const sp_grad::SpGradientVectorHandler<T>& XCurr,
            const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
            sp_grad::SpFunctionCall func);
     
     inline void
     UpdateState(const sp_grad::SpColVector<doublereal, 3>& B,
                 const sp_grad::SpColVector<doublereal, 3>& BP,
                 const sp_grad::SpColVector<doublereal, 3>& G,
                 const sp_grad::SpColVector<doublereal, 3>& GP);

     inline void
     UpdateState(const sp_grad::SpColVector<sp_grad::SpGradient, 3>& B,
                 const sp_grad::SpColVector<sp_grad::SpGradient, 3>& BP,
                 const sp_grad::SpColVector<sp_grad::SpGradient, 3>& G,
                 const sp_grad::SpColVector<sp_grad::SpGradient, 3>& GP) {
     }
#endif
     
	/* output; si assume che ogni tipo di elemento sappia, attraverso
	 * l'OutputHandler, dove scrivere il proprio output */
	void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	/* Setta i valori iniziali delle variabili (e fa altre cose)
	 * prima di iniziare l'integrazione */
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(1);
		connectedNodes[0] = pNode;
	};
	/* ************************************************ */

	/* Dati privati */
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;
};

/* AutomaticStructElem - end */

#endif /* AUTOSTR_H */

