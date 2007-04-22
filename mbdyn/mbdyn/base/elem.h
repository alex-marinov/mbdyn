/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

/*****************************************************************************
 *                                                                           *
 *                           Mbdyn - Elements                                *
 *                                                                           *
 *****************************************************************************/


#ifndef ELEM_H
#define ELEM_H

#include <vector>

#include <myassert.h>
#include <except.h>

#include <solman.h>
#include <submat.h>
#include <output.h>

#include <withlab.h>
#include <dofown.h>

#include <simentity.h>
#include <node.h>

#ifdef USE_MULTITHREAD
#include <veciter.h>
#endif /* USE_MULTITHREAD */

extern const char* psElemNames[];
extern const char* psReadControlElems[];
extern const char* psReadElemsElems[];
extern const char* psAdamsElemCode[];

/* classi dichiarate */
class ElemWithDofs;
class ElemGravityOwner;
class AerodynamicElem;
class InitialAssemblyElem;
class Rotor;

/* Elem - begin */

class Elem : public WithLabel, public SimulationEntity, public ToBeOutput
#ifdef USE_MULTITHREAD
, public InUse
#endif /* USE_MULTITHREAD */
{
	/*
	 * Tipi di Elem. Lasciare sempre UNKNOWN = -1, cosi' il primo elemento
	 * ha tipo zero, e l'ultima entry dell'enum, LAST...TYPE, e' uguale
	 * al numero di tipi definiti, quindi puo' essere usata come costante nel 
	 * dimensionamento degli arrays e come flag di fine tipi. 
	 */
public:
	enum Type {
		UNKNOWN = -1,

		/* rotors must be as early as possible */
		ROTOR = 0,
	
		AUTOMATICSTRUCTURAL,
	
		GRAVITY,
		BODY,
		JOINT,
		BEAM,
		PLATE,

		FORCE,

		ELECTRICBULK,
		ELECTRIC,

		HYDRAULIC,
	
		BULK,
		LOADABLE,
		DRIVEN,
		EXTERNAL,

		/* other aerodynamic elements must be as late as possible */
		AIRPROPERTIES,
		AEROMODAL,
		AERODYNAMIC,

		GENEL,

		SOCKETSTREAM_OUTPUT,
		RTAI_OUTPUT = SOCKETSTREAM_OUTPUT,

		LASTELEMTYPE
	};

	struct ChangedEquationStructure {};
   
public:
	Elem(unsigned int uL, flag fOut);
	virtual ~Elem(void);

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const = 0;

	/* Tipo dell'elemento (usato solo per debug ecc.) */
	virtual Elem::Type GetElemType(void) const = 0;

	/* inherited from SimulationEntity */
	virtual unsigned int iGetNumDof(void) const;
	virtual std::ostream& DescribeDof(std::ostream& out,
		char *prefix = "", bool bInitial = false, int i = -1) const;
	virtual std::ostream& DescribeEq(std::ostream& out,
		char *prefix = "", bool bInitial = false, int i = -1) const;
	virtual DofOrder::Order GetDofType(unsigned int) const;

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
	 * Non e' virtuale in quanto ritorn a NULL per tutti i tipi che non hanno
	 * dof propri.
	 * Il metodo GetDofType() ritorna, per ogni dof dell'elemento, l'ordine.
	 * E' usato per completare i singoli Dof relativi all'elemento.
	 */

	/* funzioni proprie */

	/* Dimensioni del workspace */
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const = 0;

	/* assemblaggio matrici per autovalori */
	virtual void
	AssMats(VariableSubMatrixHandler& WorkMatA,
		VariableSubMatrixHandler& WorkMatB,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr) = 0;

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef, 
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr) = 0;

	/* inverse dynamics Jacobian matrix assembly */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr);

	/* inverse dynamics residual assembly */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr, 
		const VectorHandler& XPrimeCurr,
		const VectorHandler& XPrimePrimeCurr,
		int iOrder = -1);
	
	/* returns the number of connected nodes */
	virtual inline int GetNumConnectedNodes(void) const;
	virtual inline void GetConnectedNodes(std::vector<const Node *>& connectedNodes);

	/* Adams output stuff */
	virtual inline unsigned int iGetNumDummyParts(void) const;
	virtual inline void GetDummyPartPos(unsigned int part, Vec3& x, Mat3x3& R) const;
	virtual inline void GetDummyPartVel(unsigned int part, Vec3& v, Vec3& w) const;

#ifdef USE_ADAMS
	virtual inline std::ostream& WriteAdamsDummyPartCmd(std::ostream& out,
		unsigned int part, unsigned int firstId) const;
#endif /* USE_ADAMS */
};

inline int
Elem::GetNumConnectedNodes(void) const
{
	silent_cerr(psElemNames[GetElemType()] << "(" << GetLabel() 
		<< ") cannot be used in parallel environment" << std::endl);
	throw ErrGeneric();
}
 
inline void
Elem::GetConnectedNodes(std::vector<const Node *>& connectedNodes)
{
	/* NOTE: this will cause an exception to be thrown */
	connectedNodes.resize(GetNumConnectedNodes());
}

/* Adams output stuff */
inline unsigned int
Elem::iGetNumDummyParts(void) const
{
	return 0;
}

inline void
Elem::GetDummyPartPos(unsigned int part, Vec3& x, Mat3x3& R) const
{
	throw ErrGeneric();
}

inline void
Elem::GetDummyPartVel(unsigned int part, Vec3& v, Vec3& w) const
{
	throw ErrGeneric();
}

#ifdef USE_ADAMS
inline std::ostream&
Elem::WriteAdamsDummyPartCmd(std::ostream& out, unsigned int part, unsigned int firstId) const
{
	throw ErrGeneric();
}
#endif /* USE_ADAMS */

/* Elem - end */


/* Functional object that reads elements */

/* ElemRead - begin */

struct ElemRead {
	virtual ~ElemRead( void ) { NO_OP; };
	virtual Elem *Read(const DataManager* pDM, MBDynParser& HP) = 0;
};

/* ElemRead - end */

/* register an ElemRead object */
extern bool
SetElem(const char *name, ElemRead *rf);

/* create/destroy */
extern void InitElem(void);
extern void DestroyElem(void);

/* Classe derivata da elem, relativa ad elementi che possiedono gradi di
 * liberta' propri */

/* ElemWithDofs - begin */

class ElemWithDofs : virtual public Elem, public DofOwnerOwner {
public:
	ElemWithDofs(unsigned int uL,
		const DofOwner* pDO, flag fOut);

	virtual ~ElemWithDofs(void);
};

/* ElemWithDofs - end */


/* SubjectToInitialAssembly - begin */

class SubjectToInitialAssembly {
public:
	SubjectToInitialAssembly(void);
	virtual ~SubjectToInitialAssembly(void);

	/* Numero di gradi di liberta' definiti durante l'assemblaggio iniziale
	 * e' dato dai gradi di liberta' soliti piu' le loro derivate necessarie; 
	 * tipicamente per un vincolo di posizione il numero di dof raddoppia, in
	 * quanto vengono derivate tutte le equazioni, mentre per un vincolo di 
	 * velocita' rimane inalterato. Sono possibili casi intermedi per vincoli
	 * misti di posizione e velocita' */
	virtual unsigned int iGetInitialNumDof(void) const = 0;

	/* Dimensione del workspace durante l'assemblaggio iniziale. Occorre tener
	 * conto del numero di dof che l'elemento definisce in questa fase e dei
	 * dof dei nodi che vengono utilizzati. Sono considerati dof indipendenti
	 * la posizione e la velocita' dei nodi */
	virtual void InitialWorkSpaceDim(integer* piNumRows, 
		integer* piNumCols) const = 0;

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& XCurr) = 0;

	/* Contributo al residuo durante l'assemblaggio iniziale */   
	virtual SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr) = 0;
};

/* SubjectToInitialAssembly - end */


/* InitialAssemblyElem - begin */

class InitialAssemblyElem 
: virtual public Elem, public SubjectToInitialAssembly {
public:
	InitialAssemblyElem(unsigned int uL, flag fOut);
	virtual ~InitialAssemblyElem(void);
};

/* InitialAssemblyElem - end */

#endif /* ELEM_H */

