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

#ifndef AEROELEM_H
#define AEROELEM_H

/* Elementi aerodinamici stazionari */

#include "aerodyn.h"
#include "rotor.h"
#include "beam.h"
#include "beam2.h"
#include "aerodata.h"

#include "gauss.h"
#include "aerod2.h"
#include "shape.h"

class AerodynamicOutput {
public:
	enum eOutput {
		AEROD_OUT_UNKNOWN	= -1,

		AEROD_OUT_STD		= 0x2,
		AEROD_OUT_PGAUSS	= 0x4,
		AEROD_OUT_NODE		= 0x8,

		AEROD_OUT_MASK		= (AEROD_OUT_STD | AEROD_OUT_PGAUSS | AEROD_OUT_NODE)
	};

protected:
	flag m_eOutput;

	// specific for AEROD_OUT_PGAUSS
	typedef struct {
		Vec3 f;
		doublereal alpha;
	} Aero_output;

	/* per output AERO_OUT_PGAUSS */
	std::vector<Aero_output> pOutput;
	std::vector<Aero_output>::iterator pTmpOutput;
	

public:
	AerodynamicOutput(flag f, int iNP);
	~AerodynamicOutput(void);

	void SetOutputFlag(flag f, int iNP);
	void ResetIterator(void);
	void SetData(const Vec3& v, const doublereal *pd);

	AerodynamicOutput::eOutput GetOutput(void) const;
	bool IsOutput(void) const;
	bool IsSTD(void) const;
	bool IsPGAUSS(void) const;
	bool IsNODE(void) const;
};

/* AerodynamicBody - begin */

/* Elemento aerodinamico associato ad un nodo.
 * Fisicamente corrisponde ad un rettangoloide, associato ad un nodo fisico,
 * dal quale riceve il movimento ed al quale trasmette le forze.
 * Inoltre puo' essere associato ad un elemento di rotore che gli modifica
 * la velocita' in funzione di qualche modello di velocita' indotta
 * ed al quale fornisce le forze generate per il calcolo della trazione totale
 * ecc.
 */

/* Lo modifico in modo che possieda un driver relativo allo svergolamento,
 * che consenta di inserire un contributo di svergolamento costante su
 * tutta l'apertura dell'elemento in modo da consentire la modellazione di
 * una superficie mobile equivalente */

class AerodynamicBody 
: virtual public Elem,
public AerodynamicElem,
public InitialAssemblyElem, 
public DriveOwner,
public AerodynamicOutput  {
protected:
	AeroData* aerodata;
	const StructNode* pNode;
	InducedVelocity* pIndVel;
	flag fPassiveInducedVelocity;
	
	const Vec3 f;		/* Offset del punto di riferimento */
	doublereal dHalfSpan;	/* Semiapertura del rettangoloide */
	const Mat3x3 Ra;	/* Rotaz. del sistema aerodinamico al nodo */
	const Vec3 Ra3;		/* Terza colonna della matrice Ra */
	
	const ShapeOwner Chord;		/* corda */
	const ShapeOwner ForcePoint;	/* punto di app. della forza (1/4) */
	const ShapeOwner VelocityPoint; /* punto di app. b.c. (3/4) */
	const ShapeOwner Twist;         /* svergolamento */
	
	GaussDataIterator GDI;	/* Iteratore sui punti di Gauss */
	doublereal* pdOuta;	/* Dati privati */
	doublereal** pvdOuta;	/* */
	
	Vec3 F;			/* Forza */
	Vec3 M;			/* Momento */
	
	bool bJacobian;		/* Compute Jacobian matrix contribution */
	/*
	 * overload della funzione di ToBeOutput();
	 * serve per allocare il vettore dei dati di output se il flag
	 * viene settato dopo la costruzione
	 */
	virtual void SetOutputFlag(flag f = flag(1));   
	
	/* Assemblaggio residuo */
	void AssVec(SubVectorHandler& WorkVec);
	
public:
	AerodynamicBody(unsigned int uLabel, 
			const StructNode* pN, InducedVelocity* pR,
			const Vec3& fTmp, doublereal dS,
			const Mat3x3& RaTmp,
			const Shape* pC, const Shape* pF, 
			const Shape* pV, const Shape* pT,
			integer iN, AeroData* a,
			const DriveCaller* pDC, flag fOut,
			bool bUseJacobian);
	virtual ~AerodynamicBody(void);
	
	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;
	
	/* Tipo dell'elemento (usato per debug ecc.) */
	virtual Elem::Type GetElemType(void) const {
		return Elem::AERODYNAMIC;
	};
	
	/* funzioni proprie */
	
	/* Dimensioni del workspace */
	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 6;
		*piNumCols = 6; 
	};
	
	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
	       doublereal  dCoef,
	       const VectorHandler& /* XCurr */ ,
	       const VectorHandler& /* XPrimeCurr */ ); /*{
	       	DEBUGCOUTFNAME("AerodynamicBody::AssJac");
		WorkMat.SetNullMatrix();
		return WorkMat;
	};*/
	
	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
	       doublereal dCoef,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);
	
	/*
	 * Elaborazione stato interno dopo la convergenza
	 */
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);

	/*
	 * output; si assume che ogni tipo di elemento sappia, attraverso
	 * l'OutputHandler, dove scrivere il proprio output
	 */
	virtual void Output(OutputHandler& OH) const;
	
	/* Numero di GDL iniziali */
	virtual unsigned int iGetInitialNumDof(void) const { 
		return 0;
	};
	
	/* Dimensioni del workspace */
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 6;
		*piNumCols = 1;
	};
	
	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat,	  
		      const VectorHandler& /* XCurr */) {
		DEBUGCOUTFNAME("AerodynamicBody::InitialAssJac");
		WorkMat.SetNullMatrix();
		return WorkMat;
	};
	
	/* assemblaggio residuo */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		      const VectorHandler& XCurr);
   
	/* Tipo di elemento aerodinamico */
	virtual AerodynamicElem::Type GetAerodynamicElemType(void) const {
		return AerodynamicElem::AERODYNAMICBODY;
	};
	
	/* *******PER IL SOLUTORE PARALLELO******** */        
	/*
	 * Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs
	 */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
		connectedNodes.resize(1);
		connectedNodes[0] = pNode;
	};
	
	virtual const InducedVelocity *pGetInducedVelocity(void) const {
		return pIndVel;
	};
	/* ************************************************ */
};

/* AerodynamicBody - end */


/* AerodynamicBeam - begin */

class AerodynamicBeam
: virtual public Elem,
public AerodynamicElem,
public InitialAssemblyElem,
public DriveOwner,
public AerodynamicOutput  {
protected:
	enum { NODE1 = 0, NODE2, NODE3, LASTNODE };
	enum { DELTAx1 = 0, DELTAg1, DELTAx2, DELTAg2, DELTAx3, DELTAg3 };
	
	AeroData* aerodata;
	const Beam* pBeam;
	const StructNode* pNode1;
	const StructNode* pNode2;
	const StructNode* pNode3;
	InducedVelocity* pIndVel;
	flag fPassiveInducedVelocity;
	
	const Vec3 f1;		/* Offset del punto di riferimento */
	const Vec3 f2;		/* Offset del punto di riferimento */
	const Vec3 f3;		/* Offset del punto di riferimento */
	const Mat3x3 Ra1;	/* Rotaz. del sistema aerodinamico al nodo */
	const Mat3x3 Ra2;	/* Rotaz. del sistema aerodinamico al nodo */
	const Mat3x3 Ra3;	/* Rotaz. del sistema aerodinamico al nodo */
	const Vec3 Ra1_3;	/* Terza colonna della matrice Ra */
	const Vec3 Ra2_3;	/* Terza colonna della matrice Ra */
	const Vec3 Ra3_3;	/* Terza colonna della matrice Ra */
	
	const ShapeOwner Chord;		/* corda */
	const ShapeOwner ForcePoint;    /* punto di app. della forza (1/4) */
	const ShapeOwner VelocityPoint; /* punto di applicazione b.c. (3/4) */
	const ShapeOwner Twist;         /* svergolamento */
	
	GaussDataIterator GDI;	/* Iteratore sui punti di Gauss */
	doublereal* pdOuta;	/* Dati privati */
	doublereal** pvdOuta;	/* */
	
	/*
	 * Nota: li lascio distinti perche' cosi' eventualmente ne posso fare
	 * l'output in modo agevole
	 */
	Vec3 F[LASTNODE];	/* Forza */
	Vec3 M[LASTNODE];	/* Momento */
	
	bool bJacobian;		/* Compute Jacobian matrix contribution */

	/*
	 * overload della funzione di ToBeOutput();
	 * serve per allocare il vettore dei dati di output se il flag
	 * viene settato dopo la costruzione
	 */
	virtual void SetOutputFlag(flag f = flag(1));   
	
	/* Assemblaggio residuo */
	void AssVec(SubVectorHandler& WorkVec);   
	
public:
	AerodynamicBeam(unsigned int uLabel, 
			const Beam* pB, InducedVelocity* pR,
			const Vec3& fTmp1,
			const Vec3& fTmp2,
			const Vec3& fTmp3,
			const Mat3x3& Ra1Tmp,
			const Mat3x3& Ra2Tmp,
			const Mat3x3& Ra3Tmp,
			const Shape* pC, const Shape* pF, 
			const Shape* pV, const Shape* pT,
			integer iN, AeroData* a,
			const DriveCaller* pDC, flag fOut,
			bool bUseJacobian);
	virtual ~AerodynamicBeam(void);
	
	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;
	
	/* Tipo dell'elemento (usato solo per debug ecc.) */
	virtual Elem::Type GetElemType(void) const {
		return Elem::AERODYNAMIC;
	};
	
	/* funzioni proprie */
	
	/* Dimensioni del workspace */
	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 18;
		*piNumCols = 18;
	};
	
	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
	       doublereal  dCoef  ,
	       const VectorHandler& /* XCurr */ ,
	       const VectorHandler& /* XPrimeCurr */ ); 
	       /* {
		DEBUGCOUTFNAME("AerodynamicBeam::AssJac");
		WorkMat.SetNullMatrix();
		return WorkMat;
	};*/
	
	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
	       doublereal dCoef,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);
	
	/* Numero di GDL iniziali */
	virtual unsigned int iGetInitialNumDof(void) const { 
		return 0;
	};
	
	/* Dimensioni del workspace */
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 18;
		*piNumCols = 1;
	};
	
	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat,		 
		      const VectorHandler& /* XCurr */ ) {
		DEBUGCOUTFNAME("AerodynamicBeam::AssJac");
		WorkMat.SetNullMatrix();
		return WorkMat;
	};
	
	/* assemblaggio residuo */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		      const VectorHandler& XCurr);
	
	/*
	 * Elaborazione stato interno dopo la convergenza
	 */
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);

	/*
	 * output; si assume che ogni tipo di elemento sappia, attraverso
	 * l'OutputHandler, dove scrivere il proprio output
	 */
	virtual void Output(OutputHandler& OH) const;
	
	/* Tipo di elemento aerodinamico */
	virtual AerodynamicElem::Type GetAerodynamicElemType(void) const {
		return AerodynamicElem::AERODYNAMICBEAM;
	};
	
	/*
	 *******PER IL SOLUTORE PARALLELO********
	 * Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs
	 */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
		connectedNodes.resize(3);
		connectedNodes[0] = pNode1;
		connectedNodes[1] = pNode2;
		connectedNodes[2] = pNode3;
	};
	
	virtual const InducedVelocity *pGetInducedVelocity(void) const {
		return pIndVel;
	};
   	/* ************************************************ */
};

/* AerodynamicBeam - end */

/* AerodynamicBeam2 - begin */

class AerodynamicBeam2
: virtual public Elem,
public AerodynamicElem,
public InitialAssemblyElem,
public DriveOwner,
public AerodynamicOutput  {
protected:
	enum { NODE1 = 0, NODE2, LASTNODE };
	AeroData* aerodata;
	const Beam2* pBeam;
	const StructNode* pNode1;
	const StructNode* pNode2;
	InducedVelocity* pIndVel;
	flag fPassiveInducedVelocity;
	
	const Vec3 f1;		/* Offset del punto di riferimento */
	const Vec3 f2;		/* Offset del punto di riferimento */
	const Mat3x3 Ra1;	/* Rotaz. del sistema aerodinamico al nodo */
	const Mat3x3 Ra2;	/* Rotaz. del sistema aerodinamico al nodo */
	const Vec3 Ra1_3;	/* Terza colonna della matrice Ra */
	const Vec3 Ra2_3;	/* Terza colonna della matrice Ra */
	
	const ShapeOwner Chord;		/* corda */
	const ShapeOwner ForcePoint;    /* punto di app. della forza (1/4) */
	const ShapeOwner VelocityPoint; /* punto di applicazione b.c. (3/4) */
	const ShapeOwner Twist;         /* svergolamento */
	
	GaussDataIterator GDI;	/* Iteratore sui punti di Gauss */
	doublereal* pdOuta;	/* Dati privati */
	doublereal** pvdOuta;	/* */
	
	/*
	 * Nota: li lascio distinti perche' cosi' eventualmente ne posso fare
	 * l'output in modo agevole
	 */
	Vec3 F[LASTNODE];	/* Forza */
	Vec3 M[LASTNODE];	/* Momento */
	
	/*
	 * overload della funzione di ToBeOutput();
	 * serve per allocare il vettore dei dati di output se il flag
	 * viene settato dopo la costruzione
	 */
	virtual void SetOutputFlag(flag f = flag(1));   
	
	/* Assemblaggio residuo */
	void AssVec(SubVectorHandler& WorkVec);   
	
public:
	AerodynamicBeam2(unsigned int uLabel, 
			const Beam2* pB, InducedVelocity* pR,
			const Vec3& fTmp1,
			const Vec3& fTmp2,
			const Mat3x3& Ra1Tmp,
			const Mat3x3& Ra2Tmp,
			const Shape* pC, const Shape* pF, 
			const Shape* pV, const Shape* pT,
			integer iN, AeroData* a,
			const DriveCaller* pDC, flag fOut);
	virtual ~AerodynamicBeam2(void);
	
	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;
	
	/* Tipo dell'elemento (usato solo per debug ecc.) */
	virtual Elem::Type GetElemType(void) const {
		return Elem::AERODYNAMIC;
	};
	
	/* funzioni proprie */
	
	/* Dimensioni del workspace */
	virtual void
	WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 12;
		*piNumCols = 1;
	};
	
	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
	       doublereal /* dCoef */ ,
	       const VectorHandler& /* XCurr */ ,
	       const VectorHandler& /* XPrimeCurr */ ) {
		DEBUGCOUTFNAME("AerodynamicBeam2::AssJac");
		WorkMat.SetNullMatrix();
		return WorkMat;
	};
	
	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
	       doublereal dCoef,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);
	
	/* Numero di GDL iniziali */
	virtual unsigned int iGetInitialNumDof(void) const { 
		return 0;
	};
	
	/* Dimensioni del workspace */
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 12;
		*piNumCols = 1;
	};
	
	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat,		 
		      const VectorHandler& /* XCurr */ ) {
		DEBUGCOUTFNAME("AerodynamicBeam2::AssJac");
		WorkMat.SetNullMatrix();
		return WorkMat;
	};
	
	/* assemblaggio residuo */
	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec,
		      const VectorHandler& XCurr);
	
	/*
	 * Elaborazione stato interno dopo la convergenza
	 */
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);

	/*
	 * output; si assume che ogni tipo di elemento sappia, attraverso
	 * l'OutputHandler, dove scrivere il proprio output
	 */
	virtual void Output(OutputHandler& OH) const;
	
	/* Tipo di elemento aerodinamico */
	virtual AerodynamicElem::Type GetAerodynamicElemType(void) const {
		return AerodynamicElem::AERODYNAMICBEAM;
	};
	
	/*
	 *******PER IL SOLUTORE PARALLELO********
	 * Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs
	 */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) {
		connectedNodes.resize(2);
		connectedNodes[0] = pNode1;
		connectedNodes[1] = pNode2;
	};
	
	virtual const InducedVelocity *pGetInducedVelocity(void) const {
		return pIndVel;
	};
   	/* ************************************************ */
};

/* AerodynamicBeam - end */

class DataManager;
class MBDynParser;

extern Elem *
ReadAerodynamicBody(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
extern Elem *
ReadAerodynamicBeam(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
extern Elem *
ReadAerodynamicBeam2(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

#endif /* AEROELEM_H */

