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

/* Elemento modale */

/* 
 * Copyright 1999-2017 Felice Felippone <ffelipp@tin.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 */

/* 
 * Copyright 1999-2017 Pierangelo Masarati  <masarati@aero.polimi.it>
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 *
 * Modified by Pierangelo Masarati
 */

#ifndef MODAL_H
#define MODAL_H

#include <fstream>
#include <joint.h>

#if 0
#define MODAL_USE_INV9
#endif

/* Modal - begin */

/* 
 * ATTENZIONE! 
 * per ora e' derivato da Joint; 
 * puo' darsi che venga creata una classe apposta
 */

class Modal : virtual public Elem, public Joint {
public:
	struct StrNodeData;
protected:
	const ModalNode* const pModalNode;
	const unsigned iRigidOffset;		/* 0 iff pModalNode == 0; else 12 */

	/* configuration of reference point;
	 * from ModalNode iff pModalNode == 0 */
	mutable Vec3	x;
	mutable Mat3x3	R;
	mutable Mat3x3	RT;

	const unsigned int NModes;
	const unsigned int NStrNodes;

	const unsigned int NFEMNodes; // number of FEM nodes, common
	const std::vector<std::string> IdFEMNodes; // ID of FEM nodes, common
	const Mat3xN *pXYZFEMNodes; // local position of FEM nodes, common
	const doublereal dMass; // mass, common
	const Vec3 Inv2; // undeformed static moment, common
	const Mat3x3 Inv7; // undeformed inertia moment, common

	const std::vector<unsigned int> uModeNumber;
	const MatNxN *pModalMass;
	const MatNxN *pModalStiff;
	const MatNxN *pModalDamp;

	const Mat3xN *pPHIt;
	const Mat3xN *pPHIr;
   
	const Mat3xN *pModeShapest;
	const Mat3xN *pModeShapesr;

	Mat3xN *pCurrXYZ;
	Mat3xN *pCurrXYZVel;

	const Mat3xN *pInv3;
	const Mat3xN *pInv4;
	const Mat3xN *pInv5;
	const Mat3xN *pInv8;
	const Mat3xN *pInv9;

	const Mat3xN *pInv10;
	const Mat3xN *pInv11;

	Vec3   Inv3jaj;
	Vec3   Inv3jaPj;
	Mat3x3 Inv8jaj;

	Mat3x3 Inv8jaPj;
	Mat3xN Inv5jaj;
	Mat3xN Inv5jaPj;
     
	Mat3x3 Inv9jkajak;
	Mat3x3 Inv9jkajaPk;
     
	VecN a, a0;
	VecN aPrime, aPrime0;
	VecN b;
	VecN bPrime;

public:
	struct StrNodeData {
		// constant, defined once for all at input
		const StructNode *pNode;
		std::string FEMNode;
		Vec3 OffsetFEM;
		Vec3 OffsetMB;
		Mat3x3 RotMB;

		// variable, constructed during analysis
		Vec3 d1tot;
		Mat3x3 R1tot;
		Mat3x3 R2;

		// variable, constructed during analysis
		Vec3 F;
		Vec3 M;

		bool bOut;
	};

protected:
	std::vector<StrNodeData> SND;

	/* from gravity.h */
	/* momento statico */
	Vec3 GetS_int(void) const;

	/* momento d'inerzia */
	Mat3x3 GetJ_int(void) const;
 
	Vec3 GetB_int(void) const;
	Vec3 GetG_int(void) const;

public:
	/* Costruttore non banale */
	Modal(unsigned int uL,
			const ModalNode* pModalNodeTmp, 
			const Vec3& x0,
			const Mat3x3& R0,
			const DofOwner* pDO,
			unsigned int N,
			unsigned int NS,
			unsigned int NFN,
			doublereal dMass,
			const Vec3& STmp,
			const Mat3x3& JTmp,
			const std::vector<unsigned int>& uModeNumber,
			MatNxN *pGenMass,
			MatNxN *pGenStiff,
			MatNxN *pGenDamp,
			const std::vector<std::string>& IdFEMNodes,
			Mat3xN *pN,
			const std::vector<Modal::StrNodeData>& snd,
			Mat3xN *pPHIt,
			Mat3xN *pPHIr,
			Mat3xN *pModeShapest,
			Mat3xN *pModeShapesr,
			Mat3xN *pInv3,
			Mat3xN *pInv4,
			Mat3xN *pInv5,
			Mat3xN *pInv8,
			Mat3xN *pInv9,
			Mat3xN *pInv10,
			Mat3xN *pInv11,
			VecN *a,
			VecN *aP,
			flag fOut);

	/* Distruttore */
	~Modal(void);
   
	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const;

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual unsigned int iGetNumDof(void) const;
	virtual std::ostream&
	DescribeDof(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void
	DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;
	virtual std::ostream&
	DescribeEq(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void
	DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;
	DofOrder::Order GetDofType(unsigned int i) const;
	DofOrder::Order GetEqType(unsigned int i) const;

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat, doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);

	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec, doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);

	void Output(OutputHandler& OH) const;

	/* funzioni usate nell'assemblaggio iniziale */

	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
			const VectorHandler& XCurr);   
	/* Contributo al residuo durante l'assemblaggio iniziale */   
	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
 
	/* Setta il valore iniziale delle proprie variabili */
	void SetInitialValue(VectorHandler& /* X */ );

	void SetValue(DataManager *pDM,
			VectorHandler& /* X */ , VectorHandler& /* XP */ ,
			SimulationEntity::Hints *ph = 0);

#if 0
	/* Aggiorna dati durante l'iterazione fittizia iniziale */
	virtual void DerivativesUpdate(const VectorHandler& X,
		const VectorHandler& XP);
#endif

	/* Dati privati */
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;   

	/* Funzioni che restituiscono dati che possono servire ad
	 * altri elementi (ad es. agli elementi aerodinamici modali)
	 */

	const Mat3xN& pGetPHIt(void) const {
		return *pModeShapest;
	};

	const Mat3xN& pGetPHIr(void) const {
		return *pModeShapesr;
	};

	// NOTE: not 'const' because modify internal storage
	const Mat3xN& GetCurrFEMNodesPosition(void);
	const Mat3xN& GetCurrFEMNodesVelocity(void);

	integer uGetNModes(void) const {
		return NModes;
	};

	const std::vector<unsigned int>& GetModeList(void) const {
		return uModeNumber;
	};

	const VecN& GetA(void) const {
		return a;
	};

	const VecN& GetAP(void) const {
		return aPrime;
	};

	const VecN& GetB(void) const {
		return b;
	};

	const VecN& GetBP(void) const {
		return bPrime;
	};

	integer uGetNFEMNodes(void) {
		return NFEMNodes;
	};

	integer iGetModalIndex(void) const {
		return iGetFirstIndex();
	};

	const ModalNode* pGetModalNode(void) const {
		return pModalNode;
	};

	/* from gravity.h */
	/* massa totale */
	doublereal dGetM(void) const;

	/* *******PER IL SOLUTORE PARALLELO******** */        
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(NStrNodes + (pModalNode ? 1 : 0));
		for (unsigned int j = 0; j < NStrNodes; j++) {
			connectedNodes[j] = SND[j].pNode;
		}
		if (pModalNode) {
			connectedNodes[NStrNodes] = pModalNode;
		}
	};
	/* ************************************************ */

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* Modal - end */

class DataManager;
class MBDynParser;

extern Joint *
ReadModal(DataManager* pDM, MBDynParser& HP, const DofOwner* pD0,
		unsigned int uLabel);

#endif /* MODAL_H */

