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

/* Deformable hinges */


#ifndef DRVDISP_H
#define DRVDISP_H

#include "joint.h"

/* DriveDisplacementJoint - begin */

class DriveDisplacementJoint : 
virtual public Elem, public Joint, public TplDriveOwner<Vec3> {
private:      
   
protected:
	const StructNode* pNode1;
	const StructNode* pNode2;
	mutable Vec3 f1;
	mutable Vec3 f2;

	Mat3x3 R1Ref;
	Mat3x3 RRef;

	Vec3 f2Ref;
	Vec3 dRef;

	Vec3 F;

	void AssMat(FullSubMatrixHandler& WM, doublereal dCoef);
	void AssVec(SubVectorHandler& WorkVec, doublereal dCoef);
#ifdef USE_NETCDF
	MBDynNcVar Var_d;
#endif // USE_NETCDF
public:
	/* Costruttore non banale */
	DriveDisplacementJoint(unsigned int uL,	       
			const DofOwner* pDO,
			const TplDriveCaller<Vec3>* pDC,
			const StructNode* pN1, 
			const StructNode* pN2,
			const Vec3& f1,
			const Vec3& f2, 
			flag fOut);

	/* Distruttore */
	virtual ~DriveDisplacementJoint(void);

	/* Tipo di joint */
	virtual Joint::Type GetJointType(void) const {
		return DRIVEDISP;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0);

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const;
	         
	virtual unsigned int iGetNumDof(void) const { 
		return 3;
	};

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
   
	virtual DofOrder::Order GetDofType(unsigned int i) const {
		ASSERT(i >= 0 && i <= 3);
		return DofOrder::ALGEBRAIC;
	};

	virtual void WorkSpaceDim(integer* piNumRows,
			integer* piNumCols) const { 
		*piNumRows = 15;
		*piNumCols = 15; 
	};

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
			doublereal dCoef, 
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);
   
	/* assemblaggio residuo */
	virtual SubVectorHandler& 
	AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);           
	/* funzioni usate nell'assemblaggio iniziale */
	virtual unsigned int iGetInitialNumDof(void) const { 
		return 6;
	};

	virtual void InitialWorkSpaceDim(integer* piNumRows,
			integer* piNumCols) const
	{ 
		*piNumRows = 30; 
		*piNumCols = 30;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
			const VectorHandler& XCurr);
   
	/* Contributo al residuo durante l'assemblaggio iniziale */   
	virtual SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec,
			const VectorHandler& XCurr);
   
	/* Dati privati (aggiungere magari le reazioni vincolari) */
	virtual unsigned int iGetNumPrivData(void) const;

	unsigned int iGetPrivDataIdx(const char *s) const;

	virtual doublereal dGetPrivData(unsigned int i = 0) const;

	/* *******PER IL SOLUTORE PARALLELO******** */        
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	   utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(2);
		connectedNodes[0] = pNode1;
		connectedNodes[1] = pNode2;
	};
	/* ************************************************ */

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* DriveDisplacementJoint - end */

/* DriveDisplacementPinJoint - begin */

class DriveDisplacementPinJoint : 
virtual public Elem, public Joint, public TplDriveOwner<Vec3> {
private:      
   
protected:
	const StructNode* pNode;
	mutable Vec3 f;
	mutable Vec3 x;

	Vec3 fRef;
	Vec3 dRef;

	Vec3 F;

	void AssMat(FullSubMatrixHandler& WM, doublereal dCoef);
	void AssVec(SubVectorHandler& WorkVec, doublereal dCoef);
#ifdef USE_NETCDF
	MBDynNcVar Var_d;
#endif // USE_NETCDF
public:
	/* Costruttore non banale */
	DriveDisplacementPinJoint(unsigned int uL,	       
			const DofOwner* pDO,
			const TplDriveCaller<Vec3>* pDC,
			const StructNode* pN, 
			const Vec3& f, 
			const Vec3& x,
			flag fOut);

	/* Distruttore */
	virtual ~DriveDisplacementPinJoint(void);

	/* Tipo di joint */
	virtual Joint::Type GetJointType(void) const {
		return DRIVEDISPPIN;
	};

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	void OutputPrepare(OutputHandler &OH);
	virtual void Output(OutputHandler& OH) const;

	void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0);

	virtual Hint *
	ParseHint(DataManager *pDM, const char *s) const;
	         
	virtual unsigned int iGetNumDof(void) const { 
		return 3;
	};

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
   
	virtual DofOrder::Order GetDofType(unsigned int i) const {
		ASSERT(i >= 0 && i <= 3);
		return DofOrder::ALGEBRAIC;
	};

	virtual void WorkSpaceDim(integer* piNumRows,
			integer* piNumCols) const { 
		*piNumRows = 9;
		*piNumCols = 9; 
	};

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
			doublereal dCoef, 
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr);
   
	/* assemblaggio residuo */
	virtual SubVectorHandler& 
	AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);

	/* Aggiorna le deformazioni ecc. */
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);           
	/* funzioni usate nell'assemblaggio iniziale */
	virtual unsigned int iGetInitialNumDof(void) const { 
		return 6;
	};

	virtual void InitialWorkSpaceDim(integer* piNumRows,
			integer* piNumCols) const
	{ 
		*piNumRows = 18; 
		*piNumCols = 18;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	virtual VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
			const VectorHandler& XCurr);
   
	/* Contributo al residuo durante l'assemblaggio iniziale */   
	virtual SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec,
			const VectorHandler& XCurr);
   
	/* Dati privati (aggiungere magari le reazioni vincolari) */
	virtual unsigned int iGetNumPrivData(void) const;

	unsigned int iGetPrivDataIdx(const char *s) const;

	virtual doublereal dGetPrivData(unsigned int i = 0) const;

	/* *******PER IL SOLUTORE PARALLELO******** */        
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	   utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(1);
		connectedNodes[0] = pNode;
	};
	/* ************************************************ */

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* DriveDisplacementPinJoint - end */

#endif /* DRVDISP_H */

