/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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


#ifndef SCREWJOINT_H
#define SCREWJOINT_H

#include "joint.h"
#include "constltp.h"
#include "friction.h"


/* ScrewJoint - begin */

class ScrewJoint : 
virtual public Elem, public Joint {
private:      

	integer nTheta;
	doublereal dThetaPrev, dThetaCurr, dTheta, dD;
	const doublereal dTheta0, dD0;
	doublereal dLambda;
	const doublereal dPitch;
	Vec3 Theta;
//protected:
	const StructNode* pNode1;
	const StructNode* pNode2;
	const Mat3x3 R1h;
//	const Mat3x3 R2h;
	const Vec3 f1;
	const Vec3 f2;

	Mat3x3 GammaInv;
	Mat3x3 R1, R2, R1R1h, R2R2h;
	Vec3 R1f1, R2f2, e1hz, D;
	Vec3 F1, C1, C2;

	/* friction related data */
	BasicShapeCoefficient *const Sh_c;
	BasicFriction *const fc;
	doublereal cos_pitch_angle_r;
	doublereal M3diff;
	doublereal vrel;

#ifdef USE_NETCDF
	MBDynNcVar Var_dTheta;
	MBDynNcVar Var_Theta;
	MBDynNcVar Var_vrel;
	MBDynNcVar Var_fc;
	MBDynNcVar Var_MFR;
#endif // USE_NETCDF

// 	bool bFirstRes;

	void AssMat(FullSubMatrixHandler& WM, doublereal dCoef,
			const VectorHandler& XCurr,
			const VectorHandler& XPrimeCurr
	);
	void AssVec(SubVectorHandler& WorkVec, doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr
	);

public:
	/* Costruttore non banale */
	ScrewJoint(unsigned int uL,	       
			const DofOwner* pDO,
			const StructNode* pN1, 
			const StructNode* pN2,
			const Vec3& f1Tmp, 
			const Vec3& f2Tmp,
			const Mat3x3& R1t,
//			const Mat3x3& R2t, 
			const doublereal& p, 
			flag fOut,
			BasicShapeCoefficient *const sh = 0,
			BasicFriction *const f = 0);

	/* Distruttore */
	virtual ~ScrewJoint(void);

	/* Tipo di joint */
	virtual Joint::Type GetJointType(void) const {
		return SCREWJOINT;
	};

// 	/* Contributo al file di restart */
// 	virtual std::ostream& Restart(std::ostream& out) const;

	void OutputPrepare(OutputHandler& OH);
	virtual void Output(OutputHandler& OH) const;

// 	void SetValue(DataManager *pDM,
// 			VectorHandler& X, VectorHandler& XP,
// 			SimulationEntity::Hints *ph = 0);

	virtual void AfterConvergence(const VectorHandler& X, 
		const VectorHandler& XP);

// 	virtual Hint *
// 	ParseHint(DataManager *pDM, const char *s) const;
	         
	virtual unsigned int iGetNumDof(void) const {
		if (fc) {
			return 1 + fc->iGetNumDof();
		}
		return 1;
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
		if (fc) {
			ASSERT(i >= 0 && i < 1 + fc->iGetNumDof());
		} else {
			ASSERT(i >= 0 && i < 1);
		}
		if (i < 1) {
			return DofOrder::ALGEBRAIC;
		} else {
			return fc->GetDofType(i-1);
		}
	};

	virtual DofOrder::Order GetEqType(unsigned int i) const {
		if (fc) {
			ASSERT(i >= 0 && i < 1 + fc->iGetNumDof());
		} else {
			ASSERT(i >= 0 && i < 1);
		}
		if (i < 1) {
			return DofOrder::ALGEBRAIC;
		} else {
		       return fc->GetEqType(i-1);
		}
	};

	virtual void WorkSpaceDim(integer* piNumRows,
			integer* piNumCols) const { 
		*piNumRows = 13;
		*piNumCols = 13; 
		if (fc) {
			*piNumRows += fc->iGetNumDof();
			*piNumCols += fc->iGetNumDof();
		}
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
		return 2;
	};

	virtual void InitialWorkSpaceDim(integer* piNumRows,
			integer* piNumCols) const { 
// 		*piNumRows = 26; 
// 		*piNumCols = 26;
		*piNumRows = 0; 
		*piNumCols = 0;
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

/* ScrewJoint - end */

#endif /* SCREWJOINT_H */

