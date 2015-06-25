/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

/* Trave a volumi finiti, con predisposizione per forze di inerzia consistenti
 * e legame cositutivo piezoelettico */


#ifndef BEAMSLIDER_H
#define BEAMSLIDER_H

#include <myassert.h>
#include <except.h>

#include <strnode.h>
#include <elem.h>
#include <beam.h>
#include <joint.h>

/*
 * Beam connection: 
 * struttura che contiene i nodi e gli offset di una trave a tre nodi
 *
 * BeamConn - begin
 */
class BeamConn {
protected:
	const Beam *m_pBeam;
	Vec3 m_f[3];
	Mat3x3 m_R[3];
public:
	BeamConn(const Beam *pB, 
			const Vec3& f1, const Vec3& f2, const Vec3& f3,
			const Mat3x3& R1 = Eye3, 
			const Mat3x3& R2 = Eye3, 
			const Mat3x3& R3 = Eye3);

	virtual ~BeamConn(void);
	
	const Beam *
	pGetBeam(void) const { 
		return m_pBeam; 
	};
	
	const StructNode *
	pGetNode(unsigned int i) const { 
		return m_pBeam->pGetNode(i); 
	};

	const Vec3 &
	Getf(unsigned int i) const {
		ASSERT(i >= 1 && i <= 3);
		return m_f[--i];
	};

	const Mat3x3 &
	GetR(unsigned int i) const {
		ASSERT(i >= 1 && i <= 3);
		return m_R[--i];
	};
};

/* BeamConn - end */


/* BeamSliderJoint - begin */

class BeamSliderJoint: virtual public Elem, public Joint {
public:
	enum Type { SPHERICAL, CLASSIC, SPLINE };
private:
	unsigned int nRotConstr;
	unsigned int nBeams;
	unsigned int iCurrBeam;

	enum Type iType;

	const StructNode* pNode;
	const BeamConn *const *ppBeam;
	Vec3 f;
	Mat3x3 R;
	Vec3 F;
	Vec3 m;
	Vec3 M;
	doublereal sRef;
	doublereal s;
	int activeNode;
	doublereal dL;
	doublereal dW[2];

	Vec3 xNod[Beam::NUMNODES];
	Vec3 fTmp[Beam::NUMNODES];
	Vec3 xTmp[Beam::NUMNODES];
	doublereal dN[Beam::NUMNODES];
	doublereal dNp[Beam::NUMNODES];
	doublereal dNpp[Beam::NUMNODES];
	Vec3 x;
	Vec3 l;

	Vec3 fb;
	Vec3 xc;
	Mat3x3 Rb;
   
public:
	/* Costruttore non banale */
	BeamSliderJoint(unsigned int uL, const DofOwner* pDO,
			const StructNode* pN, enum Type iT,
			unsigned int nB, const BeamConn *const *pB,
			unsigned int uIB, unsigned int uIN,
			doublereal dl,
			const Vec3& fTmp, const Mat3x3& RTmp, flag fOut);
   
	/* Distruttore */
	~BeamSliderJoint(void);

	/* Contributo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Tipo di Joint */
	virtual Joint::Type GetJointType(void) const { 
		return Joint::BEAMSLIDER;
	};
   
	virtual unsigned int iGetNumDof(void) const { 
		return 4+nRotConstr;
	};
   
	DofOrder::Order GetDofType(unsigned int i) const {
		ASSERT(i >= 0 && i < iGetNumDof());
		return DofOrder::ALGEBRAIC;
	}

	void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const { 
		*piNumRows = 6*(1+Beam::NUMNODES)+iGetNumDof();
		*piNumCols = *piNumRows;
	};
   
	VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat,
			doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);
	
	SubVectorHandler& 
	AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr);
	
	void Output(OutputHandler& OH) const;
	
	/* funzioni usate nell'assemblaggio iniziale */
	virtual unsigned int iGetInitialNumDof(void) const { 
		// return 8+2*nRotConstr;
		return 0;
	};
	virtual void 
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
		*piNumRows = 12*(1+Beam::NUMNODES)+iGetInitialNumDof(); 
		*piNumCols = *piNumRows;
	};

	/* Contributo allo jacobiano durante l'assemblaggio iniziale */
	VariableSubMatrixHandler& 
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
			const VectorHandler& XCurr);

	/* Contributo al residuo durante l'assemblaggio iniziale */
	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
 
#if 0
	/* Dati privati */
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i = 0) const;
#endif
	
	/* *******PER IL SOLUTORE PARALLELO******** */        
	/*
	 * Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs
	 */
	virtual void 
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
		connectedNodes.resize(1 + nBeams*Beam::NUMNODES);
		connectedNodes[0] = pNode;

		/* for each beam */
		for (unsigned int i = 0; i < nBeams; i++) {
			/* for each node */
			for (int j = 1; j <= Beam::NUMNODES; j++) {
				connectedNodes[Beam::NUMNODES*i + j] = ppBeam[i]->pGetNode(j);
			}
		}
	};
	/* ************************************************ */
};

#endif /* BEAMSLIDER_H */

