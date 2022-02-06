/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cfloat>

#include "dataman.h"
#include "loadable.h"
#include "point_contact.h"

/* Contatto punto - superficie */

/* STRUTTURA NEL FILE 
joint: id contatto, point_contact,
id nodo 1, id nodo ground,  posizione ground, direzione ground
*/


/* Costruttore non banale */
PointSurfaceContact::PointSurfaceContact(unsigned int uL,	       
			const DofOwner* pDO,
			const StructNode* pN1, 
			const StructNode* pNs,
			const Vec3& SDir, const doublereal Ek,
			flag fOut)
: Elem(uL,fOut),
Joint(uL, pDO, fOut), pNode1(pN1), pSup(pNs), SupDirection(SDir)		// assegnamento campi
{
		ASSERT(pNode1 != NULL);			// esecuzione solo se i nodi sono assegnati e se sono strutturali
		ASSERT(pSup != NULL);

		ASSERT(pNode1->GetNodeType() == StructNode::STRUCTURAL);
		ASSERT(pSup->GetNodeType() == StructNode::STRUCTURAL);

		ElasticStiffness = Ek;
}

/* Distruttore banale*/ 
PointSurfaceContact::~PointSurfaceContact(void)
{
		NO_OP;		
}

/* Contributo al file di restart */
std::ostream&
PointSurfaceContact::Restart(std::ostream& out) const
{
	//Joint::Restart(out) << pNode1->GetLabel() << ", "
	//		<< pSup->GetLabel() << std::endl;
	return out << "  not implemented yet: " << std::endl;

	//return out;
}

void PointSurfaceContact::Output(OutputHandler& OH)  const
{

}


/* Assemblaggio residuo */
SubVectorHandler& 
PointSurfaceContact::AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr)
{

		/* Dimensiona e resetta la matrice di lavoro */
		integer iNumRows = 0;
		integer iNumCols = 0;
		WorkSpaceDim(&iNumRows, &iNumCols);
		WorkVec.ResizeReset(iNumRows);
		
		/* Recupera gli indici globali */
		integer iNodePFirstMomIndex = pNode1->iGetFirstMomentumIndex();
		integer iNodeSupFirstMomIndex = pSup->iGetFirstMomentumIndex();
		
		/* Setta gli indici della matrice */
		for (int iCnt = 1; iCnt <=3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNodePFirstMomIndex + iCnt);
		WorkVec.PutRowIndex(iCnt+3, iNodeSupFirstMomIndex + iCnt);
		WorkVec.PutRowIndex(iCnt+6, iNodeSupFirstMomIndex + iCnt+3);
		}
		
		/* Costruisco WorkVec in AssVec */
		AssVec(WorkVec, dCoef);
		
		return WorkVec;
}
/* Costruzione matrice locale */
void PointSurfaceContact::AssVec(SubVectorHandler& WorkVec, doublereal dCoef)
{
	/* Orientazione del terreno nel sistema assoluto */
	//Vec3 
	n = pSup->GetRCurr()*SupDirection;
	
	/* Distanza punto-superficie (piano infinito) */
	Vec3 x = - pSup->GetXCurr() + pNode1->GetXCurr();	// attenzione: VERIFICARE SE TIENE 
												// CONTO DI EVENTUALI DIFFERENZE
												// DI SISTEMI DI RIFERIMENTO
 	
	/* Distanza normale */
	dDeltaL = x.Dot(n);
	
	//Vec3 
	Farm = x - n*dDeltaL;	// braccio momento di trasporto 
	
	/* se dDeltaL>0 c'è contatto */
	if (dDeltaL < 0)
	{
	/* Forze e momenti scambiati */
   	FNode1 = - n*ElasticStiffness*dDeltaL;  
    	FSup = -FNode1;
 	MSup = Farm.Cross(FSup);  
 	} 
 	else {
 	FNode1 = n*0;   
 	FSup = n*0;
 	MSup = n*0; 
 	};
 		
  	/* Costruzione vettore assemblato */
 	WorkVec.Add(1, FNode1);
  	WorkVec.Add(3+1, FSup);
  	WorkVec.Add(6+1, MSup);

 }       
 
 /*Assemblaggio Jacobiano*/
 VariableSubMatrixHandler&
PointSurfaceContact::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUT("Entering PointSurfaceContact::AssJac()" << std::endl);

	if (dDeltaL < 0) 
	{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iSupFirstPosIndex = pSup->iGetFirstPositionIndex();
	integer iSupFirstMomIndex = pSup->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iSupFirstMomIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iSupFirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iSupFirstMomIndex + iCnt + 3);
		WM.PutColIndex(6 + iCnt, iSupFirstPosIndex + iCnt + 3);
	}

	AssMat(WM, dCoef);
	}
	else
	WorkMat.SetNullMatrix();
	
	return WorkMat;
	
}


void
PointSurfaceContact::AssMat(FullSubMatrixHandler& WM, doublereal dCoef)
{
	Mat3x3 MTmp = n.Tens(n*(-dCoef*ElasticStiffness));

	WM.Add(1, 1, MTmp);
	WM.Sub(1, 3 + 1, MTmp);
	WM.Sub(3 + 1, 1, MTmp);
	WM.Add(3 + 1, 3 + 1, MTmp);

	Vec3 dFnode1_1 = -n*ElasticStiffness*n[0];
 	Vec3 dFnode1_2 = -n*ElasticStiffness*n[1];
 	Vec3 dFnode1_3 = -n*ElasticStiffness*n[2];

	// forza sul nodo pNode1	 
 	WM.Add(1,1,dFnode1_1);
	WM.Add(1,2,dFnode1_2);
	WM.Add(1,3,dFnode1_3);	
	
	WM.Sub(1,4,dFnode1_1);
	WM.Sub(1,5,dFnode1_2);
	WM.Sub(1,6,dFnode1_3);
 
 	// forza sulla superficie
 	WM.Sub(3+1,1,dFnode1_1);
	WM.Sub(3+1,2,dFnode1_2);
	WM.Sub(3+1,3,dFnode1_3);	
	
	WM.Add(3+1,4,dFnode1_1);
	WM.Add(3+1,5,dFnode1_2);
	WM.Add(3+1,6,dFnode1_3);
 
	// momento sulla superficie
	Vec3 dFarm_1 = n;//(1.,0.,0.,);//-n*n[1];
	Vec3 dFarm_2 = n;//(0.,1.,0.,);//-n*n[2];
	Vec3 dFarm_3 = n;//(0.,0.,1.,);//-n*n[3];

	WM.Add(6+1,1,dFarm_1.Cross(FSup)+Farm.Cross(-dFnode1_1));
	WM.Add(6+1,2,dFarm_2.Cross(FSup)+Farm.Cross(-dFnode1_2));
	WM.Add(6+1,3,dFarm_3.Cross(FSup)+Farm.Cross(-dFnode1_3));
 
 	WM.Sub(6+1,4,dFarm_1.Cross(FSup)+Farm.Cross(-dFnode1_1));
	WM.Sub(6+1,5,dFarm_2.Cross(FSup)+Farm.Cross(-dFnode1_2));
	WM.Sub(6+1,6,dFarm_3.Cross(FSup)+Farm.Cross(-dFnode1_3));
	

  }
 
 
 /* Contributo al residuo durante l'assemblaggio iniziale */
SubVectorHandler&
PointSurfaceContact::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering PointSurfaceContact::InitialAssRes()" << std::endl);

		/* Dimensiona e resetta la matrice di lavoro */
		integer iNumRows = 0;
		integer iNumCols = 0;
		WorkSpaceDim(&iNumRows, &iNumCols);
		WorkVec.ResizeReset(iNumRows);
		
		/* Recupera gli indici globali */
		integer iNodePFirstMomIndex = pNode1->iGetFirstMomentumIndex();
		integer iNodeSupFirstMomIndex = pSup->iGetFirstMomentumIndex();
		
		/* Setta gli indici della matrice */
		for (int iCnt = 1; iCnt <=3; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iNodePFirstMomIndex + iCnt);
		WorkVec.PutRowIndex(iCnt+3, iNodeSupFirstMomIndex + iCnt);
		WorkVec.PutRowIndex(iCnt+6, iNodeSupFirstMomIndex + iCnt+3);
		}
		
		/* Costruisco WorkVec in AssVec */
		AssVec(WorkVec, 1);

		return WorkVec;
}

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
PointSurfaceContact::InitialAssJac(VariableSubMatrixHandler& WorkMat, 
			const VectorHandler& XCurr)
{
	DEBUGCOUT("Entering PointSurfaceContact::Initial AssJac()" << std::endl);

	if (dDeltaL < 0)
	{
	FullSubMatrixHandler& WM = WorkMat.SetFull();

	/* Dimensiona e resetta la matrice di lavoro */
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	/* Recupera gli indici */
	integer iNode1FirstPosIndex = pNode1->iGetFirstPositionIndex();
	integer iNode1FirstMomIndex = pNode1->iGetFirstMomentumIndex();
	integer iSupFirstPosIndex = pSup->iGetFirstPositionIndex();
	integer iSupFirstMomIndex = pSup->iGetFirstMomentumIndex();

	/* Setta gli indici della matrice */
	for (int iCnt = 1; iCnt <= 3; iCnt++) {
		WM.PutRowIndex(iCnt, iNode1FirstMomIndex + iCnt);
		WM.PutColIndex(iCnt, iNode1FirstPosIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iSupFirstMomIndex + iCnt);
		WM.PutColIndex(3 + iCnt, iSupFirstPosIndex + iCnt);
		WM.PutRowIndex(6 + iCnt, iSupFirstMomIndex + iCnt + 3);
		WM.PutColIndex(6 + iCnt, iSupFirstPosIndex + iCnt + 3);
	}

	AssMat(WM, 1);

	}
	else
	WorkMat.SetNullMatrix();
	
	return WorkMat;
						
}

const OutputHandler::Dimensions
PointSurfaceContact::GetEquationDimension(integer index) const {
	// DOF == 0
	return OutputHandler::Dimensions::UnknownDimension;
}
