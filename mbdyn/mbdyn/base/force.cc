/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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

/* Forze */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <force.h>
#include <dataman.h>
#include <float.h>

/* Force - begin */

std::ostream& Force::Output(unsigned int NodeLabel, std::ostream& out) const
{
   return out << std::setw(8) << GetLabel() << " "
     << std::setw(8) << NodeLabel << " " 
     << (pGetDriveCaller()->dGet());
}


std::ostream& Force::Restart(std::ostream& out) const
{
   return out << "  force: " << GetLabel();
}

/* Force - end */


/* StructuralForce - begin */

/* Costruttore */
StructuralForce::StructuralForce(unsigned int uL, 
				 Force::Type T,
				 const StructNode* pN,
				 const DriveCaller* pDC, 
				 const Vec3& TmpDir,
				 flag fOut)
: Elem(uL, Elem::FORCE, fOut), 
Force(uL, T, pDC, fOut), 
pNode(pN), Dir(TmpDir)
{ 
   ASSERT(pNode != NULL);
   ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);
   ASSERT(pDC != NULL);
   ASSERT(Dir.Dot() > 0.);
}


StructuralForce::~StructuralForce(void) 
{ 
   NO_OP; 
};

/* StructuralForce - end */


/* AbstractForce - begin */

/* Costruttore non banale */

AbstractForce::AbstractForce(unsigned int uL, const Node* pN, 
			     const DriveCaller* pDC, flag fOut)
: Elem(uL, Elem::FORCE, fOut),
Force(uL, Force::ABSTRACTFORCE, pDC, fOut),
pNode(pN)
{
   NO_OP;
}


/* Contributo al file di restart */
std::ostream& AbstractForce::Restart(std::ostream& out) const
{
   Force::Restart(out) << ", abstract, " 
     << pNode->GetLabel() << ", " 
     << psReadNodesNodes[pNode->GetNodeType()] << ", " /* << iDofNumber << ", " */;
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;     
}


/* Assembla il residuo */
SubVectorHandler& AbstractForce::AssRes(SubVectorHandler& WorkVec,
					doublereal /* dCoef */ ,
					const VectorHandler& /* XCurr */ ,
					const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering AbstractForce::AssRes()" << std::endl);

   WorkVec.Resize(1);
   WorkVec.Reset(0.);

   /* Dati */   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   // integer iFirstIndex = pNode->iGetFirstRowIndex()+iDofNumber;
   integer iFirstIndex = pNode->iGetFirstRowIndex();
   WorkVec.fPutRowIndex(1, iFirstIndex+1);

   WorkVec.fPutCoef(1, dAmplitude);

   return WorkVec;
}

void AbstractForce::Output(OutputHandler& OH) const 
{
   if (fToBeOutput()) {
      Force::Output(pNode->GetLabel(), OH.Forces()) << std::endl;
   }
}

/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
AbstractForce::InitialAssRes(SubVectorHandler& WorkVec,
			     const VectorHandler& XCurr)
{
   DEBUGCOUT("Entering AbstractForce::InitialAssRes()" << std::endl);

   return AssRes(WorkVec, 1., XCurr, XCurr);
}

/* AbstractForce - end */


/* ConservativeForce - begin */

/* Costruttore non banale */

ConservativeForce::ConservativeForce(unsigned int uL, const StructNode* pN, 
				     const DriveCaller* pDC,
				     const Vec3& TmpDir, const Vec3& TmpArm,
				     flag fOut)
: Elem(uL, Elem::FORCE, fOut), 
StructuralForce(uL, Force::CONSERVATIVEFORCE, pN, pDC, TmpDir, fOut), 
Arm(TmpArm)
{ 
   NO_OP; 
};


/* Contributo al file di restart */
std::ostream& ConservativeForce::Restart(std::ostream& out) const
{
   Force::Restart(out) << ", conservative, " 
     << pNode->GetLabel() << ", reference, global, ",
     ((pNode->GetRCurr()).Transpose()*Dir).Write(out, ", ") 
       << ", reference, node, ",
     Arm.Write(out, ", ") << ", ";
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;     
}


VariableSubMatrixHandler& 
ConservativeForce::AssJac(VariableSubMatrixHandler& WorkMat,
			  doublereal dCoef,
			  const VectorHandler& /* XCurr */ ,
			  const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering ConservativeForce::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   WM.ResizeInit(3, 3, 0.);

   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex()+3;
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex()+3;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iFirstMomentumIndex+iCnt); 
      WM.fPutColIndex(iCnt, iFirstPositionIndex+iCnt);
   }      
   
   /* Dati */
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   Vec3 TmpArm(pNode->GetRRef()*Arm);
   Vec3 TmpDir = Dir*dAmplitude;
   
   /* |    F/\   |           |   F  |
    * |          | Delta_g = |      |
    * | (d/\F)/\ |           | d/\F | */
     
   WM.Add(1, 1, Mat3x3(TmpDir, TmpArm*(-dCoef)));
   
   return WorkMat;
};


/* Assembla il residuo */
SubVectorHandler& ConservativeForce::AssRes(SubVectorHandler& WorkVec,
					    doublereal /* dCoef */ ,
					    const VectorHandler& /* XCurr */ ,
					    const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering ConservativeForce::AssRes()" << std::endl);

   integer iNumRows;
   integer iNumCols;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
   }   
   
   Mat3x3 R(pNode->GetRCurr());
   Vec3 F(Dir*dAmplitude);
   Vec3 M((R*Arm).Cross(F));
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, M);

   return WorkVec;
}


void ConservativeForce::Output(OutputHandler& OH) const 
{
   if (fToBeOutput()) {      
      Force::Output(pNode->GetLabel(), OH.Forces())
	<< " " << Dir*dGet()
	  << " " << pNode->GetXCurr()+pNode->GetRCurr()*Arm << std::endl;
   }
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
ConservativeForce::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				 const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering ConservativeForce::InitialAssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   WM.ResizeInit(6, 6, 0.);

   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex()+3;
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);     
      WM.fPutRowIndex(3+iCnt, iFirstVelocityIndex+iCnt);     
      WM.fPutColIndex(iCnt, iFirstPositionIndex+iCnt);
      WM.fPutColIndex(3+iCnt, iFirstVelocityIndex+iCnt);
   }      

   /* Dati */
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   Vec3 TmpArm(pNode->GetRRef()*Arm);
   Vec3 TmpDir = Dir*dAmplitude;
   Vec3 Omega(pNode->GetWRef());
   
   /* |    F/\   |           |   F  |
    * |          | Delta_g = |      |
    * | (d/\F)/\ |           | d/\F | */
     
   WM.Add(1, 1, Mat3x3(-TmpDir, TmpArm));
   WM.Add(4, 1, Mat3x3(-TmpDir, Omega)*Mat3x3(TmpArm));
   WM.Add(4, 4, Mat3x3(-TmpDir, TmpArm));
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ConservativeForce::InitialAssRes(SubVectorHandler& WorkVec,
				 const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering ConservativeForce::InitialAssRes()" << std::endl);

   integer iNumRows;
   integer iNumCols;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
   }   
   
   Mat3x3 R(pNode->GetRCurr());
   Vec3 TmpDir(Dir*dAmplitude);
   Vec3 TmpArm(R*Arm);
   Vec3 Omega(pNode->GetWCurr());
   
   WorkVec.Add(1, TmpDir);
   WorkVec.Add(4, TmpArm.Cross(TmpDir));
   /* In 7 non c'e' nulla */
   WorkVec.Add(10, (Omega.Cross(TmpArm)).Cross(TmpDir));

   return WorkVec;
}

/* ConservativeForce - end */


/* FollowerForce - begin */

/* Costruttore non banale */

FollowerForce::FollowerForce(unsigned int uL, const StructNode* pN, 
			     const DriveCaller* pDC,
			     const Vec3& TmpDir, const Vec3& TmpArm,
			     flag fOut)
: Elem(uL, Elem::FORCE, fOut), 
StructuralForce(uL, Force::FOLLOWERFORCE, pN, pDC, TmpDir, fOut), 
Arm(TmpArm)
{ 
   NO_OP; 
};


/* Contributo al file di restart */
std::ostream& FollowerForce::Restart(std::ostream& out) const
{
   Force::Restart(out) << ", follower, "
     << pNode->GetLabel() 
     << ", reference, node, ",
     Dir.Write(out, ", ") 
     << ", reference, node, ",
     Arm.Write(out, ", ") << ", ";
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;     
}


VariableSubMatrixHandler& 
FollowerForce::AssJac(VariableSubMatrixHandler& WorkMat,
		      doublereal dCoef,
		      const VectorHandler& /* XCurr */ ,
		      const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering FollowerForce::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   integer iFirstRotationIndex = pNode->iGetFirstPositionIndex()+3;
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iFirstMomentumIndex+iCnt);     /* forza */
      WM.fPutRowIndex(3+iCnt, iFirstMomentumIndex+3+iCnt); /* coppia */
      WM.fPutColIndex(iCnt, iFirstRotationIndex+iCnt);     /* rotazione */
   }      
   
   /* Dati */
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   Mat3x3 R(pNode->GetRRef());
   Vec3 TmpDir(R*(Dir*dAmplitude));
   Vec3 TmpArm(R*Arm);
   
   /* |    F/\   |           |   F  |
    * |          | Delta_g = |      |
    * | (d/\F)/\ |           | d/\F | */
     
   WM.Add(1, 1, Mat3x3(TmpDir*dCoef));
   WM.Add(4, 1, Mat3x3(TmpArm.Cross(TmpDir*dCoef)));   
   
   return WorkMat;
};


/* Assembla il residuo */
SubVectorHandler& FollowerForce::AssRes(SubVectorHandler& WorkVec,
					doublereal /* dCoef */ ,
					const VectorHandler& /* XCurr */ ,
					const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering FollowerForce::AssRes()" << std::endl);
   
   integer iNumRows;
   integer iNumCols;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
   }   
   
   Mat3x3 R(pNode->GetRCurr());
   Vec3 TmpDir = Dir*dAmplitude;
   Vec3 F(R*TmpDir);
   Vec3 M(R*Arm.Cross(TmpDir));
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, M);

   return WorkVec;
}


void FollowerForce::Output(OutputHandler& OH) const 
{   
   if (fToBeOutput()) {
      Force::Output(pNode->GetLabel(), OH.Forces())
	<< " " << pNode->GetRCurr()*(Dir*dGet())
	  << " " << pNode->GetXCurr()+pNode->GetRCurr()*Arm << std::endl;
   }
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
FollowerForce::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			     const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering FollowerForce::InitialAssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   WM.ResizeInit(12, 6, 0.);

   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);     
      WM.fPutRowIndex(3+iCnt, iFirstPositionIndex+3+iCnt);     
      WM.fPutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);     
      WM.fPutRowIndex(9+iCnt, iFirstVelocityIndex+3+iCnt);     
      WM.fPutColIndex(iCnt, iFirstPositionIndex+3+iCnt);
      WM.fPutColIndex(3+iCnt, iFirstVelocityIndex+3+iCnt);
   }      

   /* Dati */
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   Mat3x3 R(pNode->GetRRef());
   Vec3 TmpArm(R*Arm);
   Vec3 TmpDir = R*(Dir*dAmplitude);
   Vec3 Omega(pNode->GetWRef());
   
   /* |    F/\   |           |   F  |
    * |          | Delta_g = |      |
    * | (d/\F)/\ |           | d/\F | */
     
   WM.Add(1, 1, Mat3x3(TmpDir));
   WM.Add(4, 1, Mat3x3(TmpArm.Cross(TmpDir)));
   WM.Add(7, 1, Mat3x3(Omega, TmpDir));
   WM.Add(7, 4, Mat3x3(TmpDir));
   WM.Add(10, 1, Mat3x3(Omega, TmpArm.Cross(TmpDir)));
   WM.Add(10, 4, Mat3x3(TmpArm.Cross(TmpDir)));
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
FollowerForce::InitialAssRes(SubVectorHandler& WorkVec,
			     const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering FollowerForce::InitialAssRes()" << std::endl);

   integer iNumRows;
   integer iNumCols;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);
      WorkVec.fPutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
   }   
   
   Mat3x3 R(pNode->GetRCurr());
   Vec3 TmpDir(R*(Dir*dAmplitude));
   Vec3 TmpArm(R*Arm);
   Vec3 Omega(pNode->GetWCurr());
   
   WorkVec.Add(1, TmpDir);
   WorkVec.Add(4, TmpArm.Cross(TmpDir));
   WorkVec.Add(7, Omega.Cross(TmpDir));
   WorkVec.Add(10, (Omega.Cross(TmpArm)).Cross(TmpDir)
	       +TmpArm.Cross(Omega.Cross(TmpDir)));

   return WorkVec;
}

/* FollowerForce - end */


/* ConservativeCouple - begin */

/* Costruttore non banale */

ConservativeCouple::ConservativeCouple(unsigned int uL, const StructNode* pN, 
				       const DriveCaller* pDC, 
				       const Vec3& TmpDir,
				       flag fOut)
: Elem(uL, Elem::FORCE, fOut), 
StructuralForce(uL, Force::CONSERVATIVECOUPLE, pN, pDC, TmpDir, fOut)
{ 
   NO_OP; 
};


/* Contributo al file di restart */
std::ostream& ConservativeCouple::Restart(std::ostream& out) const
{
   out << "  couple: " << GetLabel() << ", conservative, " 
     << pNode->GetLabel() << ", reference, global, ",
     Dir.Write(out, ", ")
       << ", ";
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;     
}


/* Assembla il residuo */
SubVectorHandler& ConservativeCouple::AssRes(SubVectorHandler& WorkVec,
					     doublereal /* dCoef */ ,
					     const VectorHandler& /* XCurr */ ,
					     const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering ConservativeCouple::AssRes()" << std::endl);

   integer iNumRows;
   integer iNumCols;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex()+3;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
   }   
   
   WorkVec.Add(1, Dir*dAmplitude);

   return WorkVec;
}


void ConservativeCouple::Output(OutputHandler& OH) const 
{   
   if (fToBeOutput()) {
      Force::Output(pNode->GetLabel(), OH.Forces())
	<< " " << Dir*dGet() << std::endl;
   }
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ConservativeCouple::InitialAssRes(SubVectorHandler& WorkVec,
				  const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering ConservativeCouple::InitialAssRes()" << std::endl);

   WorkVec.Resize(3);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex()+3;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);
   }   
   
   WorkVec.Add(1, Dir*dAmplitude);

   return WorkVec;
}

/* ConservativeCouple - end */


/* FollowerCouple - begin */

FollowerCouple::FollowerCouple(unsigned int uL, const StructNode* pN, 
			       const DriveCaller* pDC, const Vec3& TmpDir,
			       flag fOut)
: Elem(uL, Elem::FORCE, fOut), 
StructuralForce(uL, Force::FOLLOWERCOUPLE, pN, pDC, TmpDir, fOut)
{ 
   NO_OP; 
};


/* Contributo al file di restart */
std::ostream& FollowerCouple::Restart(std::ostream& out) const
{
   out << "  couple: " << GetLabel() << ", follower, " 
     << pNode->GetLabel() << ", reference, node, ",
     Dir.Write(out, ", ") << ", ";
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;
}


VariableSubMatrixHandler& 
FollowerCouple::AssJac(VariableSubMatrixHandler& WorkMat,
		       doublereal dCoef,
		       const VectorHandler& /* XCurr */ ,
		       const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering FollowerCouple::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   integer iFirstRotationIndex = pNode->iGetFirstPositionIndex()+3;
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex()+3;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iFirstMomentumIndex+iCnt);    /* coppia */
      WM.fPutColIndex(iCnt, iFirstRotationIndex+iCnt);    /* rotazione */
   }      
   
   /* Dati */
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   Mat3x3 R(pNode->GetRRef());
   Mat3x3 MWedge((R*Dir)*(dAmplitude*dCoef));
   
   /* | M /\| Delta_g = | M | */
     
   WM.Add(1, 1, MWedge);   
   
   return WorkMat;
};


/* Assembla il residuo */
SubVectorHandler& FollowerCouple::AssRes(SubVectorHandler& WorkVec,
					 doublereal /* dCoef */ ,
					 const VectorHandler& /* XCurr */ , 
					 const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering FollowerCouple::AssRes()" << std::endl);
   
   integer iNumRows;
   integer iNumCols;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex()+3;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
   }
       
   Mat3x3 R(pNode->GetRCurr());
   WorkVec.Add(1, (R*Dir)*dAmplitude);
   
   return WorkVec;
}


void FollowerCouple::Output(OutputHandler& OH) const 
{   
   if (fToBeOutput()) {
      Force::Output(pNode->GetLabel(), OH.Forces())
	<< " " << pNode->GetRCurr()*(Dir*dGet()) << std::endl;
   }
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
FollowerCouple::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			      const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering FollowerCouple::InitialAssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   WM.ResizeInit(6, 6, 0.);

   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex()+3;
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);     
      WM.fPutRowIndex(3+iCnt, iFirstVelocityIndex+iCnt);     
      WM.fPutColIndex(iCnt, iFirstPositionIndex+iCnt);
      WM.fPutColIndex(3+iCnt, iFirstVelocityIndex+iCnt);
   }      
   
   /* Dati */
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   Vec3 TmpDir(pNode->GetRRef()*(Dir*dAmplitude));
   Vec3 Omega(pNode->GetWRef());
   
   /* |    F/\   |           |   F  |
    * |          | Delta_g = |      |
    * | (d/\F)/\ |           | d/\F | */
     
   WM.Add(1, 1, Mat3x3(TmpDir));
   WM.Add(4, 1, Mat3x3(Omega, TmpDir));
   WM.Add(4, 4, Mat3x3(TmpDir));
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
FollowerCouple::InitialAssRes(SubVectorHandler& WorkVec,
			      const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering FollowerCouple::InitialAssRes()" << std::endl);

   WorkVec.Resize(6);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex()+3;
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   for(integer iCnt = 1; iCnt <= 3; iCnt++)
     {
	WorkVec.fPutRowIndex(iCnt, iFirstPositionIndex+iCnt);
	WorkVec.fPutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
     }   
   
   Mat3x3 R(pNode->GetRCurr());
   Vec3 TmpDir(R*(Dir*dAmplitude));
   Vec3 Omega(pNode->GetWCurr());
   
   WorkVec.Add(1, TmpDir);
   WorkVec.Add(4, Omega.Cross(TmpDir));

   return WorkVec;
}

/* FollowerCouple - end */


/* Legge una forza */

Elem* ReadForce(DataManager* pDM, 
		MBDynParser& HP, 
		unsigned int uLabel, 
		flag fCouple)
{
   const char sFuncName[] = "ReadForce()";
   DEBUGCOUT("Entering " << sFuncName << std::endl);
   
   const char* sKeyWords[] = {
#if defined(USE_STRUCT_NODES)
      "conservative",
	"follower"
#if defined(USE_ELECTRIC_NODES)
	, "abstract"
#endif
#elif defined(USE_ELECTRIC_NODES)
	"abstract"
#endif
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
#if defined(USE_STRUCT_NODES)
	CONSERVATIVE = 0,
	FOLLOWER,
#if defined(USE_ELECTRIC_NODES)
	ABSTRACT,
#endif
#elif defined(USE_ELECTRIC_NODES)
	ABSTRACT = 0,
#endif
	LASTKEYWORD 
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);
   
   /* lettura dati specifici */
   
   /* tipo di forza */
   KeyWords CurrType = KeyWords(HP.GetWord());
   if (CurrType == UNKNOWN) {
      std::cerr << std::endl << sFuncName << " at line " << HP.GetLineData() 
	<< ": unknown force type" << std::endl;      
      THROW(DataManager::ErrGeneric());
   }   
   
#ifdef DEBUG
#if defined(USE_STRUCT_NODES)
   if (CurrType == CONSERVATIVE) {      
      std::cout << "Force type: \"Conservative\"" << std::endl;
   } else if (CurrType == FOLLOWER) {      
      std::cout << "Force type: \"Follower\"" << std::endl;
   }
#if defined(USE_ELECTRIC_NODES)
   else if (CurrType == ABSTRACT) {
      std::cout << "Force type: \"Abstract\"" << std::endl;
   }
#endif
#endif
#endif // DEBUG

   Elem* pEl = NULL;
   
#if defined(USE_ELECTRIC_NODES)
   if (CurrType == ABSTRACT) {
      
      /* tabella delle parole chiave */
      KeyTable KDof((int)Node::LASTNODETYPE, psReadNodesNodes);
      HP.PutKeyTable(KDof);
      
      ScalarDof SD = ReadScalarDof(pDM, HP, 0);                          
      HP.PutKeyTable(KDof);
      
      DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
      HP.PutKeyTable(K);

      flag fOut = pDM->fReadOutput(HP, Elem::FORCE);
      
      SAFENEWWITHCONSTRUCTOR(pEl,
			     AbstractForce,
			     AbstractForce(uLabel, SD.pNode, pDC, fOut));
            
   } else
     
#endif // USE_ELECTRIC_NODES
     
#if defined(USE_STRUCT_NODES)
     { /* CurrType e' <structural> */
         
      /* nodo collegato */
      unsigned int uNode = (unsigned int)HP.GetInt();
      
      DEBUGCOUT("Linked to Node " << uNode << std::endl);
      
      /* verifica di esistenza del nodo */
      StructNode* pNode;
      if ((pNode = pDM->pFindStructNode(uNode)) == NULL) {
	 std::cerr << std::endl << sFuncName
	   << " at line " << HP.GetLineData() 
	   << ": structural node " << uNode
	   << " not defined" << std::endl;	 
	 THROW(DataManager::ErrGeneric());
      }		        
      
      /* direzione della forza */   	     
      Mat3x3 RNode(pNode->GetRCurr());     
      ReferenceFrame RF(pNode);
      Vec3 Dir(HP.GetVecRel(RF));

      /* Se la forza e' conservativa, viene passata nel sistema globale */
      if (CurrType == CONSERVATIVE) {	 
	 Dir = RNode*Dir;
      }      	           
      
      /* Normalizza la direzione */
      ASSERT(Dir.Dot() > DBL_EPSILON);
      doublereal d = Dir.Dot();
      if (d > DBL_EPSILON) {      
	 Dir /= sqrt(d);
      } else {      
	 std::cerr << "Warning, force " << uLabel 
	   << " has null direction" << std::endl;
      }        
      
      /* distanza dal nodo (vettore di 3 elementi) ( solo se e' una forza) */
      Vec3 Arm(0.);
      if (fCouple == 0) {	
	 if (!HP.IsKeyWord("null")) {	  
	    Arm = HP.GetPosRel(RF);	    
	    DEBUGCOUT("Distance is supplied" << std::endl);
	 }
      }   
            
            
#ifdef DEBUG
      if (CurrType == CONSERVATIVE) {      
	 std::cout << "Global reference frame direction: " << std::endl << Dir << std::endl;
      } else if (CurrType == FOLLOWER) {      
	 std::cout << "Node reference frame direction: " << std::endl << Dir << std::endl;
      }
      
      if (!fCouple) {      
	 std::cout << "Node reference frame arm: " << std::endl << Arm << std::endl;
      }   
#endif   
      
            
      DriveCaller* pDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
      HP.PutKeyTable(K);
      
      flag fOut = pDM->fReadOutput(HP, Elem::FORCE);
      
      /* Alloca la forza */
      if (fCouple == 0) {	
	 if (CurrType == CONSERVATIVE) {
	    SAFENEWWITHCONSTRUCTOR(pEl, 
				   ConservativeForce,
				   ConservativeForce(uLabel, pNode, pDC, 
						     Dir, Arm, fOut));
	 } else if (CurrType == FOLLOWER) {	     
	    SAFENEWWITHCONSTRUCTOR(pEl, 
				   FollowerForce,
				   FollowerForce(uLabel, pNode, pDC, 
						 Dir, Arm, fOut));
	 }	
      } else if (fCouple == 1) {
	 if (CurrType == CONSERVATIVE) {
	    SAFENEWWITHCONSTRUCTOR(pEl, 
				   ConservativeCouple,
				   ConservativeCouple(uLabel, pNode, pDC, 
						      Dir, fOut));
	 } else if (CurrType == FOLLOWER) {	  	
	    SAFENEWWITHCONSTRUCTOR(pEl, 
				   FollowerCouple,
				   FollowerCouple(uLabel, pNode, pDC, 
						  Dir, fOut));
	 }	
      }      
   }
#endif // USE_STRUCT_NODES
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      std::cerr << std::endl << sFuncName
	<< ": semicolon expected at line " << HP.GetLineData() << std::endl;      
      THROW(DataManager::ErrGeneric());
   }      
   
   return pEl;
} /* End of ReadForce() */
