/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

#include <ac/float.h>

#include <strforce.h>
#include <dataman.h>

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
      WM.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt); 
      WM.PutColIndex(iCnt, iFirstPositionIndex+iCnt);
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
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
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
      WM.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);     
      WM.PutRowIndex(3+iCnt, iFirstVelocityIndex+iCnt);     
      WM.PutColIndex(iCnt, iFirstPositionIndex+iCnt);
      WM.PutColIndex(3+iCnt, iFirstVelocityIndex+iCnt);
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
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
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
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   integer iFirstRotationIndex = pNode->iGetFirstPositionIndex()+3;
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);     /* forza */
      WM.PutRowIndex(3+iCnt, iFirstMomentumIndex+3+iCnt); /* coppia */
      WM.PutColIndex(iCnt, iFirstRotationIndex+iCnt);     /* rotazione */
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
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
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
      WM.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);     
      WM.PutRowIndex(3+iCnt, iFirstPositionIndex+3+iCnt);     
      WM.PutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);     
      WM.PutRowIndex(9+iCnt, iFirstVelocityIndex+3+iCnt);     
      WM.PutColIndex(iCnt, iFirstPositionIndex+3+iCnt);
      WM.PutColIndex(3+iCnt, iFirstVelocityIndex+3+iCnt);
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
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
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
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex()+3;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
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
      WorkVec.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
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
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   integer iFirstRotationIndex = pNode->iGetFirstPositionIndex()+3;
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex()+3;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);    /* coppia */
      WM.PutColIndex(iCnt, iFirstRotationIndex+iCnt);    /* rotazione */
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
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex()+3;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
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
      WM.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);     
      WM.PutRowIndex(3+iCnt, iFirstVelocityIndex+iCnt);     
      WM.PutColIndex(iCnt, iFirstPositionIndex+iCnt);
      WM.PutColIndex(3+iCnt, iFirstVelocityIndex+iCnt);
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
	WorkVec.PutRowIndex(iCnt, iFirstPositionIndex+iCnt);
	WorkVec.PutRowIndex(6+iCnt, iFirstVelocityIndex+iCnt);
     }   
   
   Mat3x3 R(pNode->GetRCurr());
   Vec3 TmpDir(R*(Dir*dAmplitude));
   Vec3 Omega(pNode->GetWCurr());
   
   WorkVec.Add(1, TmpDir);
   WorkVec.Add(4, Omega.Cross(TmpDir));

   return WorkVec;
}

/* FollowerCouple - end */


/* StructuralInternalForce - begin */

/* Costruttore */
StructuralInternalForce::StructuralInternalForce(unsigned int uL, 
		Force::Type T,
		const StructNode* pN1, const StructNode* pN2,
		const DriveCaller* pDC, const Vec3& TmpDir,
		flag fOut)
: Elem(uL, Elem::FORCE, fOut), 
Force(uL, T, pDC, fOut), 
pNode1(pN1), pNode2(pN2), Dir(TmpDir)
{ 
   ASSERT(pNode1 != NULL);
   ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
   ASSERT(pNode2 != NULL);
   ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
   ASSERT(pDC != NULL);
   ASSERT(Dir.Dot() > 0.);
}


StructuralInternalForce::~StructuralInternalForce(void) 
{ 
   NO_OP; 
};

/* StructuralInternalForce - end */


/* ConservativeInternalForce - begin */

/* Costruttore non banale */

ConservativeInternalForce::ConservativeInternalForce(unsigned int uL,
		const StructNode* pN1, const StructNode* pN2, 
		const DriveCaller* pDC, const Vec3& TmpDir,
		const Vec3& TmpArm1, const Vec3& TmpArm2,
		flag fOut)
: Elem(uL, Elem::FORCE, fOut), 
StructuralInternalForce(uL, Force::CONSERVATIVEFORCE, pN1, pN2, pDC,
		TmpDir, fOut), 
Arm1(TmpArm1), Arm2(TmpArm2)
{ 
   NO_OP; 
};


/* Contributo al file di restart */
std::ostream&
ConservativeInternalForce::Restart(std::ostream& out) const
{
   Force::Restart(out) << ", conservative internal, " 
     << pNode1->GetLabel() << ", reference, global, ",
     ((pNode1->GetRCurr()).Transpose()*Dir).Write(out, ", ") 
       << ", reference, node, ",
     Arm1.Write(out, ", ") << ", "
     << pNode2->GetLabel() << ", reference, global, ",
     Arm2.Write(out, ", ") << ", ";
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;     
}


VariableSubMatrixHandler& 
ConservativeInternalForce::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering ConservativeInternalForce::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   WM.ResizeInit(6, 6, 0.);

   integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex()+3;
   integer iFirstMomentumIndex1 = pNode1->iGetFirstMomentumIndex()+3;
   
   integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex()+3;
   integer iFirstMomentumIndex2 = pNode2->iGetFirstMomentumIndex()+3;
   
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iFirstMomentumIndex1+iCnt); 
      WM.PutColIndex(iCnt, iFirstPositionIndex1+iCnt);

      WM.PutRowIndex(3+iCnt, iFirstMomentumIndex2+iCnt); 
      WM.PutColIndex(3+iCnt, iFirstPositionIndex2+iCnt);
   }      
   
   /* Dati */
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   Vec3 TmpArm1(pNode1->GetRRef()*Arm1);
   Vec3 TmpArm2(pNode2->GetRRef()*Arm2);
   Vec3 TmpDir = Dir*dAmplitude;
   
   /* |    F/\   |           |   F  |
    * |          | Delta_g = |      |
    * | (d/\F)/\ |           | d/\F | */
     
   WM.Add(1, 1, Mat3x3(TmpDir, TmpArm1*(-dCoef)));
   WM.Add(4, 4, Mat3x3(TmpDir, TmpArm2*dCoef));
   
   return WorkMat;
};


/* Assembla il residuo */
SubVectorHandler&
ConservativeInternalForce::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering ConservativeInternalForce::AssRes()" << std::endl);

   integer iNumRows;
   integer iNumCols;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstMomentumIndex1 = pNode1->iGetFirstMomentumIndex();
   integer iFirstMomentumIndex2 = pNode2->iGetFirstMomentumIndex();
   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex1+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iFirstMomentumIndex2+iCnt);
   }   
   
   Vec3 F(Dir*dAmplitude);
   Vec3 M1((pNode1->GetRCurr()*Arm1).Cross(F));
   Vec3 M2(F.Cross(pNode2->GetRCurr()*Arm1));	/* - x2 /\ F */
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, M1);
   WorkVec.Sub(7, F);
   WorkVec.Sub(10, M2);

   return WorkVec;
}


void
ConservativeInternalForce::Output(OutputHandler& OH) const 
{
   if (fToBeOutput()) {      
      Force::Output(pNode1->GetLabel(), OH.Forces())
	<< " " << Dir*dGet()
	  << " " << pNode1->GetXCurr()+pNode1->GetRCurr()*Arm1
	  << " " << pNode2->GetXCurr()+pNode2->GetRCurr()*Arm1
	  << std::endl;
   }
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
ConservativeInternalForce::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering ConservativeInternalForce::InitialAssJac()"
		   << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   WM.ResizeInit(12, 12, 0.);

   integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex()+3;
   integer iFirstVelocityIndex1 = iFirstPositionIndex1+6;

   integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex()+3;
   integer iFirstVelocityIndex2 = iFirstPositionIndex2+6;

   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iFirstPositionIndex1+iCnt);     
      WM.PutRowIndex(3+iCnt, iFirstVelocityIndex1+iCnt);     

      WM.PutColIndex(iCnt, iFirstPositionIndex1+iCnt);
      WM.PutColIndex(3+iCnt, iFirstVelocityIndex1+iCnt);

      WM.PutRowIndex(6+iCnt, iFirstPositionIndex2+iCnt);     
      WM.PutRowIndex(9+iCnt, iFirstVelocityIndex2+iCnt);     

      WM.PutColIndex(6+iCnt, iFirstPositionIndex2+iCnt);
      WM.PutColIndex(9+iCnt, iFirstVelocityIndex2+iCnt);
   }      

   /* Dati */
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   Vec3 TmpArm1(pNode1->GetRRef()*Arm1);
   Vec3 TmpArm2(pNode2->GetRRef()*Arm2);
   Vec3 TmpDir = Dir*dAmplitude;
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());
   
   /* |    F/\   |           |   F  |
    * |          | Delta_g = |      |
    * | (d/\F)/\ |           | d/\F | */
     
   WM.Sub(1, 1, Mat3x3(TmpDir, TmpArm1));
   WM.Sub(4, 1, Mat3x3(TmpDir, Omega1)*Mat3x3(TmpArm1));
   WM.Sub(4, 4, Mat3x3(TmpDir, TmpArm1));
   WM.Add(7, 7, Mat3x3(TmpDir, TmpArm2));
   WM.Add(10, 7, Mat3x3(TmpDir, Omega2)*Mat3x3(TmpArm2));
   WM.Add(10, 10, Mat3x3(TmpDir, TmpArm2));
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ConservativeInternalForce::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering ConservativeInternalForce::InitialAssRes()"
		   << std::endl);

   integer iNumRows;
   integer iNumCols;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex();
   integer iFirstVelocityIndex1 = iFirstPositionIndex1+6;
   
   integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex();
   integer iFirstVelocityIndex2 = iFirstPositionIndex2+6;
   
   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstPositionIndex1+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iFirstVelocityIndex1+iCnt);

      WorkVec.PutRowIndex(12+iCnt, iFirstPositionIndex2+iCnt);
      WorkVec.PutRowIndex(18+iCnt, iFirstVelocityIndex2+iCnt);
   }   
   
   Vec3 TmpDir(Dir*dAmplitude);
   Vec3 TmpArm1(pNode1->GetRCurr()*Arm1);
   Vec3 TmpArm2(pNode2->GetRCurr()*Arm2);
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   
   WorkVec.Add(1, TmpDir);
   WorkVec.Add(4, TmpArm1.Cross(TmpDir));

   /* In 7 non c'e' nulla */
   WorkVec.Add(10, (Omega1.Cross(TmpArm1)).Cross(TmpDir));

   WorkVec.Sub(7, TmpDir);
   WorkVec.Sub(10, TmpArm2.Cross(TmpDir));

   /* In 7 non c'e' nulla */
   WorkVec.Sub(16, (Omega2.Cross(TmpArm2)).Cross(TmpDir));

   return WorkVec;
}

/* ConservativeInternalForce - end */


/* FollowerInternalForce - begin */

/* Costruttore non banale */

FollowerInternalForce::FollowerInternalForce(unsigned int uL,
		const StructNode* pN1, const StructNode* pN2, 
		const DriveCaller* pDC, const Vec3& TmpDir,
		const Vec3& TmpArm1, const Vec3& TmpArm2,
		flag fOut)
: Elem(uL, Elem::FORCE, fOut), 
StructuralInternalForce(uL, Force::FOLLOWERFORCE, pN1, pN2, pDC, TmpDir, fOut), 
Arm1(TmpArm1), Arm2(TmpArm2)
{ 
   NO_OP; 
};


/* Contributo al file di restart */
std::ostream&
FollowerInternalForce::Restart(std::ostream& out) const
{
   Force::Restart(out) << ", follower internal, "
     << pNode1->GetLabel() 
     << ", reference, node, ",
     Dir.Write(out, ", ") 
     << ", reference, node, ",
     Arm1.Write(out, ", ") << ", "
     << pNode2->GetLabel() 
     << ", reference, node, ",
     Arm2.Write(out, ", ") << ", ";
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;     
}


VariableSubMatrixHandler& 
FollowerInternalForce::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering FollowerInternalForce::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   integer iFirstRotationIndex1 = pNode1->iGetFirstPositionIndex()+3;
   integer iFirstMomentumIndex1 = pNode1->iGetFirstMomentumIndex();
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iFirstMomentumIndex1+iCnt);     /* forza */
      WM.PutRowIndex(3+iCnt, iFirstMomentumIndex1+3+iCnt); /* coppia */
      WM.PutColIndex(iCnt, iFirstRotationIndex1+iCnt);     /* rotazione */

      WM.PutRowIndex(6+iCnt, iFirstMomentumIndex1+iCnt);   /* forza */
      WM.PutRowIndex(9+iCnt, iFirstMomentumIndex1+3+iCnt); /* coppia */
      WM.PutColIndex(3+iCnt, iFirstRotationIndex1+iCnt);   /* rotazione */
   }      
   
   /* Dati */
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   Vec3 TmpDir(pNode1->GetRRef()*(Dir*(dAmplitude*dCoef)));
   Vec3 TmpArm1(pNode1->GetRRef()*Arm1);
   Vec3 TmpArm2(pNode2->GetRRef()*Arm2);
   
   /* |    F/\       0    |             |   F   |
    * |                   | Delta_g_1 = |       |
    * | (d1/\F)/\    0    |             | d1/\F |
    * |                   | Delta_g_2 = |       |
    * | -F/\d2/\  d2/\F/\ |             | d2/\F | */
     
   WM.Add(1, 1, Mat3x3(TmpDir));
   WM.Add(4, 1, Mat3x3(TmpArm1.Cross(TmpDir)));
   WM.Sub(7, 1, Mat3x3(TmpDir));
   WM.Sub(7, 1, Mat3x3(TmpArm2, TmpDir));
   WM.Add(7, 4, Mat3x3(TmpDir, TmpArm2));
   
   return WorkMat;
};


/* Assembla il residuo */
SubVectorHandler&
FollowerInternalForce::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering FollowerInternalForce::AssRes()" << std::endl);
   
   integer iNumRows;
   integer iNumCols;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstMomentumIndex1 = pNode1->iGetFirstMomentumIndex();
   integer iFirstMomentumIndex2 = pNode2->iGetFirstMomentumIndex();
   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex1+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iFirstMomentumIndex2+iCnt);
   }   
   
   Vec3 TmpDir = Dir*dAmplitude;
   Vec3 F(pNode1->GetRCurr()*TmpDir);
   Vec3 M1(pNode1->GetRCurr()*Arm1.Cross(TmpDir));
   Vec3 M2(F.Cross(pNode2->GetRCurr()*Arm2));
   
   WorkVec.Add(1, F);
   WorkVec.Add(4, M1);
   WorkVec.Sub(7, F);
   WorkVec.Add(10, M2);

   return WorkVec;
}


void FollowerInternalForce::Output(OutputHandler& OH) const 
{   
   if (fToBeOutput()) {
      Force::Output(pNode1->GetLabel(), OH.Forces())
	<< " " << pNode1->GetRCurr()*(Dir*dGet())
	<< " " << pNode1->GetXCurr()+pNode1->GetRCurr()*Arm1
	<< " " << pNode2->GetXCurr()+pNode2->GetRCurr()*Arm2
	<< std::endl;
   }
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
FollowerInternalForce::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering FollowerInternalForce::InitialAssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   integer iNumRows;
   integer iNumCols;
   WorkSpaceDim(&iNumRows, &iNumCols);
   
   /* Dimensiona e resetta la matrice di lavoro */
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex();
   integer iFirstVelocityIndex1 = iFirstPositionIndex1+6;

   integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex();
   integer iFirstVelocityIndex2 = iFirstPositionIndex2+6;

   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iFirstPositionIndex1+iCnt);     
      WM.PutRowIndex(3+iCnt, iFirstPositionIndex1+3+iCnt);     
      WM.PutRowIndex(6+iCnt, iFirstVelocityIndex1+iCnt);     
      WM.PutRowIndex(9+iCnt, iFirstVelocityIndex1+3+iCnt);     

      WM.PutColIndex(iCnt, iFirstPositionIndex1+3+iCnt);
      WM.PutColIndex(3+iCnt, iFirstVelocityIndex1+3+iCnt);
      
      WM.PutRowIndex(12+iCnt, iFirstPositionIndex2+iCnt);     
      WM.PutRowIndex(15+iCnt, iFirstPositionIndex2+3+iCnt);     
      WM.PutRowIndex(18+iCnt, iFirstVelocityIndex2+iCnt);     
      WM.PutRowIndex(21+iCnt, iFirstVelocityIndex2+3+iCnt);     

      WM.PutColIndex(12+iCnt, iFirstPositionIndex2+3+iCnt);
      WM.PutColIndex(15+iCnt, iFirstVelocityIndex2+3+iCnt);
   }      

   /* Dati */
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   Vec3 TmpArm1(pNode1->GetRRef()*Arm1);
   Vec3 TmpArm2(pNode2->GetRRef()*Arm2);
   Vec3 TmpDir = pNode1->GetRRef()*(Dir*dAmplitude);
   Vec3 Omega1(pNode1->GetWRef());
   Vec3 Omega2(pNode2->GetWRef());
   
   WM.Add(1, 1, Mat3x3(TmpDir));
   WM.Add(4, 1, Mat3x3(TmpArm1.Cross(TmpDir)));

   WM.Add(7, 1, Mat3x3(Omega1, TmpDir));
   WM.Add(7, 4, Mat3x3(TmpDir));
   WM.Add(10, 1, Mat3x3(Omega1, TmpArm1.Cross(TmpDir)));
   WM.Add(10, 4, Mat3x3(TmpArm1.Cross(TmpDir)));
   
   WM.Sub(13, 1, Mat3x3(TmpDir));
   WM.Add(16, 1, Mat3x3(TmpArm2, TmpDir));
   WM.Sub(16, 7, Mat3x3(TmpDir, TmpArm2));

   WM.Sub(19, 1, Mat3x3(Omega1, TmpDir));
   WM.Sub(19, 4, Mat3x3(TmpDir));

   WM.Add(22, 1, Mat3x3(TmpArm2, Omega1)*Mat3x3(TmpDir)
		   -Mat3x3(Omega2.Cross(TmpArm2), TmpDir));
   WM.Add(22, 4, Mat3x3(TmpArm2, TmpDir));
   WM.Add(22, 7, Mat3x3(TmpDir, Omega2)*Mat3x3(TmpArm2)
		   -Mat3x3(Omega1.Cross(TmpDir), TmpArm2));
   WM.Add(22, 10, Mat3x3(TmpDir, TmpArm2));

   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
FollowerInternalForce::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering FollowerInternalForce::InitialAssRes()" << std::endl);

   integer iNumRows;
   integer iNumCols;
   InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex();
   integer iFirstVelocityIndex1 = iFirstPositionIndex1+6;

   integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex();
   integer iFirstVelocityIndex2 = iFirstPositionIndex2+6;

   for (integer iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstPositionIndex1+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iFirstVelocityIndex1+iCnt);

      WorkVec.PutRowIndex(12+iCnt, iFirstPositionIndex2+iCnt);
      WorkVec.PutRowIndex(18+iCnt, iFirstVelocityIndex2+iCnt);
   }   
   
   Vec3 TmpDir(pNode1->GetRCurr()*(Dir*dAmplitude));
   Vec3 TmpArm1(pNode1->GetRCurr()*Arm1);
   Vec3 TmpArm2(pNode2->GetRCurr()*Arm2);
   Vec3 Omega1(pNode1->GetWCurr());
   Vec3 Omega2(pNode2->GetWCurr());
   
   WorkVec.Add(1, TmpDir);
   WorkVec.Add(4, TmpArm1.Cross(TmpDir));
   WorkVec.Add(7, Omega1.Cross(TmpDir));
   WorkVec.Add(10, (Omega1.Cross(TmpArm1)).Cross(TmpDir)
	       +TmpArm1.Cross(Omega1.Cross(TmpDir)));

   WorkVec.Sub(13, TmpDir);
   WorkVec.Sub(16, TmpArm2.Cross(TmpDir));
   WorkVec.Sub(19, Omega1.Cross(TmpDir));
   WorkVec.Sub(22, (Omega2.Cross(TmpArm2)).Cross(TmpDir)
	       +TmpArm2.Cross(Omega1.Cross(TmpDir)));

   return WorkVec;
}

/* FollowerInternalForce - end */


/* ConservativeInternalCouple - begin */

/* Costruttore non banale */

ConservativeInternalCouple::ConservativeInternalCouple(unsigned int uL,
		const StructNode* pN1, const StructNode* pN2, 
		const DriveCaller* pDC, const Vec3& TmpDir,
		flag fOut)
: Elem(uL, Elem::FORCE, fOut), 
StructuralInternalForce(uL, Force::CONSERVATIVECOUPLE, pN1, pN2, pDC,
		TmpDir, fOut)
{ 
   NO_OP; 
};


/* Contributo al file di restart */
std::ostream&
ConservativeInternalCouple::Restart(std::ostream& out) const
{
   out << "  couple: " << GetLabel() << ", conservative internal, " 
     << pNode1->GetLabel() << ", reference, global, ",
     Dir.Write(out, ", ") << ", "
     << pNode1->GetLabel() << ", ";
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;     
}


/* Assembla il residuo */
SubVectorHandler&
ConservativeInternalCouple::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering ConservativeInternalCouple::AssRes()" << std::endl);

   integer iNumRows;
   integer iNumCols;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstMomentumIndex1 = pNode1->iGetFirstMomentumIndex()+3;
   integer iFirstMomentumIndex2 = pNode2->iGetFirstMomentumIndex()+3;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex1+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iFirstMomentumIndex2+iCnt);
   }   
   
   WorkVec.Add(1, Dir*dAmplitude);
   WorkVec.Sub(4, Dir*dAmplitude);

   return WorkVec;
}


void
ConservativeInternalCouple::Output(OutputHandler& OH) const 
{   
   if (fToBeOutput()) {
      Force::Output(pNode1->GetLabel(), OH.Forces())
	<< " " << Dir*dGet() << std::endl;
   }
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
ConservativeInternalCouple::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering ConservativeInternalCouple::InitialAssRes()" << std::endl);

   WorkVec.Resize(6);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex()+3;
   integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex()+3;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstPositionIndex1+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iFirstPositionIndex2+iCnt);
   }   
   
   WorkVec.Add(1, Dir*dAmplitude);
   WorkVec.Sub(4, Dir*dAmplitude);

   return WorkVec;
}

/* ConservativeInternalCouple - end */


/* FollowerInternalCouple - begin */

FollowerInternalCouple::FollowerInternalCouple(unsigned int uL,
		const StructNode* pN1, const StructNode* pN2, 
		const DriveCaller* pDC, const Vec3& TmpDir,
		flag fOut)
: Elem(uL, Elem::FORCE, fOut), 
StructuralInternalForce(uL, Force::FOLLOWERCOUPLE, pN1, pN2, pDC, TmpDir, fOut)
{ 
   NO_OP; 
};


/* Contributo al file di restart */
std::ostream&
FollowerInternalCouple::Restart(std::ostream& out) const
{
   out << "  couple: " << GetLabel() << ", follower internal, " 
     << pNode1->GetLabel() << ", reference, node, ",
     Dir.Write(out, ", ") << ", "
     << pNode2->GetLabel() << ", ";
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;
}


VariableSubMatrixHandler& 
FollowerInternalCouple::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering FollowerInternalCouple::AssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WM.ResizeInit(iNumRows, iNumCols, 0.);

   integer iFirstRotationIndex1 = pNode1->iGetFirstPositionIndex()+3;
   integer iFirstMomentumIndex1 = pNode1->iGetFirstMomentumIndex()+3;
   integer iFirstMomentumIndex2 = pNode2->iGetFirstMomentumIndex()+3;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iFirstMomentumIndex1+iCnt);    /* coppia */
      WM.PutColIndex(iCnt, iFirstRotationIndex1+iCnt);    /* rotazione */

      WM.PutRowIndex(3+iCnt, iFirstMomentumIndex2+iCnt);    /* coppia */
   }      
   
   /* Dati */
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   Mat3x3 MWedge((pNode1->GetRRef()*Dir)*(dAmplitude*dCoef));
   
   /* | M /\| Delta_g = | M | */
     
   WM.Add(1, 1, MWedge);   
   WM.Sub(4, 1, MWedge);   
   
   return WorkMat;
};


/* Assembla il residuo */
SubVectorHandler&
FollowerInternalCouple::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ , 
		const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering FollowerInternalCouple::AssRes()" << std::endl);
   
   integer iNumRows;
   integer iNumCols;
   WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstMomentumIndex1 = pNode1->iGetFirstMomentumIndex()+3;
   integer iFirstMomentumIndex2 = pNode2->iGetFirstMomentumIndex()+3;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex1+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iFirstMomentumIndex2+iCnt);
   }
 
   Vec3 M((pNode1->GetRCurr()*Dir)*dAmplitude);
   WorkVec.Add(1, M);
   WorkVec.Sub(4, M);
   
   return WorkVec;
}


void
FollowerInternalCouple::Output(OutputHandler& OH) const 
{   
   if (fToBeOutput()) {
      Force::Output(pNode1->GetLabel(), OH.Forces())
	<< " " << pNode1->GetRCurr()*(Dir*dGet()) << std::endl;
   }
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
FollowerInternalCouple::InitialAssJac(VariableSubMatrixHandler& WorkMat,
		const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering FollowerInternalCouple::InitialAssJac()" << std::endl);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   
   /* Dimensiona e resetta la matrice di lavoro */
   WM.ResizeInit(12, 6, 0.);

   integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex()+3;
   integer iFirstVelocityIndex1 = iFirstPositionIndex1+6;
   integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex()+3;
   integer iFirstVelocityIndex2 = iFirstPositionIndex2+6;
   for (integer iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutRowIndex(iCnt, iFirstPositionIndex1+iCnt);     
      WM.PutRowIndex(3+iCnt, iFirstVelocityIndex1+iCnt);     

      WM.PutColIndex(iCnt, iFirstPositionIndex1+iCnt);
      WM.PutColIndex(3+iCnt, iFirstVelocityIndex1+iCnt);

      WM.PutRowIndex(6+iCnt, iFirstPositionIndex2+iCnt);     
      WM.PutRowIndex(9+iCnt, iFirstVelocityIndex2+iCnt);     
   }      
   
   /* Dati */
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   Vec3 TmpDir(pNode1->GetRRef()*(Dir*dAmplitude));
   Vec3 Omega1(pNode1->GetWRef());
   
   /* |    F/\   |           |   F  |
    * |          | Delta_g = |      |
    * | (d/\F)/\ |           | d/\F | */
     
   WM.Add(1, 1, Mat3x3(TmpDir));
   WM.Add(4, 1, Mat3x3(Omega1, TmpDir));
   WM.Add(4, 4, Mat3x3(TmpDir));
   
   WM.Sub(7, 1, Mat3x3(TmpDir));
   WM.Sub(10, 1, Mat3x3(Omega1, TmpDir));
   WM.Sub(10, 4, Mat3x3(TmpDir));
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
FollowerInternalCouple::InitialAssRes(SubVectorHandler& WorkVec,
			      const VectorHandler& /* XCurr */ )
{
   DEBUGCOUT("Entering FollowerInternalCouple::InitialAssRes()" << std::endl);

   WorkVec.Resize(12);
   WorkVec.Reset(0.);

   /* Dati */
   
   doublereal dAmplitude = pGetDriveCaller()->dGet();
   
   /* Indici delle incognite del nodo */
   integer iFirstPositionIndex1 = pNode1->iGetFirstPositionIndex()+3;
   integer iFirstVelocityIndex1 = iFirstPositionIndex1+6;

   integer iFirstPositionIndex2 = pNode2->iGetFirstPositionIndex()+3;
   integer iFirstVelocityIndex2 = iFirstPositionIndex2+6;

   for(integer iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iFirstPositionIndex1+iCnt);
      WorkVec.PutRowIndex(3+iCnt, iFirstVelocityIndex1+iCnt);
      
      WorkVec.PutRowIndex(6+iCnt, iFirstPositionIndex2+iCnt);
      WorkVec.PutRowIndex(9+iCnt, iFirstVelocityIndex2+iCnt);
   }   
   
   Vec3 TmpDir(pNode1->GetRCurr()*(Dir*dAmplitude));
   Vec3 Omega1(pNode1->GetWCurr());
   
   WorkVec.Add(1, TmpDir);
   WorkVec.Add(4, Omega1.Cross(TmpDir));

   WorkVec.Sub(7, TmpDir);
   WorkVec.Sub(10, Omega1.Cross(TmpDir));

   return WorkVec;
}

/* FollowerInternalCouple - end */

