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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>

#include "accj.h"

/* LinearAccelerationJoint - begin */

/* Costruttore non banale */
LinearAccelerationJoint::LinearAccelerationJoint(unsigned int uL, 
						 const DofOwner* pDO, 
						 const StructNode* pN,
						 const Vec3& TmpDir,
						 const DriveCaller* pDC,
						 flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
DriveOwner(pDC),
pNode(pN),
Dir(TmpDir),
#ifdef USE_NETCDFC   // Netcdf4 has non-pointer variables...
Var_a(0),
#endif // USE_NETCDFC
dF(0.)
{
   ASSERT(pNode != NULL);
   ASSERT(Dir.Norm() > std::numeric_limits<doublereal>::epsilon());
   Dir /= Dir.Norm();
}


/* Distruttore */
LinearAccelerationJoint::~LinearAccelerationJoint(void)
{
   NO_OP;
}


/* Tipo di Joint */
Joint::Type LinearAccelerationJoint::GetJointType(void) const
{
   return Joint::LINEARACCELERATION;
}


/* Contributo al file di restart */
std::ostream& LinearAccelerationJoint::Restart(std::ostream& out) const
{
   return out << "Not implemented yet!" << std::endl;
}


unsigned int LinearAccelerationJoint::iGetNumDof(void) const
{
   return 2;
}


DofOrder::Order LinearAccelerationJoint::GetDofType(unsigned int i) const
{
   ASSERT(i == 0 || i == 1);
   switch (i) {
    case 0:
      return DofOrder::DIFFERENTIAL;
    case 1:
      return DofOrder::ALGEBRAIC;
    default:
      silent_cerr("invalid dof number" << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
}


void LinearAccelerationJoint::WorkSpaceDim(integer* piNumRows, 
					   integer* piNumCols) const
{
   *piNumRows = 5;
   *piNumCols = 5;
}

      
VariableSubMatrixHandler& 
LinearAccelerationJoint::AssJac(VariableSubMatrixHandler& WorkMat,
				doublereal dCoef,
				const VectorHandler& /* XCurr */ , 
				const VectorHandler& /* XPrimeCurr */ )
{
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.Resize(8, 0);
   
   integer iNodeRowIndex = pNode->iGetFirstRowIndex();
   integer iNodeColIndex = pNode->iGetFirstColIndex();
   integer iIndex = iGetFirstIndex();
      
   doublereal d = Dir.dGet(1);
   WorkMat.PutItem(1, iIndex+1, iNodeColIndex+1, d);
   WorkMat.PutItem(2, iNodeRowIndex+1, iIndex+2, d);
   d = Dir.dGet(2);
   WorkMat.PutItem(3, iIndex+1, iNodeColIndex+2, d);
   WorkMat.PutItem(4, iNodeRowIndex+2, iIndex+2, d);
   d = Dir.dGet(3);
   WorkMat.PutItem(5, iIndex+1, iNodeColIndex+3, d);
   WorkMat.PutItem(6, iNodeRowIndex+3, iIndex+2, d);
   
   WorkMat.PutItem(7, iIndex+1, iIndex+1, -dCoef);
   WorkMat.PutItem(8, iIndex+2, iIndex+1, 1.);   
   
   return WorkMat;
}


SubVectorHandler& 
LinearAccelerationJoint::AssRes(SubVectorHandler& WorkVec,
				doublereal /* dCoef */ ,
				const VectorHandler& XCurr, 
				const VectorHandler& XPrimeCurr)
{
   WorkVec.Resize(5);
   
   integer iNodeRowIndex = pNode->iGetFirstRowIndex();
   integer iIndex = iGetFirstIndex();
   
   doublereal dQ = XCurr(iIndex+1);
   doublereal dQP = XPrimeCurr(iIndex+1);
   dF = XCurr(iIndex+2);
   
   Vec3 V = pNode->GetVCurr();
   
   WorkVec.PutItem(1, iNodeRowIndex+1, -dF*Dir.dGet(1));
   WorkVec.PutItem(2, iNodeRowIndex+2, -dF*Dir.dGet(2));
   WorkVec.PutItem(3, iNodeRowIndex+3, -dF*Dir.dGet(3));
   WorkVec.PutItem(4, iIndex+1, dQ-Dir.Dot(V));
   WorkVec.PutItem(5, iIndex+2, dGet()-dQP);
   
   return WorkVec;
}

void
LinearAccelerationJoint::OutputPrepare(OutputHandler &OH)
{
	if (bToBeOutput()) {
#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("Linear acceleration", OH, name);

			Var_a = OH.CreateVar<doublereal>(name + "a", "m/s^2",
				"imposed acceleration (x, y, z)");
		}
#endif // USE_NETCDF
	}
}

void LinearAccelerationJoint::Output(OutputHandler& OH) const
{
	if (bToBeOutput()) {
		if (OH.UseText(OutputHandler::JOINTS)) {
			Joint::Output(OH.Joints(), "LinearAcc", GetLabel(), 
		 		Vec3(dF, 0., 0.), Zero3, Dir*dF, Zero3) 
     			<< " " << dGet() << std::endl;
		}
	#ifdef USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
   			Joint::NetCDFOutput(OH, Vec3(dF, 0., 0.), Zero3, Dir*dF, Zero3);
			OH.WriteNcVar(Var_a, dGet());
		}
	#endif // USE_NETCDF
	}
}
 
void
LinearAccelerationJoint::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* Xp */ ,
	SimulationEntity::Hints *ph)
{
	// TODO: hints (e.g. get drive, get orientation)
	integer iIndex = iGetFirstIndex();

	// inherit initial velocity from node's
	const Vec3& V(pNode->GetVCurr());
	X(iIndex + 1) = Dir*V;
}

/* funzioni usate nell'assemblaggio iniziale */

unsigned int LinearAccelerationJoint::iGetInitialNumDof(void) const
{
   return 0;
}


void LinearAccelerationJoint::InitialWorkSpaceDim(integer* piNumRows,
						  integer* piNumCols) const
{
   *piNumRows = 0;
   *piNumCols = 0; 
}
 

/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
LinearAccelerationJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				       const VectorHandler& /* XCurr */ )
{
   WorkMat.SetNullMatrix();
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
LinearAccelerationJoint::InitialAssRes(SubVectorHandler& WorkVec,
				       const VectorHandler& /* XCurr */ )
{
   WorkVec.Resize(0);
   return WorkVec;
}

/* Dati privati */
unsigned int
LinearAccelerationJoint::iGetNumPrivData(void) const
{
   return 2;
}

/* Dati privati */
unsigned int
LinearAccelerationJoint::iGetPrivDataIdx(const char *s) const
{
   ASSERT(s != NULL);

   if (strcmp(s, "F") == 0) {
	   return 1;
   }

   if (strcmp(s, "a") == 0) {
	   return 2;
   }

   return 0;
}

doublereal
LinearAccelerationJoint::dGetPrivData(unsigned int i) const
{
   switch (i) {
    case 1:
      return dF;
    case 2:
      return dGet();
    default:
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
}

/* LinearAccelerationJoint - end */


/* AngularAccelerationJoint - begin */

/* Costruttore non banale */
AngularAccelerationJoint::AngularAccelerationJoint(unsigned int uL, 
						   const DofOwner* pDO,
						   const StructNode* pN,
						   const Vec3& TmpDir,
						   const DriveCaller* pDC,
						   flag fOut)
: Elem(uL, fOut),
Joint(uL, pDO, fOut),
DriveOwner(pDC),
pNode(pN),
Dir(TmpDir),
dM(0.)
{
   ASSERT(pNode != NULL);
   ASSERT(Dir.Norm() > std::numeric_limits<doublereal>::epsilon());
   Dir /= Dir.Norm();   
}


/* Distruttore */
AngularAccelerationJoint::~AngularAccelerationJoint(void)
{
   NO_OP;
}

   
/* Tipo di Joint */
Joint::Type AngularAccelerationJoint::GetJointType(void) const
{
   return Joint::ANGULARACCELERATION;
}


/* Contributo al file di restart */
std::ostream& AngularAccelerationJoint::Restart(std::ostream& out) const
{
   return out << "Not implemented yet!" << std::endl;
}


unsigned int AngularAccelerationJoint::iGetNumDof(void) const
{
   return 2;
}


DofOrder::Order AngularAccelerationJoint::GetDofType(unsigned int i) const
{
   ASSERT(i == 0 || i == 1);
   switch (i) {
    case 0:
      return DofOrder::DIFFERENTIAL;
    case 1:
      return DofOrder::ALGEBRAIC;
    default:
      silent_cerr("invalid dof number" << std::endl);
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
}


void AngularAccelerationJoint::WorkSpaceDim(integer* piNumRows, 
					    integer* piNumCols) const
{
   *piNumRows = 5;
   *piNumCols = 5;
}
   
      
VariableSubMatrixHandler& 
AngularAccelerationJoint::AssJac(VariableSubMatrixHandler& WorkMat,
				 doublereal dCoef,
				 const VectorHandler& /* XCurr */ , 
				 const VectorHandler& /* XPrimeCurr */ )
{   
   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.Resize(8, 0);
   
   integer iNodeRowIndex = pNode->iGetFirstRowIndex();
   integer iNodeColIndex = pNode->iGetFirstColIndex();
   integer iIndex = iGetFirstIndex();

   Vec3 TmpDir = pNode->GetRRef()*Dir;
   
   doublereal d = TmpDir.dGet(1);
   WorkMat.PutItem(1, iIndex+1, iNodeColIndex+4, d);
   WorkMat.PutItem(2, iNodeRowIndex+4, iIndex+2, d);
   d = TmpDir.dGet(2);
   WorkMat.PutItem(3, iIndex+1, iNodeColIndex+5, d);
   WorkMat.PutItem(4, iNodeRowIndex+5, iIndex+2, d);
   d = TmpDir.dGet(3);
   WorkMat.PutItem(5, iIndex+1, iNodeColIndex+6, d);
   WorkMat.PutItem(6, iNodeRowIndex+6, iIndex+2, d);     
   
   WorkMat.PutItem(7, iIndex+1, iIndex+1, -dCoef);
   WorkMat.PutItem(8, iIndex+2, iIndex+1, 1.);      
   
   return WorkMat;
}


SubVectorHandler& 
AngularAccelerationJoint::AssRes(SubVectorHandler& WorkVec,
				 doublereal /* dCoef */ ,
				 const VectorHandler& XCurr, 
				 const VectorHandler& XPrimeCurr)
{   
   WorkVec.Resize(5);
   
   integer iNodeRowIndex = pNode->iGetFirstRowIndex();
   integer iIndex = iGetFirstIndex();
   
   doublereal dQ = XCurr(iIndex+1);
   doublereal dQP = XPrimeCurr(iIndex+1);
   dM = XCurr(iIndex+2);
      
   Vec3 W = pNode->GetWCurr();
   Vec3 TmpDir = pNode->GetRCurr()*Dir;
   
   WorkVec.PutItem(1, iNodeRowIndex+4, -dM*TmpDir.dGet(1));
   WorkVec.PutItem(2, iNodeRowIndex+5, -dM*TmpDir.dGet(2));
   WorkVec.PutItem(3, iNodeRowIndex+6, -dM*TmpDir.dGet(3));
   WorkVec.PutItem(4, iIndex+1, dQ-TmpDir.Dot(W));
   WorkVec.PutItem(5, iIndex+2, dGet()-dQP);
   
   return WorkVec;
}

   
void AngularAccelerationJoint::Output(OutputHandler& OH) const
{
   Joint::Output(OH.Joints(), "AngularAcc", GetLabel(), 
		 Zero3, Vec3(dM, 0., 0.), Zero3, Dir*dM) 
     << " " << dGet() << std::endl;   
#ifdef USE_NETCDF
   Joint::NetCDFOutput(OH, Zero3, Vec3(dM, 0., 0.), Zero3, Dir*dM);
#endif // USE_NETCDF
}
 
void
AngularAccelerationJoint::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& /* Xp */ ,
	SimulationEntity::Hints *ph)
{
	// TODO: hints (e.g. get drive, get orientation)
	integer iIndex = iGetFirstIndex();

	// inherit initial angular velocity from node's
	const Vec3& W(pNode->GetWCurr());
	Vec3 TmpDir(pNode->GetRCurr()*Dir);
	X(iIndex + 1) = TmpDir*W;
}

/* funzioni usate nell'assemblaggio iniziale */
unsigned int AngularAccelerationJoint::iGetInitialNumDof(void) const
{
   return 0;
}
 

void AngularAccelerationJoint::InitialWorkSpaceDim(integer* piNumRows,
						   integer* piNumCols) const
{
   *piNumRows = 0;
   *piNumCols = 0; 
}


/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
AngularAccelerationJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
					const VectorHandler& /* XCurr */ )
{
   WorkMat.SetNullMatrix();
   return WorkMat;
}
 

/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
AngularAccelerationJoint::InitialAssRes(SubVectorHandler& WorkVec,
					const VectorHandler& /* XCurr */ )
{
   WorkVec.Resize(0);
   return WorkVec;
}

/* Dati privati */
unsigned int
AngularAccelerationJoint::iGetNumPrivData(void) const
{
   return 2;
}

/* Dati privati */
unsigned int
AngularAccelerationJoint::iGetPrivDataIdx(const char *s) const
{
   ASSERT(s != NULL);

   if (strcmp(s, "M") == 0) {
	   return 1;
   }

   if (strcmp(s, "wp") == 0) {
	   return 2;
   }

   return 0;
}

doublereal
AngularAccelerationJoint::dGetPrivData(unsigned int i) const
{
   switch (i) {
    case 1:
      return dM;
    case 2:
      return dGet();
    default:
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
}

/* AngularAccelerationJoint - end */
