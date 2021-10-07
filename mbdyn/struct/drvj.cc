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

/* Giunti di velocita' imposta */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "drvj.h"

/* LinearVelocity - begin */

/* Costruttore non banale */
LinearVelocityJoint::LinearVelocityJoint(unsigned int uL, 
					 const DofOwner* pDO,
					 const StructNode* pN,
					 const Vec3& TmpDir,
					 const DriveCaller* pDC,
					 flag fOut)
: Elem(uL, fOut), 
Joint(uL, pDO, fOut), 
DriveOwner(pDC),
pNode(pN), Dir(TmpDir), 
dF(0.)
#ifdef USE_NETCDFC // netcdfcxx4 has non-pointer vars...
,
Var_dv(0),
Var_v(0)
#endif // USE_NETCDFC
{
   NO_OP;
}


/* Distruttore */
LinearVelocityJoint::~LinearVelocityJoint(void)
{
   NO_OP;
}
   
/* Contributo al file di restart */
std::ostream& LinearVelocityJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", linear velocity, " << pNode->GetLabel() 
     << ", reference, global, ",
      Dir.Write(out, ", ") << ", ";
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;
}
 
/* dati privati */
unsigned int
LinearVelocityJoint::iGetNumPrivData(void) const
{
	return 1;
}
 
unsigned int
LinearVelocityJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	if (strcmp(s, "v") == 0) {
		return 1;
	}

	return 0;
}

doublereal
LinearVelocityJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i == 1);

	return dGet();
}

VariableSubMatrixHandler& 
LinearVelocityJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			    doublereal /* dCoef */ ,
			    const VectorHandler& /* XCurr */ ,
			    const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering LinearVelocityJoint::AssJac()" << std::endl);

   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.ResizeReset(6, 0);
   
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Dir.dGet(iCnt);
      WM.PutItem(iCnt, iFirstMomentumIndex+iCnt,
		  iFirstReactionIndex+1, d);
      WM.PutItem(3+iCnt, iFirstReactionIndex+1,
		  iFirstPositionIndex+iCnt, d);
   }
   
   return WorkMat;
}


SubVectorHandler& 
LinearVelocityJoint::AssRes(SubVectorHandler& WorkVec,
			    doublereal /* dCoef */ ,
			    const VectorHandler& XCurr,
			    const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering LinearVelocityJoint::AssRes()" << std::endl);

   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
      
   /* Indici */
   // integer iFirstPositionIndex = pNode->iGetFirstPositionIndex();
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici del nodo */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
   }
   
   WorkVec.PutRowIndex(4, iFirstReactionIndex+1);

   /* Aggiorna i dati propri */
   dF = XCurr(iFirstReactionIndex+1);
   
   /* Recupera i dati */
   Vec3 vNode(pNode->GetVCurr());
   
   
   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Add(1, -Dir*dF);
   
   
   /* Equazione di vincolo di velocita' */
   doublereal dv0 = dGet();
   WorkVec.PutCoef(4, dv0-Dir.Dot(vNode));
   
   return WorkVec;   
}


void 
LinearVelocityJoint::OutputPrepare(OutputHandler& OH)
{
	if (bToBeOutput()) {
#if USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("Linear velocity", OH, name);

			Var_dv = OH.CreateVar<doublereal>(name + "dv",
					OutputHandler::Dimensions::Dimensionless,
					"direction of imposed velocity");

			Var_v = OH.CreateVar<doublereal>(name + "v",
					OutputHandler::Dimensions::Velocity,
					"magnitude of imposed velocity");
		}
#endif // USE_NETCDF
	}
}

void LinearVelocityJoint::Output(OutputHandler& OH) const
{
   if (bToBeOutput()) {

	   Vec3 FTmp = Vec3(dF, 0., 0.);

	   if (OH.UseText(OutputHandler::JOINTS)) {
		   Joint::Output(OH.Joints(), "LinearVelocity", GetLabel(),
				   FTmp, Zero3, Dir*dF, Zero3)
			   << " " << Dir << " " << dGet() << std::endl;
	   }

#ifdef USE_NETCDF
	   if (OH.UseNetCDF(OutputHandler::JOINTS)) {
		   Joint::NetCDFOutput(OH, FTmp, Zero3, Dir*dF, Zero3);
		   OH.WriteNcVar(Var_dv, Dir);
		   OH.WriteNcVar(Var_v, dGet());
	   }
#endif // USE_NETCDF

   }   
}
 
   
/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
LinearVelocityJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				   const VectorHandler& /* XCurr */ )
{   
   DEBUGCOUT("Entering LinearVelocityJoint::InitialAssJac()" << std::endl);

   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.ResizeReset(6, 0);
   
   integer iFirstVelocityIndex = pNode->iGetFirstPositionIndex()+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = Dir.dGet(iCnt);
      WM.PutItem(iCnt, iFirstVelocityIndex+iCnt,
		  iFirstReactionIndex+1, d);
      WM.PutItem(3+iCnt, iFirstReactionIndex+1,
		  iFirstVelocityIndex+iCnt, d);
   }
   
   return WorkMat;
}


/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
LinearVelocityJoint::InitialAssRes(SubVectorHandler& WorkVec,
				   const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering LinearVelocityJoint::InitialAssRes()" << std::endl);

   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
      
   /* Indici */
   integer iFirstVelocityIndex = pNode->iGetFirstPositionIndex()+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici del nodo */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iFirstVelocityIndex+iCnt);
   }
   
   WorkVec.PutRowIndex(4, iFirstReactionIndex+1);

   /* Aggiorna i dati propri */
   dF = XCurr(iFirstReactionIndex+1);
   
   /* Recupera i dati */
   Vec3 vNode(pNode->GetVCurr());
   
   
   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Add(1, -Dir*dF);
   
   
   /* Equazione di vincolo di velocita' */
   doublereal dv0 = dGet();
   WorkVec.PutCoef(4, dv0-Dir.Dot(vNode));
   
   return WorkVec;
}

const OutputHandler::Dimensions
LinearVelocityJoint::GetEquationDimension(integer index) const {
	// DOF == 6
	OutputHandler::Dimensions dimension = OutputHandler::Dimensions::UnknownDimension;

	switch (index)
	{
		case 1:
			dimension = OutputHandler::Dimensions::Velocity;
			break;
		case 2:
			dimension = OutputHandler::Dimensions::Velocity;
			break;
		case 3:
			dimension = OutputHandler::Dimensions::Velocity;
			break;
	}

	return dimension;
}

std::ostream&
LinearVelocityJoint::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{

	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": " <<
			"velocity error" << std::endl;

	return out;
}
/* LinearVelocity - end */


/* AngularVelocity - begin */

/* Costruttore non banale */
AngularVelocityJoint::AngularVelocityJoint(unsigned int uL,
					   const DofOwner* pDO, 
					   const StructNode* pN,
					   const Vec3& TmpDir,
					   const DriveCaller* pDC,
					   flag fOut)
: Elem(uL, fOut), 
Joint(uL, pDO, fOut), 
DriveOwner(pDC), 
pNode(pN), Dir(TmpDir), 
dM(0.)
#ifdef USE_NETCDFC // netcdfcxx4 has non-pointer vars...
,
Var_dOmega(0),
Var_w(0)
#endif // USE_NETCDFC
{
   NO_OP;
}


/* Distruttore */
AngularVelocityJoint::~AngularVelocityJoint(void)
{
   NO_OP;
}
   

/* Contributo al file di restart */
std::ostream& AngularVelocityJoint::Restart(std::ostream& out) const
{
   Joint::Restart(out) << ", angular velocity, "
     << pNode->GetLabel() << ", reference, node, ",
     Dir.Write(out, ", ") << ", ";
   return pGetDriveCaller()->Restart(out) << ';' << std::endl;
}
   
/* dati privati */
unsigned int
AngularVelocityJoint::iGetNumPrivData(void) const
{
	return 1;
}
 
unsigned int
AngularVelocityJoint::iGetPrivDataIdx(const char *s) const
{
	ASSERT(s != NULL);

	if (strcmp(s, "w") == 0) {
		return 1;
	}

	return 0;
}

doublereal
AngularVelocityJoint::dGetPrivData(unsigned int i) const
{
	ASSERT(i == 1);

	return dGet();
}

VariableSubMatrixHandler& 
AngularVelocityJoint::AssJac(VariableSubMatrixHandler& WorkMat,
			    doublereal dCoef,
			    const VectorHandler& /* XCurr */ ,
			    const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering AngularVelocityJoint::AssJac()" << std::endl);

   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.ResizeReset(12, 0);
   
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex()+3;
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex()+3;
   integer iFirstReactionIndex = iGetFirstIndex();
   
   Vec3 TmpDir(pNode->GetRRef()*Dir);
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      doublereal d = TmpDir.dGet(iCnt);
      WM.PutItem(iCnt, iFirstMomentumIndex+iCnt,
		  iFirstReactionIndex+1, d);
      WM.PutItem(3+iCnt, iFirstReactionIndex+1,
		  iFirstPositionIndex+iCnt, d);
   }
   
   WM.PutCross(7, iFirstMomentumIndex, 
		iFirstPositionIndex, TmpDir*(-dM*dCoef));
   
   return WorkMat;
}


SubVectorHandler& 
AngularVelocityJoint::AssRes(SubVectorHandler& WorkVec,
			    doublereal /* dCoef */ ,
			    const VectorHandler& XCurr, 
			    const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering AngularVelocityJoint::AssRes()" << std::endl);

   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->WorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
      
   /* Indici */
   integer iFirstMomentumIndex = pNode->iGetFirstMomentumIndex()+3;
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici del nodo */
   for (int iCnt = 1; iCnt <= 3; iCnt++) {	
      WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex+iCnt);
   }
   
   WorkVec.PutRowIndex(4, iFirstReactionIndex+1);

   /* Aggiorna i dati propri */
   dM = XCurr(iFirstReactionIndex+1);
   
   /* Recupera i dati */
   Vec3 Omega(pNode->GetWCurr());
   Vec3 TmpDir(pNode->GetRCurr()*Dir);
   
   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Add(1, TmpDir*(-dM));
   
   
   /* Equazione di vincolo di velocita' */
   doublereal dw0 = dGet();
   WorkVec.PutCoef(4, dw0-TmpDir.Dot(Omega));
   
   return WorkVec;   
}

void
AngularVelocityJoint::OutputPrepare(OutputHandler &OH)
{
	if (bToBeOutput()) {
#if USE_NETCDF
		if (OH.UseNetCDF(OutputHandler::JOINTS)) {
			std::string name;
			OutputPrepare_int("Angular velocity", OH, name);

			Var_dOmega = OH.CreateVar<doublereal>(name + "dOmega",
					OutputHandler::Dimensions::AngularVelocity,
					"magnitude imposed angular velocity");

			Var_w = OH.CreateVar<Vec3>(name + "w",
					OutputHandler::Dimensions::Dimensionless,
					"direction of imposed angular velocity");
		}
#endif // USE_NETCDF
	}
}
   
void AngularVelocityJoint::Output(OutputHandler& OH) const
{
   if(bToBeOutput()) {      
	   Vec3 Tmp(pNode->GetRCurr()*Dir);
	   Vec3 MTmp = Vec3(dM, 0., 0.);


	   if (OH.UseText(OutputHandler::JOINTS)) {
		   Joint::Output(OH.Joints(), "AngularVelocity", GetLabel(),
				   Zero3, MTmp, Zero3, Tmp*dM)
			   << " " << Tmp << " " << dGet() << std::endl;     
	   }

#ifdef USE_NETCDF
	   if (OH.UseNetCDF(OutputHandler::JOINTS)) {
		   Joint::NetCDFOutput(OH, Zero3, MTmp, Zero3, Tmp*dM);
		   OH.WriteNcVar(Var_w, Tmp);
		   OH.WriteNcVar(Var_dOmega, dGet());
	   }
#endif // USE_NETCDF

   }   
}
 
   
/* Contributo allo jacobiano durante l'assemblaggio iniziale */
VariableSubMatrixHandler& 
AngularVelocityJoint::InitialAssJac(VariableSubMatrixHandler& WorkMat,
				    const VectorHandler& /* XCurr */ )
{   
   DEBUGCOUT("Entering AngularVelocityJoint::InitialAssJac()" << std::endl);

   SparseSubMatrixHandler& WM = WorkMat.SetSparse();
   WM.ResizeReset(15, 0);
   
   integer iFirstPositionIndex = pNode->iGetFirstPositionIndex()+3;
   integer iFirstVelocityIndex = iFirstPositionIndex+6;
   integer iFirstReactionIndex = iGetFirstIndex();
   
   Vec3 Tmp(Dir.Cross(pNode->GetWRef()));
   for(int iCnt = 1; iCnt <= 3; iCnt++)
     {
	doublereal d = Dir.dGet(iCnt);
	WM.PutItem(iCnt, iFirstVelocityIndex+iCnt,
		    iFirstReactionIndex+1, d);
	WM.PutItem(3+iCnt, iFirstReactionIndex+1,
		    iFirstVelocityIndex+iCnt, d);
	WM.PutItem(6+iCnt, iFirstReactionIndex+1,
		    iFirstPositionIndex+iCnt, Tmp.dGet(iCnt));
     }
   
   WM.PutCross(10, iFirstVelocityIndex, iFirstPositionIndex, Dir*(-dM));
   
   return WorkMat;
}

   
/* Contributo al residuo durante l'assemblaggio iniziale */   
SubVectorHandler& 
AngularVelocityJoint::InitialAssRes(SubVectorHandler& WorkVec,
				    const VectorHandler& XCurr)
{   
   DEBUGCOUT("Entering AngularVelocityJoint::InitialAssRes()" << std::endl);

   /* Dimensiona e resetta la matrice di lavoro */
   integer iNumRows = 0;
   integer iNumCols = 0;
   this->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   WorkVec.ResizeReset(iNumRows);
      
   /* Indici */
   integer iFirstVelocityIndex = pNode->iGetFirstPositionIndex()+9;
   integer iFirstReactionIndex = iGetFirstIndex();
   
   /* Indici del nodo */
   for(int iCnt = 1; iCnt <= 3; iCnt++)
     {	
	WorkVec.PutRowIndex(iCnt, iFirstVelocityIndex+iCnt);
     }
   
   WorkVec.PutRowIndex(4, iFirstReactionIndex+1);

   /* Aggiorna i dati propri */
   dM = XCurr(iFirstReactionIndex+1);
   
   /* Recupera i dati */
   Vec3 Omega(pNode->GetWCurr());   
   
   /* Equazioni di equilibrio, nodo 1 */
   WorkVec.Add(1, -Dir*dM);
   
   /* Equazione di vincolo di velocita' */
   doublereal dw0 = dGet();
   WorkVec.PutCoef(4, dw0-Dir.Dot(Omega));
   
   return WorkVec;
}

const OutputHandler::Dimensions
AngularVelocityJoint::GetEquationDimension(integer index) const {
	// DOF == 1
	OutputHandler::Dimensions dimension = OutputHandler::Dimensions::UnknownDimension;

	switch (index)
	{
		case 1:
			dimension = OutputHandler::Dimensions::AngularVelocity;
			break;
		case 2:
			dimension = OutputHandler::Dimensions::AngularVelocity;
			break;
		case 3:
			dimension = OutputHandler::Dimensions::AngularVelocity;
			break;
	}

	return dimension;
}

std::ostream&
AngularVelocityJoint::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{

	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << "->" << iIndex + 3 << ": " <<
			"angular velocity error" << std::endl;

	return out;
}
/* AngularVelocity - end */
