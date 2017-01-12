/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

/* 
 * vincolo che obbliga a giacere su una linea.
 * usa un set di funzioni, che devono essere consistenti,
 * e che danno la posizione, l'orientazione e le loro derivate
 * in funzione dell'ascissa curvilinea
 */

#include "loadable.h"

struct module_rollercoaster {
   // user-defined struct
   StructNode* pNode;
   
   Vec3 F;
   Vec3 M;
   doublereal s;
};

/* nota: la curva e' un cerchio nel piano x z,
 * e s e' la coordinata polare di azimuth */

static const doublereal R = 1.;

static Vec3
get_x(doublereal s)
{
   doublereal cs = cos(s);
   doublereal ss = sin(s);
   return Vec3(R*ss, 0., -R*cs);
}

static Vec3
get_xp(doublereal s)
{
   doublereal cs = cos(s);
   doublereal ss = sin(s);
   return Vec3(R*cs, 0., R*ss);
}

static Mat3x3 
get_r(doublereal s)
{
   doublereal cs = cos(s);
   doublereal ss = sin(s);
   return Mat3x3(ss, 0., -cs, 0., 1., 0., cs, 0., ss);
}

static Vec3
get_rho(doublereal s)
{
   return Vec3(0., 1., 0.);
}

/* funzioni di default */
static void*
read(LoadableElem* pEl,
	   DataManager* pDM,
	   MBDynParser& HP,
	   const DriveHandler* pDH)
{
   DEBUGCOUTFNAME("read");
   
   // allocation of user-defined struct
   module_rollercoaster* p = NULL;
   SAFENEW(p, module_rollercoaster);
   
   /* nodo */
   unsigned int uNode = (unsigned int)HP.GetInt();
       
   DEBUGCOUT("Linked to Node " << uNode << endl);
       
   /* verifica di esistenza del nodo */
   if ((p->pNode = pDM->pFindStructNode(uNode)) == NULL) {
      cerr << "line " << HP.GetLineData() 
	<< ": structural node " << uNode
	<< " not defined" << endl;    
      throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
   
   return (void *)p;
}

static unsigned int
i_get_num_dof(const LoadableElem* pEl)
{
   DEBUGCOUTFNAME("i_get_num_dof");
   return 6;
}

static DofOrder::Order
set_dof(const LoadableElem*, unsigned int i)
{
   DEBUGCOUTFNAME("set_dof");
   return DofOrder::ALGEBRAIC;
}

static void
output(const LoadableElem* pEl, OutputHandler& OH)
{
   DEBUGCOUTFNAME("output");
   
   ostream& out = OH.Loadable();
   
   module_rollercoaster* p = (module_rollercoaster *)pEl->pGetData();
   
   out << setw(8) << pEl->GetLabel() << " "
     << p->F << " "
     << p->M << " "
     << p->s << endl;
}

static ostream&
restart(const LoadableElem* pEl, ostream& out)
{
   DEBUGCOUTFNAME("restart");
   return out << "not implemented yet;" << endl;
}

static void
work_space_dim(const LoadableElem* pEl, 
		    integer* piNumRows, 
		    integer* piNumCols)
{
   DEBUGCOUTFNAME("work_space_dim");
   *piNumRows = 12;
   *piNumCols = 12;
}

static VariableSubMatrixHandler& 
ass_jac(LoadableElem* pEl, 
	VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{  
   DEBUGCOUTFNAME("ass_jac");
   integer iNumRows = 0;
   integer iNumCols = 0;
   pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeReset(iNumRows, iNumCols, 0.);
   
   module_rollercoaster* p = (module_rollercoaster *)pEl->pGetData();
   
   // set sub-matrix indices and coefs
   integer iMomIndex = p->pNode->iGetFirstMomentumIndex();
   integer iPosIndex = p->pNode->iGetFirstPositionIndex();
   integer iIndex = pEl->iGetFirstIndex();
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WM.PutRowIndex(iCnt, iMomIndex+iCnt);
      WM.PutColIndex(iCnt, iPosIndex+iCnt);
      WM.PutRowIndex(6+iCnt, iIndex+iCnt);
      WM.PutColIndex(6+iCnt, iIndex+iCnt);
   }

   Vec3 X(get_x(p->s));
   Vec3 XP(get_xp(p->s));
   Mat3x3 R(get_r(p->s));
   Vec3 rho(get_rho(p->s));
   
   Vec3 Xn(p->pNode->GetXCurr());
   Mat3x3 Rn(p->pNode->GetRCurr());
   
   Vec3 e1(R.GetVec(1));
   Vec3 e2(R.GetVec(2));
   Vec3 e3(R.GetVec(3));

   Vec3 e1n(Rn.GetVec(1));
   Vec3 e2n(Rn.GetVec(2));
   Vec3 e3n(Rn.GetVec(3));
   
   Vec3 E1(e2n.Cross(e3));
   Vec3 E2(e3n.Cross(e1));
   Vec3 E3(e1n.Cross(e2));

   Mat3x3 E(E1, E2, E3);
   Mat3x3 ET(E.Transpose());
   
   Vec3 Tmp1(rho.Cross(R*p->F));
   Vec3 Tmp2(e2n.Cross(e3.Cross(rho*p->M.dGet(1)))
	     +e3n.Cross(e1.Cross(rho*p->M.dGet(2)))
	     +e1n.Cross(e2.Cross(rho*p->M.dGet(3))));
   Vec3 Tmp3(ET*rho);
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.PutCoef(6+iCnt, iCnt, -dCoef);
      WM.PutCoef(6+iCnt, 9, XP.dGet(iCnt));
      
      WM.PutCoef(iCnt, 7, -e1.dGet(iCnt));
      WM.PutCoef(iCnt, 8, -e2.dGet(iCnt));
      
      WM.PutCoef(iCnt, 9, -Tmp1.dGet(iCnt));
      WM.PutCoef(3+iCnt, 9, -Tmp2.dGet(iCnt));
            
      WM.PutCoef(9+iCnt, 9, Tmp3.dGet(iCnt));
   }
   
   WM.Sub(4, 4, 
	  Mat3x3(e3*(p->M.dGet(1)*dCoef), e2n)
	  +Mat3x3(e1*(p->M.dGet(2)*dCoef), e3n)
	  +Mat3x3(e2*(p->M.dGet(3)*dCoef), e1n));
   
   WM.Sub(4, 10, E);
   WM.Sub(10, 4, ET*dCoef);
   
   return WorkMat;
}

static void
ass_mats(LoadableElem* pEl, 
	VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{  
   DEBUGCOUTFNAME("ass_mats");
   integer iNumRows = 0;
   integer iNumCols = 0;
   pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
   FullSubMatrixHandler& WMA = WorkMatA.SetFull();
   WMA.ResizeReset(iNumRows, iNumCols, 0.);
   
   FullSubMatrixHandler& WMB = WorkMatB.SetFull();
   WMB.ResizeReset(iNumRows, iNumCols, 0.);
   
   module_rollercoaster* p = (module_rollercoaster *)pEl->pGetData();
   
   // set sub-matrix indices and coefs
}

static SubVectorHandler& 
ass_res(LoadableElem* pEl, 
	SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
   DEBUGCOUTFNAME("ass_res");
   integer iNumRows = 0;
   integer iNumCols = 0;
   pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
   WorkVec.Resize(iNumRows);
   
   module_rollercoaster* p = (module_rollercoaster *)pEl->pGetData();

   // set sub-vector indices and coefs
   integer iMomIndex = p->pNode->iGetFirstMomentumIndex();
   integer iIndex = pEl->iGetFirstIndex();
   
   for (int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.PutRowIndex(iCnt, iMomIndex+iCnt);
      WorkVec.PutRowIndex(6+iCnt, iIndex+iCnt);
   }

   p->F = Vec3(XCurr(iIndex+1),
	       XCurr(iIndex+2),
	       0.);
   p->M = Vec3(XCurr, iIndex+4);
   p->s = XCurr(iIndex+3);
   
   Vec3 X(get_x(p->s));
   Mat3x3 R(get_r(p->s));
   
   Vec3 Xn(p->pNode->GetXCurr());
   Mat3x3 Rn(p->pNode->GetRCurr());
   
   Vec3 e1(R.GetVec(1));
   Vec3 e2(R.GetVec(2));
   Vec3 e3(R.GetVec(3));

   Vec3 e1n(Rn.GetVec(1));
   Vec3 e2n(Rn.GetVec(2));
   Vec3 e3n(Rn.GetVec(3));
   
   WorkVec.Put(1, R*p->F);
   WorkVec.Put(4, Mat3x3(e2n.Cross(e3), e3n.Cross(e1), e1n.Cross(e2))*p->M);
   
   WorkVec.Put(7, /* R.Transpose()*( */ Xn-X /* ) */ );
   WorkVec.Put(10, Vec3(e2n.Dot(e3), e3n.Dot(e1), e1n.Dot(e2)));
   
   return WorkVec;
}

static void
before_predict(const LoadableElem* pEl, 
		    VectorHandler& X,
		    VectorHandler& XP,
		    VectorHandler& XPrev,
		    VectorHandler& XPPrev)
{
   DEBUGCOUTFNAME("before_predict");
}

static void
after_predict(const LoadableElem* pEl, 
		   VectorHandler& X,
		   VectorHandler& XP)
{
   DEBUGCOUTFNAME("after_predict");
}

static void
update(LoadableElem* pEl, 
	    const VectorHandler& X,
	    const VectorHandler& XP)
{
   DEBUGCOUTFNAME("update");
}

static unsigned int
i_get_initial_num_dof(const LoadableElem* pEl)
{
   DEBUGCOUTFNAME("i_get_initial_num_dof");
   return 0;
}

static void
initial_work_space_dim(const LoadableElem* pEl, 
			    integer* piNumRows, 
			    integer* piNumCols)
{
   DEBUGCOUTFNAME("initial_work_space_dim");
   *piNumRows = 0;
   *piNumCols = 0;   
}

static VariableSubMatrixHandler& 
initial_ass_jac(LoadableElem* pEl, 
		VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("initial_ass_jac");
   integer iNumRows = 0;
   integer iNumCols = 0;
   pEl->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeReset(iNumRows, iNumCols, 0.);
   
   module_rollercoaster* p = (module_rollercoaster *)pEl->pGetData();
   
   // set sub-matrix indices and coefs
   
   return WorkMat;
}

static SubVectorHandler& 
initial_ass_res(LoadableElem* pEl, 
		SubVectorHandler& WorkVec, 
		const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("initial_ass_res");
   integer iNumRows = 0;
   integer iNumCols = 0;
   pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
   WorkVec.Resize(iNumRows);
   
   module_rollercoaster* p = (module_rollercoaster *)pEl->pGetData(); 

   // set sub-vector indices and coefs
   
   return WorkVec;
}

static void
set_value(const LoadableElem* pEl, DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph)
{
   DEBUGCOUTFNAME("set_value");
   
   integer iIndex = pEl->iGetFirstIndex();
   X.PutCoef(iIndex+3, -M_PI/2.);
}
   
static void
set_initial_value(const LoadableElem* pEl, VectorHandler& X)
{
   DEBUGCOUTFNAME("set_initial_value");
}

static unsigned int
i_get_num_priv_data(const LoadableElem* pEl)
{
   DEBUGCOUTFNAME("i_get_num_priv_data");
   return 0;
}

static doublereal
d_get_priv_data(const LoadableElem* pEl, unsigned int i)
{
   DEBUGCOUTFNAME("d_get_priv_data");
   ASSERT(pEl->iGetNumPrivData() > 0);
   if (i > pEl->iGetNumPrivData()) {
      cerr << "Module-template Elem: illegal private data index " << i << endl;      
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
   
   // return i-th priv data
   return 0.;
}

static void
destroy(LoadableElem* pEl)
{
   DEBUGCOUTFNAME("destroy");
   module_rollercoaster* p = (module_rollercoaster *)pEl->pGetData();
   SAFEDELETE(p);
}

static struct
LoadableCalls lc = {
	"rollercoaster",
	read,
	i_get_num_dof,
	set_dof,
	output,
	restart,
	work_space_dim,
	ass_jac,
	ass_mats,
	ass_res,
	before_predict,
	after_predict,
	update,
	NULL, /* after_convergence */
	i_get_initial_num_dof,
	initial_work_space_dim,
	initial_ass_jac,
	initial_ass_res,
	set_value,
	set_initial_value,
	i_get_num_priv_data,
	d_get_priv_data,
	NULL /* i_get_num_connected_nodes */ ,
	NULL /* get_connected_nodes */ ,
	destroy
};

extern "C" {
void *calls = &lc;
}

