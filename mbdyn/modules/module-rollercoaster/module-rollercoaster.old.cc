/* 
 * vincolo che obbliga a giacere su una linea.
 * usa un set di funzioni, che devono essere consistenti,
 * e che danno la posizione, l'orientazione e le loro derivate
 * in funzione dell'ascissa curvilinea
 */

#include <loadable.h>

struct module_rollercoaster {
   // user-defined struct
   StructNode* pNode;
   
   Vec3 F;
   doublereal s;
};

/* nota: la curva e' un cerchio nel piano x z,
 * e s e' la coordinata polare di azimuth */

const doublereal R = 1.;

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
void* read(LoadableElem* pEl,
	   DataManager* pDM,
	   MBDynParser& HP,
	   const DriveHandler* pDH)
{
   DEBUGCOUTFNAME("read");
   
   // allocation of user-defined struct
   module_rollercoaster* p = NULL;
   SAFENEW(p, module_rollercoaster, EMmm);
   
   /* nodo */
   unsigned int uNode = (unsigned int)HP.GetInt();
       
   DEBUGCOUT("Linked to Node " << uNode << endl);
       
   /* verifica di esistenza del nodo */
   if ((p->pNode = pDM->pFindStructNode(uNode)) == NULL) {
      cerr << "line " << HP.GetLineData() 
	<< ": structural node " << uNode
	<< " not defined" << endl;    
      THROW(DataManager::ErrGeneric());
   }
   
   return (void *)p;
}

unsigned int i_get_num_dof(const LoadableElem* pEl)
{
   DEBUGCOUTFNAME("i_get_num_dof");
   return 3;
}

DofOrder::Order set_dof(const LoadableElem*, unsigned int i)
{
   DEBUGCOUTFNAME("set_dof");
   return DofOrder::ALGEBRAIC;
}

void output(const LoadableElem* pEl, OutputHandler& OH)
{
   DEBUGCOUTFNAME("output");
   
   ostream& out = OH.Loadable();
   
   module_rollercoaster* p = (module_rollercoaster *)pEl->pGetData();
   
   out << setw(8) << pEl->GetLabel() << " "
     << p->F << " " << p->s << endl;
}

ostream& restart(const LoadableElem* pEl, ostream& out)
{
   DEBUGCOUTFNAME("restart");
   return out << "not implemented yet;" << endl;
}

void work_space_dim(const LoadableElem* pEl, 
		    integer* piNumRows, 
		    integer* piNumCols)
{
   DEBUGCOUTFNAME("work_space_dim");
   *piNumRows = 6;
   *piNumCols = 6;
}

VariableSubMatrixHandler& 
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
   WM.ResizeInit(iNumRows, iNumCols, 0.);
   
   module_rollercoaster* p = (module_rollercoaster *)pEl->pGetData();
   
   // set sub-matrix indices and coefs
   integer iMomIndex = p->pNode->iGetFirstMomentumIndex();
   integer iPosIndex = p->pNode->iGetFirstPositionIndex();
   integer iIndex = pEl->iGetFirstIndex();
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutRowIndex(iCnt, iMomIndex+iCnt);
      WM.fPutColIndex(iCnt, iPosIndex+iCnt);
      WM.fPutRowIndex(3+iCnt, iIndex+iCnt);
      WM.fPutColIndex(3+iCnt, iIndex+iCnt);
   }

   Vec3 X(get_x(p->s));
   Vec3 XP(get_xp(p->s));
   Mat3x3 R(get_r(p->s));
   Vec3 rho(get_rho(p->s));
   
   Vec3 e1(R.GetVec(1));
   Vec3 e2(R.GetVec(2));
   
   Vec3 Xn(p->pNode->GetXCurr());
   
   Vec3 Tmp(rho.Cross(R*p->F));
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WM.fPutCoef(3+iCnt, iCnt, -dCoef);
      WM.fPutCoef(3+iCnt, 6, XP.dGet(iCnt));
      
      WM.fPutCoef(iCnt, 4, -e1.dGet(iCnt));
      WM.fPutCoef(iCnt, 5, -e2.dGet(iCnt));
      
      WM.fPutCoef(iCnt, 6, -Tmp.dGet(iCnt));
   }
   
   return WorkMat;
}

void
ass_eig(LoadableElem* pEl, 
	VariableSubMatrixHandler& WorkMatA,
	VariableSubMatrixHandler& WorkMatB,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{  
   DEBUGCOUTFNAME("ass_eig");
   integer iNumRows = 0;
   integer iNumCols = 0;
   pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
   FullSubMatrixHandler& WMA = WorkMatA.SetFull();
   WMA.ResizeInit(iNumRows, iNumCols, 0.);
   
   FullSubMatrixHandler& WMB = WorkMatB.SetFull();
   WMB.ResizeInit(iNumRows, iNumCols, 0.);
   
   module_rollercoaster* p = (module_rollercoaster *)pEl->pGetData();
   
   // set sub-matrix indices and coefs
}

SubVectorHandler& 
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
   
   for (int iCnt = 1; iCnt <= 3; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iMomIndex+iCnt);
      WorkVec.fPutRowIndex(3+iCnt, iIndex+iCnt);
   }

   p->F = Vec3(XCurr.dGetCoef(iIndex+1),
	       XCurr.dGetCoef(iIndex+2),
	       0.);   
   p->s = XCurr.dGetCoef(iIndex+3);
   
   Vec3 X(get_x(p->s));
   Mat3x3 R(get_r(p->s));
   
   Vec3 Xn(p->pNode->GetXCurr());
   
   WorkVec.Put(1, R*p->F);
   WorkVec.Put(4, /* R.Transpose()* */ (Xn-X));
   
   return WorkVec;
}

void before_predict(const LoadableElem* pEl, 
		    VectorHandler& X,
		    VectorHandler& XP,
		    VectorHandler& XPrev,
		    VectorHandler& XPPrev)
{
   DEBUGCOUTFNAME("before_predict");
}

void after_predict(const LoadableElem* pEl, 
		   VectorHandler& X,
		   VectorHandler& XP)
{
   DEBUGCOUTFNAME("after_predict");
}

void update(LoadableElem* pEl, 
	    const VectorHandler& X,
	    const VectorHandler& XP)
{
   DEBUGCOUTFNAME("update");
}

unsigned int i_get_initial_num_dof(const LoadableElem* pEl)
{
   DEBUGCOUTFNAME("i_get_initial_num_dof");
   return 0;
}

void initial_work_space_dim(const LoadableElem* pEl, 
			    integer* piNumRows, 
			    integer* piNumCols)
{
   DEBUGCOUTFNAME("initial_work_space_dim");
   *piNumRows = 0;
   *piNumCols = 0;   
}

VariableSubMatrixHandler& 
initial_ass_jac(LoadableElem* pEl, 
		VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("initial_ass_jac");
   integer iNumRows = 0;
   integer iNumCols = 0;
   pEl->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeInit(iNumRows, iNumCols, 0.);
   
   module_rollercoaster* p = (module_rollercoaster *)pEl->pGetData();
   
   // set sub-matrix indices and coefs
   
   return WorkMat;
}

SubVectorHandler& 
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

void set_value(const LoadableElem* pEl, VectorHandler& X, VectorHandler& XP)
{
   DEBUGCOUTFNAME("set_value");
   
   integer iIndex = pEl->iGetFirstIndex();
   X.fPutCoef(iIndex+3, -M_PI/2.);
}
   
void set_initial_value(const LoadableElem* pEl, VectorHandler& X)
{
   DEBUGCOUTFNAME("set_initial_value");
}

unsigned int i_get_num_priv_data(const LoadableElem* pEl)
{
   DEBUGCOUTFNAME("i_get_num_priv_data");
   return 0;
}

doublereal d_get_priv_data(const LoadableElem* pEl, unsigned int i)
{
   DEBUGCOUTFNAME("d_get_priv_data");
   ASSERT(pEl->iGetNumPrivData() > 0);
   if (i > pEl->iGetNumPrivData()) {
      cerr << "Module-template Elem: illegal private data index " << i << endl;      
      THROW(ErrGeneric());
   }
   
   // return i-th priv data
   return 0.;
}

void destroy(LoadableElem* pEl)
{
   DEBUGCOUTFNAME("destroy");
   module_rollercoaster* p = (module_rollercoaster *)pEl->pGetData();
   SAFEDELETE(p, EMmm);
}
