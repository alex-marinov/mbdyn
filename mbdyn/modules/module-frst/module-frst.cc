#include <loadable.h>
#include <dataman.h>

enum {
   SLIPPING,
   STITCHING,
};

struct module_template {
   // user-defined struct
   StructNode* pNode1;
   StructNode* pNode2;
   
   Vec3 v;
   Vec3 DXref;
   doublereal mu_s;
   doublereal mu_d;
   doublereal limit;
   int status;
   
   Vec3 f;
};

/* funzioni di default */
void* read(LoadableElem* pEl,
	   DataManager* pDM,
	   MBDynParser& HP,
	   const DriveHandler* pDH)
{
   DEBUGCOUTFNAME("read");
   
   // allocation of user-defined struct
   module_template* p = NULL;
   SAFENEW(p, module_template, EMmm);

   // read data
   unsigned int label = HP.GetInt();
   p->pNode1 = (StructNode*)pDM->pFindNode(NodeType::STRUCTURAL, label);
   
   label = HP.GetInt();
   p->pNode2 = (StructNode*)pDM->pFindNode(NodeType::STRUCTURAL, label);
   
   p->mu_s = HP.GetReal();
   p->mu_d = HP.GetReal();
   p->limit = HP.GetReal();
   
   if ((p->pNode2->GetVCurr()-p->pNode1->GetVCurr()).Norm() <= p->limit) {
      p->status = STITCHING;
   } else {
      p->status = SLIPPING;
   }

   p->f = Vec3(0.);
   
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
   *piNumRows = 15;
   *piNumCols = 15;
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
   
   /*
   FullSubMatrixHandler& WM = WorkMat.SetFull();
   WM.ResizeInit(iNumRows, iNumCols, 0.);
    */
   WorkMat.SetNullMatrix();
   
   module_template* p = (module_template *)pEl->pGetData();
   
   // set sub-matrix indices and coefs
   
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
   WorkMatA.SetNullMatrix();
   WorkMatB.SetNullMatrix();
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
   
   module_template* p = (module_template *)pEl->pGetData(); 

   // set sub-vector indices and coefs
   integer iNode1FirstMomentumIndex = p->pNode1->iGetFirstMomentumIndex();
   integer iNode2FirstMomentumIndex = p->pNode2->iGetFirstMomentumIndex();
   integer iFirstReactionIndex = pEl->iGetFirstIndex();
   
   for (int i = 6; i > 0; i--) {
      WorkVec.fPutRowIndex(i, iNode1FirstMomentumIndex+i);
      WorkVec.fPutRowIndex(6+i, iNode2FirstMomentumIndex+i);
   }
   for (int i = 3; i > 0; i--) {
      WorkVec.fPutRowIndex(12+i, iFirstReactionIndex+i);
   }
   
   Vec3 vtmp(p->pNode1->GetRCurr()*p->v);
   Vec3 x(p->pNode2->GetXCurr()-p->pNode1->GetXCurr());
   Vec3 w(p->pNode2->GetVCurr()-p->pNode1->GetVCurr());
   doublereal wt = w.Norm();
     
   p->f = Vec3(XCurr, iFirstReactionIndex+1);
   doublereal fn = p->f.dGet(3);
   doublereal d = p->f.dGet(1);
   doublereal ft = d*d;
   d = p->f.dGet(2);
   ft += d*d;
   ft = sqrt(ft);

   if (p->status == STITCHING) {
      doublereal d = p->f.dGet(1);
      doublereal ft = d*d;
      d = p->f.dGet(2);
      ft += d*d;
      ft = sqrt(ft);
      if (fabs(ft) > p->mu_s*fabs(fn)) {
	 p->status = SLIPPING;
      }
   } else if (p->status == SLIPPING && wt <= p->limit) { 
      p->status = STITCHING;
      p->DXref = x;
   }
   
   if (p->status == STITCHING) {
      Vec3 ftmp = p->pNode1->GetRCurr()*p->f;
      WorkVec.Sub(1, ftmp);
      WorkVec.Sub(4, x.Cross(ftmp));
      WorkVec.Add(7, ftmp);
      WorkVec.Add(13, p->DXref-x);
   } else {
      Vec3 ftmp = vtmp*fn-w*(p->mu_d*fn/wt);
      WorkVec.Sub(1, ftmp);
      WorkVec.Sub(4, x.Cross(ftmp));
      WorkVec.Add(7, ftmp);
      WorkVec.fDecCoef(13, p->f.dGet(1));
      WorkVec.fDecCoef(14, p->f.dGet(2));
      WorkVec.fDecCoef(15, p->v.Dot(x));
   }
   
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
   WorkMat.SetNullMatrix();   
   return WorkMat;
}

SubVectorHandler& 
initial_ass_res(LoadableElem* pEl, 
		SubVectorHandler& WorkVec, 
		const VectorHandler& XCurr)
{
   DEBUGCOUTFNAME("initial_ass_res");
   WorkVec.Resize(0);
   return WorkVec;
}

void set_value(const LoadableElem* pEl, VectorHandler& X, VectorHandler& XP)
{
   DEBUGCOUTFNAME("set_value");
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
   module_template* p = (module_template *)pEl->pGetData();
   SAFEDELETE(p, EMmm);
}
