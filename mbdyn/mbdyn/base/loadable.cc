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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#include <mbconfig.h>

#include <loadable.h>
#include <dlfcn.h>
#include <dataman.h>

/* funzioni di default */
unsigned int 
__i_get_num_dof(const LoadableElem* /* pEl */ )
{
   return 0;
}

DofOrder::Order 
__set_dof(const LoadableElem*, unsigned int /* i */ )
{
   cerr << "You shouldn't be here!" << endl;
   THROW(ErrGeneric());
#ifndef USE_EXCEPTIONS
   return DofOrder::UNKNOWN;
#endif /* USE_EXCEPTIONS */
}


void 
__output(const LoadableElem* /* pEl */ , OutputHandler& /* OH */ )
{
   NO_OP;
}


ostream& 
__restart(const LoadableElem* pEl , ostream& out)
{
   return out << "loadable: " << pEl->GetLabel() << ", not implemented yet;" << endl;
}


void 
__work_space_dim(const LoadableElem* /* pEl */ ,
		 integer* piNumRows, 
		 integer* piNumCols)
{
   *piNumRows = 0;
   *piNumCols = 0;
}


VariableSubMatrixHandler& 
__ass_jac(LoadableElem* /* pEl */ ,
	  VariableSubMatrixHandler& WorkMat,
	  doublereal /* dCoef */ ,
	  const VectorHandler& /* XCurr */ ,
	  const VectorHandler& /* XPrimeCurr */ )
{
   WorkMat.SetNullMatrix();
   return WorkMat;
}


void 
__ass_eig(LoadableElem* /* pEl */ ,
	  VariableSubMatrixHandler& WorkMatA,
	  VariableSubMatrixHandler& WorkMatB,
	  const VectorHandler& /* XCurr */ ,
	  const VectorHandler& /* XPrimeCurr */ )
{
   WorkMatA.SetNullMatrix();
   WorkMatB.SetNullMatrix(); 
}


SubVectorHandler& 
__ass_res(LoadableElem* /* pEl */ ,
	  SubVectorHandler& WorkVec,
	  doublereal /* dCoef */ ,
	  const VectorHandler& /* XCurr */ ,
	  const VectorHandler& /* XPrimeCurr */ )
{
   WorkVec.Resize(0);
   return WorkVec;
}


void 
__before_predict(const LoadableElem* /* pEl */ ,
		 VectorHandler& /* X */ ,
		 VectorHandler& /* XP */ ,
		 VectorHandler& /* XPrev */ ,
		 VectorHandler& /* XPPrev */ )
{
   NO_OP;
}


void 
__after_predict(const LoadableElem* /* pEl */ ,
		VectorHandler& /* X */ ,
		VectorHandler& /* XP */ )
{
   NO_OP;
}


void 
__update(LoadableElem* /* pEl */ ,
	 const VectorHandler& /* X */ ,
	 const VectorHandler& /* XP */ )
{
   NO_OP;
}


unsigned int 
__i_get_initial_num_dof(const LoadableElem* /* pEl */ )
{
   return 0;
}


void 
__initial_work_space_dim(const LoadableElem* /* pEl */ ,
			 integer* piNumRows, 
			 integer* piNumCols)
{
   *piNumRows = 0;
   *piNumCols = 0;   
}


VariableSubMatrixHandler& 
__initial_ass_jac(LoadableElem* /* pEl */ ,
		  VariableSubMatrixHandler& WorkMat, 
		  const VectorHandler& /* XCurr */ )
{
   WorkMat.SetNullMatrix();
   return WorkMat;
}


SubVectorHandler& 
__initial_ass_res(LoadableElem* /* pEl */ ,
		  SubVectorHandler& WorkVec, 
		  const VectorHandler& /* XCurr */ )
{  
   WorkVec.Resize(0);
   return WorkVec;
}


void 
__set_value(const LoadableElem* /* pEl */ , VectorHandler& /* X */ , VectorHandler& /* XP */ )
{
   NO_OP;
}


void 
__set_initial_value(const LoadableElem* /* pEl */ , VectorHandler& /* X */ )
{
   NO_OP;
}


unsigned int 
__i_get_num_priv_data(const LoadableElem* /* pEl */ )
{
   return 0;
}


doublereal 
__d_get_priv_data(const LoadableElem* /* pEl */ , unsigned int /* i */ )
{
   cerr << "You shouldn't be here!" << endl;
   THROW(ErrGeneric());
}


void 
__destroy(LoadableElem* /* pEl */ )
{
   NO_OP;
}


const char* const 
demangled_func_names[] = {
#if 0
   "read",   
   "i_get_num_dof",
   "set_dof",
   "output",
   "restart",
   "work_space_dim",
   "ass_jac",
   "ass_eig",
   "ass_res",
   "before_predict",
   "after_predict",
   "update",
   "i_get_initial_num_dof",
   "initial_work_space_dim",
   "initial_ass_jac",
   "initial_ass_res",
   "set_value",
   "set_initial_value",
   "i_get_num_priv_data",
   "d_get_priv_data",
   "destroy",
   NULL
#else
# include "demangled.h"
#endif
};


const char* const 
mangled_func_names[] = {
#if 0
   "read__FP12LoadableElemP11DataManagerR11MBDynParserPC12DriveHandler",
   "i_get_num_dof__FPC12LoadableElem",
   "set_dof__FPC12LoadableElemUi",
   "output__FPC12LoadableElemR13OutputHandler",
   "restart__FPC12LoadableElemR7ostream",
   "work_space_dim__FPC12LoadableElemPlT1",
   "ass_jac__FP12LoadableElemR24VariableSubMatrixHandlerdRC13VectorHandlerT3",
   "ass_eig__FP12LoadableElemR24VariableSubMatrixHandlerT1RC13VectorHandlerT3",
   "ass_res__FP12LoadableElemR16SubVectorHandlerdRC13VectorHandlerT3",
   "before_predict__FPC12LoadableElemR13VectorHandlerN31",
   "after_predict__FPC12LoadableElemR13VectorHandlerT1",
   "update__FP12LoadableElemRC13VectorHandlerT1",
   "i_get_initial_num_dof__FPC12LoadableElem",
   "initial_work_space_dim__FPC12LoadableElemPlT1",
   "initial_ass_jac__FP12LoadableElemR24VariableSubMatrixHandlerRC13VectorHandler",
   "initial_ass_res__FP12LoadableElemR16SubVectorHandlerRC13VectorHandler",
   "set_value__FPC12LoadableElemR13VectorHandlerT1",
   "set_initial_value__FPC12LoadableElemR13VectorHandler",
   "i_get_num_priv_data__FPC12LoadableElem",
   "d_get_priv_data__FPC12LoadableElemUi",
   "destroy__FP12LoadableElem",
   NULL
#else
# include "mangled.h"
#endif
};


const char* const *func_names = mangled_func_names;


void* 
default_funcs[] = {
   NULL,
   
   (void*)__i_get_num_dof,
   (void*)__set_dof,
   
   (void*)__output,
   (void*)__restart,
   
   (void*)__work_space_dim,
   (void*)__ass_jac,
   (void*)__ass_eig,
   (void*)__ass_res,
   
   (void*)__before_predict,
   (void*)__after_predict,
   (void*)__update,
   
   (void*)__i_get_initial_num_dof,
   (void*)__initial_work_space_dim,
   (void*)__initial_ass_jac,
   (void*)__initial_ass_res,
   
   (void*)__set_value,
   (void*)__set_initial_value,
   
   (void*)__i_get_num_priv_data,
   (void*)__d_get_priv_data,
   
   (void*)__destroy
};


typedef void * (* p_read)(LoadableElem*,
			  DataManager*,
			  MBDynParser&,
			  const DriveHandler*);

typedef unsigned int (* p_i_get_num_dof)(const LoadableElem*);

typedef DofOrder::Order (* p_set_dof)(const LoadableElem*, unsigned int);

typedef void (* p_output)(const LoadableElem*, OutputHandler&);

typedef ostream& (* p_restart)(const LoadableElem*, ostream&);

typedef void (* p_work_space_dim)(const LoadableElem*, integer*, integer*);

typedef VariableSubMatrixHandler& (* p_ass_jac)(LoadableElem*,
						VariableSubMatrixHandler&,
						doublereal,
						const VectorHandler&,
						const VectorHandler&);

typedef void (* p_ass_eig)(LoadableElem*,
			   VariableSubMatrixHandler&,
			   VariableSubMatrixHandler&,
			   const VectorHandler&,
			   const VectorHandler&);

typedef SubVectorHandler& (* p_ass_res)(LoadableElem*,
					SubVectorHandler&,
					doublereal,
					const VectorHandler&,
					const VectorHandler&);

typedef void (* p_before_predict)(const LoadableElem* pEl, 
				  VectorHandler& X,
				  VectorHandler& XP,
				  VectorHandler& XPrev,
				  VectorHandler& XPPrev);

typedef void (* p_after_predict)(const LoadableElem* pEl, 
				 VectorHandler& X,
				 VectorHandler& XP);

typedef void (* p_update)(LoadableElem* pEl, 
			  const VectorHandler& X,
			  const VectorHandler& XP);

typedef unsigned int (* p_i_get_initial_num_dof)(const LoadableElem*);

typedef void (* p_initial_work_space_dim)(const LoadableElem*, integer*, integer*);

typedef VariableSubMatrixHandler& (* p_initial_ass_jac)(LoadableElem*,
							VariableSubMatrixHandler&, 
							const VectorHandler&);

typedef SubVectorHandler& (* p_initial_ass_res)(LoadableElem*,
						SubVectorHandler&,
						const VectorHandler&);

typedef void (* p_set_value)(const LoadableElem*, VectorHandler&, VectorHandler&);

typedef void (* p_set_initial_value)(const LoadableElem*, VectorHandler&);

typedef unsigned int (* p_i_get_num_priv_data)(const LoadableElem* pEl);

typedef doublereal (* p_d_get_priv_data)(const LoadableElem* pEl, 
					 unsigned int i);

typedef void (* p_destroy)(LoadableElem*);

		    

/* metodi effettivi */

LoadableElem::LoadableElem(unsigned int uLabel, 
			   const DofOwner* pDO, 
			   DataManager* pDM, 
			   MBDynParser& HP)
: Elem(uLabel, ElemType::LOADABLE, flag(0)),
#ifdef USE_STRUCT_NODES
InitialAssemblyElem(uLabel, ElemType::LOADABLE, flag(0)),
ElemGravityOwner(uLabel, ElemType::LOADABLE, flag(0)),
#endif /* USE_STRUCT_NODES */
ElemWithDofs(uLabel, ElemType::LOADABLE, pDO, flag(0)),
priv_data(NULL),
module_name(NULL),
handle(NULL)

{
   const char sFuncName[] = "LoadableElem::LoadableElem";
   
   ASSERT(pDM != NULL);
   
   /* nome del modulo */
   const char* s = HP.GetFileName();
   SAFENEWARR(module_name, char, strlen(s)+1, DMmm);
   strcpy(module_name, s);
   
   handle = dlopen(module_name, RTLD_NOW /* RTLD_LAZY */ );
   if (handle == NULL) {    
      cerr << sFuncName << "(" << uLabel 
	<< "): unable to open module <" << module_name 
	<< "> (dlopen returns \"" << dlerror() << "\")" << endl;
      THROW(ErrGeneric());
   }

   DEBUGCOUT("binding to function \"" << func_names[READ] 
	     << "\" (must be def'd!)" << endl);
   fsym[READ] = dlsym(handle, func_names[READ]);
   if (fsym[READ] == NULL) {
      const char* err = dlerror();
      if (err == NULL) {
	 cerr << sFuncName << "(" << uLabel 
	   << "): function \"" << func_names[READ]
	   << "\" must be defined in module <" << module_name << ">" << endl;
      } else {
	 cerr << sFuncName << "(" << uLabel 
	   << "): error in binding to function \"" << func_names[READ]
	   << "\" in module <" << module_name
	   << "> (dlsym returns \"" << err << "\")" << endl;
      }
      THROW(ErrGeneric());
   }
   
   for (int i = 1; i < LASTFUNC; i++) {
      DEBUGCOUT("binding to function \"" << func_names[i]
		<< "\" (default function available)" << endl);
      fsym[i] = dlsym(handle, func_names[i]);
      if (fsym[i] == NULL) {
	 const char* err = dlerror();
	 if (err != NULL && !fSilent) {
	    cerr << sFuncName << "(" << uLabel 
	      << "): bind error (dlsym returns \"" << err 
	      << "\" for \"" << func_names[i] 
	      << "\"); using default function" << endl;
	 }
	 fsym[i] = default_funcs[i];	 
      }
   }
   
   priv_data = (*((p_read)fsym[READ]))(this, pDM, HP, pDM->pGetDrvHdl());
   
   /* Mettere LOADABLE quando sara' definito */
   SetOutputFlag(pDM->fReadOutput(HP, ElemType::LOADABLE)); 
}


LoadableElem::~LoadableElem(void)
{
   ASSERT(fsym[DESTROY] != NULL);
   (*((p_destroy)fsym[DESTROY]))(this);
   
   ASSERT(handle != NULL);
   dlclose(handle);
   SAFEDELETEARR(module_name, DMmm);
}


ElemType::Type 
LoadableElem::GetElemType(void) const
{
   return ElemType::LOADABLE;
}


unsigned int 
LoadableElem::iGetNumDof(void) const
{
   ASSERT(fsym[IGETNUMDOF] != NULL);
   return (*((p_i_get_num_dof)fsym[IGETNUMDOF]))(this);
}


DofOrder::Order 
LoadableElem::SetDof(unsigned int i) const
{
   ASSERT(fsym[SETDOF] != NULL);
   ASSERT(i < iGetNumDof());
   return (*((p_set_dof)fsym[SETDOF]))(this, i);
}


void 
LoadableElem::Output(OutputHandler& OH) const
{
   ASSERT(fsym[OUTPUT] != NULL);
   (*((p_output)fsym[OUTPUT]))(this, OH);
}


ostream& 
LoadableElem::Restart(ostream& out) const
{
   out << "    loadable: " << GetLabel() << ", \"" << module_name << "\", ";
   ASSERT(fsym[RESTART] != NULL);
   return (*((p_restart)fsym[RESTART]))(this, out) << ';' << endl;
}


void 
LoadableElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
   ASSERT(fsym[WORKSPACEDIM] != NULL);
   (*((p_work_space_dim)fsym[WORKSPACEDIM]))(this, piNumRows, piNumCols);      
}


VariableSubMatrixHandler& 
LoadableElem::AssJac(VariableSubMatrixHandler& WorkMat,
		     doublereal dCoef, 
		     const VectorHandler& XCurr,
		     const VectorHandler& XPCurr)
{
   ASSERT(fsym[ASSJAC] != NULL);
   return (*((p_ass_jac)fsym[ASSJAC]))(this, WorkMat, dCoef, XCurr, XPCurr);
}


void
LoadableElem::AssEig(VariableSubMatrixHandler& WorkMatA,
		     VariableSubMatrixHandler& WorkMatB,
		     const VectorHandler& XCurr,
		     const VectorHandler& XPCurr)
{
   ASSERT(fsym[ASSEIG] != NULL);
   (*((p_ass_eig)fsym[ASSEIG]))(this, WorkMatA, WorkMatB, XCurr, XPCurr);
}


SubVectorHandler& 
LoadableElem::AssRes(SubVectorHandler& WorkVec,
		     doublereal dCoef,
		     const VectorHandler& XCurr, 
		     const VectorHandler& XPCurr)
{
   ASSERT(fsym[ASSRES] != NULL);
   return (*((p_ass_res)fsym[ASSRES]))(this, WorkVec, dCoef, XCurr, XPCurr);
}
   
void 
LoadableElem::BeforePredict(VectorHandler& X,
			    VectorHandler& XP,
			    VectorHandler& XPrev,
			    VectorHandler& XPPrev) const
{
   ASSERT(fsym[BEFOREPREDICT] != NULL);
   (*((p_before_predict)fsym[BEFOREPREDICT]))(this, X, XP, XPrev, XPPrev);
}


void 
LoadableElem::AfterPredict(VectorHandler& X,
			   VectorHandler& XP)
{
   ASSERT(fsym[AFTERPREDICT] != NULL);
   (*((p_after_predict)fsym[AFTERPREDICT]))(this, X, XP);
}


void 
LoadableElem::Update(const VectorHandler& XCurr, 
		     const VectorHandler& XPrimeCurr)
{
   ASSERT(fsym[UPDATE] != NULL);
   (*((p_update)fsym[UPDATE]))(this, XCurr, XPrimeCurr);
}


#ifdef USE_STRUCT_NODES
unsigned int 
LoadableElem::iGetInitialNumDof(void) const
{
   ASSERT(fsym[IGETINITIALNUMDOF] != NULL);
   return (*((p_i_get_initial_num_dof)fsym[IGETINITIALNUMDOF]))(this);
}


void 
LoadableElem::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
   ASSERT(fsym[INITIALWORKSPACEDIM] != NULL);
   (*((p_initial_work_space_dim)fsym[INITIALWORKSPACEDIM]))(this,  
							    piNumRows, 
							    piNumCols);
}


VariableSubMatrixHandler& 
LoadableElem::InitialAssJac(VariableSubMatrixHandler& WorkMat,
			    const VectorHandler& XCurr)
{
   ASSERT(fsym[INITIALASSJAC] != NULL);
   return (*((p_initial_ass_jac)fsym[INITIALASSJAC]))(this, WorkMat, XCurr);
}


SubVectorHandler& 
LoadableElem::InitialAssRes(SubVectorHandler& WorkVec, 
			    const VectorHandler& XCurr)
{
   ASSERT(fsym[INITIALASSRES] != NULL);
   return (*((p_initial_ass_res)fsym[INITIALASSRES]))(this, WorkVec, XCurr);
}


void 
LoadableElem::SetInitialValue(VectorHandler& X) const
{   
   ASSERT(fsym[SETINITIALVALUE] != NULL);
   (*((p_set_initial_value)fsym[SETINITIALVALUE]))(this, X);
}
#endif /* USE_STRUCT_NODES */


void 
LoadableElem::SetValue(VectorHandler& X, VectorHandler& XP) const
{
   ASSERT(fsym[SETVALUE] != NULL);
   (*((p_set_value)fsym[SETVALUE]))(this, X, XP);
}


unsigned int 
LoadableElem::iGetNumPrivData(void) const
{
   ASSERT(fsym[IGETNUMPRIVDATA] != NULL);
   return (*((p_i_get_num_priv_data)fsym[IGETNUMPRIVDATA]))(this);
}


doublereal 
LoadableElem::dGetPrivData(unsigned int i) const
{
   ASSERT(fsym[DGETPRIVDATA] != NULL);
   return (*((p_d_get_priv_data)fsym[DGETPRIVDATA]))(this, i);
}


Elem* 
ReadLoadable(DataManager* pDM, 
	     MBDynParser& HP, 
	     const DofOwner* pDO, 
	     unsigned int uLabel)
{
   Elem* pEl = NULL;
   
   SAFENEWWITHCONSTRUCTOR(pEl,
			  LoadableElem,
			  LoadableElem(uLabel, pDO, pDM, HP),
			  DMmm);
   return pEl;
}
