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

#ifndef LOADABLE_H
#define LOADABLE_H

#include <elem.h>
#include <dataman.h>
#ifdef USE_STRUCT_NODES
#include <gravity.h>
#ifdef USE_AERODYNAMIC_ELEMS
#include <aerodyn.h>
#endif /* USE_AERODYNAMIC_ELEMS */
#endif /* USE_STRUCT_NODES */
/* #include <aerodyn.h> */

/* Elemento caricabile dinamicamente;
 * fornisce l'interfaccia ad un modulo a parte e compilato come shared library
 * che deve provvedere le funzioni desiderate dell'elemento, altrimenti
 * vengono usate funzioni di default.
 * Deve provvedere altresi' una routine di input che legge da MBDynParser i
 * dati richiesti dall'interno del costruttore e li mette in una struttura 
 * di dati privati che l'elemento vede come un (void *) .
 * Le funzioni da fornire seguono uno schema di denominazione standard.
 * Al momento: un modulo contiene un solo elemento, le funzioni
 * hanno il nome in "stile C", ovvero AssRes diventa ass_res.
 * Gli argomenti sono gli stessi della funzione della classe Elem,
 * ma sono preceduti da un puntatore all'elemento stesso:
 * 
 *     LoadableElem::AssRes(SubVectorHandler&, 
 *              	    doublereal, 
 *      		    const VectorHandler&,
 *      		    const VectorHandler&)
 * 
 * diventa
 * 
 *     ass_res(LoadableElem*,
 *             SubVectorHandler&, 
 *             doublereal, 
 *             const VectorHandler&,
 *             const VectorHandler&)
 * 
 * La funzione di ingresso dati e':
 * 
 *     void *read(LoadableElem*,
 *                DataManager*,
 *                MBDynParser&,
 *                DriveHandler*)
 * 
 * in questo modo, tramite il puntatore all'elemento, si ha accesso 
 * a tutti i dati dell'elemento;
 * rimane da studiare la questione dei diritti di accesso.
 */

class LoadableElem;

typedef void *
(* p_read)(LoadableElem*,
	   DataManager*,
	   MBDynParser&,
	   const DriveHandler*);
typedef unsigned int (* p_i_get_num_dof)(const LoadableElem*);
typedef DofOrder::Order (* p_set_dof)(const LoadableElem*, unsigned int);
typedef void (* p_output)(const LoadableElem*, OutputHandler&);
typedef ostream& (* p_restart)(const LoadableElem*, ostream&);
typedef void (* p_work_space_dim)(const LoadableElem*, integer*, integer*);
typedef VariableSubMatrixHandler& 
(* p_ass_jac)(LoadableElem*,
	      VariableSubMatrixHandler&,
	      doublereal,
	      const VectorHandler&,
	      const VectorHandler&);
typedef void
(* p_ass_eig)(LoadableElem*,
	      VariableSubMatrixHandler&,
	      VariableSubMatrixHandler&,
	      const VectorHandler&,
	      const VectorHandler&);
typedef SubVectorHandler&
(* p_ass_res)(LoadableElem*,
	      SubVectorHandler&,
	      doublereal,
	      const VectorHandler&,
	      const VectorHandler&);
typedef void
(* p_before_predict)(const LoadableElem* pEl, 
		     VectorHandler& X,
		     VectorHandler& XP,
		     VectorHandler& XPrev,
		     VectorHandler& XPPrev);
typedef void
(* p_after_predict)(const LoadableElem* pEl, 
		    VectorHandler& X,
		    VectorHandler& XP);
typedef void
(* p_update)(LoadableElem* pEl, 
	     const VectorHandler& X,
	     const VectorHandler& XP);
typedef unsigned int (* p_i_get_initial_num_dof)(const LoadableElem*);
typedef void
(* p_initial_work_space_dim)(const LoadableElem*, 
			     integer*, 
			     integer*);
typedef VariableSubMatrixHandler&
(* p_initial_ass_jac)(LoadableElem*,
		      VariableSubMatrixHandler &, 
		      const VectorHandler&);
typedef SubVectorHandler&
(* p_initial_ass_res)(LoadableElem*,
		      SubVectorHandler&,
		      const VectorHandler&);
typedef void
(* p_set_value)(const LoadableElem*,
		VectorHandler&,
		VectorHandler&);
typedef void (* p_set_initial_value)(const LoadableElem*, VectorHandler&);
typedef unsigned int (* p_i_get_num_priv_data)(const LoadableElem* pEl);
typedef doublereal
(* p_d_get_priv_data)(const LoadableElem* pEl, 
		      unsigned int i);
typedef void (* p_destroy)(LoadableElem*);

/*
 * Struttura che contiene le chiamate alle funzioni del modulo;
 * e' l'unico contatto tra modulo ed elemento caricato runtime.
 * Se un puntatore e' vuoto, viene riempito con il metodo di default.
 */
struct LoadableCalls {
	p_read				read;
	p_i_get_num_dof			i_get_num_dof;
	p_set_dof 			set_dof;
	p_output 			output;
	p_restart 			restart;
	p_work_space_dim 		work_space_dim;
	p_ass_jac			ass_jac;
	p_ass_eig			ass_eig;
	p_ass_res			ass_res;
	p_before_predict		before_predict;
	p_after_predict			after_predict;
	p_update			update;
	p_i_get_initial_num_dof		i_get_initial_num_dof;
	p_initial_work_space_dim	initial_work_space_dim;
	p_initial_ass_jac		initial_ass_jac;
	p_initial_ass_res		initial_ass_res;
	p_set_value			set_value;
	p_set_initial_value		set_initial_value;
	p_i_get_num_priv_data		i_get_num_priv_data;
	p_d_get_priv_data		d_get_priv_data;
	p_destroy			destroy;
};
 
class LoadableElem
: virtual public Elem,
#ifdef USE_STRUCT_NODES
public InitialAssemblyElem,
#ifdef USE_AERODYNAMIC_ELEMS
public AerodynamicElem,
#endif /* USE_AERODYNAMIC_ELEMS */
public ElemGravityOwner,
#endif /* USE_STRUCT_NODES */
public ElemWithDofs {
protected:
   	void* priv_data;	/* Dati privati passati alle funzioni */
   	char* module_name;	/* Nome del modulo */
   	void* handle;		/* Handle del modulo (usato per chiusura) */
	LoadableCalls *calls;	/* Simboli delle funzioni attese */
   
public:
   	LoadableElem(unsigned int uLabel, const DofOwner* pDO,
		     DataManager* pDM, MBDynParser& HP);
   	~LoadableElem(void); 
   	virtual inline void* pGet(void) const;
   
   	inline void* pGetData(void) const;
   
   	virtual Elem::Type GetElemType(void) const;

   	virtual unsigned int iGetNumDof(void) const;   
   	virtual DofOrder::Order SetDof(unsigned int i) const;
   
   	virtual void Output(OutputHandler& OH) const;
   	virtual ostream& Restart(ostream& out) const;
   
   	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	virtual VariableSubMatrixHandler& 
     	AssJac(VariableSubMatrixHandler& WorkMat,
	       doublereal dCoef,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);
	virtual void
     	AssEig(VariableSubMatrixHandler& WorkMatA,
	       VariableSubMatrixHandler& WorkMatB,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);
   	virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
					 doublereal dCoef,
					 const VectorHandler& XCurr, 
					 const VectorHandler& XPrimeCurr);

   	virtual void BeforePredict(VectorHandler& X,
				   VectorHandler& XP,
				   VectorHandler& XPrev,
				   VectorHandler& XPPrev) const;   
   	virtual void AfterPredict(VectorHandler& X,
				  VectorHandler& XP);   
	virtual void Update(const VectorHandler& XCurr, 
                       	    const VectorHandler& XPrimeCurr);

#ifdef USE_STRUCT_NODES
   	virtual unsigned int iGetInitialNumDof(void) const;
   	virtual void 
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		      const VectorHandler& XCurr);
   	SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
   	virtual void SetInitialValue(VectorHandler& X) const;   
#endif /* USE_STRUCT_NODES */

   	virtual void SetValue(VectorHandler& X, VectorHandler& XP) const;

   	virtual unsigned int iGetNumPrivData(void) const;
   	virtual doublereal dGetPrivData(unsigned int i) const;
};

inline void* 
LoadableElem::pGet(void) const 
{ 
   	return (void*)this; 
}

inline void* 
LoadableElem::pGetData(void) const
{
   	return priv_data;
}

class DataManager;
class MBDynParser;

extern Elem* 
ReadLoadable(DataManager* pDM, MBDynParser& HP, 
	     const DofOwner* pDO, unsigned int uLabel);

#endif /* LOADABLE_H */

