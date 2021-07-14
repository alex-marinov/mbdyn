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

#ifndef LOADABLE_H
#define LOADABLE_H

/*
 * Policy for LOADABLE_VERSION changes:
 * <maj>	is increased when major changes in element API occur
 * <min>	is increased when loadable elements data structure
 * 		changes (e.g. a new call is added)
 * <fix>	is increased when minor fixes are added to loadable
 * 		elements handling, which do not involve data
 * 		structure changes.
 */

#define LOADABLE_VERSION_SET(maj, min, fix)	\
	(((maj) << 24) | ((min) << 16) | (fix))
#define LOADABLE_VERSION	LOADABLE_VERSION_SET(1, 5, 0)
#define LOADABLE_VERSION_OUT(v) \
	((v & 0xFF000000U) >> 24) << '.' << ((v & 0x00FF0000U) >> 16) << '.' << (v & 0x0000FFFFU)
/*
 * CHANGELOG:
 * 2006-02-08: 1.5.0 changed connectivity API
 * 2005-10-11: 1.4.0 changed SetValue API (hints and so...)
 * 2004-06-13: 1.3.0 allow modules to be statically linked
 * 2004-05-21: 1.2.0 added hooks for dummy part output calls
 * 2003-02-25: 1.1.0 added hook for iGetPrivDataIdx()
 * 2002-XX-XX: 1.0.0 added versioning system to detect structure conflicts
 */

#ifdef USE_RUNTIME_LOADING
#include <ltdl.h>
#endif // USE_RUNTIME_LOADING

#include "userelem.h"

/*
 * Dynamically loadable element
 *
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
 *                MBDynParser&)
 *
 * in questo modo, tramite il puntatore all'elemento, si ha accesso
 * a tutti i dati dell'elemento;
 * rimane da studiare la questione dei diritti di accesso.
 */

class LoadableElem;

typedef void *
(* p_read)(LoadableElem*, DataManager*, MBDynParser&);
typedef unsigned int (* p_i_get_num_dof)(const LoadableElem*);
typedef DofOrder::Order (* p_set_dof)(const LoadableElem*, unsigned int);
typedef void (* p_output)(const LoadableElem*, OutputHandler&);
typedef std::ostream& (* p_restart)(const LoadableElem*, std::ostream&);
typedef void (* p_work_space_dim)(const LoadableElem*, integer*, integer*);
typedef VariableSubMatrixHandler&
(* p_ass_jac)(LoadableElem*, VariableSubMatrixHandler&, doublereal,
		const VectorHandler&, const VectorHandler&);
typedef void
(* p_ass_mats)(LoadableElem*, VariableSubMatrixHandler&,
		VariableSubMatrixHandler&,
		const VectorHandler&, const VectorHandler&);
typedef SubVectorHandler&
(* p_ass_res)(LoadableElem*, SubVectorHandler&, doublereal,
		const VectorHandler&, const VectorHandler&);
typedef void
(* p_before_predict)(const LoadableElem* pEl,
		VectorHandler& X, VectorHandler& XP,
		VectorHandler& XPrev, VectorHandler& XPPrev);
typedef void
(* p_after_predict)(const LoadableElem* pEl,
		VectorHandler& X, VectorHandler& XP);
typedef void
(* p_update)(LoadableElem* pEl,
		const VectorHandler& X, const VectorHandler& XP);
typedef void
(* p_after_convergence)(const LoadableElem* pEl,
    		const VectorHandler& X, const VectorHandler& XP);
typedef unsigned int (* p_i_get_initial_num_dof)(const LoadableElem*);
typedef void
(* p_initial_work_space_dim)(const LoadableElem*, integer*, integer*);
typedef VariableSubMatrixHandler&
(* p_initial_ass_jac)(LoadableElem*, VariableSubMatrixHandler &,
		const VectorHandler&);
typedef SubVectorHandler&
(* p_initial_ass_res)(LoadableElem*, SubVectorHandler&, const VectorHandler&);
typedef void
(* p_set_value)(const LoadableElem*, DataManager *pDM,
		VectorHandler&, VectorHandler&,
		SimulationEntity::Hints *ph);
typedef void (* p_set_initial_value)(const LoadableElem*, VectorHandler&);
typedef unsigned int (* p_i_get_num_priv_data)(const LoadableElem* pEl);
typedef unsigned int
(* p_i_get_priv_data_idx)(const LoadableElem* pEl, const char *);
typedef doublereal
(* p_d_get_priv_data)(const LoadableElem* pEl, unsigned int i);
typedef int (* p_i_get_num_connected_nodes)(const LoadableElem*);
typedef void
(* p_get_connected_nodes)(const LoadableElem*,
		std::vector<const Node *>& connectedNodes);
typedef void (* p_destroy)(LoadableElem*);

/*
 * Struttura che contiene le chiamate alle funzioni del modulo;
 * e' l'unico contatto tra modulo ed elemento caricato runtime.
 * Se un puntatore e' vuoto, viene riempito con il metodo di default.
 */
struct LoadableCalls {
	/*
	 * This is mandatory; use the macro

	   LOADABLE_VERSION_SET(major, minor, fix)

	 * to initialize to the API version your module is using.
	 */
	unsigned long			loadable_version;

	/*
	 * These are optional; they are printed at module loading
	 * for informational purposes.
	 */
	const char *			name;
	const char *			version;
	const char *			vendor;
	const char *			description;

	/*
	 * This is mandatory; use it to initialize and read any data
	 * from the input stream.
	 */
	p_read				read;

	/*
	 * These are optional; fill-in as required by your module.
	 */
	p_i_get_num_dof			i_get_num_dof;
	p_set_dof 			set_dof;
	p_output 			output;
	p_restart 			restart;
	p_work_space_dim 		work_space_dim;
	p_ass_jac			ass_jac;
	p_ass_mats			ass_mats;
	p_ass_res			ass_res;
	p_before_predict		before_predict;
	p_after_predict			after_predict;
	p_update			update;
	p_after_convergence		after_convergence;
	p_i_get_initial_num_dof		i_get_initial_num_dof;
	p_initial_work_space_dim	initial_work_space_dim;
	p_initial_ass_jac		initial_ass_jac;
	p_initial_ass_res		initial_ass_res;
	p_set_value			set_value;
	p_set_initial_value		set_initial_value;
	p_i_get_num_priv_data		i_get_num_priv_data;
	p_i_get_priv_data_idx		i_get_priv_data_idx;
	p_d_get_priv_data		d_get_priv_data;
	p_i_get_num_connected_nodes	i_get_num_connected_nodes;
	p_get_connected_nodes		get_connected_nodes;

	p_destroy			destroy;

};

class LoadableElem :
	virtual public Elem,
	public UserDefinedElem
{
protected:
   	void* priv_data;	/* Dati privati passati alle funzioni */
   	char* module_name;	/* Nome del modulo */
#ifdef USE_RUNTIME_LOADING
	lt_dlhandle handle;
#endif // USE_RUNTIME_LOADING
	LoadableCalls *calls;	/* Simboli delle funzioni attese */

	void GetCalls(MBDynParser& HP);
	void BindCalls(DataManager* pDM, MBDynParser& HP);

public:
   	LoadableElem(unsigned int uLabel, const DofOwner* pDO,
		     DataManager* pDM, MBDynParser& HP);
   	LoadableElem(unsigned int uLabel, const DofOwner* pDO,
			const LoadableCalls *c,
			DataManager* pDM, MBDynParser& HP);
   	~LoadableElem(void);

   	inline void* pGetData(void) const;

   	virtual unsigned int iGetNumDof(void) const;
   	virtual DofOrder::Order GetDofType(unsigned int i) const;

   	virtual void Output(OutputHandler& OH) const;
   	virtual std::ostream& Restart(std::ostream& out) const;

   	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	virtual VariableSubMatrixHandler&
     	AssJac(VariableSubMatrixHandler& WorkMat,
	       doublereal dCoef,
	       const VectorHandler& XCurr,
	       const VectorHandler& XPrimeCurr);
	virtual void
     	AssMats(VariableSubMatrixHandler& WorkMatA,
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
   	virtual void AfterConvergence(const VectorHandler& X,
			const VectorHandler& XP);

   	virtual unsigned int iGetInitialNumDof(void) const;
   	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		      const VectorHandler& XCurr);
   	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
   	virtual void SetInitialValue(VectorHandler& X);

   	virtual void SetValue(DataManager *pDM,
			VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph = 0);

   	virtual unsigned int iGetNumPrivData(void) const;
   	virtual unsigned int iGetPrivDataIdx(const char *s) const;
   	virtual doublereal dGetPrivData(unsigned int i) const;

	/* *******PER IL SOLUTORE BLOCK JACOBI-BROYDEN******** */
     	/*
	 * Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs
	 */
     	virtual int GetNumConnectedNodes(void) const;
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
     	/* ************************************************ */

};

inline void*
LoadableElem::pGetData(void) const
{
   	return priv_data;
}

class DataManager;
class MBDynParser;

// base class for user-defined element parsing
struct LoadableElemRead : public UserDefinedElemRead {
	virtual ~LoadableElemRead( void ) { NO_OP; };
	virtual UserDefinedElem *
	Read(unsigned uLabel, const DofOwner* pDO,
		DataManager* const pDM, MBDynParser& HP) const;
};

extern Elem*
ReadLoadable(DataManager* pDM, MBDynParser& HP,
     const DofOwner* pDO, unsigned int uLabel);

#endif // LOADABLE_H
