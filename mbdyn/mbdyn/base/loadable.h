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

#ifndef LOADABLE_H
# define LOADABLE_H

# include <elem.h>
# include <dataman.h>
# ifdef USE_STRUCT_NODES
#  include <gravity.h>
# endif /* USE_STRUCT_NODES */
// #include <aerodyn.h>

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

class LoadableElem
: virtual public Elem,
# ifdef USE_STRUCT_NODES
public InitialAssemblyElem,
// public AerodynamicElem, 
public ElemGravityOwner,
# endif /* USE_STRUCT_NODES */
public ElemWithDofs {
 protected:
   void* priv_data;
   char* module_name;
   void* handle;

   /* Funzioni attese */
   enum Funcs {     
	READ = 0,
      
	IGETNUMDOF,
	SETDOF,
      
	OUTPUT,
	RESTART,
      
	WORKSPACEDIM,
	ASSJAC,
        ASSEIG,
	ASSRES,
      
        BEFOREPREDICT,
        AFTERPREDICT,
        UPDATE,
      
	IGETINITIALNUMDOF,
	INITIALWORKSPACEDIM,
	INITIALASSJAC,
	INITIALASSRES,

	SETVALUE,
	SETINITIALVALUE,

        IGETNUMPRIVDATA,
        DGETPRIVDATA,            

	// ...

	DESTROY,
	LASTFUNC
   };

   /* Simboli delle funzioni attese */
   void * fsym[LASTFUNC];
   
 public:
   LoadableElem(unsigned int uLabel, const DofOwner* pDO,
		DataManager* pDM, MBDynParser& HP);
   ~LoadableElem(void); 
   virtual inline void* pGet(void) const;
   
   inline void* pGetData(void) const;
   
   virtual ElemType::Type GetElemType(void) const;

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

# ifdef USE_STRUCT_NODES
   virtual unsigned int iGetInitialNumDof(void) const;
   virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat, const VectorHandler& XCurr);
   SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
   virtual void SetInitialValue(VectorHandler& X) const;   
# endif /* USE_STRUCT_NODES */

   virtual void SetValue(VectorHandler& X, VectorHandler& XP) const;

   virtual unsigned int iGetNumPrivData(void) const;
   virtual doublereal dGetPrivData(unsigned int i) const;
};

inline void* LoadableElem::pGet(void) const 
{ 
   return (void*)this; 
}

inline void* LoadableElem::pGetData(void) const
{
   return priv_data;
}


class DataManager;
class MBDynParser;

extern Elem* ReadLoadable(DataManager* pDM, 
			  MBDynParser& HP, 
			  const DofOwner* pDO, 
			  unsigned int uLabel);
#endif /* LOADABLE_H */
