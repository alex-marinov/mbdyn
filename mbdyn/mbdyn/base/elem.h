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

/*****************************************************************************
 *                                                                           *
 *                           Mbdyn - Elements                                *
 *                                                                           *
 *****************************************************************************/


#ifndef ELEM_H
#define ELEM_H

#include "myassert.h"
#include "except.h"

#include "solman.h"
#include "submat.h"
#include "output.h"

#include "withlab.h"
#include "dofown.h"

#include "memmans.h"

#include "node.h"

/* Tipi di Elem. Lasciare sempre UNKNOWN = -1, cosi' il primo elemento
 * ha tipo zero, e l'ultima entry dell'enum, LAST...TYPE, e' uguale
 * al numero di tipi definiti, quindi puo' essere usata come costante nel 
 * dimensionamento degli arrays e come flag di fine tipi. */
class ElemType {
 public:
   enum Type {
      UNKNOWN = -1,
	FORCE = 0,

	AUTOMATICSTRUCTURAL,
	GRAVITY,
	BODY,
	JOINT,
	BEAM,
	PLATE,

	AIRPROPERTIES,
	ROTOR,
	AERODYNAMIC,

	ELECTRICBULK,
	ELECTRIC,
	GENEL,

	HYDRAULIC,
	
	BULK,
	LOADABLE,
	DRIVEN,
	
	LASTELEMTYPE
   };      
};

extern const char* psElemNames[];
extern const char* psReadControlElems[];
extern const char* psAdamsElemCode[];


/* classi dichiarate */
class ElemWithDofs;
class ElemGravityOwner;
class AerodynamicElem;
class InitialAssemblyElem;
class Rotor;

/* Elem - begin */

class Elem : public WithLabel, public ToBeOutput {
 private:
   ElemType::Type ElemT;
   
 public:
   Elem(unsigned int uL, ElemType::Type T, flag fOut);
   virtual ~Elem(void);

   /* Scrive il contributo dell'elemento al file di restart */
   virtual ostream& Restart(ostream& out) const = 0;
   
   /* Tipo dell'elemento (usato solo per debug ecc.) */
   virtual ElemType::Type GetElemType(void) const = 0;
   
   
   /* funzioni di servizio */

   /* Il metodo iGetNumDof() serve a ritornare il numero di gradi di liberta'
    * propri che l'elemento definisce. Non e' virtuale in quanto serve a 
    * ritornare 0 per gli elementi che non possiedono gradi di liberta'.
    * Viene usato nella costruzione dei DofOwner e quindi deve essere 
    * indipendente da essi. In genere non comporta overhead in quanto il 
    * numero di dof aggiunti da un tipo e' una costante e non richede dati 
    * propri.
    * Il metodo pGetDofOwner() ritorna il puntatore al DofOwner dell'oggetto.
    * E' usato da tutti quelli che agiscono direttamente sui DofOwner.
    * Non e' virtuale in quanto ritorn a NULL per tutti i tipi che non hanno
    * dof propri.
    * Il metodo SetDof() ritorna, per ogni dof dell'elemento, l'ordine.
    * E' usato per completare i singoli Dof relativi all'elemento.
    */
   
   /* ritorna il numero di Dofs per gli elementi che sono anche DofOwners */
   virtual unsigned int iGetNumDof(void) const;
      
   /* esegue operazioni sui dof di proprieta' dell'elemento */
   virtual DofOrder::Order SetDof(unsigned int i) const;

   
   /* funzioni proprie */
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const = 0;
   
   /* assemblaggio matrici per autovalori */
   virtual void AssEig(VariableSubMatrixHandler& WorkMatA,
		       VariableSubMatrixHandler& WorkMatB,
		       const VectorHandler& XCurr,
		       const VectorHandler& XPrimeCurr);

   /* Setta i valori iniziali delle variabili (e fa altre cose) 
    * prima di iniziare l'integrazione */
   virtual void SetValue(VectorHandler& X, VectorHandler& XP) const;
      
   /* Elaborazione vettori e dati prima e dopo la predizione
    * per MultiStepIntegrator */
   virtual void BeforePredict(VectorHandler& X,
      			      VectorHandler& XP,
			      VectorHandler& XPrev,
			      VectorHandler& XPPrev) const;
   
   virtual void AfterPredict(VectorHandler& X,
			     VectorHandler& XP);
      
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr) = 0;

   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr) = 0;


   /* *******PER IL SOLUTORE BLOCK JACOBI-BROYDEN******** */
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual int GetNumConnectedNodes(void) const {
      /* FIXME: modificare in modo che questa funzione ritorni il numero
       * di nodi connessi; quindi viene chiamata la successiva con gli array
       * dimensionati opportunamente
       */
      return 0;
   };
   
   virtual void GetConnectedNodes(int& NumNodes, NodeType::Type* /* NdTyps */ , unsigned int* /* NdLabels */ ) {
#ifdef DEBUG
      cerr << "Warning: probably function Elem::GetConnectedNodes"
        " is not defined for an element Type" << endl;
#endif /* DEBUG */
      NumNodes = 0;
   };
   /* ************************************************ */

   /* Aggiorna dati in base alla soluzione */
   virtual void Update(const VectorHandler& XCurr, 
		       const VectorHandler& XPrimeCurr);
	
   
   /* Metodi per l'estrazione di dati "privati".
    * Si suppone che l'estrattore li sappia interpretare.
    * Come default non ci sono dati privati estraibili */
   virtual unsigned int iGetNumPrivData(void) const;
   virtual doublereal dGetPrivData(unsigned int i) const;
   
   /* Funzioni di casting sicuro verso elementi derivati */
   virtual inline void* pGet(void) const = 0;
   
   virtual Elem* pGetElem(void) const;
   virtual ElemWithDofs* pGetElemWithDofs(void) const;
   virtual ElemGravityOwner* pGetElemGravityOwner(void) const;
   virtual AerodynamicElem* pGetAerodynamicElem(void) const;
   virtual InitialAssemblyElem* pGetInitialAssemblyElem(void) const;
   virtual integer GetRotor(void) const {
     return -1;
   };
   
   
   /* Adams output stuff */
   virtual unsigned int iGetNumAdamsDummyParts(void) const {
      return 0;
   };
   virtual void GetAdamsDummyPart(unsigned int part, Vec3& x, Mat3x3& R) const {
      THROW(ErrGeneric());
   };
   virtual ostream& WriteAdamsDummyPartCmd(ostream& out, unsigned int part, unsigned int firstId) const {
      THROW(ErrGeneric());
   };
};

/* Elem - end */


/* Classe derivata da elem, relativa ad elementi che possiedono gradi di
 * liberta' propri */

/* ElemWithDofs - begin */

class ElemWithDofs : virtual public Elem, public DofOwnerOwner {

 public:
   ElemWithDofs(unsigned int uL, ElemType::Type T, 
		const DofOwner* pDO, flag fOut);
   
   virtual ~ElemWithDofs(void);

   /* Consente di effettuare un casting sicuro da Elem* a ElemWithDofs* */
   virtual ElemWithDofs* pGetElemWithDofs(void) const;
};

/* ElemWithDofs - end */


/* SubjectToInitialAssembly - begin */

class SubjectToInitialAssembly {
 public:
   SubjectToInitialAssembly(void);
   virtual ~SubjectToInitialAssembly(void);
      
   /* Numero di gradi di liberta' definiti durante l'assemblaggio iniziale
    * e' dato dai gradi di liberta' soliti piu' le loro derivate necessarie; 
    * tipicamente per un vincolo di posizione il numero di dof raddoppia, in
    * quanto vengono derivate tutte le equazioni, mentre per un vincolo di 
    * velocita' rimane inalterato. Sono possibili casi intermedi per vincoli
    * misti di posizione e velocita' */
   virtual unsigned int iGetInitialNumDof(void) const = 0;
   
   /* Dimensione del workspace durante l'assemblaggio iniziale. Occorre tener
    * conto del numero di dof che l'elemento definisce in questa fase e dei
    * dof dei nodi che vengono utilizzati. Sono considerati dof indipendenti
    * la posizione e la velocita' dei nodi */
   virtual void InitialWorkSpaceDim(integer* piNumRows, 
				    integer* piNumCols) const = 0;
   
   /* Contributo allo jacobiano durante l'assemblaggio iniziale */
   virtual VariableSubMatrixHandler& 
     InitialAssJac(VariableSubMatrixHandler& WorkMat,
		   const VectorHandler& XCurr) = 0;
   
   /* Contributo al residuo durante l'assemblaggio iniziale */   
   virtual SubVectorHandler& 
     InitialAssRes(SubVectorHandler& WorkVec,
		   const VectorHandler& XCurr) = 0;
};

/* SubjectToInitialAssembly - end */


/* InitialAssemblyElem - begin */

class InitialAssemblyElem 
: virtual public Elem, public SubjectToInitialAssembly {
 public:
   InitialAssemblyElem(unsigned int uL, ElemType::Type T, flag fOut);
   virtual ~InitialAssemblyElem(void);

   /* Consente di effettuare un casting sicuro da Elem* a InitialAssemblyElem* */
   virtual InitialAssemblyElem* pGetInitialAssemblyElem(void) const;
};
   
/* InitialAssemblyElem - end */

#endif
