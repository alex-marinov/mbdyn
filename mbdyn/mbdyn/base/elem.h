/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

/*****************************************************************************
 *                                                                           *
 *                           Mbdyn - Elements                                *
 *                                                                           *
 *****************************************************************************/


#ifndef ELEM_H
#define ELEM_H

#include <myassert.h>
#include <except.h>

#include <solman.h>
#include <submat.h>
#include <output.h>

#include <withlab.h>
#include <dofown.h>

#include <simentity.h>
#include <node.h>

#ifdef USE_MULTITHREAD
#include <veciter.h>
#endif /* USE_MULTITHREAD */

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

class Elem : public WithLabel, public SimulationEntity, public ToBeOutput
#ifdef USE_MULTITHREAD
, public InUse
#endif /* USE_MULTITHREAD */
{
   /*
    * Tipi di Elem. Lasciare sempre UNKNOWN = -1, cosi' il primo elemento
    * ha tipo zero, e l'ultima entry dell'enum, LAST...TYPE, e' uguale
    * al numero di tipi definiti, quindi puo' essere usata come costante nel 
    * dimensionamento degli arrays e come flag di fine tipi. 
    */
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
        AEROMODAL,
	AERODYNAMIC,

	ELECTRICBULK,
	ELECTRIC,
	GENEL,

	HYDRAULIC,
	
	BULK,
	LOADABLE,
	DRIVEN,
	EXTERNAL,

	RTAI_OUTPUT,
	
	LASTELEMTYPE
   };

 private:
   Elem::Type ElemT;

 public:
   Elem(unsigned int uL, Elem::Type T, flag fOut);
   virtual ~Elem(void);

   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const = 0;
   
   /* Tipo dell'elemento (usato solo per debug ecc.) */
   virtual Elem::Type GetElemType(void) const = 0;
  
   /* inherited from SimulationEntity */
   virtual unsigned int iGetNumDof(void) const;
   virtual DofOrder::Order GetDofType(unsigned int) const;
   
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
    * Il metodo GetDofType() ritorna, per ogni dof dell'elemento, l'ordine.
    * E' usato per completare i singoli Dof relativi all'elemento.
    */
   
   /* funzioni proprie */
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const = 0;
   
   /* assemblaggio matrici per autovalori */
   virtual void AssMats(VariableSubMatrixHandler& WorkMatA,
		       VariableSubMatrixHandler& WorkMatB,
		       const VectorHandler& XCurr,
		       const VectorHandler& XPrimeCurr);

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
      std::cerr << psElemNames[GetElemType()] << "(" << GetLabel() 
         << ") cannot be used in parallel environment" << std::endl;
      THROW(ErrGeneric());
   };
   
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* /* NdTyps */ , unsigned int* /* NdLabels */ ) {
      NumNodes = GetNumConnectedNodes();
   };
   /* ************************************************ */

   /* Funzioni di casting sicuro verso elementi derivati */
   virtual void* pGet(void) const = 0;
   
   virtual Elem* pGetElem(void) const;
   virtual ElemWithDofs* pGetElemWithDofs(void) const;
   virtual ElemGravityOwner* pGetElemGravityOwner(void) const;
   virtual AerodynamicElem* pGetAerodynamicElem(void) const;
   virtual InitialAssemblyElem* pGetInitialAssemblyElem(void) const;

#ifdef USE_ADAMS 
   /* Adams output stuff */
   virtual unsigned int iGetNumAdamsDummyParts(void) const {
      return 0;
   };
   virtual void GetAdamsDummyPart(unsigned int part, Vec3& x, Mat3x3& R) const {
      THROW(ErrGeneric());
   };
   virtual std::ostream& WriteAdamsDummyPartCmd(std::ostream& out, unsigned int part, unsigned int firstId) const {
      THROW(ErrGeneric());
   };
#endif /* USE_ADAMS */
};

/* Elem - end */


/* Classe derivata da elem, relativa ad elementi che possiedono gradi di
 * liberta' propri */

/* ElemWithDofs - begin */

class ElemWithDofs : virtual public Elem, public DofOwnerOwner {

 public:
   ElemWithDofs(unsigned int uL, Elem::Type T, 
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
   InitialAssemblyElem(unsigned int uL, Elem::Type T, flag fOut);
   virtual ~InitialAssemblyElem(void);

   /* Consente di effettuare un casting sicuro da Elem* a InitialAssemblyElem* */
   virtual InitialAssemblyElem* pGetInitialAssemblyElem(void) const;
};
   
/* InitialAssemblyElem - end */

#endif /* ELEM_H */

