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

/* Elementi di rotore */

#ifndef ROTOR_H
#define ROTOR_H

extern "C" {
#include <ac/float.h>
}

#ifdef USE_MPI
#include <mpi++.h>
#endif /* USE_MPI */

#include <aerodyn.h>
#include <strnode.h>
#include <resforces.h>

extern const char* psRotorNames[];

/* Rotor - begin */

class Rotor 
: virtual public Elem, public AerodynamicElem, public ElemWithDofs {
/* Tipi di rotori */
 public:
   enum Type {
      UNKNOWN = -1,
	NO = 0,
	UNIFORM,
	GLAUERT,
	MANGLER,
	DYNAMICINFLOW,
	
	LASTROTORTYPE
   };

 public:
   class ErrInfiniteMeanInducedVelocity {};
   
 protected:
#ifdef USE_MPI
   /* Communicator per il calcolo della trazione totale */
   bool is_parallel;
   MPI::Intracomm RotorComm;  /* Gruppo di macchine su cui si scambiano dati relativi al rotore */  
   int* pBlockLenght;    
   MPI::Aint* pDispl;        /* vettore di indirizzi di dati */
   MPI::Request ReqV;        /* request per le comunicazioni Send Receive */
   integer iForcesVecDim;    /* dimensioni del vettore forze scambiato fra i processi */ 
   doublereal*  pTmpVecR;     /* vettori temporanei per lo scambio delle forze */
   doublereal*  pTmpVecS;
   MPI::Datatype* pRotDataType; /* datatype che contiene le posizioni dei dati da scambiare 
				 * ad ogni passo */
#endif /* USE_MPI */             
   const StructNode* pCraft;
   const StructNode* pRotor;

   doublereal dOmegaRef;     /* Velocita' di rotazione di riferimento */
   
   doublereal dRadius;       /* Raggio del rotore */
   doublereal dArea;         /* Area del disco */
   doublereal dUMean;        /* Velocita' indotta media */
   doublereal dUMeanPrev;    /* Vel. indotta media al passo prec. */
   doublereal dWeight;       /* Peso della velocita' indotta media 
			      * (peso della V al passo precedente, def = 0.) */
   doublereal dCorrection;   /* Correzione (scala la velocita' indotta) */
   
   ExternResForces Res;	     /* force, couple and pole for resultants */
   SetResForces **ppRes;     /* extra forces */

   Mat3x3 RRotTranspose;     /* Trasposta della matrice di rotazione rotore */
   Vec3 RRot3;               /* Direzione dell'asse del rotore */
   Vec3 VCraft;              /* Velocita' di traslazione del velivolo */
   doublereal dPsi0;         /* Angolo di azimuth del rotore */
   doublereal dSinAlphad;    /* Angolo di influsso */
   doublereal dCosAlphad;    /* */
   doublereal dMu;           /* */
   doublereal dLambda;       /* Parametri di influsso */
   doublereal dChi;          /* */   
   
   doublereal dVelocity;     /* Velocita' di riferimento */
   doublereal dOmega;        /* Velocita' di rotazione di riferimento */

   

   /* temporaneo */
   mutable int iNumSteps;
   
   
   /* Questa funzione non viene usata, da Rotor, ma da altre classi derivate 
    * con modellazioni piu' sofisticate della velocita' indotta */
   virtual doublereal dGetPsi(const Vec3& X) const;
   
   /* Calcola la distanza di un punto dall'asse di rotazione in coordinate 
    * adimensionali */
   virtual doublereal dGetPos(const Vec3& X) const;
   
   /* Calcola la velocita' di traslazione del rotore */
   virtual void InitParam(void);
   
   /* Calcola la velocita' indotta media (uniforme) */
   virtual void MeanInducedVelocity(void);   

 public:
   Rotor(unsigned int uL, const DofOwner* pDO, 
	 const StructNode* pC, const StructNode* pR, 
	 SetResForces **ppres, flag fOut);
   virtual ~Rotor(void);      
   
   /* funzioni di servizio */

   /* Tipo dell'elemento (usato per debug ecc.) */
   virtual Elem::Type GetElemType(void) const;
   
   /* Il metodo iGetNumDof() serve a ritornare il numero di gradi di liberta'
    * propri che l'elemento definisce. Non e' virtuale in quanto serve a 
    * ritornare 0 per gli elementi che non possiedono gradi di liberta'.
    * Viene usato nella costruzione dei DofOwner e quindi deve essere 
    * indipendente da essi. In genere non comporta overhead in quanto il 
    * numero di dof aggiunti da un tipo e' una costante e non richede dati 
    * propri.
    * Il metodo pGetDofOwner() ritorna il puntatore al DofOwner dell'oggetto.
    * E' usato da tutti quelli che agiscono direttamente sui DofOwner.
    * Non e' virtuale in quanto ritorna NULL per tutti i tipi che non hanno
    * dof propri.
    * Il metodo SetDof() ritorna, per ogni dof dell'elemento, l'ordine.
    * E' usato per completare i singoli Dof relativi all'elemento.
    */   
   
   /* ritorna il numero di Dofs per gli elementi che sono anche DofOwners */
   virtual unsigned int iGetNumDof(void) const { 
      return 0;
   };
      
   /* esegue operazioni sui dof di proprieta' dell'elemento */
#ifdef DEBUG   
   virtual DofOrder::Order SetDof(unsigned int i) const
#else     
   virtual DofOrder::Order SetDof(unsigned int /* i */ ) const
#endif
    { 
       ASSERT(i >= 0 && i < this->iGetNumDof());
       return DofOrder::DIFFERENTIAL; 
    };

   
   /* funzioni proprie */
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 0;
      *piNumCols = 0;
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);
   
   /* output; si assume che ogni tipo di elemento sappia, attraverso
    * l'OutputHandler, dove scrivere il proprio output */
   virtual void Output(OutputHandler& OH) const;

   /* Contributo al file di Restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   /* Relativo ai ...WithDofs */
   virtual void SetInitialValue(VectorHandler& /* X */ ) const {
      NO_OP;
   };
   
   /* Relativo ai ...WithDofs */
   virtual void SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const {
      NO_OP;
   };   

   /* Tipo di rotore */
   virtual Rotor::Type GetRotorType(void) const = 0;

   virtual inline const Vec3& GetXCurr(void) const {
      return pRotor->GetXCurr();
   };
   
   /* accesso a dati */
   virtual inline doublereal dGetOmega(void) const {
      return dOmega;
   };
   
   virtual inline doublereal dGetRadius(void) const {
      return dRadius;
   };
   
   virtual inline doublereal dGetMu(void) const {
      return dMu;
   };
   
   virtual inline const Vec3& GetForces(void) const {
      return Res.Force();
   };
   
   virtual inline const Vec3& GetMoments(void) const {
      return Res.Couple();
   };
   
   /* Metodi per l'estrazione di dati "privati".
    * Si suppone che l'estrattore li sappia interpretare.
    * Come default non ci sono dati privati estraibili */
   virtual unsigned int iGetNumPrivData(void) const {
      return 6;
   };
   
   virtual doublereal dGetPrivData(unsigned int i) const {
      if (i >= 1 && i <= 3) {
	 return Res.Force().dGet(i);
      } else if (i >= 4 && i <= 6) {      
	 return Res.Couple().dGet(i-3);
      } else {
	 THROW(ErrGeneric());
      }
#ifndef USE_EXCEPTIONS
      return 0.;
#endif /* USE_EXCEPTIONS */
   };
   
   
   /* Somma alla trazione il contributo di un elemento */
   virtual void AddForce(unsigned int uL, 
		   const Vec3& F, const Vec3& M, const Vec3& X);
   virtual void ResetForce(void);
   
   /* Restituisce ad un elemento la velocita' indotta in base alla posizione
    * azimuthale */
   virtual Vec3 GetInducedVelocity(const Vec3& X) const = 0;
   
#ifdef USE_MPI
   void ExchangeTraction(flag fWhat);
#endif /* USE_MPI */

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps, unsigned int* NdLabels) {
     NumNodes = 2;
     NdTyps[0] = pCraft->GetNodeType();
     NdLabels[0] = pCraft->GetLabel();
     NdTyps[1] = pRotor->GetNodeType();
     NdLabels[1] = pRotor->GetLabel();
   };
   /* ************************************************ */
#ifdef USE_MPI
   void InitializeRotorComm(MPI::Intracomm* Rot);
   void  ExchangeVelocity(void); 
#endif /* USE_MPI */
};

/* Rotor - end */


/* NoRotor - begin */

class NoRotor : virtual public Elem, public Rotor {
 protected:
   
 public:
   NoRotor(unsigned int uLabel,
	   const DofOwner* pDO,
	   const StructNode* pCraft, 
	   const StructNode* pRotor, 
	   SetResForces **ppres, 
	   doublereal dR, 
	   flag fOut);
   virtual ~NoRotor(void);
   virtual inline void* pGet(void) const { return (void*)this; };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);   

   /* Contributo al file di Restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   virtual Rotor::Type GetRotorType(void) const {
     return Rotor::NO;
   };
   
   /* Somma alla trazione il contributo di un elemento */
   virtual void AddForce(unsigned int uL,
		   const Vec3& F, const Vec3& M, const Vec3& X);
   
   /* Restituisce ad un elemento la velocita' indotta in base alla posizione
    * azimuthale */
   virtual Vec3 GetInducedVelocity(const Vec3& X) const;
 


};

/* NoRotor - end */


/* UniformRotor - begin */

class UniformRotor : virtual public Elem, public Rotor {
 protected:
 public:
   UniformRotor(unsigned int uLabel, 
		const DofOwner* pDO,
		const StructNode* pCraft, 
		const StructNode* pRotor, 
		SetResForces **ppres, 
		doublereal dOR,
		doublereal dR,
		doublereal dW,
		doublereal dC,
		flag fOut);
   virtual ~UniformRotor(void);
   virtual inline void* pGet(void) const { return (void*)this; };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);   

   /* Contributo al file di Restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   Rotor::Type GetRotorType(void) const {
      return Rotor::UNIFORM;
   };
   
   /* Somma alla trazione il contributo di un elemento */
   virtual void AddForce(unsigned int uL,
		   const Vec3& F, const Vec3& M, const Vec3& X);
   
   /* Restituisce ad un elemento la velocita' indotta in base alla posizione
    * azimuthale */
   virtual Vec3 GetInducedVelocity(const Vec3& X) const;

};

/* UniformRotor - end */


/* GlauertRotor - begin */

class GlauertRotor : virtual public Elem, public Rotor {
 protected:
 public:
   GlauertRotor(unsigned int uLabel, 
		const DofOwner* pDO,
		const StructNode* pCraft,
		const StructNode* pRotor,
		SetResForces **ppres, 
		doublereal dOR,
		doublereal dR,
		doublereal dW,
		doublereal dC,
		flag fOut);
   virtual ~GlauertRotor(void);
   virtual inline void* pGet(void) const { return (void*)this; };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);   

   /* Contributo al file di Restart */
   virtual std::ostream& Restart(std::ostream& out) const;   
   
   Rotor::Type GetRotorType(void) const {
      return Rotor::GLAUERT;
   };
   
   /* Somma alla trazione il contributo di un elemento */
   virtual void AddForce(unsigned int uL,
		   const Vec3& F, const Vec3& M, const Vec3& X);
   
   /* Restituisce ad un elemento la velocita' indotta in base alla posizione
    * azimuthale */
   virtual Vec3 GetInducedVelocity(const Vec3& X) const;
   
};

/* GlauertRotor - end */


/* ManglerRotor - begin */

class ManglerRotor : virtual public Elem, public Rotor {
 protected:
 public:
   ManglerRotor(unsigned int uLabel, 
		const DofOwner* pDO,
		const StructNode* pCraft,
		const StructNode* pRotor,
		SetResForces **ppres,
		doublereal dOR,
		doublereal dR,
		doublereal dW,
		doublereal dC,
		flag fOut);
   virtual ~ManglerRotor(void);
   virtual inline void* pGet(void) const { return (void*)this; };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);   

   /* Contributo al file di Restart */
   virtual std::ostream& Restart(std::ostream& out) const;   

   Rotor::Type GetRotorType(void) const {
      return Rotor::MANGLER;
   };
   
   /* Somma alla trazione il contributo di un elemento */
   virtual void AddForce(unsigned int uL,
		   const Vec3& F, const Vec3& M, const Vec3& X);
   
   /* Restituisce ad un elemento la velocita' indotta in base alla posizione
    * azimuthale */
   virtual Vec3 GetInducedVelocity(const Vec3& X) const;

};

/* ManglerRotor - end */


/* DynamicInflowRotor - begin */

class DynamicInflowRotor : virtual public Elem, public Rotor {
 protected:
   doublereal dVConst;
   doublereal dVCosine;
   doublereal dVSine;

   doublereal dL11;
   doublereal dL13;
   doublereal dL22;
   doublereal dL31;
   doublereal dL33;
   
 public:
   DynamicInflowRotor(unsigned int uLabel,
		      const DofOwner* pDO,
		      const StructNode* pCraft, 
		      const StructNode* pRotor,
		      SetResForces **ppres, 
		      doublereal dOR,
		      doublereal dR,
		      doublereal dVConstTmp,
		      doublereal dVCosineTmp,
		      doublereal dVSineTmp,
		      flag fOut);
   virtual ~DynamicInflowRotor(void);
   virtual inline void* pGet(void) const { return (void*)this; };

   /* ritorna il numero di Dofs per gli elementi che sono anche DofOwners */
   virtual unsigned int iGetNumDof(void) const { 
      return 3;
   };
   
   /* output; si assume che ogni tipo di elemento sappia, attraverso
    * l'OutputHandler, dove scrivere il proprio output */
   virtual void Output(OutputHandler& OH) const;
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 3;
      *piNumCols = 3;
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal dCoef, 
	    const VectorHandler& XCurr,
	    const VectorHandler& XPrimeCurr);
            
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);   
   
   /* Contributo al file di Restart */
   virtual std::ostream& Restart(std::ostream& out) const;   
   
   /* Relativo ai ...WithDofs */
   virtual void SetInitialValue(VectorHandler& X) const;
   
   /* Relativo ai ...WithDofs */
   virtual void SetValue(VectorHandler& X, VectorHandler& XP) const;   

   Rotor::Type GetRotorType(void) const {
      return Rotor::DYNAMICINFLOW;
   };
   
   /* Somma alla trazione il contributo di un elemento */
   virtual void AddForce(unsigned int uL,
		   const Vec3& F, const Vec3& M, const Vec3& X);
   
   /* Restituisce ad un elemento la velocita' indotta in base alla posizione
    * azimuthale */
   virtual Vec3 GetInducedVelocity(const Vec3& X) const;

};

/* DynamicInflowRotor - end */

class DataManager;
class MBDynParser;

extern Elem* ReadRotor(DataManager* pDM,
		       MBDynParser& HP, 
		       const DofOwner* pDO, 
		       unsigned int uLabel);

#endif
