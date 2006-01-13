/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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
#include "ac/pthread.h"

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
   MPI::Intracomm RotorComm;  	/* Gruppo di macchine su cui si scambiano
				   dati relativi al rotore */  
   int* pBlockLenght;    
   MPI::Aint* pDispl;        	/* vettore di indirizzi di dati */
   MPI::Request ReqV;        	/* request per le comunicazioni Send Receive */
   integer iForcesVecDim;    	/* dimensioni vettore forze scambiato
				   fra i processi */ 
   doublereal*  pTmpVecR;	/* vettori temporanei per scambio forze */
   doublereal*  pTmpVecS;
   MPI::Datatype* pRotDataType; /* datatype che contiene le posizioni 
				   dei dati da scambiare ad ogni passo */
#endif /* USE_MPI */

#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
   mutable pthread_mutex_t forces_mutex;

   mutable pthread_mutex_t induced_velocity_mutex;
   mutable pthread_cond_t induced_velocity_cond;
   mutable bool bDone;

   void Wait(void) const;
   void Done(void) const;
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES*/

   const StructNode* pCraft;
   const StructNode* pRotor;
   const StructNode* pGround;

   doublereal dOmegaRef;     /* Velocita' di rotazione di riferimento */
   
   doublereal dRadius;       /* Raggio del rotore */
   doublereal dArea;         /* Area del disco */
   doublereal dUMean;        /* Velocita' indotta media */
   doublereal dUMeanRef;     /* Velocita' indotta media (nominale) */
   mutable doublereal dUMeanPrev;    /* Vel. indotta media al passo prec. */

   /* iterations for dUMeanRef */
   unsigned int iMaxIter;
   unsigned int iCurrIter;
   doublereal dTolerance;
   doublereal dEta;
   bool bUMeanRefConverged;

   DriveOwner Weight;
   doublereal dWeight;       /* Peso della velocita' indotta media 
			      * (peso della V al passo precedente, def = 0.) */
   doublereal dHoverCorrection; /* Correzione H (scala la velocita' indotta) */
   doublereal dForwardFlightCorrection; /* FF */
   
   ExternResForces Res;	     /* force, couple and pole for resultants */
   ResForceSet **ppRes;     /* extra forces */

   Mat3x3 RRotTranspose;     /* Trasposta della matrice di rotazione rotore */
   Mat3x3 RRot;
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

   /* Combina i due ... */
   virtual void GetPos(const Vec3& X, doublereal& dr, doublereal& dp) const;
   
   /* Calcola la velocita' di traslazione del rotore */
   virtual void InitParam(bool bComputeMeanInducedVelocity = true);
   
 public:
   Rotor(unsigned int uL, const DofOwner* pDO, 
	 const StructNode* pC, const Mat3x3& rrot, const StructNode* pR, 
	 const StructNode* pG, 
	 unsigned int iMaxIt, doublereal dTol, doublereal dE,
	 ResForceSet **ppres, flag fOut);
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
    * Il metodo GetDofType() ritorna, per ogni dof dell'elemento, l'ordine.
    * E' usato per completare i singoli Dof relativi all'elemento.
    */   
   
   /* ritorna il numero di Dofs per gli elementi che sono anche DofOwners */
   virtual unsigned int iGetNumDof(void) const { 
      return 0;
   };
      
   /* esegue operazioni sui dof di proprieta' dell'elemento */
#ifdef DEBUG   
   virtual DofOrder::Order GetDofType(unsigned int i) const
#else     
   virtual DofOrder::Order GetDofType(unsigned int /* i */ ) const
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
   
   /*
    * Elaborazione stato interno dopo la convergenza
    */
   virtual void AfterConvergence(const VectorHandler& X, 
		   const VectorHandler& XP);

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
   virtual void SetValue(DataManager *pDM,
		   VectorHandler& /* X */ , VectorHandler& /* XP */ ,
		   SimulationEntity::Hints *ph = 0)
   {
      NO_OP;
   };   

   /* Tipo di rotore */
   virtual Rotor::Type GetRotorType(void) const = 0;

   virtual inline const Vec3& GetXCurr(void) const {
      return pRotor->GetXCurr();
   };
   
   /* accesso a dati */
   virtual inline doublereal dGetOmega(void) const {
#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
      Wait();
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES*/
      return dOmega;
   };
   
   virtual inline doublereal dGetRadius(void) const {
#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
      Wait();
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */
      return dRadius;
   };
   
   virtual inline doublereal dGetMu(void) const {
#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
      Wait();
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */
      return dMu;
   };
   
   virtual inline const Vec3& GetForces(void) const {
#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
      Wait();
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */
      return Res.Force();
   };
   
   virtual inline const Vec3& GetMoments(void) const {
#if defined(USE_MULTITHREAD) && defined(MBDYN_X_MT_ASSRES)
      Wait();
#endif /* USE_MULTITHREAD && MBDYN_X_MT_ASSRES */
      return Res.Couple();
   };
   
   /* Metodi per l'estrazione di dati "privati".
    * Si suppone che l'estrattore li sappia interpretare.
    * Come default non ci sono dati privati estraibili */
   virtual unsigned int iGetNumPrivData(void) const;
   virtual unsigned int iGetPrivDataIdx(const char *s) const;
   virtual doublereal dGetPrivData(unsigned int i) const;
   
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
   virtual void GetConnectedNodes(int& NumNodes, Node::Type* NdTyps,
		   unsigned int* NdLabels) {
	   NumNodes = 2;
	   NdTyps[0] = pCraft->GetNodeType();
	   NdLabels[0] = pCraft->GetLabel();
	   NdTyps[1] = pRotor->GetNodeType();
	   NdLabels[1] = pRotor->GetLabel();
	   if (pGround != NULL) {
	      	   NumNodes++;
		   NdTyps[2] = pGround->GetNodeType();
		   NdLabels[2] = pGround->GetLabel();
	   }

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
	   const Mat3x3& rrot,
	   const StructNode* pRotor, 
	   ResForceSet **ppres, 
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
	   	const Mat3x3& rrot,
		const StructNode* pRotor, 
		const StructNode* pGround, 
		ResForceSet **ppres, 
		doublereal dOR,
		doublereal dR,
		DriveCaller *pdW,
		unsigned int iMaxIt,
		doublereal dTol,
		doublereal dE,
		doublereal dCH,
		doublereal dCFF,
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
		const Mat3x3& rrot,
		const StructNode* pRotor,
		const StructNode* pGround, 
		ResForceSet **ppres, 
		doublereal dOR,
		doublereal dR,
		DriveCaller *pdW,
		unsigned int iMaxIt,
		doublereal dTol,
		doublereal dE,
		doublereal dCH,
		doublereal dCFF,
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
		const Mat3x3& rrot,
		const StructNode* pRotor,
		const StructNode* pGround, 
		ResForceSet **ppres,
		doublereal dOR,
		doublereal dR,
		DriveCaller *pdW,
		unsigned int iMaxIt,
		doublereal dTol,
		doublereal dE,
		doublereal dCH,
		doublereal dCFF,
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

/*
 * Based on the 3 state dynamic inflow by Pitt-Peters:
 * 
 * D. M. Pitt,  D. A. Peters,
 * "Theoretical Prediction of Dynamic Inflow Derivatives",
 * Vertica, Vol. 5, pp.21-34, 1981
 *
 * as discussed by Chen in:
 *
 * R. T. N. Chen,
 * "A Survey of Nonuniform Inflow Models for Rotorcraft
 * Flight Dynamics and Control Applications"
 * Vertica, Vol 14, No. 2, pp.147-184, 1990
 */

class DynamicInflowRotor : virtual public Elem, public Rotor {
 protected:
   doublereal dVConst;
   doublereal dVSine;
   doublereal dVCosine;

   doublereal dL11;
   doublereal dL13;
   doublereal dL22;
   doublereal dL31;
   doublereal dL33;
   
 public:
   DynamicInflowRotor(unsigned int uLabel,
		      const DofOwner* pDO,
		      const StructNode* pCraft, 
		      const Mat3x3& rrot,
		      const StructNode* pRotor,
      		      const StructNode* pGround, 
		      ResForceSet **ppres, 
		      doublereal dOR,
		      doublereal dR,
		      unsigned int iMaxIt,
		      doublereal dTol,
		      doublereal dE,
		      doublereal dCH,
		      doublereal dCFF,
		      doublereal dVConstTmp,
		      doublereal dVSineTmp,
		      doublereal dVCosineTmp,
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
   virtual void SetValue(DataManager *pDM,
		   VectorHandler& X, VectorHandler& XP,
		   SimulationEntity::Hints *ph = 0);

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

#endif /* ROTOR_H */

