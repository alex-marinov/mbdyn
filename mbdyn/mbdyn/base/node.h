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

/* nodi */


#ifndef NODE_H
#define NODE_H


#include <myassert.h>

/* include del programma */
#include <output.h>
#include <withlab.h>
#include <dofown.h>

/** 
 Array dei nomi dei nodi. 
 Usato per output
 @see Node::Type
 */
extern const char* psNodeNames[];
/** 
 Array delle stringhe di identificazione dei tipi di nodi.
 Usato per input di controllo
 @see Node::Type
 */
extern const char* psReadControlNodes[];
/** 
 Array delle stringhe di identificazione dei tipi di nodi.
 Usato per input dei nodi
 @see Node::Type
 */
extern const char* psReadNodesNodes[];


/* Node - begin */

/// Nodi
class Node : public WithLabel, public DofOwnerOwner, public ToBeOutput {
 public:
   /** Enumerazione dei tipi di nodi */
   enum Type {
      UNKNOWN = -1,

	ABSTRACT = 0,
	STRUCTURAL,
	ELECTRIC,
	PARAMETER,
	HYDRAULIC,

	LASTNODETYPE
   };

 public:
   /* Errori: */
   class ErrGeneric {};
   
   /**@name Costruttori */
   //@{
   
   /** 
    Costruttore.
    @param uL label
    @param pDO puntatore al DofOwner relativo, gestito da DofOwnerOwner
    @param fOut flag di output, gestito da ToBeOutput
    @see DofOwner
    @see DofOwnerOwner
    @see ToBeOutput
    */
   Node(unsigned int uL, const DofOwner* pDO, flag fOut);
      
   /** Distruttore */
   virtual ~Node(void);
   //@}
      
   /**@name Funzioni di servizio */
   //@{

   /** Tipo del nodo (usato per debug ecc.) */
   virtual Node::Type GetNodeType(void) const = 0;

   /** Contributo del nodo al file di restart */
   virtual ostream& Restart(ostream& out) const = 0;
   //@}

   /**@name Metodi che operano sui DoF */
   //@{
      
   /** 
    Ritorna il numero di DoFs.
    Non usa il DofOwner in quanto viene usata per generale il DofOwner
    stesso (per compatibilita' con gli elementi che generano gradi di 
    liberta' ed in previsione di nodi con un numero variabile di DoF) 
    */
   virtual inline unsigned int iGetNumDof(void) const = 0;
    
   /**
    Test di validita' di un indice. 
    Nota: gli indici vanno da 1 a iGetNumDofs()
    */
   virtual flag fIsValidIndex(unsigned int i) const;
   
   /**
    Esegue operazioni sui DoF di proprieta' dell'elemento.
    In particolare ritorna il tipo di DoF in base all'indice i.
    Di default i DoF dei nodi sono assunti differenziali.
    Il tipo e' preso dall'enum DofOrder.
    Nota: gli indici sono in base 0, ovvero deve essere
    0 < i < iGetNumDof()
    @see DofOrder
    */   
   virtual DofOrder::Order SetDof(unsigned int i) const;

   /** 
    Ritorna il primo indice di riga dei DoF del nodo, in base 0.
    Ovvero, l'indice del primo DoF del nodo in un vettore a base zero.
    Per avere gli indici in un vettore a base 1 (stile Fortran),
    occorre sommare al risultato il numero del DoF.
    Vedi iGetFirstColIndex()
    */
   virtual integer iGetFirstRowIndex(void) const;
      
   /** 
    Ritorna gli indici di colonna dei DoF.
    Per la numerazione degli indici vedi iGetFirstRowIndex().
    Tipicamente gli indici di riga e di colonna sono gli stessi,
    tranne in alcuni casi notevoli.
    */
   virtual integer iGetFirstColIndex(void) const;
   
   /** 
    Restituisce il valore del DoF iDof.
    Se il nodo e' differenziale, iOrder puo' essere = 1 per avere la derivata
    */
   virtual const doublereal& dGetDofValue(int iDof, int iOrder = 0) const = 0;
   
   /**
    Setta il valore del DoF iDof a dValue.
    Se il nodo e' differenziale, iOrder puo' essere = 1 
    per operare sulla derivata
    */
   virtual void SetDofValue(const doublereal& dValue,
			    unsigned int iDof, 
			    unsigned int iOrder = 0) = 0;
   //@}
   
      
   /**@name Metodi legati all'integrazione */
   //@{
   /**
    Output.  
    Di default non fa nulla per nodi che non generano output
    */
   virtual void Output(OutputHandler& OH) const;
   virtual void Output(OutputHandler& OH,
		   const VectorHandler& X, const VectorHandler& XP) const;
   
   /**
    Setta i valori iniziali dei DoF.
    Puo' essere usata per altre inizializzazioni prima di 
    iniziare l'integrazione 
    */
   virtual void SetValue(VectorHandler& X, VectorHandler& XP) const;
            
   /**
    Elaborazione vettori e dati prima della predizione.
    Per MultiStepIntegrator
    */
   virtual void BeforePredict(VectorHandler& /* X */ ,
      			      VectorHandler& /* XP */ ,
			      VectorHandler& /* XPrev */ ,
			      VectorHandler& /* XPPrev */ ) const;
   
   /**
    Elaborazione vettori e dati dopo la predizione.
    Per MultiStepIntegrator
    */
   virtual void AfterPredict(VectorHandler& X, VectorHandler& XP) = 0;

   /**
    Aggiorna dati in base alla soluzione. 
    Usata per operazioni aggiuntive al semplice aggiornamento additivo,
    effettuato gia' dall'integratore.
    */
   virtual void Update(const VectorHandler& XCurr, 
		       const VectorHandler& XPrimeCurr) = 0;
   //@}
};

/* Node - end */


/* ScalarNode - begin */

/// Nodo Scalare
class ScalarNode : public Node {
 public:
   /**@name Costruttori */
   
   //@{
   /** Costruttore */
   ScalarNode(unsigned int uL, const DofOwner* pDO, flag fOut);
      
   /** Distruttore */
   virtual ~ScalarNode(void);
   //@}

   /**@name Funzioni di servizio */
   //@{
   /**
    Ritorna il numero di DoFs ( == 1).
    Non usa il DofOwner in quanto viene usata per generale il DofOwner stesso 
    (per compatibilita' con gli elementi che generano gradi di 
    liberta' ed in previsione di nodi con un numero variabile di DoF
    */
   virtual inline unsigned int iGetNumDof(void) const;
   //@}

   /**@name Metodi che operano sui valori del DoF.
    Funzioni che consentono l'accesso diretto ai dati privati.
    sono state definite perche' i nodi scalari sono usati nei
    modi piu' strani e quindi puo' essere necessario l'accesso diretto.
    */
   //@{
   /** Setta il valore del DoF */
   virtual void SetX(const doublereal& d) = 0;
      
   /** Ottiene il valore del DoF */
   virtual inline const doublereal& dGetX(void) const = 0;
   
   /** 
    Setta il valore della derivata.
    Definito solo per nodi differenziali 
    */
   virtual void SetXPrime(const doublereal& d) = 0;
      
   /** 
    Ottiene il valore della derivata.
    Definito solo per nodi differenziali 
    */
   virtual inline const doublereal& dGetXPrime(void) const = 0;
   //@}
};


inline unsigned int ScalarNode::iGetNumDof(void) const
{
   return 1;
}

/* ScalarNode - end */


/* ScalarDifferentialNode - begin */

/// Nodo scalare differenziale
class ScalarDifferentialNode : public ScalarNode {
 protected:
   /** Valode del DoF */
   doublereal dX;
   /** Valore della derivata del DoF */
   doublereal dXP;
   
 public:
   /**@name Costruttori */
   //@{
   /** 
    Costruttore.
    @param uL label
    @param pDO DofOwner
    @param fOut flag di output
    @param dx valore iniziale
    @param dxp valore iniziale della derivata
    */
   ScalarDifferentialNode(unsigned int uL, const DofOwner* pDO, 
			  const doublereal& dx, const doublereal& dxp, 
			  flag fOut);
   /** Distruttore */
   virtual ~ScalarDifferentialNode(void);
   //@}

   /**@name Metodi di servizio */
   //@{
   /** 
    Esegue operazioni sui DoF di proprieta' dell'elemento.
    In particolare ritorna il tipo di DoF in base all'indice i.
    */
   virtual DofOrder::Order SetDof(unsigned int i) const;
   //@}
   
   /** Metodi sui DoF */
   //@{
      
   /** 
    Restituisce il valore del DoF iDof.
    Se differenziale, iOrder puo' essere = 1 per ottenere la derivata
    */
   virtual const doublereal& dGetDofValue(int iDof, int iOrder = 0) const;

   /**
    Setta il valore del DoF iDof a dValue.
    Se differenziale, iOrder puo' essere = 1 per ottenere la derivata
    */
   virtual void SetDofValue(const doublereal& dValue, 
			    unsigned int iDof, 
			    unsigned int iOrder = 0);
   
   
   /**
    Funzione che consente l'accesso diretto ai dati privati.
    Sono state definite perche' i nodi scalari sono usati nei
    modi piu' strani e quindi puo' essere necessario l'accesso diretto 
    */
   virtual void SetX(const doublereal& d);
   
   /** Ottiene il valore del DoF. Vedi SetX() */
   virtual inline const doublereal& dGetX(void) const;
   
   /** Setta la derivata del DoF. Vedi SetX() */
   virtual void SetXPrime(const doublereal& d);
   
   /** Ottiene la derivata del DoF. Vedi GetX() */
   virtual inline const doublereal& dGetXPrime(void) const;

   /** Consente di settare il valore iniziale nel vettore della soluzione*/
   virtual void SetValue(VectorHandler& X, VectorHandler& XP) const;
      
   //@}
};


inline const doublereal& ScalarDifferentialNode::dGetX(void) const
{
   return dX;
}


inline const doublereal& ScalarDifferentialNode::dGetXPrime(void) const
{
   return dXP;
}

/* ScalarDifferentialNode - end */

/* ScalarAlgebraicNode - begin */

/// Nodo scalare algebrico. Non ha derivata.
class ScalarAlgebraicNode : public ScalarNode {
 protected:
   /** Valore del DoF */
   doublereal dX;
   
 public:
   /**@name Costruttori */
   //@{
   
   /** Costruttore */
   ScalarAlgebraicNode(unsigned int uL, const DofOwner* pDO, 
		       doublereal dx, flag fOut);
   /** Distruttore */
   virtual ~ScalarAlgebraicNode(void);
   //@}

   /**@name Metodi di servizio */
   //@{
      
   /**
    Esegue operazioni sui DoF di proprieta' dell'elemento.
    In particolare ritorna il tipo di DoF in base all'indice i. 
    */
   virtual DofOrder::Order SetDof(unsigned int i) const;   
   //@}
   
   
   /**@name Metodi che operano sul valore del DoF */
   //@{
        
   /** 
    Restituisce il valore del DoF iDof.
    Se differenziale, iOrder puo' essere = 1 per ottenere la derivata 
    */
   virtual const doublereal& dGetDofValue(int iDof, int iOrder = 0) const;
   
   /**
    Setta il valore del DoF iDof a dValue.
    Se differenziale, iOrder puo' essere = 1 per operare sulla derivata 
    */
   virtual void SetDofValue(const doublereal& dValue, 
			    unsigned int iDof, 
			    unsigned int iOrder = 0);
      
   /** 
    Funzione che consente l'accesso diretto ai dati privati.
    Sono state definite perche' i nodi astratti sono usati nei
    modi piu' strani e quindi puo' essere necessario l'accesso diretto
    */
   virtual void SetX(const doublereal& d);
   
   /** Ottiene il valore del DoF. Vedi SetX() */
   virtual inline const doublereal& dGetX(void) const;
   
   /** Non definito per nodi algebrici */
   virtual void SetXPrime(const doublereal& d);
   
   /** Non definito per nodi algebrici */
   virtual inline const doublereal& dGetXPrime(void) const;
      
   /** Consente di settare il valore iniziale nel vettore della soluzione*/
   virtual void SetValue(VectorHandler& X, VectorHandler& XP) const;
      
   //@}
};


inline const doublereal& ScalarAlgebraicNode::dGetX(void) const
{
   return dX;
}


inline const doublereal& ScalarAlgebraicNode::dGetXPrime(void) const
{
   DEBUGCERR("Error, getting derivative from algebraic dof!" << endl);
   THROW(Node::ErrGeneric());
#ifndef USE_EXCEPTIONS
   return dX; /* inutile, ma evita warnings se non si usano eccezioni */
#endif /* USE_EXCEPTIONS */
}

/* ScalarAlgebraicNode - end */


/* ParameterNode - begin */

/** Parametri.
 I nodi di tipo parametro sono derivati dai nodi scalari algebrici,
 ma non sono veri nodi. In realta' sono entita' che possiedono un valore,
 ma non generano DoFs ed equazioni. Sono usati per consentire di dare in
 modo trasparente un valore in ingresso, sotto forma di nodo, a tutti quegli
 elementi elettrici e Genel che normalmente usano un DoF scalare senza farlo
 partecipare allo jacobiano.
 */
class ParameterNode : public ScalarAlgebraicNode {   
 public:
   /**@name Costruttori */
   //@{
   /** Costruttore */
   ParameterNode(unsigned int uL, const DofOwner* pDO,
		 doublereal dx, flag fOut);
   /** Distruttore */
   virtual ~ParameterNode(void);
   //@}
     
   /**@name Metodi di servizio */
   //@{
   /** Tipo del nodo. Usato solo per debug ecc. */
   virtual Node::Type GetNodeType(void) const;
   
   /** Contributo del nodo al file di restart */
   virtual ostream& Restart(ostream& out) const;
         
   /** 
    Ritorna il numero di dofs.
    non usa il DofOwner in quanto viene usato per generale il DofOwner stesso.
    Ritorna 0 perche' il parametro non ha DoFs
    */
   virtual inline unsigned int iGetNumDof(void) const;
      
   /** 
    Verifica di validita' di un indice.
    Deve essere 0 perche' il parametro non ha DoFs 
    */
   virtual flag fIsValidIndex(unsigned int i) const;
   //@}
      
   /**@name Metodi che agiscono sul valore */
   //@{
   /** 
    Restituisce il valore del DoF iDof.
    Se differenziale, iOrder puo' essere = 1 per la derivata.
    Il parametro e' algebrico.
    */
   virtual const doublereal& dGetDofValue(int iDof, int iOrder = 0) const;
   
   /**
    Setta il valore del DoF iDof a dValue.
    Se differenziale, iOrder puo' essere = 1 per la derivata.
    Il parametro e' algebrico.
    */
   virtual void SetDofValue(const doublereal& dValue, 
			    unsigned int iDof, 
			    unsigned int iOrder = 0);
   //@}
      
   /**@name Metodi relativi al metodo di intergazione */
   //@{
   /** Output di default per nodi di cui non si desidera output */
   virtual void Output(OutputHandler& OH) const;
   
   /** Inizializzazione del valore */
   void SetValue(VectorHandler& X, VectorHandler& XP) const;
   
   /** Aggiorna dati in base alla soluzione */
   virtual void Update(const VectorHandler& XCurr,
		       const VectorHandler& XPrimeCurr);
   
   /** Elaborazione dati dopo la predizione */
   virtual void AfterPredict(VectorHandler& X, 
			     VectorHandler& XP);
   //@}
};


inline unsigned int ParameterNode::iGetNumDof(void) const
{
   return 0;
}

/* ParameterNode - end */


/* Node2Scalar - begin */

/**
 Struttura di conversione da nodo generico a nodo scalare.
 Questa struttura consente di usare un grado di liberta' di un nodo generico 
 come se fosse un nodo scalare
 */
struct NodeDof {
   /** Label del nodo */
   unsigned int uNode;
   /** DoF del nodo */
   int iDofNumber;     /* Dof of the node */
   /** Puntatore al nodo */
   Node* pNode;        /* Pointer to the node */

   /**@name Costruttori */
   //@{
   /** Costruttore di default */
   NodeDof(void);
   /** Costruttore */
   NodeDof(unsigned int u, int id, Node* p);
   /** Distruttore */
   virtual ~NodeDof(void);
   //@}
};

/** Classe di conversione da nodo generico a nodo scalare. 
 @see NodeDof */
class Node2Scalar : public ScalarNode {
 protected:
   /** Struttura che punta ad un DoF di un nodo */
   NodeDof ND;
   
 public:
   /**@name Costruttori */
   //@{
   /** Costruttore */
   Node2Scalar(const NodeDof& nd);
   /** Distruttore */
   virtual ~Node2Scalar(void);
   //@}

   /** Metodi di servizio */
   //@{
   /** Tipo del nodo. Uusato per debug ecc. */
   virtual Node::Type GetNodeType(void) const;

   /** Contributo del nodo al file di restart */
   virtual ostream& Restart(ostream& out) const;
   
   /** 
    Ritorna il numero di dofs.
    Non usa il DofOwner in quanto viene usata per generare il DofOwner stesso 
    */
   virtual inline unsigned int iGetNumDof(void) const;
      
   /** Verifica la validita' dell'indice i */
   virtual flag fIsValidIndex(unsigned int i) const;
   //@}
   
   /**@name Metodi che operano sui valori del DoF */
   //@{
   /** 
    Esegue operazioni sui DoF di proprieta' dell'elemento.
    In particolare ritorna il tipo di DoF in base all'indice i.
    */
   virtual DofOrder::Order SetDof(unsigned int i) const;

   /** 
    Ritorna gli indici di riga. 
    Tipicamente sono gli stessi di quelli di colonna 
    */
   virtual integer iGetFirstRowIndex(void) const;
      
   /** 
    Ritorna gli indici di colonna.
    Tipicamente sono gli stessi di quelli di riga. 
    @see iGetFirstRowIndex()
    */
   virtual integer iGetFirstColIndex(void) const;
   
   /**
    Restituisce il valore del DoF iDof.
    Se differenziale, iOrder puo' essere = 1 per ottenere la derivata 
    */
   virtual const doublereal& dGetDofValue(int iDof, int iOrder = 0) const;
   
   /** 
    Setta il valore del DoF iDof a dValue.
    Se differenziale, iOrder puo' essere = 1 per operare sulla derivata 
    */
   virtual void SetDofValue(const doublereal& dValue, 
			    unsigned int iDof,
			    unsigned int iOrder = 0);
      
   /**
    Funzione che consente l'accesso diretto ai dati privati.
    Sono state definite perche' i nodi astratti sono usati nei
    modi piu' strani e quindi puo' essere necessario l'accesso diretto 
    */
   virtual void SetX(const doublereal& d);
   
   /** Ottiene il valore del DoF */
   virtual inline const doublereal& dGetX(void) const;
   
   /** Setta il valore della derivata del DoF */
   virtual void SetXPrime(const doublereal& d);
   
   /** Setta il valore della derivata del DoF */
   virtual inline const doublereal& dGetXPrime(void) const;
   //@}
      
   /** Metodi relativi al metodo di integrazione */
   //@{
   /** Output di default per nodi di cui non si desidera output */
   virtual void Output(OutputHandler& /* OH */ ) const;
   
   /** Inizializza i DoF del nodo */
   void SetValue(VectorHandler& X, VectorHandler& XP) const;
   
   /** Aggiorna dati in base alla soluzione */
   virtual void Update(const VectorHandler& XCurr,
		       const VectorHandler& XPrimeCurr);

   /** Elaborazione vettori e dati prima della predizione */
   virtual void BeforePredict(VectorHandler& X,
      			      VectorHandler& XP,
			      VectorHandler& XPrev,
			      VectorHandler& XPPrev) const;
   
   /** Elaborazione vettori e dati dopo la predizione */
   virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
   //@}
};


inline unsigned int Node2Scalar::iGetNumDof(void) const
{
   return 1; 
}


inline const doublereal& Node2Scalar::dGetX(void) const
{
   return dGetDofValue(1, 0);
}


inline const doublereal& Node2Scalar::dGetXPrime(void) const
{
   return dGetDofValue(1, 1);
}

/* Node2Scalar - end */


/* ScalarDof - begin */

/** 
 Struttura che trasforma un nodo scalare in un grado di liberta' scalare.
 In pratica consente di accedere ad un DoF scalare o alla derivata di un
 nodo scalare in modo trasparente
 */
struct ScalarDof {
   /** Puntatore al nodo scalare */
   ScalarNode* pNode;
   /** Ordine del grado di liberta' */
   int iOrder;
   
   /**@name Costruttori */
   //@{
   /** Costruttore di default */
   ScalarDof(void);
   /** Costruttore */
   ScalarDof(ScalarNode* p, int i);
   /** Distruttore */   
   ~ScalarDof(void);
   //@}
      
   /**@name Funzioni che operano sui valori del DoF */
   //@{
   /** Ottiene il valore del DoF */
   doublereal dGetValue(void) const;
   //@}
};

/* ScalarDof - end */

#endif
