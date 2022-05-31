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

/* nodi */

#ifndef NODE_H
#define NODE_H

#include "myassert.h"

#include "units.h"
#include "output.h"
#include "withlab.h"
#include "dofown.h"
#include "simentity.h"

/*
 * Array dei nomi dei nodi.
 * Usato per output
 * @see Node::Type
 */
extern const char* psNodeNames[];

/*
 * Array delle stringhe di identificazione dei tipi di nodi.
 * Usato per input di controllo
 * @see Node::Type
 */
extern const char* psReadControlNodes[];

/*
 * Array delle stringhe di identificazione dei tipi di nodi.
 * Usato per input dei nodi
 * @see Node::Type
 */
extern const char* psReadNodesNodes[];

/* Node - begin */ 

class Node : public WithLabel, public SimulationEntity,
public DofOwnerOwner, public ToBeOutput
{
public:
	/* Enumerazione dei tipi di nodi */
	enum Type {
		UNKNOWN = -1,

		/* Should be Node::SCALAR; keep using Node::ABSTRACT
		 * for backward compatibility */
		ABSTRACT = 0,

		STRUCTURAL,
		ELECTRIC,
		THERMAL,
		PARAMETER,
		HYDRAULIC,

		LASTNODETYPE
	};

public:
	/* Errori: */
	class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
 	};

	/* Costruttori */

	/*
	 * Costruttore.
	 * @param uL label
	 * @param pDO puntatore al DofOwner relativo, gestito da DofOwnerOwner
	 * @param fOut flag di output, gestito da ToBeOutput
	 * @see DofOwner
	 * @see DofOwnerOwner
	 * @see ToBeOutput
	 */
	Node(unsigned int uL, const DofOwner* pDO, flag fOut);

	/* Distruttore */
	virtual ~Node(void);

	/* Funzioni di servizio */
	const Node *GetNode(void) const { return this; };

	/* Tipo del nodo (usato per debug ecc.) */
	virtual Node::Type GetNodeType(void) const = 0;

	/* Contributo del nodo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const = 0;

	/*
	 * Describe the degrees of freedom
	 */
	virtual std::ostream& DescribeDof(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void DescribeDof(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;

	virtual std::ostream& DescribeEq(std::ostream& out,
		const char *prefix = "",
		bool bInitial = false) const;
	virtual void DescribeEq(std::vector<std::string>& desc,
		bool bInitial = false,
		int i = -1) const;

	/* Metodi che operano sui DoF */

	/*
	 * Ritorna il primo indice di riga dei DoF del nodo, in base 0.
	 * Ovvero, l'indice del primo DoF del nodo in un vettore a base zero.
	 * Per avere gli indici in un vettore a base 1 (stile Fortran),
	 * occorre sommare al risultato il numero del DoF.
	 * Vedi iGetFirstColIndex()
	 */
	virtual integer iGetFirstRowIndex(void) const;

	/*
	 * Ritorna gli indici di colonna dei DoF.
	 * Per la numerazione degli indici vedi iGetFirstRowIndex().
	 * Tipicamente gli indici di riga e di colonna sono gli stessi,
	 * tranne in alcuni casi notevoli.
	 */
	virtual integer iGetFirstColIndex(void) const;

	/*
	 * Restituisce il valore del DoF iDof.
	 * Se il nodo e' differenziale, iOrder puo' essere = 1 per avere la derivata
	 */
	virtual const doublereal& dGetDofValue(int iDof, int iOrder = 0) const = 0;

	/*
	 * Restituisce il valore del DoF iDof al passo precedente.
	 * Se il nodo e' differenziale, iOrder puo' essere = 1 per avere la derivata
	 */
	virtual const doublereal& dGetDofValuePrev(int iDof, int iOrder = 0) const = 0;

	/*
	 * Setta il valore del DoF iDof a dValue.
	 * Se il nodo e' differenziale, iOrder puo' essere = 1
	 * per operare sulla derivata
	 */
	virtual void SetDofValue(const doublereal& dValue,
		unsigned int iDof,
		unsigned int iOrder = 0) = 0;

	/*
	 * priv data and dofs coincide for nodes, unless overridden */

	/*
	 * Metodi per l'estrazione di dati "privati".
	 * Si suppone che l'estrattore li sappia interpretare.
	 * Come default non ci sono dati privati estraibili
	 */
	virtual unsigned int iGetNumPrivData(void) const {
		return iGetNumDof();
	};

	/*
	 * Maps a string (possibly with substrings) to a private data;
	 * returns a valid index ( > 0 && <= iGetNumPrivData()) or 0 
	 * in case of unrecognized data; error must be handled by caller
	 */
	virtual unsigned int iGetPrivDataIdx(const char *s) const {
		return 0;
	};

	/*
	 * Returns the current value of a private data
	 * with 0 < i <= iGetNumPrivData()
	 */
	virtual doublereal dGetPrivData(unsigned int i) const {
		return dGetDofValue(i);
	};

        /* 
         * Automatic differentiation support:
         * Called on each node before assembly of the Jacobian 
         */
        virtual void UpdateJac(doublereal dCoef);
        /* 
         * Automatic differentiation support:
         * Called on each node before assembly of the Jacobian vector product Jac * Y 
         */
     virtual void UpdateJac(const VectorHandler& Y, doublereal dCoef);
};

Node::Type str2nodetype(const char *const s);

/* Node - end */

/* ScalarNode - begin */

class ScalarNode : public Node {
protected:
	/* scrive l'output */
	virtual std::ostream& Output(std::ostream& out) const;

public:
	/* Costruttori */

	/* Costruttore */
	ScalarNode(unsigned int uL, const DofOwner* pDO, flag fOut);

	/* Distruttore */
	virtual ~ScalarNode(void);

	/* Funzioni di servizio */

	/* Scrive l'output come abstract */
	virtual void Output(OutputHandler& OH) const;

	/*
	 * Ritorna il numero di DoFs ( == 1).
	 * Non usa il DofOwner in quanto viene usata per generale il DofOwner stesso
	 * (per compatibilita' con gli elementi che generano gradi di
	 * liberta' ed in previsione di nodi con un numero variabile di DoF
	 */
	virtual unsigned int iGetNumDof(void) const;

	/* Metodi che operano sui valori del DoF.
	 * Funzioni che consentono l'accesso diretto ai dati privati.
	 * sono state definite perche' i nodi scalari sono usati nei
	 * modi piu' strani e quindi puo' essere necessario l'accesso diretto.
	 */

	/* Setta il valore del DoF */
	virtual void SetX(const doublereal& d) = 0;

	/* Ottiene il valore del DoF */
	virtual const doublereal& dGetX(void) const = 0;

	/*
	 * Setta il valore della derivata.
	 * Definito solo per nodi differenziali
	 */
	virtual void SetXPrime(const doublereal& d) = 0;

	/*
	 * Ottiene il valore della derivata.
	 * Definito solo per nodi differenziali
	 */
	virtual const doublereal& dGetXPrime(void) const = 0;

	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);

	 /* returns the dimension of the component */
	const virtual MBUnits::Dimensions GetEquationDimension(integer index) const;

	/* describes the dimension of components of equation */
    virtual std::ostream& DescribeEq(std::ostream& out,
		  const char *prefix = "",
		  bool bInitial = false) const;
};

/* ScalarNode - end */


/* ScalarDifferentialNode - begin */

class ScalarDifferentialNode : virtual public ScalarNode {
protected:
	/* Valore del DoF */
	mutable doublereal dX;
	/* Valore della derivata del DoF */
	mutable doublereal dXP;

	/* Valore del DoF al passo precedente */
	doublereal dXPrev;
	/* Valore della derivata del DoF al passo precedente */
	doublereal dXPPrev;

	/* scrive l'output */
	virtual std::ostream& Output(std::ostream& out) const;

public:
	/* Costruttori */

	/*
	 * Costruttore.
	 * @param uL label
	 * @param pDO DofOwner
	 * @param fOut flag di output
	 * @param dx valore iniziale
	 * @param dxp valore iniziale della derivata
	 */
	ScalarDifferentialNode(unsigned int uL, const DofOwner* pDO,
		const doublereal& dx, const doublereal& dxp, flag fOut);
	/* Distruttore */
	virtual ~ScalarDifferentialNode(void);


	/* Metodi di servizio */

	/* Tipo di nodo */
	virtual Node::Type GetNodeType(void) const;

	/*
	 * Esegue operazioni sui DoF di proprieta' dell'elemento.
	 * In particolare ritorna il tipo di DoF in base all'indice i.
	 */
	virtual DofOrder::Order GetDofType(unsigned int i) const;

	/* Metodi sui DoF */

	/*
	 * Restituisce il valore del DoF iDof.
	 * Se differenziale, iOrder puo' essere = 1 per ottenere la derivata
	 */
	virtual const doublereal& dGetDofValue(int iDof, int iOrder = 0) const;

	/*
	 * Restituisce il valore del DoF iDof al passo precedente.
	 * Se differenziale, iOrder puo' essere = 1 per ottenere la derivata
	 */
	virtual const doublereal& dGetDofValuePrev(int iDof, int iOrder = 0) const;

	/*
	 * Setta il valore del DoF iDof a dValue.
	 * Se differenziale, iOrder puo' essere = 1 per ottenere la derivata
	 */
	virtual void SetDofValue(const doublereal& dValue,
		unsigned int iDof,
		unsigned int iOrder = 0);

	/*
	 * Funzione che consente l'accesso diretto ai dati privati.
	 * Sono state definite perche' i nodi scalari sono usati nei
	 * modi piu' strani e quindi puo' essere necessario l'accesso diretto
	 */
	virtual void SetX(const doublereal& d);

	/* Ottiene il valore del DoF. Vedi SetX() */
	virtual inline const doublereal& dGetX(void) const;

	/* Setta la derivata del DoF. Vedi SetX() */
	virtual void SetXPrime(const doublereal& d);

	/* Ottiene la derivata del DoF. Vedi GetX() */
	virtual inline const doublereal& dGetXPrime(void) const;

	/* Consente di settare il valore iniziale nel vettore della soluzione */
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	/* Aggiorna i valori interni */
	virtual void Update(const class VectorHandler&, const class VectorHandler&);

	/* restart */
	std::ostream& Restart(std::ostream& out) const;

	/*
	 * Metodi per l'estrazione di dati "privati".
	 * Si suppone che l'estrattore li sappia interpretare.
	 * Come default non ci sono dati privati estraibili
	 */
	virtual unsigned int iGetNumPrivData(void) const;

	/*
	 * Maps a string (possibly with substrings) to a private data;
	 * returns a valid index ( > 0 && <= iGetNumPrivData()) or 0 
	 * in case of unrecognized data; error must be handled by caller
	 */
	virtual unsigned int iGetPrivDataIdx(const char *s) const;

	/*
	 * Returns the current value of a private data
	 * with 0 < i <= iGetNumPrivData()
	 */
	virtual doublereal dGetPrivData(unsigned int i) const;

	 /* returns the dimension of the component */
	const virtual MBUnits::Dimensions GetEquationDimension(integer index) const;

	/* describes the dimension of components of equation */
    virtual std::ostream& DescribeEq(std::ostream& out,
		  const char *prefix = "",
		  bool bInitial = false) const;
};

inline const doublereal&
ScalarDifferentialNode::dGetX(void) const
{
	return dX;
}


inline const doublereal&
ScalarDifferentialNode::dGetXPrime(void) const
{
	return dXP;
}
/* ScalarDifferentialNode - end */

/* ScalarAlgebraicNode - begin */

class ScalarAlgebraicNode: virtual public ScalarNode {
protected:
	/* Valore del DoF */
	mutable doublereal dX;

	/* Valore del DoF al passo precedente */
	doublereal dXPrev;

	/* scrive l'output */
	virtual std::ostream& Output(std::ostream& out) const;

public:
	/* Costruttori */

	/* Costruttore */
	ScalarAlgebraicNode(unsigned int uL, const DofOwner* pDO,
		doublereal dx, flag fOut);

	/* Distruttore */
	virtual ~ScalarAlgebraicNode(void);

	/* Metodi di servizio */

	/* Tipo di nodo */
	virtual Node::Type GetNodeType(void) const;

	/*
	 * Esegue operazioni sui DoF di proprieta' dell'elemento.
	 * In particolare ritorna il tipo di DoF in base all'indice i.
	 */
	virtual DofOrder::Order GetDofType(unsigned int i) const;

	/* Metodi che operano sul valore del DoF */

	/*
	 * Restituisce il valore del DoF iDof.
	 * Se differenziale, iOrder puo' essere = 1 per ottenere la derivata
	 */
	virtual const doublereal& dGetDofValue(int iDof, int iOrder = 0) const;

	/*
	 * Restituisce il valore del DoF iDof al passo precedente.
	 * Se differenziale, iOrder puo' essere = 1 per ottenere la derivata
	 */
	virtual const doublereal& dGetDofValuePrev(int iDof, int iOrder = 0) const;

	/*
	 * Setta il valore del DoF iDof a dValue.
	 * Se differenziale, iOrder puo' essere = 1 per operare sulla derivata
	 */
	virtual void SetDofValue(const doublereal& dValue,
		unsigned int iDof,
		unsigned int iOrder = 0);

	/*
	 * Funzione che consente l'accesso diretto ai dati privati.
	 * Sono state definite perche' i nodi astratti sono usati nei
	 * modi piu' strani e quindi puo' essere necessario l'accesso diretto
	 */
	virtual void SetX(const doublereal& d);

	/* Ottiene il valore del DoF. Vedi SetX() */
	virtual inline const doublereal& dGetX(void) const;

	/* Non definito per nodi algebrici */
	virtual void SetXPrime(const doublereal& d);

	/* Non definito per nodi algebrici */
	virtual inline const doublereal& dGetXPrime(void) const;

	/* Consente di settare il valore iniziale nel vettore della soluzione */
	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	/* Aggiorna i valori interni */
	virtual void Update(const class VectorHandler&, const class VectorHandler&);

	/* restart */
	std::ostream& Restart(std::ostream& out) const;

	/*
	 * Maps a string (possibly with substrings) to a private data;
	 * returns a valid index ( > 0 && <= iGetNumPrivData()) or 0 
	 * in case of unrecognized data; error must be handled by caller
	 */
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
};

inline const doublereal&
ScalarAlgebraicNode::dGetX(void) const
{
	return dX;
}

inline const doublereal&
ScalarAlgebraicNode::dGetXPrime(void) const
{
	DEBUGCERR("Error, getting derivative from algebraic dof!" << std::endl);
	throw Node::ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/* ScalarAlgebraicNode - end */

/* ParameterNode - begin */

/* Parametri.
 * I nodi di tipo parametro sono derivati dai nodi scalari algebrici,
 * ma non sono veri nodi. In realta' sono entita' che possiedono un valore,
 * ma non generano DoFs ed equazioni. Sono usati per consentire di dare in
 * modo trasparente un valore in ingresso, sotto forma di nodo, a tutti quegli
 * elementi elettrici e Genel che normalmente usano un DoF scalare senza farlo
 * partecipare allo jacobiano.
 */
class ParameterNode : public ScalarAlgebraicNode {
public:
	/* Costruttori */

	/* Costruttore */
	ParameterNode(unsigned int uL, const DofOwner* pDO,
		doublereal dx, flag fOut);
	/* Distruttore */
	virtual ~ParameterNode(void);

	/* Metodi di servizio */

	/* Tipo del nodo. Usato solo per debug ecc. */
	virtual Node::Type GetNodeType(void) const;

	/*
	 * Ritorna il numero di dofs.
	 * non usa il DofOwner in quanto viene usato per generale il DofOwner stesso.
	 * Ritorna 0 perche' il parametro non ha DoFs
	 */
	virtual unsigned int iGetNumDof(void) const;
	virtual DofOrder::Order GetDofType(unsigned int i) const;

	/* Metodi che agiscono sul valore */

	/* Metodi relativi al metodo di intergazione */

	/* Output di default per nodi di cui non si desidera output */
	virtual void Output(OutputHandler& OH) const;

	/* Inizializzazione del valore */
	void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0);

	/* Aggiorna dati in base alla soluzione */
	virtual void Update(const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Inverse Dynamics */
	virtual void Update(const VectorHandler& XCurr, const int iOrder);

	/* Elaborazione dati dopo la predizione */
	virtual void AfterPredict(VectorHandler& X,
		VectorHandler& XP);

	/* the value! */
	virtual unsigned int iGetNumPrivData(void) const {
		return 1;
	};
};

/* ParameterNode - end */


/* Node2Scalar - begin */

/*
 * Struttura di conversione da nodo generico a nodo scalare.
 * Questa struttura consente di usare un grado di liberta'
 * di un nodo generico come se fosse un nodo scalare
 */
struct NodeDof {
	/* DoF del nodo */
	int iDofNumber;     /* Dof of the node */
	/* Puntatore al nodo */
	Node* pNode;        /* Pointer to the node */

	/* Costruttori */

	/* Costruttore di default */
	NodeDof(void);
	/* Costruttore */
	NodeDof(int id, Node* p);
	/* Distruttore */
	virtual ~NodeDof(void);
};

/*
 * Classe di conversione da nodo generico a nodo scalare.
 * @see NodeDof
 */
class Node2Scalar : public ScalarNode {
protected:
	/* Struttura che punta ad un DoF di un nodo */
	NodeDof ND;

public:
	/* Costruttori */

	/* Costruttore */
	Node2Scalar(const NodeDof& nd);
	/* Distruttore */
	virtual ~Node2Scalar(void);

	/* Metodi di servizio */
	const Node *GetNode(void) const { return ND.pNode; };

	/* Tipo del nodo. Uusato per debug ecc. */
	virtual Node::Type GetNodeType(void) const;

	/* Contributo del nodo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/*
	 * Ritorna il numero di dofs.
	 * Non usa il DofOwner in quanto viene usata per generare il DofOwner stesso
	 */
	virtual unsigned int iGetNumDof(void) const;

	/* Metodi che operano sui valori del DoF */

	/*
	 * Esegue operazioni sui DoF di proprieta' dell'elemento.
	 * In particolare ritorna il tipo di DoF in base all'indice i.
	 */
	virtual DofOrder::Order GetDofType(unsigned int i) const;

	/*
	 * Ritorna gli indici di riga.
	 * Tipicamente sono gli stessi di quelli di colonna
	 */
	virtual integer iGetFirstRowIndex(void) const;

	/*
	 * Ritorna gli indici di colonna.
	 * Tipicamente sono gli stessi di quelli di riga.
	 * @see iGetFirstRowIndex()
	 */
	virtual integer iGetFirstColIndex(void) const;

	/*
	 * Restituisce il valore del DoF iDof.
	 * Se differenziale, iOrder puo' essere = 1 per ottenere la derivata
	 */
	virtual const doublereal& dGetDofValue(int iDof, int iOrder = 0) const;

	/*
	 * Restituisce il valore del DoF iDof al passo precedente.
	 * Se differenziale, iOrder puo' essere = 1 per ottenere la derivata
	 */
	virtual const doublereal& dGetDofValuePrev(int iDof, int iOrder = 0) const;

	/*
	 * Setta il valore del DoF iDof a dValue.
	 * Se differenziale, iOrder puo' essere = 1 per operare sulla derivata
	 */
	virtual void SetDofValue(const doublereal& dValue,
		unsigned int iDof,
		unsigned int iOrder = 0);

	/*
	 * Funzione che consente l'accesso diretto ai dati privati.
	 * Sono state definite perche' i nodi astratti sono usati nei
	 * modi piu' strani e quindi puo' essere necessario l'accesso diretto
	 */
	virtual void SetX(const doublereal& d);

	/* Ottiene il valore del DoF */
	virtual inline const doublereal& dGetX(void) const;

	/* Setta il valore della derivata del DoF */
	virtual void SetXPrime(const doublereal& d);

	/* Setta il valore della derivata del DoF */
	virtual inline const doublereal& dGetXPrime(void) const;
};

inline const doublereal&
Node2Scalar::dGetX(void) const
{
	return dGetDofValue(1, 0);
}

inline const doublereal&
Node2Scalar::dGetXPrime(void) const
{
	return dGetDofValue(1, 1);
}

/* Node2Scalar - end */


/* ScalarDof - begin */

/*
 * Struttura che trasforma un nodo scalare in un grado di liberta' scalare.
 * In pratica consente di accedere ad un DoF scalare o alla derivata di un
 * nodo scalare in modo trasparente
 */
struct ScalarDof {
	/* Puntatore al nodo scalare */
	ScalarNode* pNode;
	/* Ordine del grado di liberta' */
	int iOrder;
	int iIndex;
	/* Costruttori */

	/* Costruttore di default */
	ScalarDof(void);
	ScalarDof(const ScalarDof& sd);
	/* Costruttore */
	ScalarDof(ScalarNode* p, int iO, int iI);
	/* Distruttore */
	~ScalarDof(void);

	/* Funzioni che operano sui valori del DoF */

	/* Ottiene il valore del DoF */
	const doublereal & dGetValue(void) const;

	/* Ottiene il valore del DoF al passo precedente */
	const doublereal & dGetValuePrev(void) const;

	/* Scrive nel file di Restart il contributo per la chiamata allo ScalarDof */
	std::ostream& RestartScalarDofCaller(std::ostream& out) const ;
};

/* ScalarDof - end */

#endif /* NODE_H */

