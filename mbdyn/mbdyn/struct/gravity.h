/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

/* Accelerazione di gravita'
 *
 * Elemento Gravity: contiene direzione e modulo, espresso mediante un driver,
 * dell'accelerazione di gravita'. E' un elemento unico (ne puo' essere
 * dichiarato uno solo) ed e' puntato da tutti gli elementi della classe
 * ElemGravityOwner, ovvero elementi che generano forze di inerzia
 * (per ora: Body, Beam).
 *
 * Vi e' poi la classe GravityOwner, che contiene il puntatore all'elemento
 * Gravity. Da essa e' derivata la classe ElemGravityOwner. Quando l'elemento
 * viene costruito il puntatore e' nullo. Al termine della generazione
 * degli elementi, se e' definito l'elemento Gravity, tutti gli elementi
 * ElemGravityOwner vengono inizializzati con il puntatore all'elemento
 * Gravity.
 *
 * Si e' scelta la soluzione di un elemento per contenere questi dati
 * perche' in questo modo si acquista in generalita'. Infatti e' possibile
 * dare una dinamica all'accelerazione (in vista della generalizzazione del
 * tipo di elemento) mediante l'aggiunta di gradi di liberta', ecc.
 *
 * L'accelerazione e' ottenuta mediante la chiamata della funzione propria
 * bool bGetGravity(Vec3&) da parte degli elementi ElemGravityOwner.
 * Il flag dice se e' definita l'accelerazione.
 * In caso positivo, viene copiata nel vettore passato per reference.
 */

#ifndef GRAVITY_H
#define GRAVITY_H

#include "elem.h"
#include "tpldrive.h"

/* Gravity - begin */

class Gravity : virtual public Elem {
public:
	Gravity(flag fOut);

	virtual ~Gravity(void);

	/* Tipo dell'elemento (usato solo per debug ecc.) */
	virtual Elem::Type GetElemType(void) const {
		return Elem::GRAVITY;
	};

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
	 * Non e' virtuale in quanto ritorna NULL per tutti i tipi che non hanno
	 * dof propri.
	 * Il metodo GetDofType() ritorna, per ogni dof dell'elemento, l'ordine.
	 * E' usato per completare i singoli Dof relativi all'elemento.
	 */

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

	/* assemblaggio residuo */
	virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual Vec3 GetAcceleration(const Vec3& /* X */ ) const = 0;

	virtual inline int GetNumConnectedNodes(void) const {
		return 0;
	};
};

/* Gravity - end */


/* UniformGravity - begin */

class UniformGravity : virtual public Elem, public Gravity, public TplDriveOwner<Vec3> {
protected:
	Vec3 Acc;

public:
	UniformGravity(const TplDriveCaller<Vec3>* pDC, flag fOut);

	virtual ~UniformGravity(void);

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

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
	 * Non e' virtuale in quanto ritorna NULL per tutti i tipi che non hanno
	 * dof propri.
	 * Il metodo GetDofType() ritorna, per ogni dof dell'elemento, l'ordine.
	 * E' usato per completare i singoli Dof relativi all'elemento.
	 */

	/* funzioni proprie */

	/* assemblaggio residuo */
	virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	virtual void Output(OutputHandler& OH) const;

	virtual Vec3 GetAcceleration(const Vec3& /* X */ ) const {
		return Acc;
	};
};

/* UniformGravity - end */


/* CentralGravity - begin */

class CentralGravity : virtual public Elem, public Gravity {
protected:
	Vec3 m_X0;
	doublereal m_dM;
	doublereal m_dG;
	doublereal m_dThreshold;

public:
	CentralGravity(const Vec3& X0, doublereal dM, doublereal dG, flag fOut);

	virtual ~CentralGravity(void);

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* funzioni proprie */

	virtual void Output(OutputHandler& OH) const;

	virtual Vec3 GetAcceleration(const Vec3& X) const;
};

/* CentralGravity - end */


/* GravityOwner - begin */

/* Classe base di elementi che generano forze di inerzia */

class GravityOwner {
protected:
	mutable Gravity* pGravity;

public:
	GravityOwner(void);
	virtual ~GravityOwner(void);

	void PutGravity(const Gravity* pG);
	virtual bool bGetGravity(const Vec3& X, Vec3& Acc) const;
};

/* GravityOwner - end */


/* ElemGravityOwner - begin */

class ElemGravityOwner : virtual public Elem, public GravityOwner {
	friend class NestedElem;

protected:

	/*
	 * momento statico e momento di inerzia nel sistema globale
	 */
	virtual Vec3 GetS_int(void) const {
		pedantic_cerr("ElemGravityOwner(" << GetLabel() << "): "
			"warning, using default GetS_int()" << std::endl);
		return ::Zero3;
	};

	virtual Mat3x3 GetJ_int(void) const {
		pedantic_cerr("ElemGravityOwner(" << GetLabel() << "): "
			"warning, using default GetJ_int()" << std::endl);
		return ::Zero3x3;
	};

	virtual Vec3 GetB_int(void) const {
		pedantic_cerr("ElemGravityOwner(" << GetLabel() << "): "
			"warning, using default GetB_int()" << std::endl);
		return ::Zero3;
	};

	// NOTE: gravity owners must provide the momenta moment
	// with respect to the origin of the global reference frame!
	virtual Vec3 GetG_int(void) const {
		pedantic_cerr("ElemGravityOwner(" << GetLabel() << "): "
			"warning, using default GetG_int()" << std::endl);
		return ::Zero3;
	};

public:
	ElemGravityOwner(unsigned int uL, flag fOut);
	virtual ~ElemGravityOwner(void);

	/* Usata per inizializzare la quantita' di moto */
	virtual void
	SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph = 0) = 0;

	/*
	 * massa
	 */
	virtual doublereal dGetM(void) const {
		return 0.;
	};

	Vec3 GetS(void) const {
		return GetS_int();
	};

	Mat3x3 GetJ(void) const {
		return GetJ_int();
	};

	Vec3 GetB(void) const {
		return GetB_int();
	};

	Vec3 GetG(void) const {
		return GetG_int();
	};

#ifdef DEBUG
	virtual flag fIsElemGravityOwner(void) const {
		return flag(1);
	};
#endif // DEBUG
};

/* ElemGravityOwner - end */

extern Elem *
ReadGravity(DataManager* pDM, MBDynParser& HP);

#endif // GRAVITY_H
