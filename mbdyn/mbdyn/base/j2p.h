/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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

#ifndef J2P_H
#define J2P_H

#include "node.h"
#include "elem.h"

/* Elem2Param - begin */

class Elem2Param : public ParameterNode {
protected:   
	const Elem* pElem;  
	unsigned int iNum;

public:
	Elem2Param(unsigned int uL, const DofOwner* pDO, flag fOut);
	virtual ~Elem2Param(void);
   
	virtual void Bind(const Elem* pEl, unsigned int i);

	/* Contributo del nodo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;
	std::ostream& RestartBind(std::ostream& out) const;

	/* Restituisce il valore del dof iDof;
	 * se differenziale, iOrder puo' essere = 1 per la derivata */
	virtual inline const doublereal& 
	dGetDofValue(int iDof, int iOrder = 0) const;

	/* virtual void SetX(const doublereal& d); */
	virtual inline const doublereal& dGetX(void) const;

	/* Setta il valore del dof iDof a dValue;
	 * se differenziale, iOrder puo' essere = 1 per la derivata */
	virtual void SetDofValue(const doublereal& dValue,
			unsigned int iDof, 
			unsigned int iOrder = 0);

	virtual void SetValue(DataManager *pDM,
			VectorHandler&, VectorHandler&,
			SimulationEntity::Hints *ph = 0);
};


/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
inline const doublereal&
Elem2Param::dGetDofValue(int iDof, int iOrder) const
{
	ASSERT(iDof == 1);
	ASSERT(iOrder == 0);
	return dGetX();
}

/* Restituisce il valore del dof */
inline const doublereal&
Elem2Param::dGetX(void) const
{
	/* element could be undefined (yet) */
	if (pElem != 0) {
		dX = pElem->dGetPrivData(iNum);
	}
	return dX;
}

/* Elem2Param - end */


/* StrainGageParam - begin */

class StrainGageParam : public Elem2Param {
protected:  
	doublereal dY;
	doublereal dZ;
   
public:
	StrainGageParam(unsigned int uL, const DofOwner* pDO,
     			doublereal dy, doublereal dz, flag fOut);
	virtual ~StrainGageParam(void);

	virtual void Bind(const Elem* pEl, unsigned int i);

	/* Contributo del nodo al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	/* Restituisce il valore del dof iDof;
	 * se differenziale, iOrder puo' essere = 1 per la derivata */
	virtual inline const doublereal& 
	dGetDofValue(int iDof, int iOrder = 0) const;

	/* virtual void SetX(const doublereal& d); */
	virtual inline const doublereal& dGetX(void) const;   
};


/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
inline const doublereal&
StrainGageParam::dGetDofValue(int iDof, int iOrder) const
{
	ASSERT(iDof == 1);
	ASSERT(iOrder == 0);
	return dGetX();
}


/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
inline const doublereal&
StrainGageParam::dGetX(void) const
{  
	unsigned int i = 6 * (iNum - 1);

	dX = pElem->dGetPrivData(i + 1)
		+ dZ*pElem->dGetPrivData(i + 5)
		- dY*pElem->dGetPrivData(i + 6);

	return dX;
}

/* StrainGageParam - end */

#endif /* J2P_H */

