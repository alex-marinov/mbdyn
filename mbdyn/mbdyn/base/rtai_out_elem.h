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

/* socket driver */

#ifndef RTAI_MBOX_ELEM_H
#define RTAI_MBOX_ELEM_H

#ifdef USE_RTAI

/* include del programma */

#include <elem.h>
#include <node.h>

/* RTAIMailboxElem - begin */

class RTAIMailboxElem : virtual public Elem {
protected:
	unsigned int NumChannels;
   	ScalarDof* pNodes;
   
public:
   	RTAIMailboxElem(unsigned int uL, unsigned int nmb, ScalarDof *& pn);
   	virtual ~RTAIMailboxElem(void);

	virtual inline void* pGet(void) const;

	virtual std::ostream& Restart(std::ostream& out) const;
	virtual Elem::Type GetElemType(void) const;
	virtual void WorkSpaceDim(integer* piRows, integer* piCols) const;
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WrokVec, double dCoef,
			const VectorHandler& X, const VectorHandler& XP);
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat, double dCoef,
			const VectorHandler& X, const VectorHandler& XP);

	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);
};

inline void*
RTAIMailboxElem::pGet(void) const
{ 
	return (void *)this;
}

/* RTAIMailboxElem - end */

#endif /* USE_RTAI */

#endif /* RTAI_MBOX_ELEM_H */

