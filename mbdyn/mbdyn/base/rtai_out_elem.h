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

#ifndef RTAI_OUT_ELEM_H
#define RTAI_OUT_ELEM_H

#ifdef USE_RTAI

/* include del programma */

#include <elem.h>
#include <node.h>

/* RTAIOutElem - begin */

class RTAIOutElem : virtual public Elem {
protected:
	unsigned int NumChannels;
   	ScalarDof* pNodes;

	/* MBox buffer */
	int size;
	char *buf;

	/* FIXME: store restart info as well */
	const char *host;
	unsigned long node;
	const char *name;
	bool create;
	int port;
	void *mbx;
   
public:
   	RTAIOutElem(unsigned int uL, unsigned int nmb, ScalarDof *& pn,
			const char *host, const char *m, unsigned long n,
			bool c);
   	virtual ~RTAIOutElem(void);

	virtual inline void* pGet(void) const;

	virtual std::ostream& Restart(std::ostream& out) const;
	virtual Elem::Type GetElemType(void) const;
	virtual void WorkSpaceDim(integer* piRows, integer* piCols) const;
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WrokVec, doublereal dCoef,
			const VectorHandler& X, const VectorHandler& XP);
	virtual VariableSubMatrixHandler& 
	AssJac(VariableSubMatrixHandler& WorkMat, doublereal dCoef,
			const VectorHandler& X, const VectorHandler& XP);

	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);
};

inline void*
RTAIOutElem::pGet(void) const
{ 
	return (void *)this;
}

class DataManager;
class MBDynParser;

extern Elem *
ReadRTAIOutElem(DataManager *pDM, MBDynParser& HP, unsigned int uLabel);

/* RTAIOutElem - end */

#endif /* USE_RTAI */

#endif /* RTAI_OUT_ELEM_H */

