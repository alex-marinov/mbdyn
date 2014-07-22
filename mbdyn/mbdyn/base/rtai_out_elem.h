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

/* socket driver */

#ifndef RTAI_OUT_ELEM_H
#define RTAI_OUT_ELEM_H

/* include del programma */

#include <streamoutelem.h>

/* RTMBDynOutElem - begin */

class RTMBDynOutElem : public StreamOutElem, virtual public Elem {
protected:
	/* FIXME: store restart info as well */
	std::string host;
	unsigned long node;
	bool create;
	int port;
	bool bNonBlocking;

	StreamContent *pSC;

	void *mbx;
	int (*f_send)(unsigned long node, int port, void *v_mbx,
		void *msg, int msg_size);
   
public:
   	RTMBDynOutElem(unsigned int uL, const std::string& m,
			const std::string& host, unsigned long n, bool c,
			StreamContent *pSC, bool bNonBlocking);
   	virtual ~RTMBDynOutElem(void);

	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP);

	/* Inverse Dynamics: */
	virtual void AfterConvergence(const VectorHandler& X, 
			const VectorHandler& XP,
			const VectorHandler& XPP);
};

#endif /* RTAI_OUT_ELEM_H */

