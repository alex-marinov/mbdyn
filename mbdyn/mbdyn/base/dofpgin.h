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

#ifndef DOFPGIN_H
#define DOFPGIN_H

#include "mathp.h"
#include "dataman.h"

class DofPlugIn : public MathParser::PlugIn {
protected:
	ScalarDof dof;
	DataManager *pDM;
	bool bPrev;

public:
	DofPlugIn(MathParser& mp, DataManager *pDM);
	~DofPlugIn(void);
	const char *sName(void) const;
	int Read(int argc, char *argv[]);
	TypedValue::Type GetType(void) const;
	TypedValue GetVal(void) const;

private:
	unsigned int ReadLabel(const char* s);
	Node *ReadNode(unsigned int uLabel, const char *s);
	unsigned int ReadIndex(Node *pNode, unsigned int iMaxIndex,
			const char *s);
	int ReadDofOrder(Node *pNode, unsigned int iIndex, const char *s);
};

extern MathParser::PlugIn *
dof_plugin(MathParser& mp, void *arg);

#endif /* DOFPGIN_H */

