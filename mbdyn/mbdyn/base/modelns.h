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

#ifndef MODELNS_H
#define MODELNS_H

#include "mathp.h"
#include "dataman.h"

class ModelNameSpace: public MathParser::NameSpace {
protected:
	const DataManager *pDM;

public:
	typedef MathParser::MathArgPriv_t<const DataManager *> MathArgDM;
	typedef MathParser::MathArgPriv_t<const BasicScalarFunction *> MathArgSF;
	typedef MathParser::MathArgPriv_t<Node::Type> MathArgNode;
	typedef MathParser::MathArgPriv_t<Elem::Type> MathArgElem;

	typedef std::map<std::string, MathParser::MathFunc_t *> funcType;

protected:
	funcType func;

	MathParser::MathFunc_t sf_func;
	MathParser::MathFunc_t node_func;
	MathParser::MathFunc_t elem_func;
	MathParser::MathFunc_t unique_elem_func;

public:
	ModelNameSpace(const DataManager *pDM);
	~ModelNameSpace(void);

	bool IsFunc(const char* const s) const;
	MathParser::MathFunc_t* GetFunc(const char* const s) const;
	TypedValue EvalFunc(MathParser::MathFunc_t *f,
		const MathParser::MathArgs& args) const;
};

#endif /* MODELNS_H */

