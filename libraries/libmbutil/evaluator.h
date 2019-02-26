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

/*
 * With the contribution of Ankit Aggarwal <ankit.ankit.aggarwal@gmail.com>
 * during Google Summer of Code 2015
 */

#ifndef EVALUATOR_H
#define EVALUATOR_H

#include "mathtyp.h"

class ExpressionElement {
public:
	enum EEFlags {
		EE_NONE = 0x0U,

		EE_CONSTIFY = 0x1U,

		EE_OPTIMIZE = EE_CONSTIFY
	};

protected:
	static unsigned m_uEEFlags;

public:
	virtual ~ExpressionElement(void) {};
#if 0 // TODO: check correctness whenever possible
	virtual bool Check(void) const = { return true; };
#endif
	virtual TypedValue Eval(void) const = 0;
	virtual std::ostream& Output(std::ostream& out) const = 0;

	static unsigned GetFlags(void) { return m_uEEFlags; };
	static void SetFlag(EEFlags f) { m_uEEFlags |= unsigned(f); };
	static void ClearFlag(EEFlags f) { m_uEEFlags &= !unsigned(f); };
	static bool IsFlag(EEFlags f) { return m_uEEFlags & unsigned(f); };
};

extern std::string EEStrOut(const ExpressionElement *e);

extern bool EE_Eval(TypedValue& dst, const ExpressionElement *ee);
extern bool EE_Eval(bool& dst, const ExpressionElement *ee);
extern bool EE_Eval(Int& dst, const ExpressionElement *ee);
extern bool EE_Eval(Real& dst, const ExpressionElement *ee);
extern bool EE_Eval(std::string& dst, const ExpressionElement *ee);

template <class T> bool EE_Eval(T& dst, const ExpressionElement *ee) { ASSERT(ee == 0); return false; };

#endif // EVALUATOR_H
