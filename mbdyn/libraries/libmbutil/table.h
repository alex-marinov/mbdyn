/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

/* Tabella dei simboli; richiamata da table.cc, parser.h */

#ifndef TABLE_H
#define TABLE_H

#include <cstring>
#include <map>

#include "except.h"
#include "mathtyp.h"

class Table {
	friend std::ostream& operator << (std::ostream& out, Table& T);

public:
	class ErrNameAlreadyDefined : public MBDynErrBase {
	public:
		ErrNameAlreadyDefined(MBDYN_EXCEPT_ARGS_DECL)
			: MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

	typedef std::map<std::string, NamedValue *> VM;

private:
	VM vm;

public:
	Table(bool bSetConstants);
	virtual ~Table(void);
	Var* Put(const char* const name, const TypedValue& v);
	NamedValue* Put(NamedValue* p);
	NamedValue* Get(const char* const name) const;
	VM::const_iterator begin(void) const { return vm.begin(); };
	VM::const_iterator end(void) const { return vm.end(); };
};

extern std::ostream& operator << (std::ostream& out, Table& T);

#endif // TABLE_H

