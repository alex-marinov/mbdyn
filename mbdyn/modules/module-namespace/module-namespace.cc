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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "mathp.h"
#include "dataman.h"

class TableNameSpace: public MathParser::NameSpace {
protected:
	Table *m_pTable;

public:
	TableNameSpace(const std::string& sName);
	~TableNameSpace(void);

	bool IsFunc(const std::string& fname) const;
	MathParser::MathFunc_t* GetFunc(const std::string& fname) const;
	TypedValue EvalFunc(MathParser::MathFunc_t *f,
		const MathParser::MathArgs& args) const;
	virtual Table* GetTable(void);
};

TableNameSpace::TableNameSpace(const std::string& sName)
: MathParser::NameSpace(sName)
{
	m_pTable = new Table(false);
}

TableNameSpace::~TableNameSpace(void)
{
	delete m_pTable;
}

bool
TableNameSpace::IsFunc(const std::string& fname) const
{
	return false;
}

MathParser::MathFunc_t*
TableNameSpace::GetFunc(const std::string& fname) const
{
	return 0;
}

TypedValue 
TableNameSpace::EvalFunc(MathParser::MathFunc_t *f,
	const MathParser::MathArgs& args) const
{
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

Table*
TableNameSpace::GetTable(void)
{
	return m_pTable;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	MBDynParser	*pHP = (MBDynParser *)php;

	int rc = 0;
	while (pHP->IsArg()) {
		const char *sName = pHP->GetString();

		/* registers namespace */
		MathParser::NameSpace *pNS = new TableNameSpace(sName);
		rc = pHP->GetMathParser().RegisterNameSpace(pNS);
		if (rc != 0) {
			delete pNS;
		}
	}

	return rc;
}

