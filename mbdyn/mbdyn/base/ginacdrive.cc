/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

#include "dataman.h"
#include "ginacdrive.h"

#include <sstream>

GiNaCDriveCaller::GiNaCDriveCaller(const DriveHandler* pDH, 
	const std::string& var, const std::string& expression)
: DriveCaller(pDH)
{
	GiNaC::lst l; 

	gVar = new GiNaC::symbol(var);
	l.append(*gVar);

	try {
		gExpr = GiNaC::ex(expression, l);

	} catch (std::exception e) {
		silent_cerr("GiNaCDriveCaller: "
			"expression \"" << expression << "\" parsing failed: " << e.what()
			<< std::endl);
		throw e;
	}

	try {
		gExprDVar = gExpr.diff(*gVar);

	} catch (std::exception e) {
		silent_cerr("GiNaCDriveCaller: "
			"expression \"" << expression << "\" differetiation wrt/ \"Var\" failed: " << e.what()
			<< std::endl);
		throw e;
	}

	pedantic_cout("GiNacDriveCaller: symbol=\"" << *gVar << "\" expression=\"" << gExpr << "\" derivative=\"" << gExprDVar << "\"" << std::endl);
}

GiNaCDriveCaller::~GiNaCDriveCaller(void)
{
	NO_OP;
}

/* Copia */
DriveCaller*
GiNaCDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = NULL;

	std::ostringstream var, expr;
	var << gVar;
	expr << gExpr;

	SAFENEWWITHCONSTRUCTOR(pDC, 
		GiNaCDriveCaller,
		GiNaCDriveCaller(pDrvHdl, var.str(), expr.str()));
 
	return pDC;
}

/* Restart */
std::ostream&
GiNaCDriveCaller::Restart(std::ostream& out) const
{
	out << " ginac, symbol, \"" << *gVar << "\", \"" << gExpr << "\"";

	return out;
}

bool
GiNaCDriveCaller::bIsDifferentiable(void) const
{
	return true;
}

doublereal
GiNaCDriveCaller::dGetP(const doublereal& dVar) const
{
	GiNaC::lst l;

	l.append(*gVar == dVar);

	GiNaC::ex f_expr = gExprDVar.subs(l);

	return GiNaC::ex_to<GiNaC::numeric>(f_expr).to_double();
}

DriveCaller *
GiNaCDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	const DriveHandler* pDrvHdl = 0;
	if (pDM != 0) {
		pDrvHdl = pDM->pGetDrvHdl();
	}

	std::string var;
	if (HP.IsKeyWord("symbol")) {
		var = HP.GetStringWithDelims();
		if (var.empty()) {
			silent_cerr("unable to read ginac drive caller symbol at line "
				<< HP.GetLineData() << std::endl);
			throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

	} else {
		var = "Var";
	}

	std::string expression = HP.GetStringWithDelims();
	if (expression.empty()) {
		silent_cerr("unable to read ginac drive caller expression at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	DriveCaller *pDC = 0;
	SAFENEWWITHCONSTRUCTOR(pDC,
		GiNaCDriveCaller,
		GiNaCDriveCaller(pDrvHdl, var, expression));

	return pDC;
}

