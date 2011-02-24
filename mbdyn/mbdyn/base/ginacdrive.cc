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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "dataman.h"
#include "ginacdrive.h"


GiNaCDriveCaller::GiNaCDriveCaller(const DriveHandler* pDH, 
	const std::string& expression)
: DriveCaller(pDH)
{
	GiNaC::lst l; 

	gVar = new GiNaC::symbol("Var");
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

	std::ostringstream expr;
	expr << gExpr;

	SAFENEWWITHCONSTRUCTOR(pDC, 
		GiNaCDriveCaller,
		GiNaCDriveCaller(pDrvHdl, expr.str()));
 
	return pDC;
}

/* Restart */
std::ostream&
GiNaCDriveCaller::Restart(std::ostream& out) const
{
	out << " ginac, \"" << gExpr << "\"";

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

	DriveCaller *pDC = 0;

	const char *expression = HP.GetStringWithDelims();
	if (expression == 0) {
		silent_cerr("unable read ginac drive caller expression at line "
			<< HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	SAFENEWWITHCONSTRUCTOR(pDC,
		GiNaCDriveCaller,
		GiNaCDriveCaller(pDrvHdl, expression));

	return pDC;
}

