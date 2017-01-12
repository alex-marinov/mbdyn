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

/* Classe di drivers che usano gradi di liberta' nodali
 * 
 * Il valore del drive viene fatto dipendere da un grado di liberta' nodale
 * in modo del tutto generale e trasparente. Ovvero in ingresso dati si
 * associa il drive al grado di liberta' j-esimo del nodo i-esimo, con 
 * verifiche di consistenza. In esecuzione, la funzione propria dGet() del 
 * drive restituisce il valore del grado di liberta' associato 
 * (o della derivata, se il grado di liberta' e' differenziale).
 */

#ifndef GINACDRIVE_H
#define GINACDRIVE_H

#include "drive.h"
#include <ginac/ginac.h>

class GiNaCDriveCaller : public DriveCaller
{
private:
	// parameter symbols
	GiNaC::symbol * gVar;

	// expression
	GiNaC::ex gExpr;

	// derivative
	GiNaC::ex gExprDVar;

public:
	GiNaCDriveCaller(const DriveHandler* pDH,
		const std::string& var, const std::string& expression);
	virtual ~GiNaCDriveCaller(void);

	/* Copia */
	virtual DriveCaller* pCopy(void) const;
   
	virtual std::ostream& Restart(std::ostream& out) const;
   
	inline doublereal dGet(const doublereal& dVar) const;

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const;
	virtual doublereal dGetP(const doublereal& dVar) const;
};

inline doublereal
GiNaCDriveCaller::dGet(const doublereal& dVar) const
{
	GiNaC::lst l;

	l.append(*gVar == dVar);

	GiNaC::ex f_expr = gExpr.subs(l);

	return GiNaC::ex_to<GiNaC::numeric>(f_expr).to_double();
}

struct GiNaCDCR : public DriveCallerRead {
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

#endif // GINACDRIVE_H

