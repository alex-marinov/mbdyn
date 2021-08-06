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

#ifndef USERELEM_H
#define USERELEM_H

#include "elem.h"
#include "dataman.h"
#include "gravity.h"
#include "aerodyn.h"

// base class for user-defined elements
class UserDefinedElem :
	virtual public Elem,
	public InitialAssemblyElem,
	public AerodynamicElem,
	public ElemGravityOwner
{
protected:
	bool needsAirProperties;

public:
   	UserDefinedElem(unsigned uLabel, const DofOwner* pDO);
   	virtual ~UserDefinedElem(void); 

	bool NeedsAirProperties(void) const;
	void NeedsAirProperties(bool yesno);

   	virtual Elem::Type GetElemType(void) const;
   	virtual AerodynamicElem::Type GetAerodynamicElemType(void) const;

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

class MBDynParser;
class DataManager;

// base class for user-defined element parsing
struct UserDefinedElemRead {
	virtual ~UserDefinedElemRead( void ) { NO_OP; };
	virtual UserDefinedElem *
	Read(unsigned uLabel, const DofOwner* pDO,
		DataManager* const pDM, MBDynParser& HP) const = 0;
};

// base class for user-defined element parsing
template <class UDE>
struct UDERead : public UserDefinedElemRead {
	virtual ~UDERead( void ) { NO_OP; };
	virtual UserDefinedElem *
	Read(unsigned uLabel, const DofOwner* pDO,
		DataManager* const pDM, MBDynParser& HP) const;
};

template <class UDE>
UserDefinedElem *
UDERead<UDE>::Read(unsigned uLabel, const DofOwner* pDO,
	DataManager* const pDM, MBDynParser& HP) const
{
	return new UDE(uLabel, pDO, pDM, HP);
}

// "public"
extern bool
SetUDE(const std::string &s, UserDefinedElemRead *rude);

// "private"
extern void InitUDE(void);
extern void DestroyUDE(void);
extern UserDefinedElem *
ParseUserDefinedElem(unsigned uLabel, DofOwner* pDO,
	DataManager* const pDM, MBDynParser& HP);

#endif // USERELEM_H

