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

/*
 AUTHOR: Reinhard Resch <reinhard.resch@accomp.it>
        Copyright (C) 2011(-2012) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include <cmath>
#include <cfloat>
#include <string>

#define real mbdyn_real_type
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#include <dataman.h>
#include <drive.h>
#undef real

#include "module-octave.h"

#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/toplev.h>
#include <octave/octave.h>

class OctaveDriveCaller : public DriveCaller {
public:
	explicit OctaveDriveCaller(const std::string& expr, const DataManager* pDM, MBDynParser& HP);
	virtual ~OctaveDriveCaller(void);
 
	/* Copia */
	virtual DriveCaller* pCopy(void) const;
 
	/* Scrive il contributo del DriveCaller al file di restart */   
	virtual std::ostream& Restart(std::ostream& out) const;
 
	inline doublereal dGet(const doublereal& dVar) const;
	inline doublereal dGet(void) const;

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const;
	virtual doublereal dGetP(const doublereal& dVar) const;
	virtual inline doublereal dGetP(void) const;
private:
    inline doublereal EvalExpr(doublereal dVar)const;
private:
    const std::string expr;
    const DataManager* pDM;
    MBDynParser& HP;
};

OctaveDriveCaller::OctaveDriveCaller(const std::string& expr, const DataManager* pDM, MBDynParser& HP)
: DriveCaller(0),
  expr(expr),
  pDM(pDM),
  HP(HP)
{

}

OctaveDriveCaller::~OctaveDriveCaller(void)
{

}

DriveCaller *
OctaveDriveCaller::pCopy(void) const
{
	return new OctaveDriveCaller(expr, pDM, HP);
}

std::ostream&
OctaveDriveCaller::Restart(std::ostream& out) const
{
	return out << "octave, \"" << expr << "\"";
}

inline doublereal 
OctaveDriveCaller::dGet(const doublereal& dVar) const 
{
    return EvalExpr(dVar);
}

inline doublereal
OctaveDriveCaller::dGet(void) const
{
	return EvalExpr(pDM->dGetTime());
}

inline bool
OctaveDriveCaller::bIsDifferentiable(void) const
{
	return false;
}

inline doublereal 
OctaveDriveCaller::dGetP(const doublereal&) const
{
	return 0.;
}

inline doublereal 
OctaveDriveCaller::dGetP(void) const
{
	return 0.;
}

inline doublereal
OctaveDriveCaller::EvalExpr(doublereal dVar)const
{
    octave_value_list args;

    args.append(octave_value(dVar));

	octave_value_list ans = feval(expr, args, 1);

	if ( error_state )
	{
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

    if ( ans.length() < 1 )
    {
        silent_cerr("octave error: expression \"" << expr << "\" returned no value" << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

	if ( !ans(0).is_defined() )
	{
        silent_cerr("octave error: result of expression \"" << expr << "\" is undefined" << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

    if ( !ans(0).is_scalar_type() )
    {
        silent_cerr("octave error: result of expression \"" << expr << "\" is not a scalar value" << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

	return ans(0).scalar_value(); 
}

/* prototype of the functional object: reads a drive caller */
struct OctaveDCR : public DriveCallerRead {
    OctaveDCR()
        :bFirstTime(true)
    {

    }

    virtual ~OctaveDCR()
    {
        if ( !bFirstTime )
            do_octave_atexit();
    }

	virtual DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
    {
        if ( bFirstTime )
        {
            bFirstTime = false;

            int argc = 1;
            char cmd[] = "octave";
            char* argv[] = { cmd, 0 };

            octave_main(argc, argv, 1);

            const Table& symbolTable = pDM->GetMathParser().GetSymbolTable();

            for (Table::VM::const_iterator it = symbolTable.begin(); it != symbolTable.end(); ++it)
            {
            	const std::string& mbName(it->first);
            	const NamedValue*const namedValue = it->second;
            	const TypedValue mbValue(namedValue->GetVal());
            	octave_value octValue;

            	switch ( mbValue.GetType() )
            	{
					case TypedValue::VAR_BOOL:
						octValue = mbValue.GetBool();
						break;
					case TypedValue::VAR_INT:
						octValue = mbValue.GetInt();
						break;
					case TypedValue::VAR_REAL:
						octValue = mbValue.GetReal();
						break;
					case TypedValue::VAR_STRING:
						octValue = mbValue.GetString();
						break;
					default:
						silent_cerr("octave error: data type of variable \"" << mbName << "\": not handled in switch statement " << mbValue.GetType() << std::endl);
						ASSERT(0);
            	}

            	if ( octValue.is_defined() )
            		set_global_value(mbName, octValue);
            }
        }

		const std::string expr(HP.GetStringWithDelims(HighParser::DOUBLEQUOTE));

		if ( HP.IsKeyWord("octave" "search" "path") )
		{
			while ( HP.IsStringWithDelims(HighParser::DOUBLEQUOTE) )
			{
				const std::string path(HP.GetStringWithDelims(HighParser::DOUBLEQUOTE));
				feval("addpath", octave_value(path));

				if ( error_state )
				{
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}
			}
		}

		return new OctaveDriveCaller(expr, pDM, HP);
	};
private:
    bool bFirstTime;
};


bool
octave_set(void)
{
	DriveCallerRead	*rf = new OctaveDCR;

	if (!SetDriveData("octave", rf)) {
		delete rf;

        return false;
	}

	return true;
}

#ifndef STATIC_MODULES

extern "C" 
{

int
module_init(const char *module_name, void *pdm, void *php)
{
	if (!octave_set()) {
		silent_cerr("octave: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}

} // extern "C"

#endif // ! STATIC_MODULES
