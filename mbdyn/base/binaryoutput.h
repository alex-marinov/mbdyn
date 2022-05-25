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

#ifndef BINARY_OUTPUT_H
#define BINARY_OUTPUT_H

//#include "output.h"
#include<string>
#include<cstddef>
#include<vector>

#include "myassert.h"
#include "matvec3.h"
#include "units.h"

enum class MBDynOutType {
	OutInt=1,
	OutDouble=2,
	OutChar=3
};

class BinaryOutput {
protected:
	size_t DimTime_id, DimV1_id, DimV3_id;
public:
	struct AttrVal {
		std::string attr;
		std::string val;
		AttrVal(void) { NO_OP; };
		AttrVal(const std::string& attr, const std::string& val) : attr(attr), val(val) { NO_OP; };
	};
	typedef std::vector<AttrVal> AttrValVec;

	BinaryOutput(void) :
		Start1(1,0),  // must initialize vectors otherwise can't assign
		Count1(1,1),
		Start1x3(2,0),
		Count1x3(2,1),
		Start1x3x3(3,0),
		Count1x3x3(3,1) {
	};
	
	virtual ~BinaryOutput(void);
	virtual void Open(const int format) = 0;
// 	virtual const int CreateDim(const std::string& name, integer size = -1) = 0;
	
	virtual bool isOpen(void) = 0;
	virtual void close(void) = 0;
	
	virtual void sync(void) = 0;

// 	virtual const int GetDim(const std::string& name) const = 0;

// 	inline int DimTime(void) const;
// 	inline int DimV1(void) const;
// 	inline int DimV3(void) const;

	std::vector<size_t> Start1;

	std::vector<size_t> Count1;   
	std::vector<size_t> Start1x3;
	std::vector<size_t> Count1x3;
	std::vector<size_t> Start1x3x3;
	std::vector<size_t> Count1x3x3;

	virtual void SetCurrentStep(long step) = 0;

	virtual const size_t CreateDim(const std::string& name, integer size = -1) = 0;

	const size_t DimTime() {return DimTime_id;};
	const size_t DimV1() {return DimV1_id;};
	const size_t DimV3() {return DimV3_id;};
	
	virtual void
	WriteVar(const size_t, const doublereal&) = 0;

	virtual void
	WriteVar(const size_t, const long&) = 0;

	virtual void
	WriteVar(const size_t, const int&) = 0;

	virtual void
	WriteVar(const size_t, const Vec3&) = 0;

	virtual void
	WriteVar(const size_t, const Vec3&, const size_t&) = 0;

	virtual void
	WriteVar(const size_t, const Mat3x3&) = 0;

	virtual void
	WriteVar(const size_t, const Mat3x3&, const size_t&) = 0;
	
	virtual void
	WriteVar(const size_t, const int&, const std::vector<size_t>& start, const std::vector<size_t>& count) = 0;
	
	virtual void
	WriteVar(const size_t, const long int&, const std::vector<size_t>& start, const std::vector<size_t>& count) = 0;
	
	virtual void
	WriteVar(const size_t, const unsigned&, const std::vector<size_t>& start, const std::vector<size_t>& count) = 0;
	
	virtual void
	WriteVar(const size_t, const long unsigned&, const std::vector<size_t>& start, const std::vector<size_t>& count) = 0;
	
	virtual void
	WriteVar(const size_t, const double&, const std::vector<size_t>& start, const std::vector<size_t>& count) = 0;
	
	virtual void
	PutVar(const size_t var, const doublereal* data) = 0;

// 	template <class Tvar>
// 	void
// 	WriteVar(const size_t, const Tvar&);
// 	
// 	template <class Tvar, class Tstart>
// 	virtual void
// 	WriteVar(const size_t, const Tvar&, const Tstart&);
// 
// 	template <class Tvar, class Tstart>
// 	void
// 	WriteVar(const size_t, const Tvar&, 
// 			const std::vector<Tstart>&, 
// 			const std::vector<size_t>& = std::vector<size_t>(1,1));
	
	size_t
	CreateVar(const std::string& name, const std::string& type);

	template <class T>
	size_t
	CreateVar(const std::string& name,
		const MBUnits::Dimensions phys_dim, const std::string& description);

	size_t
	CreateVar(const std::string& name, const MBDynOutType& type,
		const AttrValVec& attrs, const std::vector<size_t>& dims);

	size_t
	CreateRotationVar(const std::string& name_prefix,
		const std::string& name_postfix,
		OrientationDescription od,
		const std::string& description);
};


#endif //BINARY_OUTPUT_H
